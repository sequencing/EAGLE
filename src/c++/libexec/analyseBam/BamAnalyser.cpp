/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <fstream>
#include <boost/iostreams/device/file.hpp>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "genome/SharedFastaReference.hh"
#include "genome/BamMetadata.hh"
#include "io/BamIndexer.hh"
#include "model/Nucleotides.hh"
#include "BamAnalyser.hh"


using namespace std;
//#define DEBUG


namespace eagle
{
namespace main
{

BamAnalyser::BamAnalyser (const BamAnalyserOptions &options )
    : options_( options )
{
    genome::SharedFastaReference::init( options_.referenceGenome );
}

void BamAnalyser::run()
{
    // Create a boost iostream flow that takes a BAM file as input, sends it to our BAM parser, which will use a callback for each read
    bool bamIndexer  = 0;
    bool bamAnalyser = 0;
    bool mismatchTableGenerator = 1;

    if (bamIndexer)
    {
        // Add BAM Indexer
        const boost::filesystem::path baiPath( "out.bam.bai" );
        boost::iostreams::file_sink baiSink(baiPath.string());
        bgzfStream_.push( eagle::io::bam::BamIndexer<boost::iostreams::file_sink>(baiSink) );
    }

    if (bamAnalyser)
    {
        bgzfStream_.push( BamReadDumper() );
    }

    if (mismatchTableGenerator)
    {
        vector<bool> requestedMetrics( 64 ); // max 64 metrics
        vector<bool> requestedTables ( 64 ); // max 64 tables
        for (int i=std::min(requestedMetrics.size(),options_.requestedMetrics.size())-1; i>=0; --i)
        {
            requestedMetrics[i] = (options_.requestedMetrics[i] == '1');
        }
        for (int i=std::min(requestedTables.size(),options_.requestedTables.size())-1; i>=0; --i)
        {
            requestedTables[i] = (options_.requestedTables[i] == '1');
        }
        bgzfStream_.push( MetricsComputer( requestedMetrics, requestedTables ) );
    }

    // Terminate the iostreams chain
    bgzfStream_.push( boost::iostreams::null_sink() );

    // Send the BAM data to the stream
    ifstream inp( options_.bamFile.string().c_str() );
    bgzfStream_.exceptions( ios_base::badbit ); // Required to rethrow any error thrown by our iostreams filters (leading to better error messages)
    while (inp.good())
    {
        char buf[128];
        inp.read( buf, sizeof(buf) );
        bgzfStream_.write( buf, inp.gcount() );
    }
}


MetricsComputer::MetricsComputer( vector<bool>& requestedMetrics, vector<bool>& requestedTables )
    : requestedMetrics_( requestedMetrics )
    , requestedTables_ ( requestedTables )
    , mismatches_           ( 256, vector< unsigned int >( 256, 0 ) )
    , homopolymerDeletions_ ( 102, vector< unsigned int >( 102, 0 ) )
    , homopolymerInsertions_( 102, vector< unsigned int >( 102, 0 ) )
    , homopolymerCount_     ( 102, 0 )
    , mismatchCountPerRead_ ( 102, 0 )
    , mismatchCount_ ( 0 )
    , insertionCount_( 0 )
    , deletionCount_ ( 0 )
    , baseCount_     ( 0 )
    , mismatchCountPer10kWindow_( 10000 )
    , baseCountPer10kWindow_    ( 10000 )
{
}

void MetricsComputer::parsedAlignment( const io::bam::BamAlignment& alignment, const io::bam::VirtualOffset& virtualOffset, const io::bam::VirtualOffset& virtualEndOffset )
{
    model::IUPAC baseConverter;
    unsigned long mismatchCountForThisRead = 0;

    static unsigned int lastRefId = 1337;
    if (lastRefId != alignment.refId)
    {
        cout << "Now reading refId " << alignment.refId << endl;
        lastRefId = alignment.refId;
    }
//    if (alignment.refId != 0) { return; }

#ifdef DEBUG
    cout << "Read \"" << alignment.getReadName() << endl;
    cout << " Pos   = (" << alignment.refId << "):" << alignment.pos << endl;

    // Print CIGAR
    cout << " CIGAR = ";
    const unsigned int *cigar = alignment.getCigar();
    const unsigned int cigarLength = alignment.getNCigarOp();
    for (unsigned int i=0; i<cigarLength; ++i)
    {
        const string letters( "MIDNSHP=X" );
        const unsigned int opLen = cigar[i] >> 4;
        unsigned int letterNum = cigar[i] & 0xF;
        assert( letterNum < 9 );
        const char letter = letters[letterNum];
        cout << opLen << letter << ",";
    }
    cout << endl;

    // Print read seq
    cout << " Seq   = ";
    const unsigned char *seq = alignment.getSeq();
    unsigned int seqLength = alignment.lSeq;
    for (unsigned int i=0; i<seqLength; ++i)
    {
        char binBase = (i%2==0)?( seq[i/2] >> 4 ):( seq[i/2] & 0xF );
        char base = baseConverter.binToIupac( binBase );
        cout << base;
    }
    cout << endl;

    // Print ref
    cout << " Ref   = ";
    unsigned long globalPos = alignment.pos;
    bool overlapContigBoundary;
    for (unsigned int i=0; i<seqLength; ++i)
    {
        char base = genome::SharedFastaReference::get()->get( globalPos, i, overlapContigBoundary );
        cout << base;
    }
    cout << endl;

    // Print ruler
    cout << " Ruler : ";
    for (unsigned int i=0; i<seqLength; ++i)
    {
        cout << (i%10);
    }
    cout << endl;
    cout << "         ";
    for (unsigned int i=0; i<seqLength; ++i)
    {
        if (i%10 == 0)
        {
            cout << (i/10);
        }
        else
        {
            cout << " ";
        }
    }
    cout << endl;

#else //ifdef DEBUG

    const unsigned int *cigar = alignment.getCigar();
    const unsigned int cigarLength = alignment.getNCigarOp();
    const unsigned char *seq = alignment.getSeq();
    unsigned int seqLength = alignment.lSeq;
    unsigned long globalPos = alignment.pos;
    bool overlapContigBoundary;

#endif //ifdef DEBUG

    if (requestedMetrics_[2])
    {
        baseCount_ += seqLength;
    }
    if (requestedMetrics_[3] || requestedMetrics_[4] || requestedTables_[3])
    {
        baseCountPer10kWindow_.incRegion( globalPos, seqLength );
    }
    if (requestedTables_[6])
    {
        unsigned int mapQ = (alignment.binMqNl>>8) & 0xFF;
        if (alignment.refId == alignment.nextRefId && mapQ>200)
        {
            unsigned long insertSize = abs(static_cast<int>(alignment.pos) - static_cast<int>(alignment.nextPos)) + seqLength;
            if(insertSize >= insertSizes_.size())
            {
                insertSizes_.resize( insertSize+1 );
            }
            if(insertSize < insertSizes_.size())
            {
                insertSizes_[insertSize]++;
            }
        }
    }
    if (requestedTables_[5])
    {
        observedCoverage_.forgetBefore( globalPos );
        readsAddedForTable5_.forgetBefore( globalPos );

        observedCoverage_.incRegion( globalPos, seqLength );
        unsigned int observedCov = observedCoverage_.get( globalPos );
        unsigned int addedCov = readsAddedForTable5_.get( globalPos );

        // Discard all previous segments that finish too early to be connected to any new segment
        for (unsigned int i=0; i<table5Data_.size(); ++i)
        {
            if (table5Data_[i].second < globalPos - 5)
            {
                cout << "Table 5:\t" << table5Data_[i].first << "\t" << table5Data_[i].second << "\t" << (table5Data_[i].second - table5Data_[i].first + 1) << endl;
//                cout << "  (globalPos=" << globalPos << ", observedCov=" << observedCov << ", addedCov=" << addedCov << ")" << endl;
                table5Data_.erase(table5Data_.begin()+i);
                --i;
            }
        }

//        while ( (observedCov > 100) && (observedCov > 80 + addedCov) )
        if (observedCov > 80 + addedCov)
        {
            // Find first item that can be connected to the new segment
            bool found = false;
            for (unsigned int i=0; i<table5Data_.size(); ++i)
            {
                if (table5Data_[i].second < globalPos + 5)
                {
//                    cout << "Adding to table 5: " << globalPos << "-" << (globalPos+seqLength-1) << endl;
//                    cout << "  Extending " << table5Data_[i].first << "-" << table5Data_[i].second << " to " << (globalPos+seqLength-1) << endl;
                    table5Data_[i].second = globalPos+seqLength-1;
                    readsAddedForTable5_.incRegion( globalPos, seqLength );
                    found = true;
                    break;
                }
            }

            if (!found && (observedCov > 100+addedCov))
            {
//                cout << "Adding to table 5: " << globalPos << "-" << (globalPos+seqLength-1) << endl;
//                cout << "  observedCov=" << observedCov << " , addedCov=" << addedCov << endl;
                table5Data_.push_back( std::make_pair( globalPos, globalPos+seqLength-1 ) );
                readsAddedForTable5_.incRegion( globalPos, seqLength );
            }
        }
    }

    // Mismatches analysis
    if (requestedTables_[0] || requestedTables_[1] || requestedMetrics_[1] || requestedTables_[2] || requestedMetrics_[2] || requestedMetrics_[3] || requestedMetrics_[4] || requestedTables_[3])
    {
        unsigned int posInReadSeq = 0;
        unsigned int posInRefSeq = 0;
        for (unsigned int i=0; i<cigarLength; ++i)
        {
            const string letters( "MIDNSHP=X" );
            const unsigned int opLen = cigar[i] >> 4;
            unsigned int letterNum = cigar[i] & 0xF;
            assert( letterNum < 9 );
            const char letter = letters[letterNum];
#ifdef DEBUG
            cout << "Analysing: " << opLen << letter << endl;
#endif //ifdef DEBUG
            switch (letter)
            {
            case 'M':
                if (requestedTables_[0] || requestedMetrics_[1] || requestedTables_[2] || requestedMetrics_[2] || requestedMetrics_[3] || requestedMetrics_[4] || requestedTables_[3])
                {
                    for (unsigned int j=0; j<opLen; ++j)
                    {
                        char seqBinBase = (posInReadSeq%2==0)?( seq[posInReadSeq/2] >> 4 ):( seq[posInReadSeq/2] & 0xF );
                        char seqBase = baseConverter.binToIupac( seqBinBase );
                        char refBase = toupper( genome::SharedFastaReference::get()->get( globalPos, posInRefSeq, overlapContigBoundary ) );
                        if (seqBase != refBase)
                        {
#ifdef DEBUG
                            cout << (boost::format(" Mismatch detected!!!: %c->%c (@pos %d in ref / %d in read)") % refBase % seqBase % posInRefSeq % posInReadSeq).str() << endl;
#endif //ifdef DEBUG
                             if (requestedTables_[0])
                             {
                                 ++mismatches_[(unsigned int)refBase][(unsigned int)seqBase];
                             }
                             if (requestedMetrics_[1] || requestedTables_[2] || requestedMetrics_[2] || requestedMetrics_[3] || requestedMetrics_[4] || requestedTables_[3])
                             {
                                 ++mismatchCountForThisRead;
                                 ++mismatchCount_;
                                 mismatchCountPer10kWindow_.inc( globalPos + posInRefSeq );
                             }

                        }
                        ++posInReadSeq;
                        ++posInRefSeq;
                    }
                }
                else
                {
                    posInReadSeq += opLen;
                    posInRefSeq += opLen;
                }
                break;
            case 'D':
                if (requestedMetrics_[2])
                {
                    ++deletionCount_;
                }
                if (requestedTables_[1])
                {
                    unsigned int homopolymerLength = 1;
                    unsigned int deletionLength = opLen;

                    char refBase = toupper( genome::SharedFastaReference::get()->get( globalPos, posInRefSeq, overlapContigBoundary ) );
                    while (toupper( genome::SharedFastaReference::get()->get( globalPos, posInRefSeq+homopolymerLength, overlapContigBoundary ) )
                           == refBase)
                    {
                        ++homopolymerLength;
                    }
                    if( homopolymerLength >= homopolymerDeletions_.size() )
                    {
                        homopolymerDeletions_.resize( homopolymerLength+1 );
                    }
                    if( deletionLength < homopolymerDeletions_[0].size() )
                    {
                        if( deletionLength >= homopolymerDeletions_[homopolymerLength].size() )
                        {
                            homopolymerDeletions_[homopolymerLength].resize( deletionLength+1 );
                        }
                        ++homopolymerDeletions_[homopolymerLength][deletionLength];
                    }
#ifdef DEBUG
                    cout << (boost::format(" Homopolymer deletion detected!!!: poly%c-Length=%d, del=%d (@pos %d in ref / %d in read)") % refBase % homopolymerLength % deletionLength % posInRefSeq % posInReadSeq).str() << endl;
#endif //ifdef DEBUG
                }
                posInRefSeq += opLen;
                break;
            case 'I':
                if (requestedMetrics_[2])
                {
                    ++insertionCount_;
                }
                if (requestedTables_[1])
                {
                    unsigned int homopolymerLength = 1;
                    unsigned int insertionLength = opLen;

                    char refBase = toupper( genome::SharedFastaReference::get()->get( globalPos, posInRefSeq, overlapContigBoundary ) );
                    while (toupper( genome::SharedFastaReference::get()->get( globalPos, posInRefSeq+homopolymerLength, overlapContigBoundary ) )
                           == refBase)
                    {
                        ++homopolymerLength;
                    }
                    if( homopolymerLength >= homopolymerInsertions_.size() )
                    {
                        homopolymerInsertions_.resize( homopolymerLength+1 );
                    }
                    if( insertionLength < homopolymerInsertions_[0].size() )
                    {
                        if( insertionLength >= homopolymerInsertions_[homopolymerLength].size() )
                        {
                            homopolymerInsertions_[homopolymerLength].resize( insertionLength+1 );
                        }
                        ++homopolymerInsertions_[homopolymerLength][insertionLength];
                    }
#ifdef DEBUG
                    cout << (boost::format(" Homopolymer insertion detected!!!: poly%c-Length=%d, ins=%d (@pos %d in ref / %d in read)") % refBase % homopolymerLength % insertionLength % posInRefSeq % posInReadSeq).str() << endl;
#endif //ifdef DEBUG
                }
                posInReadSeq += opLen;
                break;
            case 'S':
                posInReadSeq += opLen;
                posInRefSeq += opLen;
                break;
            default:
                EAGLE_ERROR( (boost::format("Unexpected CIGAR letter '%c'") % letter).str() );
            }
        }
    }

    // Count homopolymers
    if (requestedTables_[1])
    {
        unsigned int homopolymerLength = 0;

        char prevBase = 'x';
        for (unsigned int i=0; i<=seqLength; ++i) {
            char base = (i<seqLength)?toupper( genome::SharedFastaReference::get()->get( globalPos, i, overlapContigBoundary ) ):'x';
            if (base == prevBase) {
                ++homopolymerLength;
            }
            else
            {
                if (homopolymerLength > 1)
                {
#ifdef DEBUG
                    cout << (boost::format(" Homopolymer size %d @%d") % homopolymerLength % i).str() << endl;
#endif //ifdef DEBUG
                    if( homopolymerLength >= homopolymerCount_.size() )
                    {
                        homopolymerCount_.resize( homopolymerLength + 1 );
                    }
                    ++homopolymerCount_[homopolymerLength];
                }
                prevBase = base;
                homopolymerLength = 1;
            }
        }
    }

    // Update mismatch count per read for histogram
    if (requestedMetrics_[1] || requestedTables_[2])
    {
        if ( mismatchCountForThisRead >= mismatchCountPerRead_.size() )
        {
            mismatchCountPerRead_.resize( mismatchCountForThisRead + 1 );
        }
        ++mismatchCountPerRead_[mismatchCountForThisRead];
    }

    // QScores count
    if (requestedTables_[7])
    {
        const char *quals = alignment.getQual();
        unsigned int seqLength = alignment.lSeq;

        int qSum = 0;
        for (unsigned int i=0; i<seqLength; ++i) {
            int qual = quals[i];
            qSum += qual;
        }
        int averageQ = qSum / seqLength;

#define MAX_QSCORE 50
        unsigned int templateNum = max(0, min(MAX_QSCORE, averageQ));
        if (qualityTable_.size() <= templateNum)
        {
            qualityTable_.resize( templateNum+1 );
        }
        if (qualityTable_[templateNum].size() < seqLength)
        {
            qualityTable_[templateNum].resize( seqLength );
            for (unsigned int i=0; i<seqLength; ++i)
                qualityTable_[templateNum][i].resize( MAX_QSCORE+1 );
        }
        for (unsigned int i=0; i<seqLength; ++i) {
            int qual = quals[i];
            qual = max(0, min(MAX_QSCORE, qual));
            qualityTable_[templateNum][i][qual]++;
        }
    }
}

void MetricsComputer::finishedParsing()
{
    if (requestedTables_[0])
    { // Print mismatches table
        for (unsigned int i=0; i<mismatches_.size(); ++i)
        {
            for (unsigned int j=0; j<mismatches_[i].size(); ++j)
            {
                if (mismatches_[i][j])
                {
                    cout << (boost::format("%c->%c\t%d") % (char)i % (char)j % mismatches_[i][j]).str() << endl;
                }
            }
        }
    }

    if (requestedTables_[1])
    { // Print homopolymer table
        // Print deletions
        for (unsigned int i=0; i<homopolymerDeletions_.size(); ++i)
        {
            for (unsigned int j=0; j<homopolymerDeletions_[i].size(); ++j)
            {
                if (homopolymerDeletions_[i][j])
                {
                    cout << (boost::format("homo %d\tdel %d\t%d\tout of %d\t= %f%%") % i % j % homopolymerDeletions_[i][j] % homopolymerCount_[i] % (homopolymerDeletions_[i][j]*100.0/homopolymerCount_[i])).str() << endl;
                }
            }
        }

        // Print insertions
        for (unsigned int i=0; i<homopolymerInsertions_.size(); ++i)
        {
            for (unsigned int j=0; j<homopolymerInsertions_[i].size(); ++j)
            {
                if (homopolymerInsertions_[i][j])
                {
                    cout << (boost::format("homo %d\tins %d\t%d\tout of %d\t= %f%%") % i % j % homopolymerInsertions_[i][j] % homopolymerCount_[i] % (homopolymerInsertions_[i][j]*100.0/homopolymerCount_[i])).str() << endl;
                }
            }
        }
    }

    if (requestedMetrics_[1] || requestedTables_[2])
    { // Print mismatch counts for histogram
        unsigned long sum = 0;
        for (unsigned int i=0; i<mismatchCountPerRead_.size(); ++i)
        {
            if (requestedTables_[2])
            {
            cout << (boost::format("mismatchCountPerRead[%d]\t%d") % i % mismatchCountPerRead_[i]).str() << endl;
            }
            sum += mismatchCountPerRead_[i];
        }
        if (requestedMetrics_[1])
        {
            unsigned long threshold1a = sum - (sum / 100000);
            cout << (boost::format("Metric 1a threshold=%d") % threshold1a).str() << endl;
            unsigned int i;
            for (i=0; i<mismatchCountPerRead_.size(); ++i)
            {
                if (threshold1a < mismatchCountPerRead_[i])
                    break;
                else
                    threshold1a -= mismatchCountPerRead_[i];
            }
            cout << (boost::format("Metric 1a: %d") % i).str() << endl;

            unsigned long threshold1b = sum / 100000;
            cout << (boost::format("Metric 1b threshold=%d") % threshold1b).str() << endl;
            for (i=0; i<mismatchCountPerRead_.size(); ++i)
            {
                if (mismatchCountPerRead_[i] < threshold1b)
                    break;
            }
            cout << (boost::format("Metric 1b: %d") % i).str() << endl;
        }
    }

    if (requestedMetrics_[2])
    { // Print mismatch rate
        double mismatchRate = (double)mismatchCount_ / (double)baseCount_;
        double insertionRate = (double)insertionCount_ / (double)baseCount_;
        double deletionRate = (double)deletionCount_ / (double)baseCount_;
        cout << "Metric 2a (mismatch rate): \t" << 100*mismatchRate  << "%" << endl;
        cout << "Metric 2b (insertion rate):\t" << 100*insertionRate << "%" << endl;
        cout << "Metric 2c (deletion rate): \t" << 100*deletionRate  << "%" << endl;
    }

    if (requestedMetrics_[3] || requestedMetrics_[4] || requestedTables_[3])
    {
#define MIN_COVERAGE 5
        unsigned int totalMismatchCount = 0;
        unsigned int totalBaseCount = 0;
        unsigned int totalPassedWindows = 0;
        for (unsigned int i = 0; i < baseCountPer10kWindow_.data_.size(); ++i)
        {
            unsigned int mismatchCount = mismatchCountPer10kWindow_.get( i );
            unsigned int baseCount = baseCountPer10kWindow_.get( i );
            double mismatchRate = 100.0 * (double)mismatchCount / (double)baseCount;
            if (requestedTables_[3])
            {
                cout << (boost::format("Table 3\t%d\t%d\t%f") % mismatchCount % baseCount % mismatchRate).str();
            }
            if (baseCount > MIN_COVERAGE * baseCountPer10kWindow_.windowSize_)
            {
                totalMismatchCount += mismatchCount;
                totalBaseCount += baseCount;
                ++totalPassedWindows;
                if (requestedTables_[3])
                {
                    cout << "\tpassed";
                }
            }
            else
            {
                if (requestedTables_[3])
                {
                    cout << "\tignored";
                }
            }
            if (requestedTables_[3])
            {
                cout << endl;
            }
        }
        double averageMismatchRate = 100.0 * (double)totalMismatchCount / (double)totalBaseCount;
        double averageCoverage = (double)totalBaseCount / ((double)totalPassedWindows * (double)baseCountPer10kWindow_.windowSize_);
        if (requestedTables_[3])
        {
            cout << (boost::format("Table3 Summary:\t%d\t%d\t%f\t%f") % totalMismatchCount % totalBaseCount % averageMismatchRate % averageCoverage).str() << endl;
        }

        if (requestedMetrics_[3] || requestedMetrics_[4])
        {
            // Calculate metrics 3&4
            double coverageVariance = 0;
            double mismatchRateVariance = 0;
            for (unsigned int i = 0; i < baseCountPer10kWindow_.data_.size(); ++i)
            {
                unsigned int mismatchCount = mismatchCountPer10kWindow_.get( i );
                unsigned int baseCount = baseCountPer10kWindow_.get( i );
                double mismatchRate = 100.0 * (double)mismatchCount / (double)baseCount;
                if (baseCount > MIN_COVERAGE * baseCountPer10kWindow_.windowSize_)
                {
                    coverageVariance += pow( ((double)baseCount / (double)baseCountPer10kWindow_.windowSize_) - averageCoverage, 2 );
                    mismatchRateVariance += pow( mismatchRate - averageMismatchRate, 2 );
                }
            }
            mismatchRateVariance /= (double)totalPassedWindows;
            coverageVariance /= (double)totalPassedWindows;
            if (requestedMetrics_[3])
            {
                cout << "Metric 3:\t" << sqrt( mismatchRateVariance ) << endl;
            }
            if (requestedMetrics_[4])
            {
                cout << "Metric 4:\t" << sqrt( coverageVariance ) << endl;
            }
        }
    }

    if (requestedTables_[6])
    {
        for (unsigned int i=0; i<insertSizes_.size(); ++i)
        {
            if (insertSizes_[i])
            {
                cout << (boost::format("insertSize\t%d\t%d") % i % insertSizes_[i]).str() << endl;
            }
        }
    }

    if (requestedTables_[7])
    {
        cout << " *** Quality table ***" << endl;
        cout << "Output to Qualitytable.xx" << endl;
        ofstream ofs( "QualityTable.xx" );
        ofs << "#templateNum\tcycle\tQ1:#Q1\tQ2:#Q2\t..." << endl;

        ofs << "0\t0";
        for (unsigned int templateNum = 1; templateNum < qualityTable_.size(); ++templateNum)
        {
            unsigned long long count = 0;
            //for (unsigned int cycle = 0; cycle < qualityTable_[templateNum].size(); ++cycle)
            if (qualityTable_[templateNum].size() > 1)
            {
                unsigned int cycle = 1; // any cycle would do
                for (unsigned int qscore = 0; qscore < qualityTable_[templateNum][cycle].size(); ++qscore)
                {
                    count += qualityTable_[templateNum][cycle][qscore];
                }
            }

            if (count)
            {
                ofs << '\t' << templateNum << ':' << count;
            }
        }
        ofs << endl;

        for (unsigned int templateNum = 1; templateNum < qualityTable_.size(); ++templateNum)
        {
            for (unsigned int cycle = 0; cycle < qualityTable_[templateNum].size(); ++cycle)
            {
                ofs << templateNum << '\t' << (cycle + 1);
                for (unsigned int qscore = 0; qscore < qualityTable_[templateNum][cycle].size(); ++qscore)
                {
                    unsigned long long count = qualityTable_[templateNum][cycle][qscore];
                    if (count)
                        ofs << '\t' << qscore << ':' << count;
                }
                ofs << endl;
            }
        }
        ofs << endl;
    }
}


} // namespace main
} // namespace eagle
