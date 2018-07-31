/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <boost/algorithm/string/predicate.hpp>
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/format.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lexical_cast.hpp>
#include "common/Exceptions.hh"
#include "io/Text.hh"
#include "model/Phred.hh"
#include "genome/QualityModel.hh"
#include "libzoo/util/AutoGrowVector.hh"

using namespace std;

//#define DEBUG_QUALITIES


namespace eagle
{
namespace genome
{

class MyDiscreteDist
{
public:
    MyDiscreteDist() {}

    MyDiscreteDist( const vector<double> w )
    {
        setWeights( w );
    }

    void setWeights( const vector<double> w )
    {
        weights_.reserve( w.size() );
        double sum = std::accumulate( w.begin(), w.end(), 0.0 );
        double cur = 0;
        for (unsigned int i=0; i<w.size(); ++i)
        {
            cur += w[i];
            weights_.push_back( cur?cur/sum:-1 );
        }
    }

    unsigned int operator()( const double val ) {
        vector<double>::iterator it = lower_bound( weights_.begin(), weights_.end(), val );
        if ( it != weights_.end() )
            return (&*it - &weights_[0]);
        else
            return weights_.size() - 1;
    }

private:
    vector<double> weights_;
};


QualityModel::QualityModel( const vector<boost::filesystem::path>& qualityTableFiles )
{ 
    useNewStuff_ = false;
    if (useNewStuff_)
    {
        parseBigQualityTableFile( qualityTableFiles.front() );
        return;
    }
    unsigned int lastCycle = 0;
    BOOST_FOREACH( const boost::filesystem::path& file, qualityTableFiles)
    {
        lastCycle = parseQualityTableFile( file, lastCycle );
    }
}

unsigned int QualityModel::parseBigQualityTableFile( const boost::filesystem::path& filename )
{
    assert (useNewStuff_);
    unsigned long totalSize = boost::filesystem::file_size( filename );
    bigTable_.resize( totalSize / sizeof(unsigned int) );
    ifstream is( filename.string().c_str(), ios_base::binary );
    is.read( reinterpret_cast<char*>(&bigTable_[0]), totalSize );
    return 0;
}

unsigned int QualityModel::parseQualityTableFile( const boost::filesystem::path& filename, const int cycleOffset )
// Returns last cycle processed + 1 (= first cycle of next file to process)
{
    assert (!useNewStuff_);

    io::DsvReader tsvReader( filename );
    vector<string> tokens;
    unsigned int cycle = cycleOffset;
/*
    double overallExpectedErrors = 0.0;
    double overallOutOf = 0.0;
    double overallExpectedErrorsForCurrentRead = 0.0;
    double overallOutOfForCurrentRead = 0.0;
*/

    if ( boost::algorithm::ends_with( filename.string(), ".qtable2" ) )
    {
        // Parsing of new file format

        AutoGrowVector< AutoGrowVector< AutoGrowVector< double > > > qualityTable;
        while( tsvReader.getNextLineFields<'\t'>(tokens) )
        {
            assert( tokens.size() >= 3 && "Quality table should have at least 3 entries per line (profileId, cycle, and 1 or more quality:count)" );

            unsigned int profileId = boost::lexical_cast<unsigned int>( tokens[0] );
            cycle = boost::lexical_cast<unsigned int>( tokens[1] );
            if (cycle == 0 && profileId == 0)
            {
                cycle = 1; // just a convention that profileId 0 and cycle 1 determines the proportion of profileIds
            }
            assert( cycle > 0 ); // cycles are 1-based by convention
            cycle += cycleOffset;

            for (unsigned int i=2; i<tokens.size(); ++i)
            {
                if ( tokens[i].empty() )
                    continue;

                unsigned int quality;
                double count = 0;
                istringstream ss( tokens[i] );
                char colon;
                ss >> quality >> colon >> count;
                assert( colon == ':' && count > 0 && "quality:count pairs must be separated by ':' in quality table files" );
                qualityTable[profileId][cycle][quality] = count;
            }

            // Extend quality vector as appropriate
            if (qualityDistPerCyclePerLastQuality.size() <= cycle)
            {
                qualityDistPerCyclePerLastQuality.resize( cycle+1 );
            }
            if (qualityDistPerCyclePerLastQuality[cycle].size() <= profileId)
            {
                qualityDistPerCyclePerLastQuality[cycle].resize( profileId+1 );
            }
            qualityDistPerCyclePerLastQuality[cycle][profileId] = boost::random::discrete_distribution<>(qualityTable[profileId][cycle]);
        }
        return cycle;
    }

    // Parsing of old file format
    while( tsvReader.getNextLineFields<'\t'>(tokens) )
    {
        assert( (tokens.size() == 43 || tokens.size() == 53) && "Quality table should have 43 or 53 entries per line" );

        cycle = boost::lexical_cast<unsigned int>( tokens[0] );
        assert( cycle > 0 ); // cycles are 1-based by convention
        cycle += cycleOffset;
        unsigned int lastQ = boost::lexical_cast<unsigned int>( tokens[1] );

        vector<double> values;
        try
        {
            std::transform( tokens.begin()+2, tokens.end(), std::back_inserter(values), boost::bind( &boost::lexical_cast<double,std::string>, _1) );
        }
        catch (const boost::bad_lexical_cast &e)
        {
            EAGLE_ERROR("Error while reading quality table: a numerical field seems to contain non-numerical characters");
        }
        assert( values.size() == 41 || values.size() == 51 );

        // Ensure that the distribution to determine the quality level cannot generate a value of zero
        if (lastQ==0 && values[0]!=0)
        {
            EAGLE_ERROR("Error while reading quality table: Column 3 must be zero if column 2 is zero (the distribution to determine the quality level cannot generate a value of zero)");
        }

        // Extend quality vector as appropriate
        if (qualityDistPerCyclePerLastQuality.size() <= cycle)
        {
            qualityDistPerCyclePerLastQuality.resize( cycle+1 );
        }
        if (qualityDistPerCyclePerLastQuality[cycle].size() <= lastQ)
        {
            qualityDistPerCyclePerLastQuality[cycle].resize( lastQ+1 );
        }
        qualityDistPerCyclePerLastQuality[cycle][lastQ] = boost::random::discrete_distribution<>(values);
/*
        // Expected errors statistics
        double expectedErrors = 0.0;
        double outOf = 0.0;

        if (cycle==102 && lastQ==0)
        {
            cout << (boost::format("Quality stats for read 1: expected mismatch rate = %f%%") % (100.0*overallExpectedErrorsForCurrentRead/overallOutOfForCurrentRead) ).str() << endl;
            overallExpectedErrorsForCurrentRead = 0.0;
            overallOutOfForCurrentRead = 0.0;
        }

        for (unsigned int qValue=0; qValue<=40; ++qValue)
        {
            double qCount = values[qValue];
            expectedErrors += qValue?(model::Phred::qualToProb(qValue)*qCount):0; // don't include Q0 in error count
            outOf += qCount;
#ifdef DEBUG_QUALITIES
            cout << (boost::format("Q=%d count=%d => error-rate=%f => %f expected mismatches out of %f bases = %f%%") % qValue % qCount % qualToProb(qValue) % expectedErrors % outOf % (100.0*expectedErrors/outOf) ).str() << endl;
#endif // DEBUG_QUALITIES
        }

        if (outOf > 0)
        {
            overallExpectedErrors += expectedErrors;
            overallOutOf += outOf;
            overallExpectedErrorsForCurrentRead += expectedErrors;
            overallOutOfForCurrentRead += outOf;
        }
#ifdef DEBUG_QUALITIES
        cout << (boost::format("Quality stats for cycle %d: %f expected mismatches out of %f bases = %f%%") % cycle % expectedErrors % outOf % (100.0*expectedErrors/outOf) ).str() << endl;
        cout << (boost::format("Quality stats over all cycles : %f expected mismatches out of %f bases = %f%%") % overallExpectedErrors % overallOutOf % (100.0*overallExpectedErrors/overallOutOf) ).str() << endl;
#endif // DEBUG_QUALITIES
*/
    }
/*
    cout << (boost::format("Quality stats for read 2: expected mismatch rate = %f%%") % (100.0*overallExpectedErrorsForCurrentRead/overallOutOfForCurrentRead) ).str() << endl;
    cout << (boost::format("Quality stats overall   : expected mismatch rate = %f%%") % (100.0*overallExpectedErrors/overallOutOf) ).str() << endl;
*/
    return cycle;
}

/*
static int Qbins[41] = {
        0, 0, 0, 0, 1, 1, 1, 1, 1, 1,
        2, 2, 2, 2, 2, 2, 2, 2, 2, 3,
        3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
        4, 5, 5, 5, 5, 5, 6, 6, 6, 7, 7 };
static int bin2Q[8] = { 0, 6, 15, 22, 27, 33, 37, 40 };
*/

/*
unsigned int QualityModel::getQuality( boost::mt19937& randomGen, const unsigned int cycle, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
    unsigned int result = bin2Q[ 5 ];
    assert (useNewStuff_);

    unsigned int cycleNum = cycle;
    unsigned int precedingKmer = clusterErrorModelContext.qualityModelContext.kmer & 0x3FF;
    unsigned int newBaseNum = bclBase;
    unsigned int prevQualityBin = Qbins[ clusterErrorModelContext.qualityModelContext.profileNumber ];
    vector<unsigned int> counts(8);
    unsigned int sum = 0;
    for (unsigned int newQualityBin = 0; newQualityBin<8; ++newQualityBin)
    {
//        assert( cycleNum >= 1 && cycleNum <= cycleCount );
        assert( precedingKmer < 1024 );
        assert( newBaseNum < 4 );
        assert( prevQualityBin < 8 );
        assert( newQualityBin < 8 );
        unsigned int entryNum =
            (cycleNum-1) * (1024*4*8*8) + 
            precedingKmer * (4*8*8) +
            newBaseNum * (8*8) +
            prevQualityBin * (8) +
            newQualityBin;
        assert( entryNum < bigTable_.size() );
        unsigned int val = bigTable_[entryNum];
        counts[newQualityBin] = val;
        sum += val;
    }

    if (sum != 0)
    {
        unsigned int rnd = (unsigned int)( (double)randomGen() / randomGen.max() * sum );
        for (unsigned int newQualityBin = 0; newQualityBin<8; ++newQualityBin)
        {
            if (rnd >= counts[newQualityBin])
            {
                rnd -= counts[newQualityBin];
            }
            else
            {
                result = bin2Q[ newQualityBin ];
                break;
            }
        }
    }
    else
    {
        result = bin2Q[ 5 ]; //clusterErrorModelContext.qualityModelContext.profileNumber;
    }

    clusterErrorModelContext.qualityModelContext.kmer = (clusterErrorModelContext.qualityModelContext.kmer << 2 ) | (bclBase & 3);
    clusterErrorModelContext.qualityModelContext.profileNumber = result;

    return result;
}
*/

unsigned int QualityModel::getQuality( boost::mt19937& randomGen, const unsigned int cycle, ClusterErrorModelContext& clusterErrorModelContext )
{
    assert (!useNewStuff_);

    if (cycle >= qualityDistPerCyclePerLastQuality.size())
    {
        EAGLE_ERROR( "The quality table doesn't model as many cycles as necessary for this simulation" );
    }

    if (clusterErrorModelContext.qualityModelContext.profileNumber == 0)
    {
        // No profile number assigned to this read => find the last cycle containing a profile spec and use it
        unsigned int cycleForProfileNumber = cycle;
        while (qualityDistPerCyclePerLastQuality[cycleForProfileNumber][0].max() == 0)
        {
            if (cycleForProfileNumber == 0)
            {
                EAGLE_ERROR("Cannot find quality level distribution in quality tables (there should at least be an entry for {cycle=0, profile=0})");
            }
            --cycleForProfileNumber;
        }
        clusterErrorModelContext.qualityModelContext.profileNumber = qualityDistPerCyclePerLastQuality[cycleForProfileNumber][0](randomGen);
        assert( clusterErrorModelContext.qualityModelContext.profileNumber > 0 );
    }
    unsigned int profileNumber = clusterErrorModelContext.qualityModelContext.profileNumber;
    assert( profileNumber > 0 );
    if (profileNumber >= qualityDistPerCyclePerLastQuality[cycle].size())
    {
        EAGLE_ERROR( (boost::format("The quality table doesn't contain the required entry for {cycle=%d, profileNumber=%d}") % cycle % profileNumber).str() );
    }
    unsigned int quality = qualityDistPerCyclePerLastQuality[cycle][profileNumber](randomGen);

#ifdef DEBUG_QUALITIES
    cout << (boost::format("Q=%d => error-rate=%f") % quality % model::Phred::qualToProb(quality)).str() << endl;
#endif // DEBUG_QUALITIES
    return quality;
}


SequencingMismatchModel::SequencingMismatchModel( const boost::filesystem::path& mismatchTableFilename )
{
    // Mismatch model for base 'x': x->A, x->C, x->G, x->T, del, x->insertedA, x->insertedC, x->insertedG, x->insertedT (dup=='x->insertedx')
    // Parse mismatch table file: each line based on the model, for x={A,C,G,T}
    if (mismatchTableFilename == "")
    { // Use default values
        vector< double > errorA = boost::assign::list_of(0.0)(1.0)(1.0)(1.0)(0.0)(0.0)(0.0)(0.0)(0.0);
        vector< double > errorC = boost::assign::list_of(1.0)(0.0)(1.0)(1.0)(0.0)(0.0)(0.0)(0.0)(0.0);
        vector< double > errorG = boost::assign::list_of(1.0)(1.0)(0.0)(1.0)(0.0)(0.0)(0.0)(0.0)(0.0);
        vector< double > errorT = boost::assign::list_of(1.0)(1.0)(1.0)(0.0)(0.0)(0.0)(0.0)(0.0)(0.0);

        errorDistPerBase.push_back( boost::random::discrete_distribution<>(errorA) );
        errorDistPerBase.push_back( boost::random::discrete_distribution<>(errorC) );
        errorDistPerBase.push_back( boost::random::discrete_distribution<>(errorG) );
        errorDistPerBase.push_back( boost::random::discrete_distribution<>(errorT) );
    }
    else
    { // Parse mismatch table file
        io::DsvReader tsvReader( mismatchTableFilename );
        vector<string> tokens;
        const vector<string> expectedRowHeaders = boost::assign::list_of("A")("C")("G")("T");
        BOOST_FOREACH( const string& expectedRowHeader, expectedRowHeaders)
        {
            (void)expectedRowHeader; // prevents "unused variable" warning when asserts are not compiled
            tsvReader.getNextLineFields<'\t'>(tokens);
            assert( tokens.size() == 10 && "There should be 10 entries per line" );
            assert( tokens[0] == expectedRowHeader && "Unexpected value in first column (The first 4 entries - ignoring the lines starting with '#' - should be A,C,G,T)" );
            vector<double> values;
            try
            {
                std::transform( tokens.begin()+1, tokens.end(), std::back_inserter(values), boost::bind( &boost::lexical_cast<double,std::string>, _1) );
            }
            catch (const boost::bad_lexical_cast &e)
            {
                EAGLE_ERROR("Error while reading mismatch table: a numerical field seems to contain non-numerical characters");
            }
            assert( values.size() == 9 );
            errorDistPerBase.push_back( boost::random::discrete_distribution<>(values) );
        }
    }
}

void SequencingMismatchModel::apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
    if (randomGen() > errorRate * randomGen.max())
    {
        randomErrorType = ErrorModel::NoError;
    }
    else
    {
        unsigned char errorType = errorDistPerBase[bclBase](randomGen);
        switch (errorType)
        {
        case 0: // x->A
        case 1: // x->C
        case 2: // x->G
        case 3: // x->T
            randomErrorType = ErrorModel::BaseSubstitution;
            bclBase = errorType;
            break;

        case 4: // del
            randomErrorType = ErrorModel::BaseDeletion;
            break;

        case 5: // x->insertedA
        case 6: // x->insertedC
        case 7: // x->insertedG
        case 8: // x->insertedT
            randomErrorType = ErrorModel::BaseInsertion;
            bclBase = errorType-5;
            break;

        default:
            assert( false );
        }
    }

#ifdef DEBUG_QUALITIES
    // Debugging stats
    static unsigned long mismatches = 0;
    static unsigned long noMismatches = 0;
    if (randomErrorType == ErrorModel::NoError)
    {
        noMismatches++;
    }
    else
    {
        mismatches++;
    }
    cout << (boost::format("mismatch rate = %d mismatches and %d noMismatch out of %d => %f%%") % mismatches % noMismatches % (mismatches + noMismatches) % (100.0*(float)mismatches/(float)(mismatches + noMismatches))).str() << endl;
#endif // DEBUG_QUALITIES
}


HomopolymerIndelModel::HomopolymerIndelModel( const boost::filesystem::path& homopolymerIndelTableFilename )
{
    if (homopolymerIndelTableFilename == "")
    { // Use default values: no indels
        homoDeletionTable_ = homoInsertionTable_ = boost::assign::list_of
            (0.0)     // homopolymer of size 0 or more
            ;
    }
    else
    { // Parse homopolymer indel table file
        io::DsvReader tsvReader( homopolymerIndelTableFilename );
        vector<string> tokens;
        while ( tsvReader.getNextLineFields<'\t'>(tokens) )
        {
            if ( tokens.size() == 0 ) { continue; }
            assert ( tokens.size() == 3 && "There should be 3 entries per line" );
            vector<double> values;
            try
            {
                std::transform( tokens.begin(), tokens.end(), std::back_inserter(values), boost::bind( &boost::lexical_cast<double,std::string>, _1) );
            }
            catch (const boost::bad_lexical_cast &e)
            {
                EAGLE_ERROR("Error while reading homopolymer indel table: a numerical field seems to contain non-numerical characters");
            }
            assert( values.size() == 3 );
            assert( values[0] == homoDeletionTable_.size() && values[0] == homoInsertionTable_.size() && "First column should contain consecutive numbers starting from 0" );
            homoDeletionTable_.push_back ( values[1] );
            homoInsertionTable_.push_back( values[2] );
        }
    }
}

void HomopolymerIndelModel::apply( boost::mt19937& randomGen, const double errorRate, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
    if (bclBase != clusterErrorModelContext.homopolymerModelContext.lastBase)
    {
        clusterErrorModelContext.homopolymerModelContext.lastBase = bclBase;
        clusterErrorModelContext.homopolymerModelContext.homopolymerLength = 1;
    }
    else
    {
        ++(clusterErrorModelContext.homopolymerModelContext.homopolymerLength);
        unsigned int length = clusterErrorModelContext.homopolymerModelContext.homopolymerLength;
        assert( length >= 1 );
        unsigned int tableEntry = std::min<unsigned int>( length, homoDeletionTable_.size()-1 );
        double delErrorRate = homoDeletionTable_ [tableEntry];
        double insErrorRate = homoInsertionTable_[tableEntry];
        double randomValue = (double)randomGen() / (double)randomGen.max();
        switch (clusterErrorModelContext.homopolymerModelContext.errorDirection)
        {
        case 0: // This homopolymer didn't get any insertion or deletion yet
            if (randomValue < delErrorRate )
            {
                randomErrorType = ErrorModel::BaseDeletion;
                clusterErrorModelContext.homopolymerModelContext.errorDirection = -1;
            }
            else if (randomValue < delErrorRate + insErrorRate )
            {
                randomErrorType = ErrorModel::BaseInsertion;
                clusterErrorModelContext.homopolymerModelContext.errorDirection = +1;
            }
            break;
        case +1: // This homopolymer already got one or more insertions
            if (randomValue < insErrorRate )
            {
                randomErrorType = ErrorModel::BaseInsertion;
            }
            break;
        case -1: // This homopolymer already got one or more deletions
            if (randomValue < delErrorRate )
            {
                randomErrorType = ErrorModel::BaseDeletion;
            }
            break;
        }
    }
}


class MotifRepeatQualityDropInfo
{
public:
   MotifRepeatQualityDropInfo()
       : meanQualityDrop( 0 )
        {}

    float meanQualityDrop;
    MyDiscreteDist distribution;
};

eagle::model::IUPAC baseConverter_;

uint64_t kmerStringToInt64( const string &s )
{
    uint64_t result = 0;
    for (unsigned int i=0; i<s.size(); ++i)
    {
        unsigned int binValue;
        if (s[i]>='0' && s[i] <='3')
            binValue = s[i] - '0';
        else
            binValue = baseConverter_.bin(s[i]);
        result = (result << 2) | binValue;
    }
    return result;
}

#define MAX_MOTIF_KMER_LENGTH 10
#define AVERAGE_QUALITY 34
MotifQualityDropModel::MotifQualityDropModel( const boost::filesystem::path& tableFilename )
    : active_( false )
{
    if (tableFilename != "")
    {
        tableData_.resize( MAX_MOTIF_KMER_LENGTH+1 );

        // Parse motif quality drop table file
        io::DsvReader tsvReader( tableFilename );
        vector<string> tokens;
        while ( tsvReader.getNextLineFields<'\t'>(tokens) )
        {
            if ( tokens.size() == 0 ) { continue; }
            assert ( tokens.size() > 3 && "There should be more than 4 entries per line: {kmer, repeatCount, coveredBaseCount, meanQ, (quality:count)+}" );
            uint64_t kmer = kmerStringToInt64( tokens[0] );
            unsigned int repeatCount = boost::lexical_cast<unsigned int, string>( tokens[1] );
            float meanQualityDrop = boost::lexical_cast<double, string>( tokens[2] );
            unsigned int kmerLength = tokens[0].size();
            unsigned int kmerLengthInBits = 2 * kmerLength;
            uint64_t kmerMask = (1ull << kmerLengthInBits) - 1;

            AutoGrowVector<double> distribution;
            for (unsigned int i=3; i<tokens.size(); ++i)
            {
                unsigned int quality;
                double count = 0;
                istringstream ss( tokens[i] );
                char colon;
                ss >> quality >> colon >> count;
                assert( colon == ':' && count > 0 && "quality:count pairs must be separated by ':' in motif Qdrop table" );
                distribution[quality] = count;
            }

            active_ = true;
            boost::shared_ptr< std::vector< MotifRepeatQualityDropInfo > > &mapData = tableData_[kmerLength][ kmer ];
            if (mapData.get() == NULL)
            {
                mapData.reset( new std::vector< MotifRepeatQualityDropInfo > );

                uint64_t kmerPermutation = kmer;
                for (unsigned int permutation=1; permutation < kmerLength; ++permutation)
                {
                    // Permute kmer
                    uint64_t leftMostBase = kmerPermutation >> (kmerLengthInBits - 2);
                    kmerPermutation = ((kmerPermutation << 2) & kmerMask) | leftMostBase;

                    //            cout << "MotifQualityDropModel constructor: Adding " << kmer << "\t" << repeatCount << "\t" << roundf(37 - meanQualityDrop) << endl;
                    boost::shared_ptr< std::vector< MotifRepeatQualityDropInfo > > &mapDataPermutation = tableData_[kmerLength][ kmerPermutation ];
                    assert (mapDataPermutation.get() == NULL);
                    mapDataPermutation = mapData;
                }
            }

            if (repeatCount >= mapData->size())
            {
                mapData->resize( repeatCount+1 );
            }
            assert( (*mapData)[repeatCount].meanQualityDrop == 0 && "Same motif repeat is described in 2 lines of motif Qdrop table" );
            (*mapData)[repeatCount].meanQualityDrop = AVERAGE_QUALITY - meanQualityDrop;
            (*mapData)[repeatCount].distribution.setWeights( distribution );
        }
    }
}

MotifRepeatQualityDropInfo* MotifQualityDropModel::getMotifRepeatQualityDrop( const uint64_t kmer1, const unsigned int repeatKmerLength, const unsigned int repeatCount )
{
//    cout << "getMotifRepeatQualityDrop:" << endl;
//    cout << "  kmer1=" << kmer1 << endl;
//    cout << "  repeatKmerLength=" << repeatKmerLength << endl;
//    cout << "  repeatCount=" << repeatCount << endl;

//    int result = ( repeatCount - 1 ) * repeatKmerLength;
//    cout << " => result=" << result << endl;

    std::map< uint64_t, boost::shared_ptr< std::vector< MotifRepeatQualityDropInfo > > >::iterator it = tableData_[repeatKmerLength].find( kmer1 );
    if (it != tableData_[repeatKmerLength].end())
    {
//        cout << "kmer found up to repeat " << it->second->size() << endl;
        unsigned int repeatCount2 = min<unsigned int>( repeatCount, it->second->size() );
        if (repeatCount2 > 0)
        {
            if ( (*it->second).size() <= repeatCount2 )
                repeatCount2 = (*it->second).size() - 1;
            MotifRepeatQualityDropInfo *info = &((*it->second)[repeatCount2]);
//            int result2 = (int)(info.meanQualityDrop);
//            cout << " => result2=" << result2 << endl;
            return info;
        }
    }

    return 0;
}

void MotifQualityDropModel::applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext, const unsigned int cycle, boost::mt19937& randomGen )
{
    if (!active_) { return; }

    uint64_t &kmer = clusterErrorModelContext.motifQualityDropModelContext.kmer;
    unsigned int &kmerLength = clusterErrorModelContext.motifQualityDropModelContext.kmerLength;
    double &qualityDropLevel = clusterErrorModelContext.motifQualityDropModelContext.qualityDropLevel;
    MotifRepeatQualityDropInfo* &shortTermEffect = clusterErrorModelContext.motifQualityDropModelContext.shortTermEffect;
    float &shortTermQualityDrop = clusterErrorModelContext.motifQualityDropModelContext.shortTermQualityDrop;
    int &qualityDropDueToPhasing = clusterErrorModelContext.phasingContext.qualityDrop;

    bool repeatDetected = false;
    MotifRepeatQualityDropInfo *strongestRepeat_effect = 0;
    unsigned int strongestRepeat_repeatLengthExcludingFirst = 0;

    if (kmerLength >= 5) // No point looking for short repeats
    {
        // detect repeated kmer
        const unsigned int maxKmerLength = min<unsigned int>( MAX_MOTIF_KMER_LENGTH, kmerLength );

        for ( unsigned int repeatKmerLength = 1; repeatKmerLength <= maxKmerLength; ++repeatKmerLength )
        {
            unsigned int repeatKmerLengthInBits = 2 * repeatKmerLength;
            uint64_t kmerMask = (1ull << repeatKmerLengthInBits) - 1;
            uint64_t kmer0 = kmer;
            uint64_t kmer1 = kmer & kmerMask;
            uint64_t kmer2;
            unsigned int repeatCount = 0;
            do
            {
                kmer0 >>= repeatKmerLengthInBits;
                ++repeatCount;
                kmer2 = kmer0 & kmerMask;
            } while (kmer1 == kmer2 &&
                     repeatKmerLength * (repeatCount + 1) <= kmerLength &&
                     repeatKmerLengthInBits * (repeatCount + 1) <= 64
                );

            unsigned int repeatLengthExcludingFirst = ( repeatCount - 1 ) * repeatKmerLength;
            const unsigned int repeatLengthThreshold = 4;
            if ( repeatLengthExcludingFirst >= repeatLengthThreshold
                 && repeatLengthExcludingFirst > strongestRepeat_repeatLengthExcludingFirst )
            {
                MotifRepeatQualityDropInfo *effect = getMotifRepeatQualityDrop( kmer1, repeatKmerLength, repeatCount );
                if ( ( effect && strongestRepeat_effect && (int)effect->meanQualityDrop > (int)strongestRepeat_effect->meanQualityDrop )
                     || ( effect && !strongestRepeat_effect && (int)effect->meanQualityDrop > 0 ) )
                {
                    repeatDetected = true;
                    strongestRepeat_effect = effect;
                    strongestRepeat_repeatLengthExcludingFirst = repeatLengthExcludingFirst;
                }
            }
        }

        if ( repeatDetected
             && strongestRepeat_effect->meanQualityDrop > (shortTermEffect?shortTermEffect->meanQualityDrop:0) )
        {
            if (qualityDropLevel == 0)
            {
                qualityDropLevel = (double)randomGen() / randomGen.max();
            }
            int newQuality = strongestRepeat_effect->distribution( qualityDropLevel );
            assert( newQuality != 0 );
            int newQualityDrop = max<int>(AVERAGE_QUALITY - newQuality, 0);

            if (shortTermQualityDrop < newQualityDrop)
            {
                qualityDropDueToPhasing -= (int)shortTermQualityDrop;
                qualityDropDueToPhasing += newQualityDrop;
                shortTermQualityDrop = newQualityDrop;
            }
            shortTermEffect = strongestRepeat_effect;
        }
        else
        {
            // Pulls down the "quality drop effect to cancel"
            float reducedMeanQualityDrop = strongestRepeat_effect?strongestRepeat_effect->meanQualityDrop:0;
            double newShortTermQualityDrop;
            if (shortTermEffect && shortTermEffect->meanQualityDrop)
            {
                assert( reducedMeanQualityDrop <= shortTermEffect->meanQualityDrop );
                newShortTermQualityDrop = shortTermQualityDrop * reducedMeanQualityDrop / shortTermEffect->meanQualityDrop;
            }
            else
            {
                assert( reducedMeanQualityDrop <= 0 );
                newShortTermQualityDrop = reducedMeanQualityDrop;
            }
            const double attenuation = 1.0;
            shortTermQualityDrop = (shortTermQualityDrop * (1.0 - attenuation)) + (newShortTermQualityDrop * attenuation);
            shortTermEffect = strongestRepeat_effect;
        }
    }

    kmer = (kmer << 2) | (bclBase & 3);
    ++kmerLength;
}


RandomQualityDropModel::RandomQualityDropModel( /*const boost::filesystem::path& tableFilename*/ )
{
}

void RandomQualityDropModel::applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
}


QualityGlitchModel::QualityGlitchModel( /*const boost::filesystem::path& tableFilename*/ )
{
}

void QualityGlitchModel::applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
}


HappyPhasingModel::HappyPhasingModel( /*const boost::filesystem::path& tableFilename*/ )
{
}

void HappyPhasingModel::applyQualityDrop( unsigned int& quality, const char bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
    // Use the values from the MotifQualityDrop plugin, if available
    const int qualityDropDueToPhasing = clusterErrorModelContext.phasingContext.qualityDrop;
    if (qualityDropDueToPhasing < 5)
        return; // Such a small quality drop is not worth correcting

    const unsigned int kmerLength = clusterErrorModelContext.motifQualityDropModelContext.kmerLength;
    if (kmerLength > 3)
    {
        uint64_t kmer = clusterErrorModelContext.motifQualityDropModelContext.kmer;
        float concordance = 0.0;
        float maxConcordance = 0.0;
        float concordanceWeight = 1.0;
        float decay = 0.7;
        for (int i=0; i<3; ++i)
        {
            kmer >>= 2;
            if ((unsigned char)bclBase == (kmer & 3))
                concordance += concordanceWeight;
            maxConcordance += concordanceWeight;
            concordanceWeight *= decay;
        }
        concordance /= maxConcordance;
        concordance *= concordance; // squared seems nice

        // restore quality in happy phasing areas
        quality += (int)(qualityDropDueToPhasing * concordance);
    }
}


QQTable::QQTable( const boost::filesystem::path& qqTableFilename )
{
    if (qqTableFilename == "")
    { // Use default Phred values
        for (unsigned int q=0; q<model::Phred::QUALITY_MAX; ++q)
        {
            qualityToProbability_.push_back( model::Phred::qualToProb( q ) );
        }
        qualityToProbability_.push_back( 0.0 ); // QUALITY_MAX has this special property
    }
    else
    { // Parse QQ table file
        io::DsvReader tsvReader( qqTableFilename );
        vector<string> tokens;
        while ( tsvReader.getNextLineFields<'\t'>(tokens) )
        {
            if ( tokens.size() == 0 ) { continue; }
            assert ( tokens.size() == 2 && "There should be 2 entries per line" );
            vector<double> values;
            try
            {
                std::transform( tokens.begin(), tokens.end(), std::back_inserter(values), boost::bind( &boost::lexical_cast<double,std::string>, _1) );
            }
            catch (const boost::bad_lexical_cast &e)
            {
                EAGLE_ERROR("Error while reading QQ table: a numerical field seems to contain non-numerical characters");
            }
            assert( values.size() == 2 );
            assert( values[0] == qualityToProbability_.size() && "First column should contain consecutive numbers starting from 0" );
            qualityToProbability_.push_back ( values[1] );
        }
    }

    if (qualityToProbability_.size() < model::Phred::QUALITY_MAX+1)
    {
        EAGLE_ERROR("QQ table doesn't contain enough values");
    }
}

double QQTable::qualToErrorRate( const unsigned int qual )
{
    if (qual >= qualityToProbability_.size())
    {
        BOOST_THROW_EXCEPTION( eagle::common::EagleException( 0, "Requested quality is higher than allowed max") );
    }
    double prob = qualityToProbability_[qual];
    return prob;
}


ErrorModel::ErrorModel( const vector<boost::filesystem::path>& qualityTableFilenames, const boost::filesystem::path& mismatchTableFilename, const boost::filesystem::path& homopolymerIndelTableFilename, const boost::filesystem::path& motifQualityDropTableFilename, const boost::filesystem::path& qqTableFilename, const std::vector< std::string >& errorModelOptions )
    : qualityModel_           ( qualityTableFilenames )
    , sequencingMismatchModel_( mismatchTableFilename )
    , homopolymerIndelModel_  ( homopolymerIndelTableFilename )
    , motifQualityDropModel_  ( motifQualityDropTableFilename )
    , randomQualityDropModel_ ()
    , qualityGlitchModel_     ()
    , happyPhasingModel_      ()
    , longreadBaseDuplicationModel_( errorModelOptions )
    , longreadDeletionModel_       ( errorModelOptions )
    , qqTable_                ( qqTableFilename )
{
}

void ErrorModel::getQualityAndRandomError( boost::mt19937& randomGen, const unsigned int cycle, const char base, unsigned int& quality, unsigned int& randomErrorType, char& bclBase, ClusterErrorModelContext& clusterErrorModelContext )
{
    bclBase = baseConverter_.normalizedBcl( base );
    if (bclBase==4)
    {
        // base is N
        randomErrorType = NoError;
        bclBase = 0;
        quality = 0;
        return;
    }
    assert( bclBase >= 0 && bclBase < 4 );

//    quality = qualityModel_.getQuality( randomGen, cycle, bclBase, clusterErrorModelContext );
    quality = qualityModel_.getQuality( randomGen, cycle, clusterErrorModelContext );
    motifQualityDropModel_.applyQualityDrop( quality, bclBase, clusterErrorModelContext, cycle, randomGen );
    randomQualityDropModel_.applyQualityDrop( quality, bclBase, clusterErrorModelContext );
    qualityGlitchModel_.applyQualityDrop( quality, bclBase, clusterErrorModelContext );
    happyPhasingModel_.applyQualityDrop( quality, bclBase, clusterErrorModelContext );

    // Apply quality drop due to phasing, using an additive strategy
    // This quality drop was calculated as part of the previous "applyQualityDrop" methods
    if ((int)quality > clusterErrorModelContext.phasingContext.qualityDrop)
    {
      quality -= clusterErrorModelContext.phasingContext.qualityDrop;
    }
    else
    {
      quality = 0;
    }

    // Make sure quality scores stay above 2
    if ( quality < 2 )
        quality = 2;

    double errorRate = qqTable_.qualToErrorRate( quality );

//#define REPORT_ERROR_RATE
#ifdef REPORT_ERROR_RATE
    static double totalErrors = 0;
    static double totalBases  = 0;
    static double threshold   = 1;
    totalErrors += errorRate;
    totalBases++;
    if (totalBases >= threshold)
    {
        clog << (boost::format("%f @ %f => %f") % totalErrors % totalBases % (totalErrors/totalBases)).str() << endl;
        threshold *= 2;
    }
#endif //ifdef REPORT_ERROR_RATE

    sequencingMismatchModel_.apply( randomGen, errorRate, randomErrorType, bclBase, clusterErrorModelContext );
    homopolymerIndelModel_.apply( randomGen, errorRate, randomErrorType, bclBase, clusterErrorModelContext );
    longreadBaseDuplicationModel_.apply( randomGen, errorRate, randomErrorType, bclBase, clusterErrorModelContext );
    longreadDeletionModel_.apply( randomGen, errorRate, randomErrorType, bclBase, clusterErrorModelContext );
}

} // namespace genome
} // namespace eagle
