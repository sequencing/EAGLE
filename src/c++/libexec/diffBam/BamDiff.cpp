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
#include "genome/BamAdapters.hh"
#include "io/Bam.hh"
#include "io/BamIndexer.hh"
#include "io/BgzfCompressor.hh"
#include "io/StorableBamAlignment.hh"
#include "model/Fragment.hh"
#include "BamDiff.hh"


using namespace std;
#define DEBUG 2
bool skipNotPF = true;
bool doRomanStuff = false;


namespace eagle
{
namespace main
{


vector< io::bam::StorableBamAlignment > DiffComputer::masterAlignments_;
vector< io::bam::StorableBamAlignment > DiffComputer::slaveAlignments_;
vector< unsigned int > DiffComputer::slave2MasterRefId_;
model::FragmentPosResolver DiffComputer::fragmentPosResolver_;
vector<unsigned int> matchCountForMapQ(500);
vector<unsigned int> mismatchCountForMapQ(500);
ofstream outCorrect("diffCorrect.bam.txt");
ofstream outFP     ("diffFP.bam.txt");
ofstream outFN     ("diffFN.bam.txt");

// Output streams
boost::iostreams::filtering_ostream bamStreamFP_;
boost::iostreams::filtering_ostream bamStreamFN_;
boost::iostreams::filtering_ostream bamStreamTP_;
boost::iostreams::filtering_ostream bgzfStreamFP_;
boost::iostreams::filtering_ostream bgzfStreamFN_;
boost::iostreams::filtering_ostream bgzfStreamTP_;



BamDiff::BamDiff (const BamDiffOptions &options )
    : options_( options )
{
//    genome::SharedFastaReference::init( options_.referenceGenome );
}

BamDiff::~BamDiff ()
{
    bgzfStreamFP_.pop();
    bgzfStreamFN_.pop();
    bgzfStreamTP_.pop();
    io::serializeBgzfFooter( bamStreamFP_ );
    io::serializeBgzfFooter( bamStreamFN_ );
    io::serializeBgzfFooter( bamStreamTP_ );
}

void BamDiff::run()
{
    // Create BAM input stream
    {
        // Create a boost iostream flow that takes a BAM file as input, sends it to our BAM parser, which will use a callback for each read
        masterBgzfStream_.push( DiffComputer( 0 ) );
        slaveBgzfStream_.push( DiffComputer( 1 ) );

        // Terminate the iostreams chain
        masterBgzfStream_.push( boost::iostreams::null_sink() );
        slaveBgzfStream_.push( boost::iostreams::null_sink() );
    }

    // Create BAM output streams
    CreateBamOutputStream( "outFP.bam", bamStreamFP_, bgzfStreamFP_ );
    CreateBamOutputStream( "outFN.bam", bamStreamFN_, bgzfStreamFN_ );
    CreateBamOutputStream( "outTP.bam", bamStreamTP_, bgzfStreamTP_ );

    // Send the BAM data to the stream
    ifstream inp1( options_.masterBamFile.string().c_str() );
    ifstream inp2( options_.slaveBamFile.string().c_str() );
    masterBgzfStream_.exceptions( ios_base::badbit ); // Required to rethrow any error thrown by our iostreams filters (leading to better error messages)
    slaveBgzfStream_.exceptions( ios_base::badbit ); // Required to rethrow any error thrown by our iostreams filters (leading to better error messages)
    while (inp1.good() || inp2.good())
    {
        char buf[128];
        inp1.read( buf, sizeof(buf) );
        masterBgzfStream_.write( buf, inp1.gcount() );

        inp2.read( buf, sizeof(buf) );
        slaveBgzfStream_.write( buf, inp2.gcount() );
    }

    // Ensure all the data is flushed through the streams
    masterBgzfStream_.pop();
    slaveBgzfStream_.pop();
}

void BamDiff::CreateBamOutputStream( const boost::filesystem::path& outFilename, boost::iostreams::filtering_ostream& bamStream_, boost::iostreams::filtering_ostream& bgzfStream_ )
{
    const boost::filesystem::path bamPath( outFilename );
    const int bamGzipLevel = boost::iostreams::gzip::best_speed;
    const vector<string> argv_;

    boost::iostreams::file_sink bamSink(bamPath.string());
    if (!bamSink.is_open()) {
        BOOST_THROW_EXCEPTION(common::IoException(errno, "Failed to open output BAM file " + bamPath.string()));
    }
    clog << "Creating BAM file: " << bamPath << endl;

    bamStream_.push(bamSink);

    bgzfStream_.push(eagle::io::bam::BgzfCompressor(bamGzipLevel));

    { // Add BAM Index
        boost::filesystem::path baiPath( outFilename.string() + ".bai" );
        boost::iostreams::file_sink baiSink(baiPath.string());
        bgzfStream_.push(eagle::io::bam::BamIndexer<boost::iostreams::file_sink>(baiSink));
    }

    bgzfStream_.push(bamStream_);
}




DiffComputer::DiffComputer( unsigned int fileNum )
    : fileNum_( fileNum )
{
//#define TEST_POS_RESOLVER
#ifdef TEST_POS_RESOLVER
    //e.g. "unknown-flowcell_0:5:2208:775123:0"
    unsigned int lane = 5;
    unsigned int tile = 2208;
    unsigned long cluster = 775123;
    unsigned long globalSamplePos = fragmentPosResolver_.GetSimulatedPosInSampleGenome( lane, tile, cluster );
    clog << (boost::format("lane=%d, tile=%d, cluster=%d => globalSamplePos=%d") % lane % tile % cluster % globalSamplePos).str() << endl;
    exit(0);
#endif // TEST_POS_RESOLVER
}

DiffComputer::~DiffComputer()
{
}

void DiffComputer::parsedRefSeqInfo( const vector< BamRefInfoItemType >& bamRefInfo )
{
    if (fileNum_ == 0)
    { // We use the master reference for the output
        const vector<string> argv;
        io::serializeHeader<genome::EagleBamHeaderAdapter>( bgzfStreamFP_, argv, genome::EagleBamHeaderAdapter( bamRefInfo ) );
        io::serializeHeader<genome::EagleBamHeaderAdapter>( bgzfStreamFN_, argv, genome::EagleBamHeaderAdapter( bamRefInfo ) );
        io::serializeHeader<genome::EagleBamHeaderAdapter>( bgzfStreamTP_, argv, genome::EagleBamHeaderAdapter( bamRefInfo ) );
    }

    // Fill slave2MasterRefId_ in
    static vector< vector< BamRefInfoItemType > > allrefChroms;
    if (allrefChroms.size() <= fileNum_)
    {
        allrefChroms.resize( fileNum_ + 1 );
    }
    allrefChroms[fileNum_] = bamRefInfo;

    if (allrefChroms.size() == 2
        && !allrefChroms[0].empty()
        && !allrefChroms[1].empty()
        )
    {
        slave2MasterRefId_.resize( allrefChroms[1].size() );
        for (unsigned int i=0; i<allrefChroms[1].size(); ++i)
        {
            bool found = false;
            for (unsigned int j=0; j<allrefChroms[0].size(); ++j)
            {
                if (allrefChroms[1][i].first == allrefChroms[0][j].first)
                {
                    clog << (boost::format("Slave chromosome %d \"%s\" mapped to master chromosome %d") % i % allrefChroms[1][i].first % j).str() << endl;
                    slave2MasterRefId_[i] = j;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                EAGLE_WARNING( "Slave chromosome " << allrefChroms[1][i].first << " not found in master" );
            }
        }
    }
}

void DiffComputer::parsedAlignment( const io::bam::BamAlignment& alignment, const io::bam::VirtualOffset& virtualOffset, const io::bam::VirtualOffset& virtualEndOffset )
{
    if (DEBUG>=3) cout << "received: " << fileNum_ << " " << io::bam::StorableBamAlignment(alignment) << endl;

    if (alignment.getFlag() & 0x4) { return; } // Test that segment is mapped
    if (skipNotPF && alignment.getFlag() & 0x200) { return; } // Test that segment passed filter, if required

    {
        switch (fileNum_)
        {
        case 0:
            masterAlignments_.push_back( alignment );
            break;
        case 1:
            slaveAlignments_.push_back( alignment );
            break;
        }
    }
}

void DiffComputer::finishedParsing()
{
    static int allDone = 0;
    if (++allDone < 2) { return; }

    if (DEBUG>=1) cout << "finishedParsing" << endl;
    const unsigned int window1 = 10;
    const unsigned int window2 = 128;
    bool processMaster = true;
    bool processSlave  = true;
    vector< io::bam::StorableBamAlignment > unmatchedMasterAlignments;
    vector< io::bam::StorableBamAlignment > unmatchedSlaveAlignments;
    vector< io::bam::StorableBamAlignment >::iterator masterItr = masterAlignments_.begin();
    vector< io::bam::StorableBamAlignment >::iterator slaveItr = slaveAlignments_.begin();

    do
    {
        processMaster = true;
        processSlave  = true;
        if (masterItr == masterAlignments_.end()) { processMaster = false; }
        if (slaveItr  == slaveAlignments_.end())  { processSlave  = false; }
        if (processMaster && processSlave)
        {
            if (masterItr->pos <= slaveItr->pos + window1)
            {
                processSlave = false;
            }
            else
            {
                processMaster = false;
            }
        }

        // This is the part that behaves "as if" we were receiving alignments one by one from PUMA, ordered by position
        // Everything above this was just to simulate PUMA's work until we get the real one
        if (processMaster)
        {
            if (DEBUG>=2) cout << "Master: " << *masterItr << endl;
            unmatchedMasterAlignments.push_back( *masterItr );
            ++masterItr;
        }
        else if (processSlave)
        {
            if (DEBUG>=2) cout << "Slave: " << *slaveItr << endl;

            // Trying to match the alignment
            bool matchFound = false;
            for (vector< io::bam::StorableBamAlignment >::reverse_iterator unmatchedMasterItr = unmatchedMasterAlignments.rbegin();
                 unmatchedMasterItr != unmatchedMasterAlignments.rend();
                 ++unmatchedMasterItr)
            {
                if (unmatchedMasterItr->pos < slaveItr->pos - window2) { break; }
                if (io::bam::StorableBamAlignment::SeqCompare( *unmatchedMasterItr, *slaveItr ) == 0)
                {
                    ProcessMatch( *unmatchedMasterItr, *slaveItr );
                    ProcessCorrect( *unmatchedMasterItr, *slaveItr );
                    unmatchedMasterAlignments.erase( (unmatchedMasterItr+1).base() );
                    matchFound = true;
                    break;
                }
            }
            if (!matchFound)
            {
                if (DEBUG>=2) cout << "No match found" << endl;
                ProcessUnmatchedSlave( *slaveItr ); // For Roman's stuff
                unmatchedSlaveAlignments.push_back( *slaveItr );
                ProcessFP( *slaveItr );
            }

            ++slaveItr;
        }
     }
    while (processMaster || processSlave);

    FinaliseFP();

    // Report all the unmatched master alignments
    BOOST_FOREACH( const io::bam::StorableBamAlignment& masterAlignment, unmatchedMasterAlignments )
    {
        if (masterAlignment.refId == 0xFFFFFFFF) { break; }
        ProcessFN( masterAlignment );
    }

    // End of BAM separation into FP/FN/TP. Below is for Roman's MapQ histogram
    if (!doRomanStuff) { return; }

    // Sort the 2 lists of items that we didn't manage to match yet, for matching by sequence
    if (DEBUG>=1) cout << (boost::format("Sorting %d unmatched master alignments by sequence") % unmatchedMasterAlignments.size()).str() << endl;
    sort( unmatchedMasterAlignments.begin(), unmatchedMasterAlignments.end(), io::bam::StorableBamAlignment::SeqCompareLt );
    if (DEBUG>=1) cout << (boost::format("Sorting %d unmatched slave alignments by sequence") % unmatchedSlaveAlignments.size()).str() << endl;
    sort( unmatchedSlaveAlignments.begin(), unmatchedSlaveAlignments.end(), io::bam::StorableBamAlignment::SeqCompareLt );

    // Go through the 2 lists to find the missing matches by comparing their sequences
    masterItr = unmatchedMasterAlignments.begin();
    slaveItr = unmatchedSlaveAlignments.begin();
    while (masterItr != unmatchedMasterAlignments.end() &&
           slaveItr  != unmatchedSlaveAlignments.end() )
    {
        if (DEBUG>=2) cout << "comparing " << *masterItr << "\n     with " << *slaveItr << endl;;
        switch (io::bam::StorableBamAlignment::SeqCompare( *masterItr, *slaveItr ))
        {
        case -1:
            if (DEBUG>=2) cout << " => -1" << endl;
            ++masterItr;
            break;
        case 1:
            if (DEBUG>=2) cout << " => 1" << endl;
            ++slaveItr;
            break;
        case 0:
            if (DEBUG>=2) cout << " => 0" << endl;
            // Check how many matches there are:
            //   If 1 slave seq matches 1 master seq, it's perfect
            //   If 1 slave seq matches multiple master seqs, we can't conclude anything so we skip them all
            //   If multiple slave seqs match 1 master seq, it's a bug - it should never happen
            //   If multiple slave seqs match multiple master seqs, we can't conclude anything so we skip them all
            if ( (masterItr+1 != unmatchedMasterAlignments.end() && io::bam::StorableBamAlignment::SeqCompare( *masterItr, *(masterItr+1) ) == 0) ||
                 (slaveItr+1 != unmatchedSlaveAlignments.end() && io::bam::StorableBamAlignment::SeqCompare( *slaveItr, *(slaveItr+1) ) == 0) )
            {
                // multiple matches => skip all
                if (DEBUG>=2) cout << "Multiple matches" << endl;
                while (masterItr+1 != unmatchedMasterAlignments.end() && io::bam::StorableBamAlignment::SeqCompare( *masterItr, *(masterItr+1) ) == 0)
                {
                    ++masterItr;
                }
                while (slaveItr+1 != unmatchedSlaveAlignments.end() && io::bam::StorableBamAlignment::SeqCompare( *slaveItr, *(slaveItr+1) ) == 0)
                {
                    ++slaveItr;
                }
            }
            else
            {
                UnprocessUnmatchedSlave( *slaveItr ); // For Roman's stuff
                ProcessMatch( *masterItr, *slaveItr );
                // Mark those alignments to be ignored later
                masterItr->refId = 0xFFFFFFFF;
                slaveItr->refId = 0xFFFFFFFF;
            }
            ++masterItr;
            ++slaveItr;
            break;
        default:
            assert( false && "should never reach here" );
        }
    }

    // Sort the 2 lists of remaining items that we didn't manage to match
    if (DEBUG>=1) cout << (boost::format("Sorting %d unmatched master alignments by position") % unmatchedMasterAlignments.size()).str() << endl;
    sort( unmatchedMasterAlignments.begin(), unmatchedMasterAlignments.end(), io::bam::StorableBamAlignment::PosCompareLt );
    if (DEBUG>=1) cout << (boost::format("Sorting %d unmatched slave alignments by position") % unmatchedSlaveAlignments.size()).str() << endl;
    sort( unmatchedSlaveAlignments.begin(), unmatchedSlaveAlignments.end(), io::bam::StorableBamAlignment::PosCompareLt );

/*
    // Report all the unmatched alignments
    BOOST_FOREACH( const io::bam::StorableBamAlignment& masterAlignment, unmatchedMasterAlignments )
    {
        if (masterAlignment.refId == 0xFFFFFFFF) { break; }
        ProcessFN( masterAlignment );
    }
    BOOST_FOREACH( io::bam::StorableBamAlignment& slaveAlignment, unmatchedSlaveAlignments )
    {
        if (slaveAlignment.refId == 0xFFFFFFFF) { break; }
        ProcessFP( slaveAlignment );
    }
*/

    // Final report
    Finalise();
}


void DiffComputer::ProcessMatch( const io::bam::StorableBamAlignment& masterAlignment, const io::bam::StorableBamAlignment& slaveAlignment )
{
    if (DEBUG>=2) cout << "Match found between master " << masterAlignment << " and slave " << slaveAlignment << endl;
    unsigned int mapQ = slaveAlignment.getMapQ();
    assert( mapQ < matchCountForMapQ.size() );
    ++(matchCountForMapQ[mapQ]);
}

void DiffComputer::ProcessUnmatchedMaster( const io::bam::StorableBamAlignment& alignment )
{
}

void DiffComputer::ProcessUnmatchedSlave( const io::bam::StorableBamAlignment& alignment )
{
    unsigned int mapQ = alignment.getMapQ();
    assert( mapQ < mismatchCountForMapQ.size() );
    ++(mismatchCountForMapQ[mapQ]);
}

void DiffComputer::UnprocessUnmatchedSlave( const io::bam::StorableBamAlignment& alignment )
{
    unsigned int mapQ = alignment.getMapQ();
    assert( mapQ < mismatchCountForMapQ.size() );
    --(mismatchCountForMapQ[mapQ]);
}


void DiffComputer::ProcessCorrect( const io::bam::StorableBamAlignment& masterAlignment, io::bam::StorableBamAlignment& slaveAlignment )
{
//    outCorrect << masterAlignment << "\n" << slaveAlignment << "\n" << endl;
    slaveAlignment.refId = slave2MasterRefId_[ slaveAlignment.refId ];
    slaveAlignment.nextRefId = slave2MasterRefId_[ slaveAlignment.nextRefId ];
    genome::EagleBamAlignmentAdapter eagleBamAlignmentAdapter( slaveAlignment );
    io::serializeAlignment( bgzfStreamTP_, eagleBamAlignmentAdapter );
}



void DiffComputer::ParseReadName( io::bam::StorableBamAlignment& alignment, unsigned int& lane, unsigned int& tile, unsigned long& cluster)
{
    //e.g. "unknown-flowcell_0:5:2208:775123:0"
    string readName = alignment.getReadNameAsString();
    clog << "ParseReadNameAndAddSimulatedPosInfo: readName=" << readName << endl;
    static bool firstTime = true;
    if (firstTime)
    {
        // Check that the format is what we expect
        assert( readName.substr( 0, 19 ) == "unknown-flowcell_0:" );
        firstTime = false;
    }
    char* ptr = (char*)readName.c_str()+19;
    lane = atoi(ptr); //boost::lexical_cast<unsigned int>( ptr );
    while (*(++ptr) != ':') {}
    ++ptr;
    tile = atoi(ptr); //boost::lexical_cast<unsigned int>( ptr );
    while (*(++ptr) != ':') {}
    ++ptr;
    cluster = atoi(ptr); //boost::lexical_cast<unsigned int>( ptr );
    clog << (boost::format("%d %d %d") % lane % tile % cluster).str() << endl;
}

void DiffComputer::ParseReadName( io::bam::StorableBamAlignment& alignment, unsigned int& fullTileNum, unsigned long& cluster)
{
    unsigned int lane;
    unsigned int tile;
    ParseReadName( alignment, lane, tile, cluster );
    fullTileNum = (lane-1)*32 +  + (((tile/1000)-1)%10)*16 + (((tile/100)-1)%10)*8 + (tile-1)%10;
}

void DiffComputer::ParseReadNameAndAddSimulatedPosInfo( io::bam::StorableBamAlignment& alignment )
{
    unsigned int lane;
    unsigned int tile;
    unsigned long cluster;
    ParseReadName( alignment, lane, tile, cluster );

    /*unsigned long globalSamplePos = */fragmentPosResolver_.GetSimulatedPosInSampleGenome( lane, tile, cluster );
}

struct Idea1Struct
{
    unsigned int fullTileNum_;
    unsigned long cluster_;
    io::bam::StorableBamAlignment alignment_;
    unsigned long sampleGlobalPos_;

    Idea1Struct( 
        unsigned int fullTileNum,
        unsigned long cluster,
        io::bam::StorableBamAlignment alignment
        )
        : fullTileNum_(fullTileNum)
        , cluster_(cluster)
        , alignment_(alignment)
        , sampleGlobalPos_(0xFFFFFFFFFFFFFFFF)
        {}

    static bool comparisonByClusterNumber( const Idea1Struct& lhs, const Idea1Struct& rhs )
    {
        return (lhs.cluster_ < rhs.cluster_);
    }
};
vector< vector< Idea1Struct > > idea1_FPs;
void DiffComputer::ProcessFP( io::bam::StorableBamAlignment& alignment )
{
//    outFP << alignment << endl;
/*
IDEA 1:
save each cluster#+ptr to alignment per lane+tile
reorder each of them per cluster#
go through fragments.tile and fragments.pos linearly and advance in all the lane+tile vector at the same time to find the sample global pos of each alignment
reorder everything together per sample global pos
add sample local pos info to each alignment
use segments mapping info to convert sample global pos to ref global pos
add ref local pos info to each alignment

OR

IDEA 2:
    // lane+tile+fragment => global sample pos
    //
    // global sample pos => global ref pos
keep a cache of lane+tile+fragment => global sample pos
 for each non-cached entry: start from closest smaller entry and go through the fragments.tile and fragments.pos files. Update cache
add sample local pos info to each alignment
use segments mapping info to convert sample global pos to ref global pos
add ref local pos info to each alignment
*/

#define IDEA2
#ifdef IDEA1
    clog << "Idea 1, checkpoint 1" << endl;
    unsigned int fullTileNum;
    unsigned long cluster;
    ParseReadName( alignment, fullTileNum, cluster );
    if (idea1_FPs.size() != 256+1)
    {
        idea1_FPs.resize(256+1);
    }
    assert( fullTileNum <= 256 );
    idea1_FPs[fullTileNum].push_back( Idea1Struct( fullTileNum, cluster, alignment ) );

    // test thingy:
    static int nb = 0;
    if (++nb>10)
    {
        FinaliseFP();
        idea1_FPs.clear();
    }
#endif //ifdef IDEA1

#ifdef IDEA2
    // Adjust slave alignment's refId to match the master's (as we use the master's chromosomes for the output)
    alignment.refId = slave2MasterRefId_[ alignment.refId ];
    alignment.nextRefId = slave2MasterRefId_[ alignment.nextRefId ];

/*
    // Add ("real") simulated position information
    {
        ParseReadNameAndAddSimulatedPosInfo( alignment );
    }
*/

    genome::EagleBamAlignmentAdapter eagleBamAlignmentAdapter( alignment );
    io::serializeAlignment( bgzfStreamFP_, eagleBamAlignmentAdapter );
#endif //ifdef IDEA2
}

void DiffComputer::FinaliseFP()
{
#ifdef IDEA1
    clog << "Idea 1, checkpoint 2" << endl;
    // save each cluster#+ptr to alignment per lane+tile
    //  done above

    // reorder each of them per cluster#
    BOOST_FOREACH( vector< Idea1Struct >& tileFPs, idea1_FPs )
    {
        std::sort( tileFPs.begin(), tileFPs.end(), Idea1Struct::comparisonByClusterNumber );
    }

    // go through fragments.tile and fragments.pos linearly and advance in all the lane+tile vector at the same time to find the sample global pos of each alignment
    vector<unsigned int> clusterNumPerTile(256);
    model::FragmentList fragmentList;
    while (1)
    {
        unsigned int tile;
        model::Fragment fragment = fragmentList.getNext( 0, 0, &tile );
        if (!fragment.isValid()) { break; }
        if (!idea1_FPs[tile].empty())
        {
            if (++(clusterNumPerTile[tile]) == idea1_FPs[tile][0].fullTileNum_)
            {
                EAGLE_WARNING( "We have a match! " <<  " : " << fragment );
            }
        }
    }

    // reorder everything together per sample global pos
    // add sample local pos info to each alignment
    // use segments mapping info to convert sample global pos to ref global pos
    // add ref local pos info to each alignment
    clog << "Idea 1, checkpoint 9" << endl;
    exit(0);
#endif //ifdef IDEA1
}

void DiffComputer::ProcessFN( const io::bam::StorableBamAlignment& alignment )
{
//    outFN << alignment << endl;
    genome::EagleBamAlignmentAdapter eagleBamAlignmentAdapter( alignment );
    io::serializeAlignment( bgzfStreamFN_, eagleBamAlignmentAdapter );
}


void DiffComputer::Finalise()
{
    for (unsigned int i=0; i<matchCountForMapQ.size(); ++i)
    {
        unsigned int matches = matchCountForMapQ[i];
        unsigned int mismatches = mismatchCountForMapQ[i];
        unsigned int total = matches + mismatches;
        if (total > 0)
        {
            double Q = -10.0 * log10( (double)mismatches / (double)total );
            cout << (boost::format("%d: %d+%d=%d => %f") % i % matches % mismatches % total % Q).str() << endl;
        }
    }
}


} // namespace main
} // namespace eagle
