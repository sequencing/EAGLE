/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include <boost/ptr_container/ptr_deque.hpp>
#include <boost/format.hpp>
#include <boost/shared_ptr.hpp>
#include <fstream>

#include "common/Exceptions.hh"
#include "common/Logger.hh"
#include "common/Semaphore.hh"
#include "genome/SharedFastaReference.hh"
#include "genome/ReferenceToSample.hh"
#include "io/Bcl.hh"
#include "io/Fastq.hh"
#include "genome/BamMetadata.hh"
#include "model/Fragment.hh"
#include "model/PassFilter.hh"
#include "SequencerSimulator.hh"


using namespace std;
using eagle::io::BclTile;
using eagle::io::FastqTile;


namespace eagle
{
namespace main
{

SequencerSimulator::SequencerSimulator (const SequencerSimulatorOptions &options )
    : options_        ( options )
    , runInfo_        ( options.runInfo )
    , tileNum_        ( (options_.lane-1)*options_.tilesPerLane + (options_.tileNum-1) )
    , fragmentList_   ( options_.fragmentsDir /*options.binNum*(options.readCount/options.binCount)*/ )
    , readClusterFactory_( runInfo_, options_.sampleGenomeDir, options_.qualityTableFiles, options_.mismatchTableFile, options_.homopolymerIndelTableFile, options_.motifQualityDropTableFile, options_.qqTableFile, options_.randomSeed, options_.errorModelOptions )
{
}

void SequencerSimulator::run()
{
    if (options_.generateBclTile)
    {
        generateBclTile();
    }

    if (options_.generateFastqTile)
    {
        generateFastqTile();
    }

    if (options_.generateBam)
    {
        generateBam();
    }

    if (options_.generateSampleBam)
    {
        generateSampleBam();
    }
}


void SequencerSimulator::generateBclTile()
{    // for each tile in the to-be-processed set:
    unsigned long long tileReadCount = fragmentList_.getTileSize( tileNum_ );
    clog << "SequencerSimulator::generateBclTile: tile=" << tileNum_ << ", readCount=" << tileReadCount << endl;

    unsigned int tileId = options_.tileId;
    string bclFilenameTemplate = (boost::format("%s/Data/Intensities/BaseCalls/L%03g/C%%d.1/s_%d_%g.bcl") % options_.outDir.string() % options_.lane % options_.lane % tileId).str();
    string statsFilenameTemplate = (boost::format("%s/Data/Intensities/BaseCalls/L%03g/C%%d.1/s_%d_%g.stats") % options_.outDir.string() % options_.lane % options_.lane % tileId).str();
    string filterFilename = (boost::format("%s/Data/Intensities/BaseCalls/L00%d/s_%d_%04g.filter") % options_.outDir.string() % options_.lane % options_.lane % tileId).str();
    string posFilename = (boost::format("%s/Data/Intensities/s_%d_%04g_pos.txt") % options_.outDir.string() % options_.lane % tileId).str();
    string clocsFilename = (boost::format("%s/Data/Intensities/L00%d/s_%d_%04g.clocs") % options_.outDir.string() % options_.lane % options_.lane % tileId).str();
    string controlFilename = (boost::format("%s/Data/Intensities/BaseCalls/L00%d/s_%d_%04g.control") % options_.outDir.string() % options_.lane % options_.lane % tileId).str();
    unsigned int clusterLength = runInfo_.getClusterLength();
    BclTile bclTile( tileReadCount, clusterLength, bclFilenameTemplate, statsFilenameTemplate, filterFilename, clocsFilename, controlFilename );

    for (unsigned int i=0; i<tileReadCount; ++i)
    {
        // Read next paired read position for our tile(s) of interest
        eagle::model::Fragment fragment = fragmentList_.getNext( tileNum_ );
        eagle::genome::ReadClusterWithErrors readClusterWithErrors = readClusterFactory_.getReadClusterWithErrors( fragment );

        // Apply errors
        //automatically done by class//            readClusterWithErrors.generateErrors();

        // Output BCL
        const char *bclCluster = readClusterWithErrors.getBclCluster();
        bool isPassingFilter = model::PassFilter::isBclClusterPassingFilter( bclCluster, clusterLength );
        bclTile.addClusterToRandomLocation( bclCluster, isPassingFilter );
    }

    // Flush tile to disk (TODO: ...in parallel with next tile's creation)
    if (options_.maxConcurrentWriters > 0)
    {
        clog << "Ready to flush tile. Waiting for semaphore." << endl;
        eagle::common::Semaphore semaphore( "EagleSemaphore", options_.maxConcurrentWriters );
        semaphore.wait();
        bclTile.flushToDisk();
        semaphore.post();
    }
    else
    {
        bclTile.flushToDisk();
    }
}

void SequencerSimulator::generateFastqTile()
{    // for each tile in the to-be-processed set:
    unsigned long long tileReadCount = fragmentList_.getTileSize( tileNum_ );
    clog << "SequencerSimulator::generateFastqTile: tile=" << tileNum_ << ", readCount=" << tileReadCount << endl;

    //unsigned int tileId = options_.tileId;
    string fastqFilenameTemplate = (boost::format("%s/EAGLE_S%d_L%03g_R%%d_001.fastq") % options_.outDir.string() % (tileNum_ + 1) % options_.lane).str();
    string read1FastqFilename = (boost::format(fastqFilenameTemplate) % 1).str();
    string read2FastqFilename = (boost::format(fastqFilenameTemplate) % 2).str();
    unsigned int clusterLength = runInfo_.getClusterLength();
    FastqTile fastqTile( tileReadCount, clusterLength, read1FastqFilename, read2FastqFilename, runInfo_, options_.lane, options_.tileId );

    for (unsigned int i=0; i<tileReadCount; ++i)
    {
        // Read next paired read position for our tile(s) of interest
        eagle::model::Fragment fragment = fragmentList_.getNext( tileNum_ );
        eagle::genome::ReadClusterWithErrors readClusterWithErrors = readClusterFactory_.getReadClusterWithErrors( fragment );

        // Output FASTQ
        const char *bclCluster = readClusterWithErrors.getBclCluster();
        bool isPassingFilter = model::PassFilter::isBclClusterPassingFilter( bclCluster, clusterLength );

        const string read1Nucleotides = readClusterWithErrors.getNucleotideOrQualitySequenceForRead(0,true,false,true);
        const string read1Qualities = readClusterWithErrors.getNucleotideOrQualitySequenceForRead(0,false,false,true);
        const string read2Nucleotides = readClusterWithErrors.getNucleotideOrQualitySequenceForRead(1,true,false,true);
        const string read2Qualities = readClusterWithErrors.getNucleotideOrQualitySequenceForRead(1,false,false,true);

        fastqTile.addCluster( read1Nucleotides, read1Qualities, read2Nucleotides, read2Qualities, isPassingFilter );
    }

    fastqTile.finaliseAndWriteInfo();
}

/*
bool fetchNextFragment( const genome::RefToSampleSegment& segment,
                        const unsigned long currentPos,
                        const unsigned long lastPosToProcess,
                        model::Fragment& resultFragment )
{
    return false;
}
*/

void SequencerSimulator::generateBam()
{
    clock_t startTime = clock();
    const string currentChr( options_.bamRegion );
    genome::RefToSampleSegmentReader refToSampleSegmentReader( options_.sampleGenomeDir / "segmentsFromRef.tsv", currentChr );
//    refToSampleSegmentReader.jumpToChromosome( currentChr );
    vector< genome::RefToSampleSegment > segmentsToMerge;
    vector< model::Fragment > fragmentForSegment;
    vector< model::FragmentList* > fragmentListForSegment;
    vector< string > sharedFastReaderIdForSegment;
    unsigned long currentPos = 0;

    PreferredFastaReader mainReferenceGenome( options_.sampleGenomeDir / ".." / "reference_genome" );
    eagle::genome::BamOrMetadataOutput bamOrMetadataOutput( options_.outDir / options_.outFilename, runInfo_, &mainReferenceGenome /*, currentChr*/);

    // Find the global ref position of the desired chromosmes
    unsigned long chrGlobalPosInRef = 0;
    const vector<string> refContigNames = mainReferenceGenome.allContigNames();
    const vector<unsigned long> refContigLengths = mainReferenceGenome.allContigLengths();
    bool refChrFound = false;
    for (unsigned int i=0; i<refContigNames.size(); ++i)
    {
        if (refContigNames[i] == currentChr)
        {
            refChrFound = true;
            break;
        }
        else
        {
            chrGlobalPosInRef += refContigLengths[i];
        }
    }
    assert( refChrFound && "Error: specified chromosome not found in reference genome" );
    (void)refChrFound; // prevents "variable set but not used" warning when asserts are not compiled

    genome::RefToSampleSegment refToSampleSegment;
    bool segmentAvailable;
    do
    {
        segmentAvailable = refToSampleSegmentReader.getNextSegment( refToSampleSegment );
        if (segmentAvailable && ( currentPos == refToSampleSegment.refPos_ || segmentsToMerge.size() == 0 ) )
        {
            // Add this segment to the list of segments to be merged together
            segmentsToMerge.push_back( refToSampleSegment );
            fragmentForSegment.push_back( model::Fragment() );
            unsigned long firstGlobalPosInSample = refToSampleSegment.getSampleGlobalStartPos();
            unsigned long lastGlobalPosInSample = refToSampleSegment.getSampleGlobalEndPos();
            fragmentListForSegment.push_back( new model::FragmentList( options_.fragmentsDir, firstGlobalPosInSample, lastGlobalPosInSample, 500 ) ); //!!!!!!!!get this 500 dynamically
            cout << "Adding allele to merge: " << refToSampleSegment << " (sample global pos range: " << firstGlobalPosInSample << "-" << lastGlobalPosInSample << ")" << endl;
            currentPos = refToSampleSegment.refPos_;
            sharedFastReaderIdForSegment.push_back( refToSampleSegment.sampleChrAllele_ );
        }
        else
        {
            unsigned long lastPosToProcess;
            if (segmentAvailable)
            {
                lastPosToProcess = refToSampleSegment.refPos_ - 1;
                refToSampleSegmentReader.goBack( 1 );
            }
            else
            {
                lastPosToProcess = 0xFFFFFFFFFFFFFFFFul - 1;
            }

            // Start processing segments merge
            cout << (boost::format("Processing %s from %d to %d: merging %d allele(s)") % currentChr % currentPos % lastPosToProcess % segmentsToMerge.size()).str() << endl;
            do
            {
                unsigned long minFragmentPos = 0xFFFFFFFFFFFFFFFFul;
                unsigned int minFragmentPosIndex = 0;
                signed long minFragmentGlobalPosShift = 0;
                for (unsigned int i=0; i<segmentsToMerge.size(); ++i)
                {
                    if (fragmentForSegment[i].isValid()
                        || fragmentListForSegment[i]->getNext( fragmentForSegment[i] ))
                    {
                        assert( fragmentForSegment[i].isValid() );
                        unsigned long firstGlobalPosInSample = segmentsToMerge[i].getSampleGlobalStartPos();
                        unsigned long posInRef = fragmentForSegment[i].startPos_ - firstGlobalPosInSample + segmentsToMerge[i].refPos_;
                        if (posInRef < minFragmentPos)
                        {
                            minFragmentPos = posInRef;
                            minFragmentPosIndex = i;
                            minFragmentGlobalPosShift = segmentsToMerge[i].refPos_ - firstGlobalPosInSample + chrGlobalPosInRef - 1;
                        }
                    }
                    else
                    {
                        cout << "Removing allele to merge: " << segmentsToMerge[i] << endl;
                        segmentsToMerge.erase( segmentsToMerge.begin()+i );
                        fragmentForSegment.erase( fragmentForSegment.begin()+i );
                        if (fragmentListForSegment[i])
                        {
                            delete fragmentListForSegment[i];
                        }
                        fragmentListForSegment.erase( fragmentListForSegment.begin()+i );
                        sharedFastReaderIdForSegment.erase( sharedFastReaderIdForSegment.begin()+i );
                        --i;
                    }
                }

                // Add read to BAM
                if (minFragmentPos <= lastPosToProcess && minFragmentPos != 0xFFFFFFFFFFFFFFFFul)
                {
                    eagle::genome::ReadClusterWithErrors readClusterWithErrors = readClusterFactory_.getReadClusterWithErrors( fragmentForSegment[minFragmentPosIndex] );
                    static int count = 0;
//                    if (++count > 1000000) { return; }
                    clock_t currentTime = clock();
                    if (currentTime - startTime > 20000 )
                    {
                        clog << "time to add BAM line " << count  << ": " << eagle::common::displayTime(currentTime - startTime) << endl;
                    }
                    startTime = currentTime;
                    genome::SharedFastaReference::setActive( sharedFastReaderIdForSegment[minFragmentPosIndex] );

                    unsigned long firstPosToProcess = segmentsToMerge[minFragmentPosIndex].refPos_-1 + chrGlobalPosInRef;
                    unsigned long lastPosToProcess  = segmentsToMerge[minFragmentPosIndex].getRightMostRefPos()-1 + chrGlobalPosInRef;
                    bamOrMetadataOutput.addRebased( readClusterWithErrors, minFragmentGlobalPosShift, firstPosToProcess, lastPosToProcess, options_.dropLastBase, segmentsToMerge[minFragmentPosIndex] );
                    fragmentForSegment[minFragmentPosIndex].invalidate(); // to make it load the next fragment at the next loop
                    currentPos = minFragmentPos;
                }
                else
                {
                    currentPos = lastPosToProcess + 1; 
                }
            }
            while (currentPos <= lastPosToProcess);

            cout << (boost::format("Finished processing %s ") % currentChr).str();
            if (lastPosToProcess == 0xFFFFFFFFFFFFFFFFul - 1)
            {
                cout << "until the end";
            }
            else
            {
                cout << (boost::format("until %d") % lastPosToProcess).str();
            }
            cout << (boost::format(": %d allele(s) left") % segmentsToMerge.size()).str() << endl;

            currentPos = lastPosToProcess + 1;
        }
    }
    while (segmentAvailable);
}


// GenerateSampleBam: Generate BAM file aligned on the sample genome, instead of being aligned on the reference genome
void SequencerSimulator::generateSampleBam()
{
    unsigned long long readCount = fragmentList_.size(); // fragmentList_.getTileSize( tileNum_ );
    clog << "SequencerSimulator::generateBam: readCount=" << readCount << endl;

    eagle::genome::BamOrMetadataOutput bamOrMetadataOutput( options_.outDir / options_.outFilename, runInfo_ );

    for (unsigned long long i=0; i<readCount; ++i)
    {
//        if (i>100000) { break; } // TODO: remove this, it was for testing

        // Read next paired read position for any tile
        eagle::model::Fragment fragment = fragmentList_.getNext();
        eagle::genome::ReadClusterWithErrors readClusterWithErrors = readClusterFactory_.getReadClusterWithErrors( fragment );

        // Apply errors
        //automatically done by class//            readClusterWithErrors.generateErrors();

        // Output BAM
        bamOrMetadataOutput.add(readClusterWithErrors);
    }
}


} // namespace main
} // namespace eagle
