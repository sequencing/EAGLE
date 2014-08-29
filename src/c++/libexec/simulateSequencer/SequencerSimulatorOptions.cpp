/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'sequencerSimulator'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "SequencerSimulatorOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

SequencerSimulatorOptions::SequencerSimulatorOptions()
    : generateBclTile(false)
    , generateBam(false)
    , generateSampleBam(false)
    , sampleGenomeDir()
    , fragmentsDir()
    , mismatchTableFile("")
    , homopolymerIndelTableFile("")
    , motifQualityDropTableFile("")
    , qqTableFile("")
    , outDir()
    , outFilename()
    , readCount(1000000)
//    , binNum(0)
//    , binCount(1)
    , laneCount(8)
    , tilesPerLane(32)
    , lane(0)
    , tileNum(0)
    , tileId(0)
    , maxConcurrentWriters(0)
    , randomSeed(1)
    , dropLastBase(false)
{
    parameters_.add_options()
        ("run-info",            bpo::value< bfs::path >(&runInfo),
                                "[input]  \tFull path to the RunInfo.xml file")
        ("sample-genome-dir,s", bpo::value< bfs::path >(&sampleGenomeDir),
                                "[input]  \tFull path to the directory containing the sample's genome FASTA files")
        ("fragments-dir,f",  bpo::value< bfs::path >(&fragmentsDir),
                                "[input]  \tFull path to the directory containing the fragments.* files")
        ("quality-table,q",     bpo::value< std::vector< bfs::path > >(&qualityTableFiles),
                                "[input]  \tFile containing the quality table: 1 line per cycle, tab-separated pairs \"quality:occurrences\" items")
        ("mismatch-table",      bpo::value< bfs::path >(&mismatchTableFile),
                                "[input]  \tFile containing the mismatch table (default: equal probabilities for each SNP, no indel)")
        ("homopolymer-indel-table",      bpo::value< bfs::path >(&homopolymerIndelTableFile),
                                "[input]  \tFile containing the homopolymer indel table (default: no indel)")
        ("motif-quality-drop-table",   bpo::value< bfs::path >(&motifQualityDropTableFile),
                                "[input]  \tFile containing the motif quality drop table (default: no quality drop)")
        ("qq-table",            bpo::value< bfs::path >(&qqTableFile),
                                "[input]  \tFile containing the QQ table (default: Phred values: error-rate=10^(-Q/10))")
        ("output-dir,o",        bpo::value< bfs::path >(&outDir)->default_value(outDir),
                                "[output] \tFull path to the output directory")
        ("output-filename",     bpo::value< bfs::path >(&outFilename),
                                "[output] \tFile name to be used for BAM output")
        ;
    namedOptions_.add_options()
        ("generate-bcl-tile", bpo::value< bool >(&generateBclTile)->zero_tokens(), "Generates BCL tile identified by the following parameters:")
        ("read-count,n", bpo::value<unsigned int>(&readCount)->default_value(readCount), "Number of reads")
        ("lane-count,m", bpo::value<unsigned int>(&laneCount)->default_value(laneCount), "Number of lanes")
        ("tiles-per-lane,u", bpo::value<unsigned int>(&tilesPerLane)->default_value(tilesPerLane), "Number of tiles per lane")
        ("lane,l", bpo::value<unsigned int>(&lane)->default_value(lane), "Lane of the tile to be processed (1-based value)")
        ("tile-num", bpo::value<unsigned int>(&tileNum)->default_value(tileNum), "Tile number to be processed in the specified lane (1-based value)")
        ("tile-id", bpo::value<unsigned int>(&tileId)->default_value(tileId), "Tile id corresponding to the provided tile number for the desired naming scheme")
        ("max-concurrent-writers", bpo::value<unsigned int>(&maxConcurrentWriters)->default_value(maxConcurrentWriters), "Number of EAGLE processes allowed to flush their tile simultaneously (0=unlimited). This is per computer. Some disks exhibit better performance when this is set to 1.")
        ("random-seed", bpo::value<unsigned int>(&randomSeed)->default_value(randomSeed), "Multiplier used to calculate the actual seeds used for the generation of mismatches for each read")
        ("generate-bam", bpo::value< bool >(&generateBam)->zero_tokens(), "Generates BAM file aligned on the reference genome")
        ("generate-sample-bam", bpo::value< bool >(&generateSampleBam)->zero_tokens(), "Generates BAM file aligned on the sample genome")
        ("bam-region", bpo::value<std::string>(&bamRegion), "Bam region to generate (e.g. chr1 or chr1:1000-2000)")
        ("drop-last-base", bpo::value< bool >(&dropLastBase)->zero_tokens(), "Don't include the last base of each read in BAM output (e.g. read length 101 becomes 100)")
        ("error-model-options", bpo::value< std::vector< std::string > >(&errorModelOptions), "Used to initialise an error model plugin. value should be plugin-name:key=value:key=value:etc.\nDefault values:\n LONGREAD-deletion:prob=0.0:dist-file=filename\n LONGREAD-base-duplication:prob=0.0")
//        ("bin-num", bpo::value<unsigned int>(&binNum)->default_value(binNum), "Bin number to generate")
//        ("bin-count", bpo::value<unsigned int>(&binCount)->default_value(binCount), "Total number of bins")
        ;
}

void SequencerSimulatorOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    std::string mode = check.mutuallyExclusiveOptions(
        boost::assign::list_of("generate-bcl-tile")("generate-bam")("generate-sample-bam")
        );

    check.requiredOptions(boost::assign::list_of
                          ("run-info")
                          ("sample-genome-dir")
                          ("fragments-dir")
                          ("quality-table")
                          );

    if (mode == "generate-bcl-tile")
    {
        check.requiredOptions(boost::assign::list_of
                              ("lane")
                              ("tile-num")
                              ("tile-id")
                              );
    }
    else if (mode == "generate-bam")
    {
        check.requiredOptions(boost::assign::list_of
                              ("output-filename")
                              ("bam-region")
                              );
    }
    else // if (mode == "generate-sample-bam")
    {
        check.requiredOptions(boost::assign::list_of
                              ("output-filename")
                              );
    }
}

} //namespace main
} // namespace eagle
