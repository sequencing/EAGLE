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

#ifndef EAGLE_MAIN_SEQUENCER_SIMULATOR_OPTIONS_HH
#define EAGLE_MAIN_SEQUENCER_SIMULATOR_OPTIONS_HH

#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class SequencerSimulatorOptions : public eagle::common::Options
{
public:
    SequencerSimulatorOptions();
    std::map<std::string,unsigned int> exceptionPloidy() const;
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       simulateSequencer [parameters] [options]");}
    void postProcess(boost::program_options::variables_map &vm);

public:
    bool generateBclTile;
    bool generateFastqTile;
    bool generateBam;
    bool generateSampleBam;
    boost::filesystem::path runInfo;
    boost::filesystem::path sampleGenomeDir;
    boost::filesystem::path fragmentsDir;
    std::vector< boost::filesystem::path > qualityTableFiles;
    boost::filesystem::path mismatchTableFile;
    boost::filesystem::path homopolymerIndelTableFile;
    boost::filesystem::path motifQualityDropTableFile;
    boost::filesystem::path qqTableFile;
    boost::filesystem::path outDir;
    boost::filesystem::path outFilename;
    unsigned int readCount;
//    unsigned int binNum;
//    unsigned int binCount;
    unsigned int laneCount;
    unsigned int tilesPerLane;
    unsigned int lane;
    unsigned int tileNum;
    unsigned int tileId;
    unsigned int maxConcurrentWriters;
    unsigned int randomSeed;
    std::string bamRegion;
    bool dropLastBase;
    std::vector< std::string > errorModelOptions;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_SEQUENCER_SIMULATOR_OPTIONS_HH
