/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#ifndef EAGLE_MAIN_SEQUENCER_SIMULATOR_HH
#define EAGLE_MAIN_SEQUENCER_SIMULATOR_HH

#include "io/RunInfo.hh"
#include "genome/EnrichedFragment.hh"
#include "genome/ReadCluster.hh"
#include "SequencerSimulatorOptions.hh"


namespace eagle
{
namespace main
{


class SequencerSimulator
{
public:
    SequencerSimulator (const SequencerSimulatorOptions &options);
    void run();
    void generateBclTile();
    void generateFastqTile();
    void generateBam();
    void generateSampleBam();

private:
    const SequencerSimulatorOptions &options_;
    eagle::io::RunInfo runInfo_;
    unsigned int tileNum_;
    eagle::model::FragmentList fragmentList_;
    eagle::genome::ReadClusterFactory readClusterFactory_;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_SEQUENCER_SIMULATOR_HH
