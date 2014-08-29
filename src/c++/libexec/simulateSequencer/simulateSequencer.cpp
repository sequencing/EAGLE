/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "SequencerSimulatorOptions.hh"
#include "SequencerSimulator.hh"


static void sequencerSimulatorLauncher(const eagle::main::SequencerSimulatorOptions &options)
{
    eagle::main::SequencerSimulator sequencerSimulator( options );
    sequencerSimulator.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(sequencerSimulatorLauncher, argc, argv);
}
