/**
 ** Copyright (c) 2018 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "FastqAnalyserOptions.hh"
#include "FastqAnalyser.hh"


static void fastqAnalyserLauncher(const eagle::main::FastqAnalyserOptions &options)
{
    eagle::main::FastqAnalyser fastqAnalyser( options );
    fastqAnalyser.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(fastqAnalyserLauncher, argc, argv);
}
