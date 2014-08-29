/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "BamAnalyserOptions.hh"
#include "BamAnalyser.hh"


static void bamAnalyserLauncher(const eagle::main::BamAnalyserOptions &options)
{
    eagle::main::BamAnalyser bamAnalyser( options );
    bamAnalyser.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(bamAnalyserLauncher, argc, argv);
}
