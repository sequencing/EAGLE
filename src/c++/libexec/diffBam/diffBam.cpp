/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "BamDiffOptions.hh"
#include "BamDiff.hh"


static void bamDiffLauncher(const eagle::main::BamDiffOptions &options)
{
    eagle::main::BamDiff bamDiff( options );
    bamDiff.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(bamDiffLauncher, argc, argv);
}
