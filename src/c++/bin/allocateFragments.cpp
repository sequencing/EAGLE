/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "main/FragmentsAllocatorOptions.hh"
#include "main/FragmentsAllocator.hh"


static void allocateFragmentsLauncher(const eagle::main::FragmentsAllocatorOptions &options)
{
    eagle::main::FragmentsAllocator fragmentsAllocator( options );
    fragmentsAllocator.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(allocateFragmentsLauncher, argc, argv);
}
