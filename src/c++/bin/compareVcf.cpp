/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "main/VcfComparatorOptions.hh"
#include "main/VcfComparator.hh"


static void vcfComparatorLauncher(const eagle::main::VcfComparatorOptions &options)
{
    eagle::main::VcfComparator vcfComparator(options);
    vcfComparator.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(vcfComparatorLauncher, argc, argv);
}
