/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "HomopolymerIndelVcfGeneratorOptions.hh"
#include "HomopolymerIndelVcfGenerator.hh"


static void homopolymerIndelVcfGeneratorLauncher(const eagle::main::HomopolymerIndelVcfGeneratorOptions &options)
{
    eagle::main::HomopolymerIndelVcfGenerator homopolymerIndelVcfGenerator( options );
    homopolymerIndelVcfGenerator.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(homopolymerIndelVcfGeneratorLauncher, argc, argv);
}
