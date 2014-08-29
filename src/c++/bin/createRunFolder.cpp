/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "main/RunFolderGeneratorOptions.hh"
#include "main/RunFolderGenerator.hh"


static void runFolderGeneratorLauncher(const eagle::main::RunFolderGeneratorOptions &options)
{
    eagle::main::RunFolderGenerator runFolderGenerator(options);
    runFolderGenerator.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(runFolderGeneratorLauncher, argc, argv);
}
