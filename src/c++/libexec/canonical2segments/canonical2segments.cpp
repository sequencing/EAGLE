/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \author Lilian Janin
 **/

#include "CanonicalToSegmentsConverterOptions.hh"
#include "CanonicalToSegmentsConverter.hh"


static void canonicalToSegmentsConverterLauncher(const eagle::main::CanonicalToSegmentsConverterOptions &options)
{
    eagle::main::CanonicalToSegmentsConverter canonicalToSegmentsConverter( options );
    canonicalToSegmentsConverter.run();
}


int main(int argc, char *argv[])
{
    eagle::common::run(canonicalToSegmentsConverterLauncher, argc, argv);
}
