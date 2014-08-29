/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'canonicalToSegmentsConverter'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "CanonicalToSegmentsConverterOptions.hh"


namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

CanonicalToSegmentsConverterOptions::CanonicalToSegmentsConverterOptions()
    : input()
    , outputDir()
{
    parameters_.add_options()
        ("input,i",             bpo::value< bfs::path >(&input),
                                "[input]  \tFull path to the canonical.vcf file")
        ("output-dir,o",        bpo::value< bfs::path >(&outputDir)->default_value(outputDir),
                                "[output] \tFull path to the output directory")
        ;
}

void CanonicalToSegmentsConverterOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(boost::assign::list_of
                          ("input")
                          ("output-dir")
                          );
}

} //namespace main
} // namespace eagle
