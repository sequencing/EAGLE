/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'vcfComparator'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "main/VcfComparatorOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

VcfComparatorOptions::VcfComparatorOptions()
{
    parameters_.add_options()
        ("simulated-variants,s", bpo::value< std::vector< bfs::path > >(&simulatedVariants),
                                 "[input]  \tFull path to the simulated variants VCF file (multiple occurrences allowed)")
        ("called-variants,c",    bpo::value< std::vector< bfs::path > >(&calledVariants),
                                 "[input]  \tFull path to the called variants VCF file (multiple occurrences allowed)")
        ;
}

void VcfComparatorOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(
        boost::assign::list_of("simulated-variants")("called-variants")
    );
}

} //namespace main
} // namespace eagle
