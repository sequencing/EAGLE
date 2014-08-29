/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'homopolymerIndelVcfGenerator'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "HomopolymerIndelVcfGeneratorOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

HomopolymerIndelVcfGeneratorOptions::HomopolymerIndelVcfGeneratorOptions()
{
    parameters_.add_options()
        ("indel-probabilities,i", bpo::value< bfs::path >(&indelProbabilitiesFile),
                                "[input]  \tFile containing the indel probabilities")
        ("reference-genome,r"   , bpo::value< bfs::path >(&referenceGenome),
                                "[input]  \tFull path to the reference genome FASTA files")
        ;
}

std::string HomopolymerIndelVcfGeneratorOptions::usageSuffix() const
{
    return "";
}

void HomopolymerIndelVcfGeneratorOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(boost::assign::list_of
                          ("indel-probabilities")
                          ("reference-genome")
                          );
}

} //namespace main
} // namespace eagle
