/**
 ** Copyright (c) 2018 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'fastqAnalyser'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "FastqAnalyserOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

FastqAnalyserOptions::FastqAnalyserOptions()
    : requestedMetrics( "11111111111111111111111111111111" )  // all metrics activated by default
    , requestedTables ( "00000000000000000000000000000000" )  // no table activated by default
{
    parameters_.add_options()
        ("fastq,f",              bpo::value< bfs::path >(&fastqFile),
                                "[input]  \tFASTQ file to analyse")
        ("reference-genome,r", bpo::value< bfs::path >(&referenceGenome),
                                "[input]  \tFull path to the reference genome FASTA files")
        ;

    namedOptions_.add_options()
        ("requested-metrics,m", bpo::value< std::string >(&requestedMetrics)->default_value(requestedMetrics), 
         "Binary value of requested metrics\n (see --help)"
            )
        ("requested-tables,t", bpo::value< std::string >(&requestedTables)->default_value(requestedTables),
         "Binary value of requested tables\n (see --help)"
            )
        ;
}

std::string FastqAnalyserOptions::usageSuffix() const
{
    return 
        "Metrics and table requests must be specified in reverse binary order (i.e. item 0 is first char).\n"
        "Available metrics:\n"
        "\n"
        "Available tables:\n"
        "  Table 7 (\"xxxxxxx1\"): Quality scores table\n"
        ;
}

void FastqAnalyserOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(boost::assign::list_of
                          ("fastq")
//                          ("reference-genome")
                          );
}

} //namespace main
} // namespace eagle
