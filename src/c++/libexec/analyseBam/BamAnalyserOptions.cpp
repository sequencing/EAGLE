/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'bamAnalyser'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "BamAnalyserOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

BamAnalyserOptions::BamAnalyserOptions()
    : requestedMetrics( "11111111111111111111111111111111" )  // all metrics activated by default
    , requestedTables ( "00000000000000000000000000000000" )  // no table activated by default
{
    parameters_.add_options()
        ("bam,b",              bpo::value< bfs::path >(&bamFile),
                                "[input]  \tBAM file to analyse")
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

std::string BamAnalyserOptions::usageSuffix() const
{
    return 
        "Metrics and table requests must be specified in reverse binary order (i.e. item 0 is first char).\n"
        "Available metrics:\n"
        "  Metric 0 (\"1\"):       Unused\n"
        "  Metric 1 (\"x1\"):      Levels where mismatch count per read crosses thresholds (1a=area threshold, 1b=level threshold)\n"
        "  Metric 2 (\"xx1\"):     2a=mismatch rate, 2b=insertion rate, 2c=deletion rate\n"
        "  Metric 3 (\"xxx1\"):    Standard deviation of mismatch rate over 10k windows\n"
        "  Metric 4 (\"xxxx1\"):   Standard deviation of coverage over 10k windows\n"
        "\n"
        "Available tables:\n"
        "  Table 0 (\"1\"):        Mismatch table\n"
        "  Table 1 (\"x1\"):       Homopolymer indel table\n"
        "  Table 2 (\"xx1\"):      Histogram data of number of mismatches per read\n"
        "  Table 3 (\"xxx1\"):     For each 10k window: {mismatch count, base count, mismatch rate}\n"
        "  Table 4 (\"xxxx1\"):    Unused\n"
        "  Table 5 (\"xxxxx1\"):   Areas with coverage >100\n"
        "  Table 6 (\"xxxxxx1\"):  Insert sizes histogram\n"
        ;
}

void BamAnalyserOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(boost::assign::list_of
                          ("bam")
                          ("reference-genome")
                          );
}

} //namespace main
} // namespace eagle
