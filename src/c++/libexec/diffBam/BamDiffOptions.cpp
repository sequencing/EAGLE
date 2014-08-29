/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'bamDiff'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "BamDiffOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

BamDiffOptions::BamDiffOptions()
{
    parameters_.add_options()
        ("master-bam,i"      , bpo::value< bfs::path >(&masterBamFile),
                                "[input]  \tMain BAM file, assumed to contain the true alignments")
        ("slave-bam,j"       , bpo::value< bfs::path >(&slaveBamFile),
                                "[input]  \tBAM file assumed to contain some incorrect alignments")
        ;

/*
    namedOptions_.add_options()
        ("requested-metrics,m", bpo::value< std::string >(&requestedMetrics)->default_value(requestedMetrics), 
         "Binary value of requested metrics\n (see --help)"
            )
        ("requested-tables,t", bpo::value< std::string >(&requestedTables)->default_value(requestedTables),
         "Binary value of requested tables\n (see --help)"
            )
        ;
*/
}

std::string BamDiffOptions::usageSuffix() const
{
    return "";
}

void BamDiffOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(boost::assign::list_of
                          ("master-bam")
                          ("slave-bam")
                          );
}

} //namespace main
} // namespace eagle
