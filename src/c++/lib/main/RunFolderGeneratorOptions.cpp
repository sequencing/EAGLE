/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'runFolderGenerator'
 **
 ** \author Lilian Janin
 **/

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include "main/RunFolderGeneratorOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

RunFolderGeneratorOptions::RunFolderGeneratorOptions()
    : outputDir("")
    , tileIdList("")
{
    parameters_.add_options()
        ("run-info,i",   bpo::value< bfs::path >(&runInfo),   "[input]  \tFull path to the RunInfo.xml file")
        ;

    namedOptions_.add_options()
        ("tile-id,t", bpo::value< std::string >(&tileIdList), "Comma-separated list of tile Ids")
        ;

    unnamedOptions_.add_options()
        ("output-dir,o", bpo::value< bfs::path >(&outputDir), "[output] \tOutput dir")
        ;
    positionalOptions_.add("output-dir",1);
}

void RunFolderGeneratorOptions::postProcess(bpo::variables_map &vm)
{
    eagle::common::OptionsHelper check(vm);

    check.requiredOptions(
        boost::assign::list_of("output-dir")("run-info")("tile-id")
    );

    boost::split( tileId, tileIdList, boost::is_any_of(",") );
}

} //namespace main
} // namespace eagle
