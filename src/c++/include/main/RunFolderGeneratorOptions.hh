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

#ifndef EAGLE_MAIN_RUN_FOLDER_GENERATOR_OPTIONS_HH
#define EAGLE_MAIN_RUN_FOLDER_GENERATOR_OPTIONS_HH

#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class RunFolderGeneratorOptions : public eagle::common::Options
{
public:
    RunFolderGeneratorOptions();
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       createRunFolder <path/to/RunFolder> [parameters] [options]");}
    std::string usageSuffix() const {std::stringstream usage;
                                     usage << std::endl;
                                     usage << "Note:" << std::endl;
                                     usage << "  -o [ --output-dir ] arg \t[output] Full path to the RunFolder" << std::endl;
                                     usage << "is also allowed for legacy purpose. In this case, the 'output-dir' parameter" << std::endl;
                                     usage << "replaces the positional argument <path/to/RunFolder>. For example:" << std::endl;
                                     usage << "       createRunFolder -i path/to/RunInfo.xml -o /path/to/RunFolder" << std::endl;
                                     usage << "    Or:" << std::endl;
                                     usage << "       createRunFolder --run-info /path/to/RunInfo.xml --output-dir /path/to/RunFolder" << std::endl;
                                     return usage.str();}
    void postProcess(boost::program_options::variables_map &vm);

public:
    boost::filesystem::path outputDir;
    boost::filesystem::path runInfo;
    std::vector< std::string > tileId;

private:
    std::string tileIdList;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_RUN_FOLDER_GENERATOR_OPTIONS_HH
