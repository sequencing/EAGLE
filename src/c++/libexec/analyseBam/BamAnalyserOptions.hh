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

#ifndef EAGLE_MAIN_BAM_ANALYSER_OPTIONS_HH
#define EAGLE_MAIN_BAM_ANALYSER_OPTIONS_HH

#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class BamAnalyserOptions : public eagle::common::Options
{
public:
    BamAnalyserOptions();
    std::map<std::string,unsigned int> exceptionPloidy() const;
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       analyseBam [parameters] [options]");}
    std::string usageSuffix() const;
    void postProcess(boost::program_options::variables_map &vm);

public:
    boost::filesystem::path bamFile;
    boost::filesystem::path referenceGenome;
    std::string requestedMetrics;
    std::string requestedTables;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_BAM_ANALYSER_OPTIONS_HH
