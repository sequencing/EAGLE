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

class BamDiffOptions : public eagle::common::Options
{
public:
    BamDiffOptions();
    std::map<std::string,unsigned int> exceptionPloidy() const;
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       diffBam [parameters] [options]");}
    std::string usageSuffix() const;
    void postProcess(boost::program_options::variables_map &vm);

public:
    boost::filesystem::path masterBamFile;
    boost::filesystem::path slaveBamFile;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_BAM_ANALYSER_OPTIONS_HH
