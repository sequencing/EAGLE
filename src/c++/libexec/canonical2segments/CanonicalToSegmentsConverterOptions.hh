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

#ifndef EAGLE_MAIN_CANONICAL_TO_SEGMENTS_CONVERTER_OPTIONS_HH
#define EAGLE_MAIN_CANONICAL_TO_SEGMENTS_CONVERTER_OPTIONS_HH

#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class CanonicalToSegmentsConverterOptions : public eagle::common::Options
{
public:
    CanonicalToSegmentsConverterOptions();
    std::map<std::string,unsigned int> exceptionPloidy() const;
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       canonical2segments [parameters] [options]");}
    void postProcess(boost::program_options::variables_map &vm);

public:
    boost::filesystem::path input;
    boost::filesystem::path outputDir;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_CANONICAL_TO_SEGMENTS_CONVERTER_OPTIONS_HH
