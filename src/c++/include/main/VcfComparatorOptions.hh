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

#ifndef EAGLE_MAIN_VCF_COMPARATOR_OPTIONS_HH
#define EAGLE_MAIN_VCF_COMPARATOR_OPTIONS_HH

#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class VcfComparatorOptions : public eagle::common::Options
{
public:
    VcfComparatorOptions();
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       compareVcf [parameters] [options]");}
    void postProcess(boost::program_options::variables_map &vm);

public:
    std::vector< boost::filesystem::path > simulatedVariants;
    std::vector< boost::filesystem::path > calledVariants;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_VCF_COMPARATOR_OPTIONS_HH
