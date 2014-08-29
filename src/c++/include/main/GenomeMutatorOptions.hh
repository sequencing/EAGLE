/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Command line options for 'variantModelling'
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_MAIN_GENOME_MUTATOR_OPTIONS_HH
#define EAGLE_MAIN_GENOME_MUTATOR_OPTIONS_HH

#include <string>
#include <map>
#include <boost/filesystem.hpp>

#include "common/Program.hh"

namespace eagle
{
namespace main
{

class GenomeMutatorOptions : public eagle::common::Options
{
public:
    enum Modes
    {
         SAFE_MODE,
         WHOLE_DIR
    };
    GenomeMutatorOptions();
    std::map<std::string,unsigned int> exceptionPloidy() const;
private:
    std::string usagePrefix() const {return std::string("Usage:\n")
                                          + std::string("       applyVariants [parameters] [options]");}
    std::string usageSuffix() const {std::stringstream usage;
                                     usage << "Examples:" << std::endl;
                                     usage << " * Safe Mode" << std::endl;
                                     usage << "       applyVariants -v /path/to/VariantList.vcf \\" << std::endl;
                                     usage << "                     -r /path/to/ReferenceDir/reference_1.fa \\" << std::endl;
                                     usage << "                     -r /path/to/ReferenceDir/reference_2.fa \\" << std::endl;
                                     usage << "                     ... etc ... \\" << std::endl;
                                     usage << "                     [options]" << std::endl;
                                     usage << " * Whole-dir Mode" << std::endl;
                                     usage << "       applyVariants -v /path/to/VariantList.vcf \\" << std::endl;
                                     usage << "                     -R /path/to/ReferenceDir \\" << std::endl;
                                     usage << "                     [options]" << std::endl;
                                     return usage.str();}
    void postProcess(boost::program_options::variables_map &vm);

public:
    std::vector<boost::filesystem::path> referenceGenome;
    boost::filesystem::path wholeGenome;
    boost::filesystem::path sampleGenome;
    std::vector<boost::filesystem::path> variantList;
    boost::filesystem::path annotatedVariantList;
    unsigned int organismPloidy;
    std::vector<std::string> ploidyChromosome;
    std::vector<unsigned int> ploidyLevel;
    std::string prefixToAdd;
    bool noTranslocationError;
    bool onlyPrintOutputContigNames;
    //bool withGenomeSize;
    bool force;

    Modes mode;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_GENOME_MUTATOR_OPTIONS_HH
