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

#include <string>
#include <vector>
#include <boost/assign.hpp>
#include <boost/foreach.hpp>
#include <boost/assert.hpp>

#include <boost/format.hpp>

#include "main/GenomeMutatorOptions.hh"

namespace eagle
{
namespace main
{

namespace bpo = boost::program_options;
namespace bfs = boost::filesystem;

GenomeMutatorOptions::GenomeMutatorOptions()
    : referenceGenome()
    , wholeGenome()
    , sampleGenome(bfs::current_path() / "sample_genome")
    , variantList()
    , annotatedVariantList(bfs::current_path() / "canonical.vcf")
    , organismPloidy(1)
    , ploidyChromosome()
    , ploidyLevel()
    , prefixToAdd("")
    , noTranslocationError(false)
    , onlyPrintOutputContigNames(false)
    //, withGenomeSize(false)
    , force(false)
{
    parameters_.add_options()
        ("reference-genome,r",       bpo::value< std::vector<bfs::path> >(&referenceGenome),
                                     "[input]  \tFull path to the reference genome FASTA file (multiple references allowed)")
        ("whole-genome,R",           bpo::value< bfs::path >(&wholeGenome),
                                     "[input]  \tFull path to the reference genome dir (single directory containing multiple FASTA files)")
        ("variant-list,v",           bpo::value< std::vector<bfs::path> >(&variantList),
                                     "[input]  \tFull path to the file containing the list of variants in VCF format (multiple lists allowed)")
        ("sample-genome,s",          bpo::value< bfs::path >(&sampleGenome)->default_value( sampleGenome, "./sample_genome" ),
                                     "[output] \tFull path to the output directory that will contain the sample reference (may write multiple FASTA files)")
        ("annotated-variant-list,a", bpo::value< bfs::path >(&annotatedVariantList)->default_value( annotatedVariantList, "./canonical.vcf" ),
                                     "[output] \tFull path to the annotated variant list (single VCF file)")
        ;
    namedOptions_.add_options()
        ("organism-ploidy,p",        bpo::value< unsigned int >(&organismPloidy)->default_value(organismPloidy),
                                     "Default ploidy level: HAPLOID(1), DIPLOID(2), TRIPLOID(3), TETRAPLOID(4), etc.")
        ("ploidy-chromosome,c",      bpo::value< std::vector<std::string> >(&ploidyChromosome),
                                     "Name of chromosome to be forced at a ploidy 'ploidy-level'")
        ("ploidy-level,l",           bpo::value< std::vector<unsigned int> >(&ploidyLevel),
                                     "Level of ploidy for chromosome 'ploidy-chromosome'")
        ("prefix",                   bpo::value< std::string >(&prefixToAdd),
                                     "Prefix to add to the output contig (and file) names")
        ("no-translocation-error",   bpo::value< bool >(&noTranslocationError)->zero_tokens(),
                                     "Do not issue a final error when translocations didn't get applied")
        ("only-print-output-contig-names",
                                     bpo::value< bool >(&onlyPrintOutputContigNames)->zero_tokens(),
                                     "Print output contig names and exit")
        //("with-genome-size,g",       bpo::value< bool >(&withGenomeSize)->default_value(false)->zero_tokens(),
        //                             "Read FASTA metadata from (and generate) a genome_size.xml file, instead of the usual *.fai")
        ;
}

std::map<std::string,unsigned int> GenomeMutatorOptions::exceptionPloidy() const
{
    BOOST_ASSERT(ploidyChromosome.size() == ploidyLevel.size());

    std::map<std::string,unsigned int> result;
    for (unsigned int i=0; i < ploidyChromosome.size(); ++i)
    {
        result.insert(std::make_pair( ploidyChromosome[i], ploidyLevel[i] ));
    }
    return result;
}

void GenomeMutatorOptions::postProcess(bpo::variables_map &vm)
{
    using eagle::common::InvalidOptionException;
    using boost::format;
    eagle::common::OptionsHelper check(vm);
    // carry this value downstream
    force = static_cast<bool>(vm.count("force"));

    std::string reference = check.mutuallyExclusiveOptions(
                                  boost::assign::list_of("reference-genome")("whole-genome")
                            );
    mode = (reference == "whole-genome" ? WHOLE_DIR : SAFE_MODE);

    check.addPathOptions(referenceGenome,"reference-genome");
    check.addPathOptions(wholeGenome,"whole-genome");
    check.addPathOptions(variantList,"variant-list");
    check.inputPathsExist();

    if (!onlyPrintOutputContigNames)
    {
        check.clearPathOptions();
        check.addPathOptions(sampleGenome, "sample-genome");
        check.outputDirsWriteable();

        check.clearPathOptions();
        check.addPathOptions(annotatedVariantList, "annotated-variant-list");
        check.outputFilesWriteable();
    }

    // check that ploidy is sound
    check.inRange<unsigned int>(std::make_pair(organismPloidy,"organism-ploidy"),1);  // 1 <= ploidy < inf
    if(ploidyChromosome.size() != ploidyLevel.size())
    {
        const std::string message("\n   *** The number of occurrences of 'ploidy-chromosome' does not match its 'ploidy-level' counterpart ***\n");
        BOOST_THROW_EXCEPTION(InvalidOptionException(message));
    }

}

} //namespace main
} // namespace eagle
