/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description Top level component to induce variants in a reference.
 **
 ** \author Mauricio Varea
 **/

#ifndef EAGLE_MAIN_GENOME_MUTATOR_HH
#define EAGLE_MAIN_GENOME_MUTATOR_HH

#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>

#include "genome/Reference.hh"
#include "genome/VariantList.hh"
#include "model/Genotype.hh"
#ifdef DISTRIBUTED_GENOME_MUTATOR
#include "genome/SharedFastaReference.hh"
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR

#include "GenomeMutatorOptions.hh"


namespace eagle
{
namespace main
{

class GenomeMutator
{
public:
    GenomeMutator(  // SAFE_MODE
        const std::vector<boost::filesystem::path> inputReferences,
        const std::vector<boost::filesystem::path> variantFiles,
        const boost::filesystem::path outputReference,
        const boost::filesystem::path outputVariants,
        const eagle::model::Ploidy ploidy,
        const std::string prefixToAdd,
        const bool overwrite,
        const GenomeMutatorOptions &options)
    : genome_( inputReferences )
    , sample_( outputReference, overwrite )
    , variantList_( variantFiles, outputVariants, ploidy, overwrite )
    , prefixToAdd_(prefixToAdd)
    , outputReference_( outputReference )
    , options_( options )
    {
#ifdef DISTRIBUTED_GENOME_MUTATOR
        genome::SharedFastaReference::init( inputReferences );
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR
    }
    GenomeMutator(  // WHOLE_DIR
        const boost::filesystem::path inputReference,
        const std::vector<boost::filesystem::path> variantFiles,
        const boost::filesystem::path outputReference,
        const boost::filesystem::path outputVariants,
        const eagle::model::Ploidy ploidy,
        const std::string prefixToAdd,
        const bool overwrite,
        const GenomeMutatorOptions &options)
    : genome_( inputReference )
    , sample_( outputReference, overwrite )
    , variantList_( variantFiles, outputVariants, ploidy, overwrite )
    , prefixToAdd_(prefixToAdd)
    , outputReference_( outputReference )
    , options_( options )
    {
#ifdef DISTRIBUTED_GENOME_MUTATOR
        genome::SharedFastaReference::init( inputReference );
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR
    }
    void run();
private:
    unsigned int process(unsigned int i,
                         eagle::genome::EventIterator startPosition,
                         eagle::model::Contig& contigOut,
                         eagle::model::Direction direction,
                         eagle::model::Locus& finalPosition
#ifdef DISTRIBUTED_GENOME_MUTATOR
                         , std::vector<bool>& appliedStructuralVariants
#endif //ifdef DISTRIBUTED_GENOME_MUTATOR
                         );
    eagle::genome::MultiFastaReference genome_;
    eagle::genome::MultiFastaReference sample_;
    eagle::genome::VariantList         variantList_;
    const std::string prefixToAdd_;
    const boost::filesystem::path outputReference_;
    const GenomeMutatorOptions &options_;
};

} // namespace main
} // namespace eagle

#endif // EAGLE_MAIN_GENOME_MUTATOR_HH
