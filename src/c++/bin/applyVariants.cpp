/**
 ** Copyright (c) 2014 Illumina, Inc.
 **
 ** This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
 ** covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
 **
 ** \description The main for the reference mutator.
 **
 ** \author Mauricio Varea
 **/

#include "main/GenomeMutatorOptions.hh"
#include "main/GenomeMutator.hh"
#include "model/Genotype.hh"

static void genomeMutatorLauncher(const eagle::main::GenomeMutatorOptions &options)
{
    if (eagle::main::GenomeMutatorOptions::WHOLE_DIR == options.mode)
    {
        std::clog << (boost::format("Looking for a reference genome in %s ...") % options.wholeGenome ).str() << std::endl;
        eagle::main::GenomeMutator genomeMutator(
            options.wholeGenome,
            options.variantList,
            options.sampleGenome,
            options.annotatedVariantList,
            eagle::model::Ploidy( options.organismPloidy, options.exceptionPloidy() ),
            options.prefixToAdd,
            options.force,
            options
            );
        genomeMutator.run();
    } else {
        eagle::main::GenomeMutator genomeMutator(
            options.referenceGenome,
            options.variantList,
            options.sampleGenome,
            options.annotatedVariantList,
            eagle::model::Ploidy( options.organismPloidy, options.exceptionPloidy() ),
            options.prefixToAdd,
            options.force,
            options
            );
        genomeMutator.run();
    }
}

int main(int argc, char *argv[])
{
    eagle::common::run(genomeMutatorLauncher, argc, argv);
}


