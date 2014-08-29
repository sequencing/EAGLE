################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## description Configuration makefile to run the ORION tests (without noise).
##
## author Mauricio Varea
##
################################################################################

include $(VALIDATION_DIR)/testHuman1/config.mk

# TestHuman2 = testHuman1 with the following overrides
EAGLE_TEST_REF       = /illumina/scratch/iGenomes/Homo_sapiens/NCBI/build37.1/Sequence/Chromosomes
VARIANTS_VCF         = $(EAGLE_DATADIR)/Variants/Stephens2009.vcf

ALIGNMENT_EXTRAS      = --PAS_PARAMS "--max_patterns_to_store 100000"

# female. Override for male.
EAGLE_EXTRAS           = --genome-mutator-options='--organism-ploidy=2 --ploidy-chromosome=chrY --ploidy-level=0 --ploidy-chromosome=chrM --ploidy-level=100'
