################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## description Configuration makefile to configure validation tests.
##
## author Lilian Janin
##
################################################################################

include $(VALIDATION_DIR)/testHuman1/config.mk

# TestHuman3 = testHuman1 with the following overrides
EAGLE_TEST_REF       = /illumina/scratch/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes/chr21.fa
