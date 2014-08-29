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


EAGLE_TEST_REF       = /illumina/scratch/iGenomes/Enterobacteriophage_lamdba/NCBI/1993-04-28/Sequence/Chromosomes
VARIANTS_VCF         = $(EAGLE_DATADIR)/Variants/None.vcf
RUN_INFO_XML         = $(EAGLE_DATADIR)/RunInfo/RunInfo_SingleRead1x1Tiles.xml
DEPTH                = 30

CASAVA_ANALYSIS      = eland_extended
CASAVA_USE_BASES     = y*n
