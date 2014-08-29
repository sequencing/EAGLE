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


include $(VALIDATION_DIR)/test1/config.mk

# Test3 = test1 with the following overrides
RUN_INFO_XML         = $(EAGLE_DATADIR)/RunInfo/RunInfo_PairedReads1x1Tiles.xml

CASAVA_ANALYSIS      = eland_pair
CASAVA_USE_BASES     = y*n,y*n
