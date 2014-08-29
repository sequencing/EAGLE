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


include $(VALIDATION_DIR)/test8/config.mk

# Test9 = test8 with the following overrides
RUN_INFO_XML         = $(EAGLE_DATADIR)/RunInfo/RunInfo_SingleRead1x1Tiles.xml

CASAVA_ANALYSIS      = eland_extended
CASAVA_USE_BASES     = y*n
