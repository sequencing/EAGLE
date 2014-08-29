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


include $(VALIDATION_DIR)/test3/config.mk

# Test4 = test3 with the following overrides
RUN_INFO_XML         = $(EAGLE_DATADIR)/RunInfo/RunInfo_PairedReadsBarcode1x1Tiles.xml
