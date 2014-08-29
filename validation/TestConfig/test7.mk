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


include $(VALIDATION_DIR)/test5/config.mk

# Test7 = test5 with the following overrides
VARIANTS_VCF         = $(EAGLE_DATADIR)/Variants/OneAtTelomere.vcf
