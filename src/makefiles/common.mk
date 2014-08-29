################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## description Main workflow.
##
## author Mauricio Varea
##
################################################################################

SHELL := /bin/bash
EAGLE_BINDIR := @EAGLE_FULL_BINDIR@
EAGLE_LIBEXECDIR := @EAGLE_FULL_LIBEXECDIR@

NULL    :=
SPACE   :=$(NULL) $(NULL)
COMMA   :=,
PIPE    :=|
BACKSLASHED_PIPE :=\|

AND      = &&
OR       = ||
MKDIR    = mkdir -p
SLEEP    = sleep
TOUCH    = touch
WIPE_OUT = rm -rf
TIME     = /usr/bin/time -v -o $(EAGLE_OUTDIR)/time.log --append

APPLY_VARIANTS     := $(EAGLE_BINDIR)/applyVariants
ALLOCATE_FRAGMENTS := $(EAGLE_BINDIR)/allocateFragments
CREATE_RUN_FOLDER  := $(EAGLE_BINDIR)/createRunFolder
SIMULATE_SEQUENCER := $(EAGLE_LIBEXECDIR)/simulateSequencer
CANONICAL2SEGMENTS := $(EAGLE_LIBEXECDIR)/canonical2segments

#Index = $(shell echo "$(2)" | sed -r -e "s/[ \t]+/\n/g" | grep -n $(1) | cut -d ':' -f 1)

.PRECIOUS: %/.sentinel
%/.sentinel:
	$(MKDIR) $* $(OR) ( $(SLEEP) 5 $(AND) $(MKDIR) $* ); $(TOUCH) $@
