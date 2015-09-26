################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## description lane-by-lane workflow.
##
## author Mauricio Varea
##
################################################################################

ifeq (,$(DefaultList))
$(error empty list of lanes in 'default' workflow)
endif

lane := $(word 1, $(DefaultList))
fmtLane := $(shell printf "L%03i" $(lane))
DefaultList := $(wordlist 2, $(words $(DefaultList)), $(DefaultList))

assign.keys:=$(TILES:%=tile.%)
assign.values:=$(shell echo {1..$(words $(TILES))})
include $(MAKEFILES_DIR)/assign.mk

# Implicit rule that contains action:
$(EAGLE_OUTDIR)/.$(fmtLane)_%.bcl.completed:
	$(TIME) $(SIMULATE_SEQUENCER) $(EAGLE_FORCE) --generate-bcl-tile \
	        --run-info=$< \
	        --sample-genome-dir="$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)" \
	        $(QUALITY_TABLE:%=--quality-table=%) \
	        $(QQ_TABLE_OPTION) \
	        $(MISMATCH_TABLE_OPTION) \
	        $(HOMOPOLYMER_INDEL_TABLE_OPTION) \
	        $(MOTIF_QUALITY_DROP_TABLE_OPTION) \
	        $(ERROR_MODEL_OPTIONS:%=--error-model-options=%) \
	        --fragments-dir="$(dir $(word 2,$^))" \
	        --output-dir="$(dir $<)" \
	        --lane-count=$(words $(LANES)) \
	        --tiles-per-lane=$(words $(TILES)) \
	        --lane=$(@:$(EAGLE_OUTDIR)/.L00%_$(*).bcl.completed=%) \
	        --tile-num=$(tile.$*) --tile-id=$* \
	        $(RANDOM_SEED_OPTION) \
	        $(SEQUENCER_SIMULATOR_OPTIONS) \
	$(AND) $(TOUCH) $@

$(EAGLE_OUTDIR)/sge/$(fmtLane)_%.bcl.completed: 
	$(MAKE) -n $(EAGLE_OUTDIR)/.L00$(@:$(EAGLE_OUTDIR)/sge/L00%_$(*).bcl.completed=%)_$*.bcl.completed | sed -e 's/make/#make/g' \
	         > $(EAGLE_OUTDIR)/sge/$(@:$(EAGLE_OUTDIR)/sge/%_$(*).bcl.completed=%)_$*.script \
	$(AND) sync \
	$(AND) sync \
	$(AND) qsub -sync y -cwd -v PATH $(EAGLE_OUTDIR)/sge/$(@:$(EAGLE_OUTDIR)/sge/%_$(*).bcl.completed=%)_$*.script

# Explicit rule for all tiles in a lane:
$(foreach t,$(TILES), $(EAGLE_OUTDIR)/.$(fmtLane)_$(t).bcl.completed): $(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml \
                                                                       $(EAGLE_OUTDIR)/fragments/fragments.done

$(foreach t,$(TILES), $(EAGLE_OUTDIR)/sge/$(fmtLane)_$(t).bcl.completed): $(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml \
                                                                          $(EAGLE_OUTDIR)/fragments/fragments.done \
                                                                          $(EAGLE_OUTDIR)/sge/.sentinel

