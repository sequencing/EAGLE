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

include $(MAKEFILES_DIR)/common.mk

ifneq (,$(COVERAGE_DEPTH))
COVERAGE_DEPTH_OPTION = --coverage-depth=$(COVERAGE_DEPTH)
endif
ifneq (,$(RANDOM_SEED))
RANDOM_SEED_OPTION = --random-seed=$(RANDOM_SEED)
endif
ifneq (,$(QQ_TABLE))
QQ_TABLE_OPTION = --qq-table=$(QQ_TABLE)
endif
ifneq (,$(MISMATCH_TABLE))
MISMATCH_TABLE_OPTION = --mismatch-table=$(MISMATCH_TABLE)
endif
ifneq (,$(HOMOPOLYMER_INDEL_TABLE))
HOMOPOLYMER_INDEL_TABLE_OPTION = --homopolymer-indel-table=$(HOMOPOLYMER_INDEL_TABLE)
endif

.PHONY: all
all: $(foreach lane,$(LANES),$(foreach tile,$(TILES), \
     $(EAGLE_OUTDIR)/.$(shell printf "L%03i" $(lane))_$(tile).bcl.completed ))

TileCount := $(shell dc -e "$(words $(TILES)) $(words $(LANES)) *p")

sge: $(foreach lane,$(LANES),$(foreach tile,$(TILES), \
     $(EAGLE_OUTDIR)/sge/$(shell printf "L%03i" $(lane))_$(tile).bcl.completed ))

.PRECIOUS: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(APPLY_VARIANTS)).log
.PHONY: sample $(notdir $(APPLY_VARIANTS))
structure: $(notdir $(APPLY_VARIANTS))
$(notdir $(APPLY_VARIANTS)): $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml
$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml: $(REFERENCE_GENOME)
	( $(TIME) $(APPLY_VARIANTS) $(EAGLE_FORCE) \
	  $(^:%=--$(REFERENCE_MODE)=%) \
	  $(VARIANT_LIST:%=--variant-list=%) \
	  --sample-genome=$(dir $@) \
	  --annotated-variant-list=$(dir $@)/structural.vcf \
	  $(GENOME_MUTATOR_OPTIONS) \
	  $(AND) touch $@ \
	) &> $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(APPLY_VARIANTS)).log

$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/%.fa: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/%.struct
	cp $< $@

tmpRule:

.PHONY: sample
sample: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/sample_complete.txt
$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/sample_complete.txt: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml
	echo here1 $(EAGLE_OUTDIR)/$(SAMPLE_GENOME) $(MAKE) tmpRule $(patsubst %.struct,%.fa,$(wildcard $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/*.struct)) $(eval $(patsubst %.struct,%.fa,$(wildcard $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/*.struct))) \
	$(AND) echo here2 $(shell ls $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/*.struct) \
	$(AND) echo here3 $(patsubst %.struct,%.fa,$(shell ls $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/*.struct)) \
	$(AND) $(MAKE) tmpRule $(patsubst %.struct,%.fa,$(shell ls $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/*.struct)) \
	$(AND) touch $@

.PRECIOUS: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(ALLOCATE_FRAGMENTS)).log
.PHONY: fragments $(notdir $(ALLOCATE_FRAGMENTS))
fragments: $(notdir $(ALLOCATE_FRAGMENTS))
$(notdir $(ALLOCATE_FRAGMENTS)): $(EAGLE_OUTDIR)/fragments.pos
$(EAGLE_OUTDIR)/fragments.pos: sample
	( $(TIME) $(ALLOCATE_FRAGMENTS) $(EAGLE_FORCE) \
	  --sample-genome-dir=$(EAGLE_OUTDIR)/$(SAMPLE_GENOME) \
	  --output-dir=$(dir $@) \
	  $(COVERAGE_DEPTH_OPTION) \
	  $(RANDOM_SEED_OPTION) \
	  --bases-per-cluster=$(BASES_PER_CLUSTER) \
	  --tiles=$(TileCount) \
	  $(FRAGMENTS_ALLOCATOR_OPTIONS) \
	) &> $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(ALLOCATE_FRAGMENTS)).log

$(notdir $(CREATE_RUN_FOLDER)): $(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml
$(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml: $(EAGLE_OUTDIR)/$(RUN_FOLDER).xml
	$(TIME) $(CREATE_RUN_FOLDER) $(EAGLE_FORCE) $(dir $@) --run-info=$< --tile-id="$(subst $(SPACE),$(COMMA),$(TILES))" $(RUN_FOLDER_GENERATOR_OPTIONS)

DefaultList:=$(LANES)
include $(foreach ln, $(LANES), $(MAKEFILES_DIR)/defaultLane.mk)

.PHONY: bam
bam: test.bam
test.bam: $(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml $(EAGLE_OUTDIR)/fragments.pos
	$(TIME) $(SIMULATE_SEQUENCER) $(EAGLE_FORCE) --generate-bam \
	        --run-info=$< \
	        --sample-genome-dir=$(EAGLE_OUTDIR)/$(SAMPLE_GENOME) \
	        $(QUALITY_TABLE:%=--quality-table=%) \
	        $(QQ_TABLE_OPTION) \
	        $(MISMATCH_TABLE_OPTION) \
	        $(HOMOPOLYMER_INDEL_TABLE_OPTION) \
	        --insert-positions=$(word 2,$^) \
	        --output-dir=$(dir $<) \
	        --lane-count=$(words $(LANES)) \
	        --tiles-per-lane=$(words $(TILES)) \
	        $(RANDOM_SEED_OPTION) \
	        $(SEQUENCER_SIMULATOR_OPTIONS)
