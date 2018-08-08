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
ifneq (,$(TEMPLATE_LENGTH_TABLE))
TEMPLATE_LENGTH_TABLE_OPTION = --template-length-table=$(TEMPLATE_LENGTH_TABLE)
endif
ifneq (,$(MISMATCH_TABLE))
MISMATCH_TABLE_OPTION = --mismatch-table=$(MISMATCH_TABLE)
endif
ifneq (,$(HOMOPOLYMER_INDEL_TABLE))
HOMOPOLYMER_INDEL_TABLE_OPTION = --homopolymer-indel-table=$(HOMOPOLYMER_INDEL_TABLE)
endif
ifneq (,$(MOTIF_QUALITY_DROP_TABLE))
MOTIF_QUALITY_DROP_TABLE_OPTION = --motif-quality-drop-table=$(MOTIF_QUALITY_DROP_TABLE)
endif
ifneq (,$(GC_COVERAGE_FIT_TABLE))
GC_COVERAGE_FIT_TABLE_OPTION = --gc-coverage-fit-table=$(GC_COVERAGE_FIT_TABLE)
endif

CHROMOSOME_ALLELES_FORWARD_AND_REVERSE = $(CHROMOSOME_ALLELES) $(CHROMOSOME_ALLELES:%=%_rev) 


.PHONY: all
all: $(foreach lane,$(LANES),$(foreach tile,$(TILES), \
     $(EAGLE_OUTDIR)/.$(shell printf "L%03i" $(lane))_$(tile).bcl.completed ))

TileCount := $(shell dc -e "$(words $(TILES)) $(words $(LANES)) *p")

sge: $(foreach lane,$(LANES),$(foreach tile,$(TILES), \
     $(EAGLE_OUTDIR)/sge/$(shell printf "L%03i" $(lane))_$(tile).bcl.completed ))

.PHONY: print-output-contig-names
print-output-contig-names: $(REFERENCE_GENOME)
	( $(TIME) $(APPLY_VARIANTS) --only-print-output-contig-names $(EAGLE_FORCE) \
	  $(^:%=--$(REFERENCE_MODE)=%) \
	  $(VARIANT_LIST:%=--variant-list=%) \
	  --sample-genome=$(dir $@) \
	  --annotated-variant-list=$(dir $@)/canonical.vcf \
	  $(GENOME_MUTATOR_OPTIONS) \
	) 2> /dev/null | grep -P "^Chromosome allele: "

.PRECIOUS: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(APPLY_VARIANTS)).log
.PHONY: sample $(notdir $(APPLY_VARIANTS))
sample: $(notdir $(APPLY_VARIANTS))
$(notdir $(APPLY_VARIANTS)): $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml
$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml: $(REFERENCE_GENOME)
	( $(TIME) $(APPLY_VARIANTS) $(EAGLE_FORCE) \
	  $(^:%=--$(REFERENCE_MODE)=%) \
	  $(VARIANT_LIST:%=--variant-list=%) \
	  --sample-genome=$(dir $@) \
	  --annotated-variant-list=$(dir $@)/canonical.vcf \
	  $(GENOME_MUTATOR_OPTIONS) \
	) &> $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(APPLY_VARIANTS)).log

.PRECIOUS: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(ALLOCATE_FRAGMENTS)).log
.PHONY: fragments $(notdir $(ALLOCATE_FRAGMENTS))
fragments: $(notdir $(ALLOCATE_FRAGMENTS))
$(notdir $(ALLOCATE_FRAGMENTS)): $(EAGLE_OUTDIR)/fragments/fragments.done
$(EAGLE_OUTDIR)/fragments_single_process.pos: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml
	( $(TIME) $(ALLOCATE_FRAGMENTS) $(EAGLE_FORCE) \
	  --sample-genome-dir="$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)" \
	  --output-dir="$(dir $@)" \
	  $(COVERAGE_DEPTH_OPTION) \
	  $(RANDOM_SEED_OPTION) \
	  $(TEMPLATE_LENGTH_TABLE_OPTION) \
	  --bases-per-cluster=$(BASES_PER_CLUSTER) \
	  --tiles=$(TileCount) \
	  $(GC_COVERAGE_FIT_TABLE_OPTION) \
	  $(FRAGMENTS_ALLOCATOR_OPTIONS) \
	) &> $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(ALLOCATE_FRAGMENTS)).log

$(EAGLE_OUTDIR)/fragments/fragments.done: $(foreach chrAllele, $(CHROMOSOME_ALLELES_FORWARD_AND_REVERSE), $(subst $(PIPE),$(BACKSLASHED_PIPE),$(EAGLE_OUTDIR)/fragments/fragments_$(chrAllele)/fragments.done))
	( $(TIME) $(ALLOCATE_FRAGMENTS) $(EAGLE_FORCE) \
	  --merge-existing-fragments \
	  --sample-genome-dir="$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)" \
	  --output-dir="$(dir $@)" \
	  $(COVERAGE_DEPTH_OPTION) \
	  $(RANDOM_SEED_OPTION) \
	  $(TEMPLATE_LENGTH_TABLE_OPTION) \
	  --bases-per-cluster=$(BASES_PER_CLUSTER) \
	  --tiles=$(TileCount) \
	  $(GC_COVERAGE_FIT_TABLE_OPTION) \
	  $(FRAGMENTS_ALLOCATOR_OPTIONS) \
	) 2>&1 \
	$(AND) \
	  touch $@ \
	 >> $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(ALLOCATE_FRAGMENTS)).log

$(EAGLE_OUTDIR)/fragments/fragments_%/fragments.done: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml
	mkdir -p "$(EAGLE_OUTDIR)/fragments/fragments_$*" \
	$(AND) \
	( $(TIME) $(ALLOCATE_FRAGMENTS) $(EAGLE_FORCE) \
	  --contig="$*" \
	  --sample-genome-dir="$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/${*:_rev=}.fa" \
	  --output-dir="$(dir $@)" \
	  $(COVERAGE_DEPTH_OPTION) \
	  $(RANDOM_SEED_OPTION) \
	  $(TEMPLATE_LENGTH_TABLE_OPTION) \
	  --bases-per-cluster=$(BASES_PER_CLUSTER) \
	  --tiles=$(TileCount) \
	  $(GC_COVERAGE_FIT_TABLE_OPTION) \
	  $(FRAGMENTS_ALLOCATOR_OPTIONS) \
	) 2>&1 \
	$(AND) \
	  touch "$@" \
	 >> $(EAGLE_OUTDIR)/$(SAMPLE_GENOME).$(notdir $(ALLOCATE_FRAGMENTS)).log

$(notdir $(CREATE_RUN_FOLDER)): $(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml
$(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml: $(EAGLE_OUTDIR)/$(RUN_FOLDER).xml
	$(TIME) $(CREATE_RUN_FOLDER) $(EAGLE_FORCE) $(dir $@) --run-info=$< --tile-id="$(subst $(SPACE),$(COMMA),$(TILES))" $(RUN_FOLDER_GENERATOR_OPTIONS)

DefaultList:=$(LANES)
include $(foreach ln, $(LANES), $(MAKEFILES_DIR)/defaultLane.mk)

$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/segmentsFromRef.tsv: $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/genome_size.xml
	$(TIME) $(CANONICAL2SEGMENTS) -i $(EAGLE_OUTDIR)/$(SAMPLE_GENOME)/canonical.vcf -o $(EAGLE_OUTDIR)/$(SAMPLE_GENOME) &> $(EAGLE_OUTDIR)/canonical2segments.log

.PHONY: bam bams
bam: eagle.bam.bai

eagle.bam.bai: eagle.bam
	$(TIME) $(SAMTOOLS) index eagle.bam

eagle.bam: $(foreach bamChr, $(BAM_CHROMOSOMES), $(subst $(PIPE),$(BACKSLASHED_PIPE),eagle_$(bamChr).bam))
ifeq (1,$(words $(BAM_CHROMOSOMES)))
	cp $^ $@
else
	$(SAMTOOLS) view -H "eagle_$(word 1,${BAM_CHROMOSOMES}).bam" | tr -d '\0' | grep --text SN | sed 's/.*\tSN:\(.*\)\t.*/eagle_\1.bam/' | xargs ls -f 2> /dev/null | xargs $(TIME) $(SAMTOOLS) cat -o eagle.bam
endif
#	$(TIME) $(SAMTOOLS) cat -o eagle.bam $^

eagle_%.bam: $(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml $(EAGLE_OUTDIR)/fragments/fragments.done $(EAGLE_OUTDIR)/sample_genome/segmentsFromRef.tsv
	$(TIME) $(SIMULATE_SEQUENCER) $(EAGLE_FORCE) --generate-bam \
	        --run-info=$< \
	        --sample-genome-dir="$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)" \
	        $(QUALITY_TABLE:%=--quality-table=%) \
	        $(QQ_TABLE_OPTION) \
	        $(MISMATCH_TABLE_OPTION) \
	        $(HOMOPOLYMER_INDEL_TABLE_OPTION) \
	        $(MOTIF_QUALITY_DROP_TABLE_OPTION) \
	        $(ERROR_MODEL_OPTIONS:%=--error-model-options=%) \
	        --fragments-dir="$(EAGLE_OUTDIR)/fragments" \
	        --output-dir="$(EAGLE_OUTDIR)" \
	        --output-filename="eagle_$*.bam" \
	        --lane-count=$(words $(LANES)) \
	        --tiles-per-lane=$(words $(TILES)) \
	        $(RANDOM_SEED_OPTION) \
	        $(SEQUENCER_SIMULATOR_OPTIONS) \
			--bam-region="$*"

.PHONY: sample-bam
sample-bam: eagle.sample.bam
eagle.sample.bam: $(EAGLE_OUTDIR)/$(RUN_FOLDER)/RunInfo.xml $(EAGLE_OUTDIR)/fragments/fragments.done
	$(TIME) $(SIMULATE_SEQUENCER) $(EAGLE_FORCE) --generate-sample-bam \
	        --run-info=$< \
	        --sample-genome-dir="$(EAGLE_OUTDIR)/$(SAMPLE_GENOME)" \
	        $(QUALITY_TABLE:%=--quality-table=%) \
	        $(QQ_TABLE_OPTION) \
	        $(MISMATCH_TABLE_OPTION) \
	        $(HOMOPOLYMER_INDEL_TABLE_OPTION) \
	        $(MOTIF_QUALITY_DROP_TABLE_OPTION) \
	        $(ERROR_MODEL_OPTIONS:%=--error-model-options=%) \
	        --fragments-dir="$(EAGLE_OUTDIR)/fragments" \
	        --output-dir="$(EAGLE_OUTDIR)" \
	        --output-filename=eagle.sample.bam \
	        --lane-count=$(words $(LANES)) \
	        --tiles-per-lane=$(words $(TILES)) \
	        $(RANDOM_SEED_OPTION) \
	        $(SEQUENCER_SIMULATOR_OPTIONS)


.PHONY: fastq
fastq: $(foreach lane,$(LANES),$(foreach tile,$(TILES), \
       $(EAGLE_OUTDIR)/.$(shell printf "L%03i" $(lane))_$(tile).fastq.completed ))
