################################################################################
##
## Copyright (c) 2014 Illumina, Inc.
##
## This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
## covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
##
## description Tumour-normal cancer workflow.
##
## author Lilian Janin
##
################################################################################

CHROMOSOMES   ?= chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr1 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
VCF_GENERATOR  = /home/ljanin/workspace/git/EAGLE/src/analysis/perl/vcfGenerator.pl
FASTA_PATH     = /illumina/scratch/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/Chromosomes
REPEAT_MASK    = /illumina/development/EAGLE/Resources/UCSC_RepeatMasker_all.bed
RANDOM_SEED   ?= 1
DBSNP_VCF     ?= /illumina/scratch/eagle/RealisticInputs/dbSnp_set2_UCSC.vcf

TRANSLOCATIONS_GENERATOR = /home/ljanin/workspace/git/EAGLE/src/analysis/cpp/add_translocations/AddTranslocations
GAP_FILE  = /home/ljanin/workspace/git/EAGLE/src/analysis/cpp/add_translocations/HumanNCBI37_UCSCgaps.bed


all: variants
	@echo -e "VCF files generated.\n  To simulate, try:\n    configureEAGLE_Normal+Tumour.pl --shared-variants=small_variants_normal.vcf --shared-variants=translocations_normal.vcf --tumour-only-variants=translocations_tumour_only.vcf --tumour-only-variants=small_variants_tumour_only.vcf --normal-coverage=40 --tumour-overall-coverage=80 --tumour-purity=0.7 --reference-genome=genome.fa --strelka-reference-genome=genome.fa --isaac-sorted-reference=/illumina/development/iSAAC/iGenomes/Homo_sapiens/UCSC/hg19/Sequence/iSAACIndex.20121213/sorted-reference.xml\n"


sim: EAGLE/.done

variants: variants_normal variants_tumour_only

variants_normal: translocations_normal.vcf small_variants_normal.vcf

variants_tumour_only: translocations_tumour_only.vcf small_variants_tumour_only.vcf


genome.fa: $(foreach chr,$(CHROMOSOMES),$(FASTA_PATH)/$(chr).fa)
	cat $^ > $@

empty.vcf:
	echo "" > empty.vcf

small_variants_normal.vcf: $(foreach chr,$(CHROMOSOMES),small_variants_normal_$(chr)/combined.vcf)
	cat $^ | awk 'BEGIN { OFS=FS="\t" } { $$3="normal_"$$3 ; print}' > $@

small_variants_normal_%/.sentinel:
	mkdir $(dir $@) \
	&& touch $@

.PRECIOUS: small_variants_normal_%/.done small_variants_normal_%/fromDbSnp.vcf small_variants_normal_%/RepeatMasker.bed
small_variants_normal_%/RepeatMasker.bed: small_variants_normal_%/.sentinel
	grep -P "^$*\t" $(REPEAT_MASK) > $@

#      -invFrq  0.000005
#      -delLog  \
#      -insLog  \

small_variants_normal_%/.done: small_variants_normal_%/.sentinel small_variants_normal_%/RepeatMasker.bed
	( cd $(dir $@) \
	; \
	$(VCF_GENERATOR) -v \
      -outDir  . \
      -nploidy 2 \
      -fasFile $(FASTA_PATH)/$*.fa \
      -repFile $(notdir $(word 2, $^)) \
      -chrName $* \
      -homFrq  0.3333 \
      -snpFrq  0.0002 \
      -delFrq  0.00005 \
      -delMin  1 \
      -delMax  100 \
      -insFrq  0.00005 \
      -insMin  1 \
      -insMax  100 \
      -invFrq  0 \
      -invMin  500 \
      -invMax  1000 \
      -cnvFrq  0.000005 \
      -cnvMin  2000 \
      -cnvMax  20000 \
      -log-width 0.1 \
      -gapFile /illumina/development/EAGLE/Resources/UCSC_Gaps/$*_gap.txt \
	  -random-seed $(RANDOM_SEED) \
	) && \
	touch $@

small_variants_normal_%/fromDbSnp.vcf: $(DBSNP_VCF) small_variants_normal_%/.sentinel
	grep -P "^$*\t" $< | awk '{if (NR % 300 != 0) { print $0 } }' > $@

small_variants_normal_%/combined.vcf: small_variants_normal_%/.done small_variants_normal_%/fromDbSnp.vcf
	( cat $(dir $(word 1,$^))/$*-A.mutMap.vcf \
	; grep '1/0' $(dir $<)/$*-B.mutMap.vcf | sed 's:1/0:0/1:' \
	; cat $(word 2,$^) \
	) | sort -gk 2 > $@

translocations_normal_stage1.vcf: genome.fa empty.vcf
	$(TRANSLOCATIONS_GENERATOR) -g $(GAP_FILE) -r $^ -vv -c empty.vcf -t 0 --random-seed $(RANDOM_SEED)1 -o $@

translocations_tumour_only_stage1.vcf: genome.fa empty.vcf
	$(TRANSLOCATIONS_GENERATOR) -g $(GAP_FILE) -r $^ -vv -c empty.vcf -t 0 --random-seed $(RANDOM_SEED)2 -o $@

#  -t 110    -preFile /illumina/scratch/optimusprime/Repeat_Simulation/input_transformGenome/$*-preDefined-repeats.txt 

translocations_normal.vcf: translocations_normal_stage1.vcf
	awk 'BEGIN { OFS=FS="\t" } { $$3="normal_"$$3 ; print}' $^ > $@

translocations_tumour_only.vcf: translocations_tumour_only_stage1.vcf
	awk 'BEGIN { OFS=FS="\t" } { $$3="tumour_"$$3 ; print}' $^ > $@




small_variants_tumour_only.vcf: $(foreach chr,$(CHROMOSOMES),small_variants_tumour_only_$(chr)/combined.vcf)
	cat $^ | awk 'BEGIN { OFS=FS="\t" } { $$3="tumour_"$$3 ; print}' > $@

small_variants_tumour_only_%/.sentinel:
	mkdir $(dir $@) \
	&& touch $@

.PRECIOUS: small_variants_tumour_only_%/.done small_variants_tumour_only_%/fromDbSnp.vcf small_variants_tumour_only_%/RepeatMasker.bed
small_variants_tumour_only_%/RepeatMasker.bed: small_variants_tumour_only_%/.sentinel
	grep -P "^$*\t" $(REPEAT_MASK) > $@

#      -delLog  \
#      -insLog  \

small_variants_tumour_only_%/.done: small_variants_tumour_only_%/.sentinel small_variants_tumour_only_%/RepeatMasker.bed
	( cd $(dir $@) \
	; \
	$(VCF_GENERATOR) -v \
      -outDir  . \
      -nploidy 2 \
      -fasFile $(FASTA_PATH)/$*.fa \
      -repFile $(notdir $(word 2, $^)) \
      -chrName $* \
      -homFrq  0.3333 \
      -snpFrq  0.000001 \
      -delFrq  0.000002 \
      -delMin  1 \
      -delMax  100000 \
      -insFrq  0.000002 \
      -insMin  1 \
      -insMax  1000 \
      -invFrq  0.00000002 \
      -invMin  10 \
      -invMax  10000 \
      -cnvFrq  0.000000005 \
      -cnvMin  10000 \
      -cnvMax  1000000 \
      -log-width 0.1 \
      -gapFile /illumina/development/EAGLE/Resources/UCSC_Gaps/$*_gap.txt \
	  -random-seed $(RANDOM_SEED)3 \
	) && \
	touch $@

small_variants_tumour_only_%/fromDbSnp.vcf: $(DBSNP_VCF) small_variants_tumour_only_%/.sentinel
	grep -P "^$*\t" $< | awk '{if (NR % 300 == 0) { print $0 } }' > $@

small_variants_tumour_only_%/combined.vcf: small_variants_tumour_only_%/.done small_variants_tumour_only_%/fromDbSnp.vcf
	( cat $(dir $(word 1,$^))/$*-A.mutMap.vcf \
	; grep '1/0' $(dir $<)/$*-B.mutMap.vcf | sed 's:1/0:0/1:' \
	; cat $(word 2,$^) \
	) | sort -gk 2 > $@
