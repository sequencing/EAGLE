#!/usr/bin/env perl

## MODULE: transformGenomeNploidyEX.pl
## AUTHOR: Bret D. Barnes
##
## Copyright (c) 2008 Illumina Cambridge Limited. All rights reserved.
## This source file is covered by the "Solexa Public Source License"
## agreement and bound by the terms therein.

## Adapted from transformGenome_v3.0.pl
## Add SNPs, Indels, CNV's, inversions and translocations to chromosomes and allows any ploidy,
##  up to 26. Variants can be simulated or predefined.
##
##
## Program inputs:
##
## -fas, -f     Chromosoem fasta files
## -pre, -p     Pre-defined variant files
## -rep, -r     Repeat mapping files
##
## -hom, -h     Homozygous frequency
##
## -insFrq      Insertion frequency
## -insMin      Insertion min size
## -insMax      Insertion max size
## -insMdl      Insertion size distribution model [linear, logarithmic], default: linear
## -insLog      Turns insMdl from linear to logorithmic (exponetial-decay)
## -insExt      Additional singleton insert size to create with a single copy guarenteed
##

##
## Workflow:
##
##  I.   Preprocessor
##  II.  Load inputs
##       A. Pre-defined indels
##       B. Fasta file info
##       C. Repeat regions
##  III. Add variant: addVariant();
##       A. Inversions
##       B. CNV's
##       C. Deletions
##       D. Insertions
##       E. SNP's
##  IV.  Transform genome
##  V.   Add translocations
##

use warnings;
use strict;
use Getopt::Long;

use constant FALSE => 0;
use constant TRUE  => 1;

## BED Format
use constant BED_CHR    => 0;
use constant BED_START  => 1;
use constant BED_END    => 2;
use constant BED_NAME   => 3;
use constant BED_SCORE  => 4;
use constant BED_STRAND => 5;

## Pre-defined variants file format defs
use constant PRE_TYPE_IDX  => 0;
use constant PRE_CHR_IDX   => 1;
use constant PRE_START_IDX => 2;
use constant PRE_END_IDX   => 3;
use constant PRE_REP_IDX   => 4;
use constant PRE_SEQ_IDX   => 5;

## Variant data structure format
use constant VAR_TYPE_IDX     => 0;
use constant VAR_CHR_IDX      => 1;
use constant VAR_START_IDX    => 2;
use constant VAR_END_IDX      => 3;
use constant VAR_REP_IDX      => 4;
use constant VAR_ZYGO_IDX     => 5;
use constant VAR_SRC_IDX      => 6;
use constant VAR_ORG_SEQ_IDX  => 7;
use constant VAR_NEW_SEQ_IDX  => 8;
use constant VAR_UP_ONE_IDX   => 9;
use constant VAR_DN_ONE_IDX   => 10;

my $mmapHeader = "NEW_START\tVAR_TYPE\tCHR\tREF_START\tREF_END\tREPEAT_CNT\t" .
    "ZYGOSITY\tSOURCE\tREF_SEQ\tNEW_SEQ";
my $cmapHeader = "ALN_IDX\tREF_POS\tNEW_POS";

##
## Global options defaults:
##
my %options = ();
my %params  = ();

$options{'nploidy'}    = 1;

$options{'homFrq'}     = 0.00;
$options{'snpFrq'}     = 0.00;

$options{'del'}{'frq'}     = 0.00;
$options{'del'}{'min'}     = 1;
$options{'del'}{'max'}     = 1;
$options{'del'}{'mdl'}     = "linear";
#$options{'del'}{'ext'}     = ();

$options{'ins'}{'frq'}     = 0.00;
$options{'ins'}{'min'}     = 1;
$options{'ins'}{'max'}     = 1;
$options{'ins'}{'mdl'}     = "linear";
#$options{'ins'}{'ext'}     = ();

$options{'inv'}{'frq'}     = 0.00;
$options{'inv'}{'min'}     = 1000;
$options{'inv'}{'max'}     = 10000;
$options{'inv'}{'mdl'}     = "linear";
#$options{'inv'}{'ext'}     = ();

$options{'cnv'}{'frq'}     = 0.00;
$options{'cnv'}{'min'}     = 1000;
$options{'cnv'}{'max'}     = 20000;
$options{'cnv'}{'mdl'}     = "linear";
#$options{'cnv'}{'ext'}     = ();
$options{'cnv'}{'copyMax'} = 3;

$options{'bases'}      = 'ACGT';
$options{'log-width'}  = 0.25;
$options{'rep-width'}  = 3;
$options{'rep-res'}    = 100;
$options{'varSep'}     = 1;

my $skipNext;

##
## Define legal size model types
##
$options{'legMdl'}{'linear'} = 1;
$options{'legMdl'}{'log'}    = 1;

##
## Allowed variant types
##
my @allowedTypes = qw(trn inv cnv ins del snp);
foreach my $type (@allowedTypes) {
    $options{'allowed-types'}{$type} = 1;
}

## Long conversion table
my %typeMap = ("ins" => "nucleotide_insertion",
               "del" => "nucleotide_deletion",
               "inv" => "inversion",
               "cnv" => "duplication",
               "snp" => "snv",
               "trn" => "translocation",
               );

##
## VCF Score Default
##
$options{'vcf-score'}        = 100;
$options{'genotype-quality'} = 100;
$options{'read-depth'}       = 30;
$options{'sample-name'}      = "Sample_Sim_1";

my %vcfCnts = ();
foreach my $key (keys %typeMap) {
    $vcfCnts{$key} = 0;
}

my @nPloidNames = qw(A B C D E F G H I J K L M N O P Q R S T U V W X Y Z);

## Stats hash
my %stats = (
	     'hom-snp-cnt'    => 0,
	     'hom-indel-cnt'  => 0,
	     'hom-snp-rate'   => 0,
	     'hom-indel-rate' => 0,

	     'hom-inv-cnt'    => 0,
	     'hom-inv-rate'   => 0,
	     );
##
## Global data structures:
##
my %preDef = ();   ## Hash for holding pre-defined variants
my @varDat = ();   ## Holds all variant infomration
my @varMap = ();   ## Maps variants to chr and position
my @varCnt = ();   ## Keeps track of the number of variants of each type
my @useMap = ();   ## Marks used space
my @posMap = ();   ## Maps likleyhood of a variant to genomic position
my %lenMdl = ();   ## Holds size models for each variant type

##
## Global variables:
##
my $chrSeq;
my $chrId;
my $chrLen;
my $chrRange;

GetOptions(
	   "fasFile|fas|f=s@" => \$options{'fasFile'},
	   "repFile|rep|r=s@" => \$options{'repFile'},
	   "preFile|pre|p=s@" => \$options{'preFile'},
	   "gapFile|gap|g=s@" => \$options{'gapFile'},
	   "chrName|chr|c=s@" => \$options{'chrName'},

	   "outDir|out|o=s"   => \$options{'outDir'},
	   "nploidy|n=i"      => \$options{'nploidy'},

	   "homFrq=f"         => \$options{'homFrq'},
	   "snpFrq=f"         => \$options{'snp'}{'frq'},
	   "log-width=f"      => \$options{'log-width'},

	   "delFrq=f"         => \$options{'del'}{'frq'},
	   "delMin=i"         => \$options{'del'}{'min'},
	   "delMax=i"         => \$options{'del'}{'max'},
	   "delMdl=s"         => \$options{'del'}{'mdl'},
	   "delLog"           => \$options{'del'}{'log'},
	   #"delExt=i@"        => \$options{'del'}{'ext'},

	   "insFrq=f"         => \$options{'ins'}{'frq'},
	   "insMin=i"         => \$options{'ins'}{'min'},
	   "insMax=i"         => \$options{'ins'}{'max'},
	   "insMdl=s"         => \$options{'ins'}{'mdl'},
	   "insLog"           => \$options{'ins'}{'log'},
	   #"insExt=i@"        => \$options{'ins'}{'ext'},

	   "invFrq=f"         => \$options{'inv'}{'frq'},
	   "invMin=i"         => \$options{'inv'}{'min'},
	   "invMax=i"         => \$options{'inv'}{'max'},
	   "invMdl=s"         => \$options{'inv'}{'mdl'},
	   "invLog"           => \$options{'inv'}{'log'},
	   #"invExt=i@"        => \$options{'inv'}{'ext'},

	   "cnvFrq=f"         => \$options{'cnv'}{'frq'},
	   "cnvMin=i"         => \$options{'cnv'}{'min'},
	   "cnvMax=i"         => \$options{'cnv'}{'max'},
	   "cnvMdl=s"         => \$options{'cnv'}{'mdl'},
	   "cnvLog"           => \$options{'cnv'}{'log'},
           "copyMax=i"        => \$options{'cnv'}{'copyMax'},
	   #"cnvExt=i@"        => \$options{'cnv'}{'ext'},

	   "write-cmap"       => \$options{'write-cmap'},
	   "vcf-score=f"      => \$options{'vcf-score'},
	   "vcf-format|vcf"   => \$options{'vcf-format'},
	   "skipNext"         => \$skipNext,

	   "random-seed=i"   => \$options{'random-seed'},

	   "verbose|v"        => \$options{'verbose'},
	   "help|h"           => \$options{'help'},
	   );

my $programStartTime = time();

# Random seed intialization
if (defined $options{'random-seed'}) {
  srand( $options{'random-seed'} );
  print "Random seed: $options{'random-seed'}\n";
}
else {
  print "Random seed: " . srand() . "\n";
}

## Simple hack to skip the first ploid name
if ($skipNext) {
    my $tmp = shift(@nPloidNames);
}

##
## Run preprocessor: Check input parameters
##
preprocessor();

##
## Simulate each chromosome
##
foreach my $chr (keys %params) {
    ##
    ## Clear data structures and variables:
    ##
    $chrSeq   = "";
    $chrId    = "";
    $chrLen   = 0;
    $chrRange = 0;

    $options{'reference'} = "simulated-$chr";

    %preDef = ();   ## Hash for holding pre-defined variants
    @posMap = ();   ## Maps likleyhood of a variant to genomic position

    @varDat = ();   ## Holds all variant infomration
    @varMap = ();   ## Maps variants to chr and position
    @useMap = ();   ## Marks used space

    #%lenMdl = ();   ## Holds size models for each variant type

    ## Get input files
    my $fasFile = $params{$chr}{'fasFile'};
    my $preFile = $params{$chr}{'preFile'};
    my $repFile = $params{$chr}{'repFile'};
    my $gapFile = $params{$chr}{'gapFile'};

    ## Load input file: (chr seq, repeat pre-defined variants, repeat map positions)
    loadFas($fasFile); ## Load chromsome sequence
    loadPre($preFile, $chr); ## Load chromsome pre-defined variants
    loadRep($repFile, $gapFile, $chr); ## Load chromsome repeat mapping positions
    loadVar($chr);     ## Simulate all varaint types
    loadSeq();         ## Simulate/extract all variant sequence infomration
    tranChr($chr);     ## Transform each chromsome and write all outputs

    ## Write summary stats
    #writeStats();
}

$options{'run-time'} = (time() - $programStartTime) / 3600;
mssg("\# Run Time: $options{'run-time'}\n");

## Write params to params file
#writeParams();

exit(0);
## End of program execution

# ########## ########## ########## ########## ########## #

##
## Subroutines
##
sub revComp {
    my ($seq) = @_;

    $seq = uc($seq);
    $seq = reverse($seq);
    $seq =~ tr/ACTG/TGAC/;

    return($seq);
}

sub writeVcfHeader {
    my ($vcfFH) = @_;

    my $header
      = join("\n",
             '##fileformat=VCFv4.1',
             '##reference=' . $options{'reference'} . ';',
             '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">',
             '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
             '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',
             '##INFO=<ID=BKPT,Number=-1,Type=String,Description="Breakpoint type of an imprecise structural variation">',
             '##INFO=<ID=SVLEN,Number=-1,Type=Integer,Description="Difference in length between REF and ALT alleles">',
             '##INFO=<ID=MATEID,Number=1,Type=String,Description="ID of partner breakend together with which novel adjacency is formed">',
             '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
             '##INFO=<ID=CN,Number=1,Type=Integer,Description="Copy number">',
	     '##INFO=<ID=ZS,Number=1,Type=Integer,Description="Start position in simulated (transformed) genome">',
	     '##INFO=<ID=ZN,Number=1,Type=String,Description="Name of simulated chromosome copy">',
             '##ALT=<ID=DEL,Description="Deletion">',
             '##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">',
             '##ALT=<ID=INS,Description="Insertion of novel sequence">',
             '##ALT=<ID=INV,Description="Inversion">',
             '##ALT=<ID=ADJ,Description="Novel adjacency">',
	     '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
	     '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
	     '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">');

    # '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total Depth">',
    # '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
    # '##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">',
    # '##ALT=<ID=DEL:ME:ALU,Description="Deletion of ALU element">',
    # '##ALT=<ID=DEL:ME:L1,Description="Deletion of L1 element">',
    # '##ALT=<ID=DUP,Description="Duplication">',

    print($vcfFH "$header\n");
    print($vcfFH
          join("\t",
               '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
               'INFO', 'FORMAT', $options{'sample-name'}), "\n");
}

sub writeVcfLine {
    my ($fh, $dat, $newStart, $chromName) = @_;

    my @data = @{$dat};

    my $vcfSeqId  = $data[VAR_CHR_IDX];
    my $vcfPos    = $data[VAR_START_IDX];
    my $vcfEnd    = $data[VAR_END_IDX];
    my $vcfId     = $data[VAR_CHR_IDX];
    my $vcfRef    = "N";
    my $vcfAlt    = ".";
    my $vcfScore  = $options{'vcf-score'};
    my $vcfFilter = "PASS";
    my $vcfInfo   = "NS=1";
    my $upEnd     = $data[VAR_UP_ONE_IDX];
    my $dnEnd     = $data[VAR_DN_ONE_IDX];
    my $varType   = $data[VAR_TYPE_IDX];
    my $copyNum   = $data[VAR_REP_IDX];

    ##
    ## Format fields
    ##

    ## Genotype
    my $format = "GT";
    my $sampleInfo = '1/.'; # variant is present; het / hom unknown
    if ($data[VAR_ZYGO_IDX] eq "hom") {
        $sampleInfo = '1/1';
    } elsif ($data[VAR_ZYGO_IDX] eq "het") {
	  #$sampleInfo = '1/2';
      $sampleInfo = '1/0';
    }

    ## Genotype-quality
    $format .= ':GQ';
    $sampleInfo = join(':', $sampleInfo, $options{'genotype-quality'});

    ## Read depth
    $format .= ':DP';
    $sampleInfo = join(':', $sampleInfo, $options{'read-depth'});

    ## Add/fix copy number
    if ($varType eq "cnv") {
	$copyNum += 2;
    } else {
	$copyNum = 2;
    }
    $vcfInfo .= ";CN=$copyNum";

    ## Add start position in simulated genome
    $vcfInfo .= ";ZS=$newStart";
    ## Add chromosome copy name
    $vcfInfo .= ";ZN=$chromName";

    ## Build each variant case
    my ($vcfLen, $vcfType, $vcfEvent, $vcfMateId);

    if ($varType eq "snp") {
	$vcfRef    = $data[VAR_ORG_SEQ_IDX];
	$vcfAlt    = $data[VAR_NEW_SEQ_IDX];
	$vcfLen    = 0;
	$vcfType   = "SNV";
	$vcfEvent  = "SNV$vcfCnts{$varType}";
	$vcfId     = "snv_$vcfCnts{$varType}_$vcfSeqId";

	$vcfInfo .= ";SVTYPE=$vcfType;SVLEN=$vcfLen;END=$vcfEnd;EVENT=$vcfEvent";
    if ($vcfRef ne "N" and $vcfRef ne "") {
      print $fh join("\t", $vcfSeqId, $vcfPos, $vcfId, $vcfRef, $vcfAlt,
                     $vcfScore, $vcfFilter, $vcfInfo, $format, $sampleInfo) . "\n";
    }
    } elsif ($varType eq "ins") {
	$vcfRef    = $upEnd;
	$vcfAlt    = $upEnd . $data[VAR_NEW_SEQ_IDX];
	$vcfLen    = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
	$vcfType   = "INS";
	$vcfEvent  = "INS$vcfCnts{$varType}";
	$vcfPos--;
	$vcfEnd    = $vcfPos;
	$vcfId     = "ins_$vcfCnts{$varType}_$vcfSeqId";
	
	$vcfInfo .= ";SVTYPE=$vcfType;SVLEN=$vcfLen;END=$vcfEnd;EVENT=$vcfEvent";
    if ($vcfRef ne "N" and $vcfRef ne "") {
      print $fh join("\t", $vcfSeqId, $vcfPos, $vcfId, $vcfRef, $vcfAlt,
                     $vcfScore, $vcfFilter, $vcfInfo, $format, $sampleInfo) . "\n";
    }
    } elsif ($varType eq "del") {
	$vcfRef    = $upEnd . $data[VAR_ORG_SEQ_IDX];
	$vcfAlt    = $upEnd;
	$vcfLen    = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
	$vcfType   = "DEL";
	$vcfEvent  = "DEL$vcfCnts{$varType}";
	$vcfPos--;
	$vcfId     = "del_$vcfCnts{$varType}_$vcfSeqId";
	
	$vcfInfo .= ";SVTYPE=$vcfType;SVLEN=-$vcfLen;END=$vcfEnd;EVENT=$vcfEvent";
    if ($vcfAlt ne "N" and $vcfAlt ne "" and $vcfLen > 0) {
      print $fh join("\t", $vcfSeqId, $vcfPos, $vcfId, $vcfRef, $vcfAlt,
                     $vcfScore, $vcfFilter, $vcfInfo, $format, $sampleInfo) . "\n";
    }
    } elsif ($varType eq "cnv") {
	$vcfPos--;
	for (my $cn = 0; $cn < $copyNum - 2; $cn++) {
	    ## VCF 4.1 Spec: NOTE: remove for loop
	    #$vcfRef    = $upEnd;
	    #$vcfAlt    = "<DUP:TANDEM>";
	    #$vcfLen    = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
	    #$vcfType   = "DUP:TANDEM";
	    #$vcfEvent  = "DUP:TANDEM$vcfCnts{$varType}";
	    #$vcfPos--;
	    #$vcfId     = "dupTandem_$vcfCnts{$varType}_$vcfSeqId";
	    
	    ## EAGLE SPEC: Needs for loop
	    my $char    = substr($data[VAR_ORG_SEQ_IDX], 0, 1);
	    $vcfLen     = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
	    my $start   = $vcfPos + $vcfLen - 1;
	    my $matePos = $vcfPos;
	    $vcfRef     = $char;
	    $vcfAlt     = $char . "[$vcfSeqId:$matePos" . "[";
	    $vcfType    = "DUP:TANDEM";
	    $vcfEvent   = "DUP:TANDEM$vcfCnts{$varType}";
	    $vcfId      = "dupTandem_$vcfCnts{$varType}_$vcfSeqId" . "_$cn.1";
	    
        if ($vcfRef ne "N" and $vcfRef ne "") {
          $vcfInfo .= ";SVTYPE=$vcfType;SVLEN=$vcfLen;END=$vcfEnd;EVENT=$vcfEvent";
          print $fh join("\t", $vcfSeqId, $start, $vcfId, $vcfRef, $vcfAlt,
                         $vcfScore, $vcfFilter, $vcfInfo, $format, $sampleInfo) . "\n";
	    
          $start      = $vcfPos;
          $matePos    = $vcfPos + $vcfLen - 1;
          $vcfAlt     = "]$vcfSeqId:$matePos]$char";
          $vcfId      = "dupTandem_$vcfCnts{$varType}_$vcfSeqId" . "_$cn.2";
          print $fh join("\t", $vcfSeqId, $start, $vcfId, $vcfRef, $vcfAlt,
                         $vcfScore, $vcfFilter, $vcfInfo, $format, $sampleInfo) . "\n";
        }
      }
    } elsif ($varType eq "inv") {    
	## bnd_w
	my $altPos = $vcfEnd;
	$vcfPos--;
	$vcfRef    = $upEnd;
    if ($vcfRef ne "N" and $vcfRef ne "") {
      $vcfAlt    = "$upEnd]$vcfSeqId:$altPos]";
      $vcfLen    = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
      $vcfType   = "BND";
      $vcfEvent  = "INV$vcfCnts{$varType}";
      $vcfId     = "inv_bnd_W_$vcfCnts{$varType}_$vcfSeqId";
      $vcfMateId = "inv_bnd_X_$vcfCnts{$varType}_$vcfSeqId";
      $vcfInfo .= ";SVTYPE=$vcfType;SVLEN=$vcfLen;END=$vcfEnd;EVENT=$vcfEvent";
      print $fh join("\t", $vcfSeqId, $vcfPos, $vcfId, $vcfRef, $vcfAlt,
                     $vcfScore, $vcfFilter, "$vcfInfo;MATEID=$vcfMateId", $format, $sampleInfo) . "\n";
	
      ## bnd_V
      $vcfPos++;
      $altPos = $vcfEnd + 1;
      $vcfRef    = substr($data[VAR_ORG_SEQ_IDX], 0, 1);
      $vcfAlt    = "[$vcfSeqId:$altPos" . "[" . substr($data[VAR_ORG_SEQ_IDX], 0, 1);
      $vcfLen    = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
      $vcfType   = "BND";
      $vcfId     = "inv_bnd_V_$vcfCnts{$varType}_$vcfSeqId";
      $vcfMateId = "inv_bnd_U_$vcfCnts{$varType}_$vcfSeqId";
      print $fh join("\t", $vcfSeqId, $vcfPos, $vcfId, $vcfRef, $vcfAlt,
                     $vcfScore, $vcfFilter, "$vcfInfo;MATEID=$vcfMateId", $format, $sampleInfo) . "\n";
	
      ## bnd_U
      my $char   = substr($data[VAR_ORG_SEQ_IDX], length($data[VAR_ORG_SEQ_IDX]) - 1, 1);
      $vcfPos    = $data[VAR_END_IDX];
      $altPos    = $data[VAR_START_IDX] - 1;
      $vcfRef    = $char;
      $vcfAlt    = "$char]$vcfSeqId:$altPos]";
      $vcfLen    = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
      $vcfType   = "BND";
      $vcfId     = "inv_bnd_U_$vcfCnts{$varType}_$vcfSeqId";
      $vcfMateId = "inv_bnd_V_$vcfCnts{$varType}_$vcfSeqId";
      print $fh join("\t", $vcfSeqId, $vcfPos, $vcfId, $vcfRef, $vcfAlt,
                     $vcfScore, $vcfFilter, "$vcfInfo;MATEID=$vcfMateId", $format, $sampleInfo) . "\n";
	
      ## bnd_X
      $vcfPos    = $data[VAR_END_IDX] + 1;
      $altPos    = $data[VAR_START_IDX];
      $vcfRef    = $dnEnd;
      $vcfAlt    = "[$vcfSeqId:$altPos" . "[$dnEnd";
      $vcfLen    = $data[VAR_END_IDX] - $data[VAR_START_IDX] + 1;
      $vcfType   = "BND";
      $vcfId     = "inv_bnd_X_$vcfCnts{$varType}_$vcfSeqId";
      $vcfMateId = "inv_bnd_W_$vcfCnts{$varType}_$vcfSeqId";
      print $fh join("\t", $vcfSeqId, $vcfPos, $vcfId, $vcfRef, $vcfAlt,
                     $vcfScore, $vcfFilter, "$vcfInfo;MATEID=$vcfMateId", $format, $sampleInfo) . "\n";
    }
	$vcfCnts{$varType}++;
    } elsif ($varType eq "trn") {
	die("[ERROR]: Currently do NOT support inter (chimeric) chromosomal variants!");
    } else {
	die("[ERROR]: Attempting to convert unknown type: $varType!");
    }
    $vcfCnts{$varType}++;
}

sub writeBedLine {
    my ($fh, $dat, $newStart, $ploidyChr, $idx) = @_;

    my $chr    = $$dat[VAR_CHR_IDX];
    my $start  = $$dat[VAR_START_IDX];
    my $end    = $$dat[VAR_END_IDX] + 1;

    ## Update start and stop
    if ($$dat[VAR_TYPE_IDX] eq "ins") {
	$end = $start + 1;
    }    
    print $fh "$chr\t$start\t$end\t$idx\t$ploidyChr;$idx";
    foreach my $val (@{$dat}) {
	print $fh ";$val";
    }
    print $fh "\n";
}

sub tranChr {
    my ($chrName) = @_;
    ##
    ## Description: By looping backwards through each chromosome
    ##  we can add all vairants into the new chromosome and create
    ##  alignment strings and coordinate maps. This should be 
    ##  easy because the sequences have already been calculated.
    ##
    ## All output files are also written:
    ##  outDir/chrX[A-Z].mutMap.txt
    ##  outDir/chrX[A-Z].crdMap.txt
    ##  outDir/chrX[A-Z].alnMap.fas
    ##  outDir/chrX[A-Z].trnChr.fas
    ##

    for (my $nploid = 0; $nploid < $options{'nploidy'}; $nploid++) {
	#my $chrName = $options{'chrName'}[$nploid];
	my $name = $nPloidNames[$nploid];
	my $mutMapFile = "$options{'outDir'}/$chrName-$name.mutMap.txt";
	my $crdMapFile = "$options{'outDir'}/$chrName-$name.crdMap.txt";
	my $alnMapFile = "$options{'outDir'}/$chrName-$name.alnMap.fas";
	my $trnChrFile = "$options{'outDir'}/$chrName-$name.trnChr.fas";
	my $vcfMapFile = "$options{'outDir'}/$chrName-$name.mutMap.vcf";

	mssg("\# Transforming: $chrName-$name\n");

	open(MMAP, ">$mutMapFile")
	    or die "Could not open $mutMapFile for writing: $!";
	open(CMAP, ">$crdMapFile")
	    or die "Could not open $crdMapFile for writing: $!";
	open(VMAP, ">$vcfMapFile")
	    or die "Could not open $vcfMapFile for writing: $!";
	open(AGEN, ">$alnMapFile")
	    or die "Could not open $alnMapFile for writing: $!";
	open(TGEN, ">$trnChrFile")
	    or die "Could not open $trnChrFile for writing: $!";

        my $tgSeq = uc($chrSeq);  ## Transformed genome seq
        my $oaSeq = uc($chrSeq);  ## Original genome aligned seq
        my $taSeq = uc($chrSeq);  ## Transformed genome aligned seq
	my %tgMap = ();           ## Transformed genome mutation map

	my @tgPosHeap = ();       ## Array to keep track of all transformed
	                          ##  positions in order
	my @moves = ();           ## Array to keep track of all movements

        my $lenDelta  = 0;
        my $lenChange = 0;

        ##
        ## Transform the genome by inserting mutation backwards
        ##
        my $prevStart = -1;
        my $prevEnd   = -1;

	my $tenMapCnt = int(scalar(keys %{$varMap[$nploid]}) / 10);
	my $mapCount  = 0;
	my $mssgCnt   = 0;

	foreach my $varPos (sort {$b <=> $a} keys %{$varMap[$nploid]}) {
	    if ($mapCount % $tenMapCnt == 0) {
		mssg("\#\tOn $mssgCnt\t$mapCount\n");
		$mssgCnt++;
	    }
	    $mapCount++;

	    my $move = 0;

	    my $varIdx = $varMap[$nploid]{$varPos};

	    $varPos--; ## Fix to 0 based genomic positions

	    my $varType   = $varDat[$varIdx][VAR_TYPE_IDX];
	    my $varStart  = $varDat[$varIdx][VAR_START_IDX];
	    my $varEnd    = $varDat[$varIdx][VAR_END_IDX];
	    my $varRep    = $varDat[$varIdx][VAR_REP_IDX];
	    my $varZygo   = $varDat[$varIdx][VAR_ZYGO_IDX];
	    my $varOrgSeq = $varDat[$varIdx][VAR_ORG_SEQ_IDX];
	    my $varNewSeq = $varDat[$varIdx][VAR_NEW_SEQ_IDX];

	    #$prevStart = $varPos;
	    #$prevEnd   = $varPos + length();
next;
	    if ($varType eq "ins") {
                ## Update transfomed genome sequence and alignment seqs
                substr($tgSeq, $varPos, 0, $varNewSeq);
                substr($oaSeq, $varPos, 0, "-" x length($varNewSeq));
                substr($taSeq, $varPos, 0, lc($varNewSeq));

                $lenChange += length($varNewSeq);
                $lenDelta  += length($varNewSeq);
                $move      += length($varNewSeq);
	    } elsif ($varType eq "cnv") {
		## We treat CNV's as multiple insertions right next to
		##  one another
		for (my $cnt = 0; $cnt < $varRep; $cnt++) {
		    ## Update transfomed genome sequence and alignment seqs
		    substr($tgSeq, $varPos, 0, $varOrgSeq);
		    substr($oaSeq, $varPos, 0, "-" x length($varOrgSeq));
		    substr($taSeq, $varPos, 0, lc($varOrgSeq));

		    $lenChange += length($varOrgSeq);
		    $lenDelta  += length($varOrgSeq);
		    $move      += length($varOrgSeq);
		}
	    } elsif ($varType eq "inv") {
		##
                ## [SANITY CHECK]: Make sure all original nucleotides match
		##
                for (my $ii = 0; $ii < length($varOrgSeq); $ii++) {
                    my $org_nuc = substr($oaSeq, $varPos + $ii, 1);
                    my $map_nuc = substr($varOrgSeq, $ii, 1);
                    if (uc($org_nuc) ne uc($map_nuc)) {
                        my $org_seq = substr($oaSeq, $varPos, length($varOrgSeq));
                        errMssg("During genome transformation (inversion): " .
				"mapped index ($varPos + $ii, $map_nuc) " .
				"does not match genomic value: $org_nuc\n" . 
				"[ERROR]: Genome delete seq : $org_seq\n" .
				"[ERROR]: Extract delete seq: $varOrgSeq\n");
                    }
                }

		my $invSeq = reverse($varOrgSeq);
                $invSeq =~ tr/ACTG/TGAC/;

                substr($tgSeq, $varPos, length($varOrgSeq), $invSeq);
                substr($oaSeq, $varPos, length($varOrgSeq), lc($varOrgSeq));
                substr($taSeq, $varPos, length($varOrgSeq), lc($invSeq));

                #$lenChange += length($varOrgSeq);
                #$lenDelta  -= length($varOrgSeq);
                #$move      -= length($varOrgSeq);
	    } elsif ($varType eq "del") {
		##
                ## [SANITY CHECK]: Make sure all original nucleotides match
		##
                for (my $ii = 0; $ii < length($varOrgSeq); $ii++) {
                    my $org_nuc = substr($oaSeq, $varPos + $ii, 1);
                    my $map_nuc = substr($varOrgSeq, $ii, 1);
                    if (uc($org_nuc) ne uc($map_nuc)) {
                        my $org_seq = substr($oaSeq, $varPos, length($varOrgSeq));
                        errMssg("During genome transformation (deletion): " .
				"mapped index ($varPos + $ii, $map_nuc) " .
				"does not match genomic value: $org_nuc\n" . 
				"[ERROR]: Genome delete seq : $org_seq\n" .
				"[ERROR]: Extract delete seq: $varOrgSeq\n");
                    }
                }

                substr($tgSeq, $varPos, length($varOrgSeq), "");
                substr($oaSeq, $varPos, length($varOrgSeq), lc($varOrgSeq));
                substr($taSeq, $varPos, length($varOrgSeq), "-" x length($varOrgSeq));

                $lenChange += length($varOrgSeq);
                $lenDelta  -= length($varOrgSeq);
                $move      -= length($varOrgSeq);
	    } elsif ($varType eq "snp") {
		##
                ## [SANITY CHECK]: Make sure all original nucleotides match
		##
                my $orig_nuc = substr($oaSeq, $varPos, 1);
                if (uc($orig_nuc) ne uc($varOrgSeq)) {
                    errMssg("During genome transformation: mapped index ($varPos, $varOrgSeq) " .
                        "does not match genomic value: $orig_nuc\n");
                }
                substr($tgSeq, $varPos, 1, $varNewSeq);
                substr($oaSeq, $varPos, 1, lc($varOrgSeq));
                substr($taSeq, $varPos, 1, lc($varNewSeq));
	    } else {
		errMssg("[line: " . __LINE__ . "] Type ($varType) is not valid\n");
		errMssg("[line: " . __LINE__ . "] Occured on Nploidy: $nploid, chr: $chrName\n");
	    }

            ##
            ## Update transformed mutation map
            ##
	    $varPos++;

	    unshift @tgPosHeap, $varPos;
	    unshift @moves, $move;

            #$tgMap{$varPos} = 0;
            #foreach my $tgPos (keys %tgMap) {
            #    $tgMap{$tgPos} += $move;
            #}
	}

	##
	## Build the tgMap hash
	##
	if (scalar(@tgPosHeap) != scalar(@moves)) {
	    errMssg("[line: " . __LINE__ . "] " .
		    "tgPosHeap != moves: " . scalar(@tgPosHeap) . " != " . scalar(@moves) . "\n");
	}
	my $sumMove = 0;
	for (my $idx = 0; $idx < @tgPosHeap; $idx++) {
	    my $tgPos = $tgPosHeap[$idx];
	    my $move  = $moves[$idx];
	    $tgMap{$tgPos} = $sumMove;
	    $sumMove += $move;
	}

	##
        ## [SANITY CHECK]: Check that the original and tranformed 
	##  aligned genome sequence are the same length.
	##
        if (length($oaSeq) != length($taSeq)) {
            errMssg("Original and trasformed aligned genome sequences " .
		    "are different length: " . length($oaSeq) . " != " . 
		    length($taSeq) . "\n");
        }

	##
        ## Print mutation map file
	##
	## Print MMAP Header
	mssg("\# Writing mutation map\n");
	writeVcfHeader(\*VMAP);

	print MMAP "\# This file contains variant mappings. Chromosome source and name " .
	    "can be found in fasta file format, the remaining fields are described below.\n";
	print MMAP "\#! $mmapHeader\n";
	print MMAP ">$chrName-$name\n";

	my $mmCnt = 0;
        foreach my $varPos (sort {$a <=> $b} keys %{$varMap[$nploid]}) {
	    my $varIdx    = $varMap[$nploid]{$varPos};
	    my $varOrgSeq = $varDat[$varIdx][VAR_ORG_SEQ_IDX];

            if (! defined($tgMap{$varPos})) {
                #errMssg("No transformed map position for original position $varPos\n")
              $tgMap{$varPos} = 0;
            }

            ## Correct position for mutation length...
            my $nPos = $varPos + $tgMap{$varPos};

	    ##
	    ## Output format:
	    ##
	    ## 1. Transfom Genomic Position
	    ## 2. Original Genomic Position
	    ## 3. Rest of @varDat fields
	    writeVcfLine(\*VMAP, \@{$varDat[$varIdx]}, $nPos, "$chrName-$name");
	    writeBedLine(\*MMAP, \@{$varDat[$varIdx]}, $nPos, $name, $mmCnt);

	    #print MMAP "$nPos";
	    #foreach my $val (@{$varDat[$varIdx]}) {
	    #print MMAP "\t$val";
	    #}
	    #print MMAP "\n";
	    $mmCnt++;
        }

	## Print aligned reference/transformed genome
    if (0) {
	mssg("\# Writing reference to transformed alignment map\n");
        print AGEN ">reference-$chrName-$name\n$oaSeq\n";
        print AGEN ">transform-$chrName-$name\n$taSeq\n";

        ## Print transformed genome
	mssg("\# Writing transformed genome\n");
        print TGEN ">$chrName-$name\n$tgSeq\n";

        ## Write mutation coordinate map file
        my $delta1 = 0;
        my $delta2 = 0;
        my $ogIdx  = 0;
        my $tgIdx  = 0;
	if ($options{'write-cmap'}) {
	    mssg("\# Writing coordinate map\n");
	    print CMAP "\# This file contains coordinate mappings indexed for the complete " .
		"reference to transformed chromosome alignements. Chromosome source and name " .
		"can be found in fasta file format, the remaining fields are described below.\n";
	    print CMAP "\#! $cmapHeader\n";
	    print CMAP ">$chrName-$name\n";
	    for (my $ogPos = 0; $ogPos < length($oaSeq); $ogPos++) {
		if (substr($oaSeq, $ogPos, 1) eq "-") { $ogIdx--; }
		if (substr($taSeq, $ogPos, 1) eq "-") { $tgIdx--; }
		
		$ogIdx++;
		$tgIdx++;
		print CMAP "$ogPos\t$ogIdx\t$tgIdx\n";
	    }
	}
  }

	close(MMAP);
	close(CMAP);
	close(VMAP);
	close(AGEN);
	close(TGEN);
    }

}

sub loadSeq {
    ##
    ## Description: For each variant we will sequentialy substring
    ##  out the variant sequence. This is done in order so that our
    ##  substring method is not doing random cacheing. At the same
    ##  time we can pick the new variant sequecnes for SNPs, inversion,
    ##  and insertions. For deletions and CNVs the new sequence 
    ##  should be trival
    ##
    mssg("\# Loading variant sequence data\n");

    for (my $varIdx = 0; $varIdx < @varDat; $varIdx++) {
	my $type  = $varDat[$varIdx][VAR_TYPE_IDX];
	my $start = $varDat[$varIdx][VAR_START_IDX] - 1;
	my $end   = $varDat[$varIdx][VAR_END_IDX] - 1;
	my $len   = $end - $start + 1;

	## Quick santiy check
	if ($start + $len > $chrLen - 1) {
	    foreach my $val (@{$varDat[$varIdx]}) {
		print STDERR "DAT\t$val\n";
	    }
#	    errMssg("Start + len ($start + $len) > chrLen ($chrLen)\n");
        next;
	}

	## First set the original sequence
	my $orgSeq = substr($chrSeq, $start, $len);
	$varDat[$varIdx][VAR_ORG_SEQ_IDX] = $orgSeq;

	## Set the one base pair up and downstream sequences
	$varDat[$varIdx][VAR_UP_ONE_IDX] = substr($chrSeq, $start - 1, 1);
	$varDat[$varIdx][VAR_DN_ONE_IDX] = substr($chrSeq, $start + $len, 1);

	## Generate new sequence
	if ($type eq "snp") {
	    ## Ensure the SNP is of length 1
	    if (length($orgSeq) != 1) {
		errMssg("[line: " . __LINE__ . "] Type ($type) is SNP, but " .
			"length($orgSeq) != 1, during loadSeq\n");
	    }

	    ## Ensure the new SNP is different
	    my $rIdx   = rand(length($options{'bases'}));
	    my $newSeq = substr($options{'bases'}, $rIdx, 1);
	    while($newSeq eq $orgSeq) {
		$rIdx   = rand(length($options{'bases'}));
		$newSeq = substr($options{'bases'}, $rIdx, 1);
	    }

	    ## Set new sequence
	    $varDat[$varIdx][VAR_NEW_SEQ_IDX] = $newSeq;
	} elsif ($type eq "ins") {
	    ## Generate new insertion sequence from previous insertion
	    ##  sequence by doing triplicate swapping.
	    ## The idea is that this will make the insert sequence more
	    ##  realistic.
	    ##
	    ## NOTE: If the event is a real event and the sequence already
	    ##  exist we can skip this step
	    if ($varDat[$varIdx][VAR_SRC_IDX] eq "real" && defined($varDat[$varIdx][VAR_NEW_SEQ_IDX]) &&
		length($varDat[$varIdx][VAR_NEW_SEQ_IDX]) != 0) {
		## Do nothing...
	    } else {
		my $newSeq = swapSeq($orgSeq);

		## Set new sequence
		$varDat[$varIdx][VAR_NEW_SEQ_IDX] = $newSeq;
	    }
	} elsif ($type eq "cnv" ||
		 $type eq "inv" ||
		 $type eq "del") {
	    ## Do nothing, we won't save the new seqs
	    ##  because they're trivial to calculate (inv) 
	    ##  or they're the same (cnv) or they're
	    ##  nothing (del)

	    ## Set new sequence
	    $varDat[$varIdx][VAR_NEW_SEQ_IDX] = "";
	} else {
	    errMssg("[line: " . __LINE__ . "] Unknown type: $type during loadSeq\n");
	}
    }

    mssg("\# Done loading all variant sequence data\n\n");
}

sub loadVar {
    my ($chrName) = @_;

    ##
    ## Descriptions: For each variant type we will simulate new 
    ##  variants until we've reached the number of intended variants
    ##  of each type. Sequence content is not determined at this point
    ##  for random access substring limitations. That is done in the
    ##  loadSeq subroutine in genomic order.
    ##

    foreach my $varType (@allowedTypes) {
	next if ($varType eq "trn");

	if (! defined($options{$varType}{'frq'})) {
	    errMssg("No frequency defined for $varType\n");
	}
	next if ($options{$varType}{'frq'} == 0);

	mssg("\# Adding Variants: $varType...\n");

	##
	## Loop over each ploidy chromosome to be simulated
	##
	for (my $nploid = 0; $nploid < $options{'nploidy'}; $nploid++) {
	    ## Intialize counters
	    my $varIdx  = 0;
	    my $curTime = time();
	    $stats{$varType}{'target'} = int($chrLen * $options{$varType}{'frq'});

	    if (! defined($varCnt[$nploid]{$varType})) {
		$varCnt[$nploid]{$varType} = 0;
	    }

	    mssg("\#\tOn $nploid target: $stats{$varType}{'target'}, target: $stats{$varType}{'target'}\n"); 
	    
	    ##
	    ## Continue to get or generate variants until variant 
	    ##  threshold is met
	    ##
	    while ($stats{$varType}{'target'} > $varCnt[$nploid]{$varType}) {
        my $modVal = int($stats{$varType}{'target'} / 10);
        if ($modVal != 0 && $varCnt[$nploid]{$varType} % $modVal == 0) {
		    mssg("\#\t\t$varCnt[$nploid]{$varType}\n");
		}

		my ($varStart, $varEnd, $rep, $varLen, $lenIdx);
		my $chr = $chrName;
		my $src = "sim";

		##
		## Generate a new mutation at random
		##
		if ($varType ne "snp") {
		    $lenIdx = int(rand(scalar(@{$lenMdl{$varType}})));
		    $varLen = $lenMdl{$varType}[$lenIdx];
		} else {
		    $varLen = 1;
		}

		## This should generate the new random position, but with 
		##  repeat region preference
		my $posIdx = int(rand(@posMap));
		#$varStart = ($posMap[$posIdx] * $options{'rep-res'}) + int(rand($options{'rep-res'}));
		if (! defined($posMap[$posIdx])) {
		    $varStart = $posIdx;
		} else {
		    $varStart = $posMap[$posIdx] + int(rand($options{'rep-res'}));
		}
		$varEnd   = $varStart + $varLen - 1;

		## Shorten event if it extends off the end of a chromosome
		if ($varEnd >= $chrLen - ($varLen / 2)) {
            $varLen = int($varLen / 2);
		    $varEnd = $chrLen - $varLen;
		}

		## Generate copy number stats
		$rep = 1;
		if ($varType eq "cnv") {
		    $rep = int(rand($options{$varType}{'copyMax'} - 1)) + 1;
		}

		## Determine zygosity
		my $zygo = "het";
		if ($nploid == 0) {
		    my $r = rand();
		    if ($r < $options{'homFrq'}) {
			$zygo = "hom";
		    }
		}

		## Load @varDat, @varMap, @useMap
		if ($zygo eq "hom") {
		    ## Loop over all chromosomes and make sure it can be added
		    my $passed = TRUE;
		    for (my $n = 0; $n < $options{'nploidy'}; $n++) {
			for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
			    if (defined($useMap[$n]{$pos})) {
				$passed = FALSE;
				last;
			    }
			}
			last if (! $passed);
		    }
		    if ($passed) {
			## Write all map positions
			for (my $n = 0; $n < $options{'nploidy'}; $n++) {
			    for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
				$useMap[$n]{$pos} = 1;
			    }
			    ## Save individual ploidy mapping
			    $varMap[$n]{$varStart} = scalar(@varDat);
			    
			    ## Update count metric
			    if (! defined($varCnt[$n]{$varType})) {
				$varCnt[$n]{$varType} = 1;
			    } else {
				$varCnt[$n]{$varType}++;
			    }
			}
			## Save to pre-def hash
			push @{$preDef{$varType}}, [ $chr, $varStart, $varEnd, $rep, $zygo, $src ];
			## Save variant
			push @varDat, [ $varType, $chr, $varStart, $varEnd, $rep, $zygo, $src ];
		    } else {
			#mssg("[Warning]: Unable to add mut $varType, $varStart, $varEnd because it was previously " .
			#     "defined during het assignment!\n");
		    }
		} else {
		    ## Only mark chr A
		    ## Loop over all chromosomes and make sure it can be added
		    my $passed = TRUE;
		    for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
			if (defined($useMap[$nploid]{$pos})) {
			    $passed = FALSE;
			    last;
			}
		    }
		    if ($passed) {
			## Write all map positions
			for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
			    $useMap[$nploid]{$pos} = 1;
			}
			## Update count metric
			if (! defined($varCnt[$nploid]{$varType})) {
			    $varCnt[$nploid]{$varType} = 1;
			} else {
			    $varCnt[$nploid]{$varType}++;
			}
			
			## Save individual ploidy mapping
			$varMap[$nploid]{$varStart} = scalar(@varDat);
			## Save to pre-def hash
			push @{$preDef{$varType}}, [ $chr, $varStart, $varEnd, $rep, $zygo, $src ];
			## Save variant
			push @varDat, [ $varType, $chr, $varStart, $varEnd, $rep, $zygo, $src ];
		    } else {
			#mssg("[Warning]: Unable to add mut $varType, $varStart, $varEnd because it was previously " .
			#     "defined during hom assignment!\n");
		    }
		}
	    }
	    ##
	    ## Write stats message: for current ploidy chromosome
	    ##
	}

	##
	## Write stats message: for all ploidy chromsomes
	##

    }

    mssg("\# Done loading\n\#\n");
}

##
## File loading subroutines
##
sub loadFas {
    my ($file) = @_;

    $chrSeq = "";
    $chrId  = "";
    $chrLen = 0;

    open(IN, "<$file")
	or die "Could not open $file for reading: $!";

    mssg("\# Loading fasta file file: $file\n");
    while(<IN>) {
	s/^\s+//;
	s/\s+$//;
	next if /^(\s)*$/;

	my $line = $_;
	if ($line =~ m/^>(.*)$/) {
	    $chrId = $1;
	} else {
	    $chrSeq .= uc($line);
	}
    }
    close(IN);

    ## Quick fix...
    #$chrSeq = substr($chrSeq, 0, length($chrSeq) - 20000);

    $chrLen = length($chrSeq);

    mssg("\# Loaded id: $chrId with length: $chrLen\n");
}

sub loadPre {
    my ($file, $curChr) = @_;

    %preDef = ();

    if (! defined($file)) {
	return(0);
    }

    open(IN, "<$file")
	or die "Could not open $file for reading: $!";

    mssg("\# Loading pre-def file file: $file\n");
    while(<IN>) {
	s/^\s+//;
	s/\s+$//;
	next if /^(\s)*$/;

	my $line = $_;
	my @lineArray = split(/\t/, $line);

	my $varType  = lc($lineArray[PRE_TYPE_IDX]);
	my $chr      = $lineArray[PRE_CHR_IDX];
	my $varStart = $lineArray[PRE_START_IDX];
	my $varEnd   = $lineArray[PRE_END_IDX];
	my $rep      = $lineArray[PRE_REP_IDX];
	my $src      = "real";
	my $seq      = "";

	next if ($lineArray[PRE_CHR_IDX] ne $curChr);

	if (defined($lineArray[PRE_SEQ_IDX]) && length($lineArray[PRE_SEQ_IDX]) != 0) {
	    $seq = $lineArray[PRE_SEQ_IDX];
	    $varEnd = $varStart + length($seq) - 1; 
	}

	if (! defined($rep) || $rep < 0) {
	    $rep = 1;
	}

	## Determine zygosity
	my $r = rand();
	my $zygo = "het";
	if ($r < $options{'homFrq'}) {
	    $zygo = "hom";
	}

	## We will skip translocations for now...
	if ($varType eq "trn") {
	    push @{$preDef{$varType}}, [ $chr, $varStart, $varEnd, $rep, $zygo ];
	    next;
	}

	## Load @varDat, @varMap, @useMap
	if ($zygo eq "hom") {
	    ## Loop over all chromosomes and make sure it can be added
	    my $passed = TRUE;
	    for (my $nploid = 0; $nploid < $options{'nploidy'}; $nploid++) {
		for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
		    if (defined($useMap[$nploid]{$pos})) {
			$passed = FALSE;
			last;
		    }
		}
		last if (! $passed);
	    }
	    if ($passed) {
		## Write all map positions
		for (my $nploid = 0; $nploid < $options{'nploidy'}; $nploid++) {
		    for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
			$useMap[$nploid]{$pos} = 1;
		    }
		    ## Save individual ploidy mapping
		    $varMap[$nploid]{$varStart} = scalar(@varDat);

		    ## Update count metric
		    if (! defined($varCnt[$nploid]{$varType})) {
			$varCnt[$nploid]{$varType} = 1;
		    } else {
			$varCnt[$nploid]{$varType}++;
		    }
		}
		## Save to pre-def hash
		push @{$preDef{$varType}}, [ $chr, $varStart, $varEnd, $rep, $zygo, $src ];
		## Save variant
		push @varDat, [ $varType, $chr, $varStart, $varEnd, $rep, $zygo, $src, "", $seq ];
	    } else {
		mssg("[Warning]: Unable to add mut $varType, $varStart, $varEnd because it was previously " .
		     "defined during het assignment!\n");
	    }
	} else {
	    ## Only mark chr A
	    ## Loop over all chromosomes and make sure it can be added
	    my $passed = TRUE;
	    my $nploid = 0;
	    for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
		if (defined($useMap[$nploid]{$pos})) {
		    $passed = FALSE;
		    last;
		}
	    }
	    if ($passed) {
		## Write all map positions
		my $nploid = 0;
		for (my $pos = $varStart - 1 - $options{'varSep'}; $pos < $varEnd + $options{'varSep'}; $pos++) {
		    $useMap[$nploid]{$pos} = 1;
		}
		## Update count metric
		if (! defined($varCnt[$nploid]{$varType})) {
		    $varCnt[$nploid]{$varType} = 1;
		} else {
		    $varCnt[$nploid]{$varType}++;
		}

		## Save individual ploidy mapping
		$varMap[$nploid]{$varStart} = scalar(@varDat);
		## Save to pre-def hash
		push @{$preDef{$varType}}, [ $chr, $varStart, $varEnd, $rep, $zygo, $src ];
		## Save variant
		push @varDat, [ $varType, $chr, $varStart, $varEnd, $rep, $zygo, $src, "", $seq ];
	    } else {
		mssg("[Warning]: Unable to add mut $varType, $varStart, $varEnd because it was previously " .
		     "defined during hom assignment!\n");
	    }
	}
    }
    close(IN);

    foreach my $varType (sort keys %preDef) {
	mssg("\# Loaded " . scalar(@{$preDef{$varType}}) . " pre defined variants of type: $varType\n");
    }
}

sub loadRep {
    my ($repFile, $gapFile, $curChr) = @_;

    @posMap = ();

    if (! defined($repFile)) {
	for (my $ii = 0; $ii < $chrLen; $ii++) {
	    $posMap[$ii] = $ii;
	}
	return(0);
    }

    ##
    ## First load the gap regions
    ##
    open(IN, "<$gapFile")
	or die "Could not open $gapFile for reading: $!";

    mssg("\# Loading gap region file: $gapFile\n");
    my @gaps = ();
    while(<IN>) {
	s/^\s+//;
	s/\s+$//;
	next if /^(\s)*$/;
	next if /^\#.*$/;

	my $line = $_;
	my @fields = split(/\t/, $line);

	next if ($fields[BED_CHR] ne $curChr);
	#if ($fields[BED_CHR] ne $options{'chrName'}) {
	    #print STDERR "Skipping: '$fields[BED_CHR]' != '$options{'chrName'}'\n";
	    #next;
	#}

	my $gapStart = $fields[BED_START] - 1;
	my $gapEnd   = $fields[BED_END] - 1;
	my $len   = $gapEnd - $gapStart;

	if ($len < 0) {
	    errMssg("Start > end ($gapStart > $gapEnd) on line: $line\n");
	}
	if ($len < $options{'rep-res'}) {
	    mssg("[Warning]: Repeat region less repeat resolution ($len < $options{'rep-res'}), " .
		 "skipping region!\n");
	    next;
	}

	push @gaps, [ $gapStart, $gapEnd ];
    }
    close(IN);

    mssg("\# Loaded: " . @gaps . " gapped regions.\n");

    open(IN, "<$repFile")
	or die "Could not open $repFile for reading: $!";

    mssg("\# Loading repeat region file: $repFile\n");
    my $curIdx = 0;
    my $curPos = 0;
    my $regCnt = 0;
    my $skpCnt = 0;
    my $mskCnt = 0;

    while(<IN>) {
	s/^\s+//;
	s/\s+$//;
	next if /^(\s)*$/;
	next if /^\#.*$/;

	my $line = $_;
	my @fields = split(/\t/, $line);

	next if ($fields[BED_CHR] ne $curChr);

	my $repStart = $fields[BED_START] - 1;
	my $repEnd   = $fields[BED_END] - 1;
	my $len   = $repEnd - $repStart;

	if ($len < 0) {
	    errMssg("Start > end ($repStart > $repEnd) on line: $line\n");
	}
	if ($len < $options{'rep-res'}) {
	    #mssg("[Warning]: Repeat region less repeat resolution ($len < $options{'rep-res'}), " .
	    #"skipping region!\n");
	    $skpCnt++;
	    next;
	}

	while($curPos < $repStart) {
	    if (! checkInter(\@gaps, $curPos)) {
		$posMap[$curIdx] = $curPos;
		$curIdx++;
	    } else {
		$mskCnt++;
	    }
	    $curPos += $options{'rep-res'};
	}

	while($curPos <= $repEnd) {
	    if (! checkInter(\@gaps, $curPos)) {
		for (my $ii = 0; $ii < $options{'rep-width'}; $ii++) {
		    $posMap[$curIdx] = $curPos;
		    $curIdx++;
		}
	    } else {
		$mskCnt++;
	    }
	    $curPos += $options{'rep-res'};
	}
	$regCnt++;
    }
    close(IN);

    while($curPos < $chrLen) {
	if (! checkInter(\@gaps, $curPos)) {
	    $posMap[$curIdx] = $curPos;
	    $curIdx++;
	} else {
	    $mskCnt++;
	}
	$curPos += $options{'rep-res'};
    }

    ## Set the chr search representative range
    $chrRange = $curPos;

    mssg("\# Loaded $regCnt regions, $curPos positions mapped to $curIdx space\n");
    mssg("\# Skipped: $skpCnt via too small and maksed via gaps: $mskCnt\n\#\n");
}

##
## Preprocessing subroutines
##
sub printHelp
{
    print STDERR "Options:\n\n";
    print STDERR "\t-fasFile, -fas, -f        Reference genome fasta file [REQUIRED]\n";
    print STDERR "\t-repFile, -rep, -r        Repeat region file [OPTIONAL]\n";
    print STDERR "\t-preFile, -pre, -p        Pre-defined variants file [OPTIONAL]\n";
    print STDERR "\t-gapFile, -gap, -g        Contig gaps in genome assembly file [OPTIONAL]\n\n";

    print STDERR "\t-nploidy, -n              Target ploidy to simulate (DEFAULT 1)\n";
    print STDERR "\t-homFrq                   Homozygous frequency between all chrs for all variatns (DEFAULT 0)\n";
    print STDERR "\t-log-width                Width paramater for logorimic variant size simulation (DEFAULT 0.25)\n\n";

    print STDERR "\t-snpFrq                   SNP frequency (DEFAULT 0)\n\n";

    print STDERR "\t-delFrq                   Deletion frequency (DEFAULT 0)\n";
    print STDERR "\t-delMin                   Deletion minimum size (DEFAULT 1)\n";
    print STDERR "\t-delMax                   Deletion maximum size (DEFAULT 1)\n";
    print STDERR "\t-delMdl                   Deletion size model (DEFAULT 'linear', other option 'log')\n";
    print STDERR "\t-delLog                   Change Deletion size model to 'log' (DEFAULT FALSE)\n\n";

    print STDERR "\t-insFrq                   Insertion frequency (DEFAULT 0)\n";
    print STDERR "\t-insMin                   Insertion minimum size (DEFAULT 1)\n";
    print STDERR "\t-insMax                   Insertion maximum size (DEFAULT 1)\n";
    print STDERR "\t-insMdl                   Insertion size model (DEFAULT 'linear', other option 'log')\n";
    print STDERR "\t-insLog                   Change Insertion size model to 'log' (DEFAULT FALSE)\n\n";

    print STDERR "\t-invFrq                   Inversion frequency (DEFAULT 0)\n";
    print STDERR "\t-invMin                   Inversion minimum size (DEFAULT 1000)\n";
    print STDERR "\t-invMax                   Inversion maximum size (DEFAULT 10000)\n";
    print STDERR "\t-invMdl                   Inversion size model (DEFAULT 'linear', other option 'log')\n";
    print STDERR "\t-invLog                   Change Inversion size model to 'log' (DEFAULT FALSE)\n\n";

    print STDERR "\t-cnvFrq                   CNV frequency (DEFAULT 0)\n";
    print STDERR "\t-cnvMin                   CNV minimum size (DEFAULT 1000)\n";
    print STDERR "\t-cnvMax                   CNV maximum size (DEFAULT 20000)\n";
    print STDERR "\t-cnvMdl                   CNV size model (DEFAULT 'linear', other option 'log')\n";
    print STDERR "\t-cnvLog                   Change CNV size model to 'log' (DEFAULT FALSE)\n";
    print STDERR "\t-copyMax                  Max number of copies, use 1 as min with linear model (DEFAULT 3)\n\n";

    print STDERR "\t-write-cmap               Will write perfect coordinate maps between ref and transformed genome. " .
	"Warning: very slow and large files\n\n";
    print STDERR "\t-out-dir, -o              Output directory [REQUIRED]\n";

    print STDERR "\t-verbose, -v              Verbose flag\n";
    print STDERR "\t-help, -h                 Display this help message\n\n\n";
}

sub preprocessor
{
    if (defined($options{'help'})) {
	print STDERR "Usage: $0 -g refGenomeFile -o outputDir [options]\n";
	printHelp();
	exit(0);
    }

    ##
    ## Throw help if required options are not specified
    ##
    if (! defined($options{'fasFile'}) ||
	! defined($options{'outDir'})
	) {
	print STDERR "Usage: $0 -f fasFile -o outputDir [options]\n";
	print STDERR "Usage: Missing arguments! Exiting...\n";
	printHelp();
	exit(1);
    }

    ##
    ## Ensure out directory exist
    ##
    $options{'outDir'} =~ s/\/+$//;  ## Clean dir tail
    if (! -e $options{'outDir'}) {
	buildDir($options{'outDir'});
    }

    ##
    ## Ensure that there is at least 1 nploidy
    ##
    if ($options{'nploidy'} < 1) {
	errMssg("Can't have nploidy (-n) less than 1: $options{'nploidy'}\n");
    }

    ##
    ## Process all variant type frequencies, min/max sizes, size models
    ##
    foreach my $varType (@allowedTypes) {
	## Skip translocations and snp checks for now
	next if ($varType eq "trn");
	next if ($varType eq "snp");

	## Swap size models if required
	if ($options{$varType}{'log'}) { $options{$varType}{'mdl'} = "log"; }

	## Make sure size models are legit
	if (! defined($options{'legMdl'}{$options{$varType}{'mdl'}})) {
	    errMssg("$varType size model: $options{$varType}{'mdl'} not allowed!\n");
	}

	## Make sure all size values are legit
	if ($options{$varType}{'min'} < 1) {
	    mssg("\# $varType minimum size < 1: $options{$varType}{'min'}, setting to 1.");
	    $options{$varType}{'min'} = 1;
	}
	if ($options{$varType}{'max'} < $options{$varType}{'min'}) {
	    mssg("Maximum $varType size is less than minimum: $options{$varType}{'max'} < " .
		 "$options{$varType}{'min'} ,setting equal.");
	    $options{$varType}{'max'} = $options{$varType}{'min'}
	}

	## Build length model
	if ($options{$varType}{'mdl'} eq "log") {
	    my $range = $options{$varType}{'max'} - $options{$varType}{'min'} + 1;
	    for (my $size = 0; $size < $range; $size++) {
		my $numBins = exp(-1 * $size * $options{'log-width'}) * $range;
		if ($numBins < 1) { $numBins = 1; }
		for (my $bin = 0; $bin < $numBins; $bin++) {
		    push @{$lenMdl{$varType}}, $options{$varType}{'min'} + $size;
		}
	    }
	} elsif ($options{$varType}{'mdl'} eq "linear") {
	    for (my $size = $options{$varType}{'min'}; $size < $options{$varType}{'max'} + 1; $size++) {
		push @{$lenMdl{$varType}}, $size;
	    }
	} else {
	    errMssg("Impossible to get here!!!");
	}
    }

    ##
    ## Make sure each fasta file has a name
    ##
    my %chrNames = ();
    foreach my $name (@{$options{'chrName'}}) {
	$chrNames{$name} = 1;
    }

    for (my $ii = 0; $ii < @{$options{'fasFile'}}; $ii++) {
	if (! defined($options{'chrName'}[$ii])) {
	    my $jj = 1;
	    my $nextName = "chr$jj";
	    while(defined($chrNames{$nextName})) {
		$jj++;
		$nextName = "chr$jj";
	    }
	    mssg("\# No name provided for fasta index: $ii, using: $nextName\n");
	    push @{$options{'chrName'}}, $nextName;
	}
	## Assign the params defitions
	$params{$options{'chrName'}[$ii]}{'fasFile'} = $options{'fasFile'}[$ii];
	if (defined($options{'preFile'}[$ii])) {
	    $params{$options{'chrName'}[$ii]}{'preFile'} = $options{'preFile'}[$ii];
	}
	if (defined($options{'repFile'}[$ii])) {
	    $params{$options{'chrName'}[$ii]}{'repFile'} = $options{'repFile'}[$ii];
	}
	if (defined($options{'gapFile'}[$ii])) {
	    $params{$options{'chrName'}[$ii]}{'gapFile'} = $options{'gapFile'}[$ii];
	}
    }

    ##
    ## Write params
    ##
    mssg("\#\n\# *** Parameters ***\n\#\n");
    foreach my $param (sort keys %params) {
	if (ref($params{$param}) eq "SCALAR") {
	    mssg("\# Set $param to $params{$param}\n");
	} elsif (ref($params{$param}) eq "ARRAY") {
	    my $cnt = 0;
	    foreach my $val (@{$params{$param}}) {
		mssg("\# Set $param, $cnt to $val\n");
		$cnt++;
	    }
	} elsif (ref($params{$param}) eq "HASH") {
	    my $cnt = 0;
	    foreach my $key (sort keys %{$params{$param}}) {
		mssg("\# Set $param, $key to $params{$param}{$key}\n");
		$cnt++;
	    }
	} else {
	    if (defined($params{$param})) {
		mssg("\# Set $param to $params{$param}\n");
	    } else {
		mssg("\# No value defined for $param\n");
	    }
	}
    }
    mssg("\#\n");

    ##
    ## Write options
    ##
    mssg("\#\n\# *** Options ***\n\#\n");
    foreach my $param (sort keys %options) {
	if (ref($options{$param}) eq "SCALAR") {
	    mssg("\# Set $param to $options{$param}\n");
	} elsif (ref($options{$param}) eq "ARRAY") {
	    my $cnt = 0;
	    foreach my $val (@{$options{$param}}) {
		mssg("\# Set $param, $cnt to $val\n");
		$cnt++;
	    }
	} elsif (ref($options{$param}) eq "HASH") {
	    my $cnt = 0;
	    foreach my $key (sort keys %{$options{$param}}) {
		next if (! defined($options{$param}{$key}));
		mssg("\# Set $param, $key to $options{$param}{$key}\n");
		$cnt++;
	    }
	} else {
	    if (defined($options{$param})) {
		mssg("\# Set $param to $options{$param}\n");
	    } else {
		mssg("\# No value defined for $param\n");
	    }
	}
    }
    mssg("\#\n");

    mssg("\# Done with preprocessing\n\#\n");
}

##
## Sequence subroutines
##
sub checkInter {
    my ($dat, $pos) = @_;

    foreach my $row (@{$dat}) {
	my $start = $$row[0];
	my $end   = $$row[1];

	if ($start <= $pos && $pos <= $end) {
	    return(TRUE);
	}
    }
    return(FALSE);
}

sub swapSeq {
    my ($orgSeq) = @_;
    $options{'swap-size'} = 3;
    $options{'swap-rep'}  = 3;
    $options{'swap-min'}  = 10;

    my $len = length($orgSeq);

    ##
    ## First we need to make sure the sequence does not 
    ##  contain any N's or non-ACTG chars
    for (my $ii = 0; $ii < $len; $ii++) {
	if (substr($orgSeq, $ii, 1) !~ m/^[$options{'bases'}]$/) {
	    substr($orgSeq, $ii, 1, substr($options{'bases'}, rand(length($options{'bases'})), 1));
	}
    }

    ## Make sure insert is greater than min length
    ##  otherwise generate a random sequence
    if ($len <= $options{'swap-min'}) {
	my $newSeq = "";
	for (my $cnt = 0; $cnt < $len; $cnt++) {
	    $newSeq .= substr($options{'bases'}, rand(length($options{'bases'})), 1);
	}
	return($newSeq);
    }

    my $trys = $options{'swap-rep'} * $len;
    for (my $cnt = 0; $cnt < $trys; $cnt++) {
	my $pos1 = rand($len - $options{'swap-size'});
	my $pos2 = rand($len - $options{'swap-size'});

	my $seq1 = substr($orgSeq, $pos1, $options{'swap-size'});
	my $seq2 = substr($orgSeq, $pos2, $options{'swap-size'});

	#$orgSeq = substr($orgSeq, $pos1, $options{'swap-size'}, $seq2);
	#$orgSeq = substr($orgSeq, $pos2, $options{'swap-size'}, $seq1);
	substr($orgSeq, $pos1, $options{'swap-size'}, $seq2);
	substr($orgSeq, $pos2, $options{'swap-size'}, $seq1);
    }
    return($orgSeq);
}

##
## General subroutines
##
sub buildDir
{
    my ($dir) = @_;

    if (-e $dir) {
        my $clean_cmd = "rm -rf $dir/*";
	mssg("\# Cleaning '$dir'\n");
        system($clean_cmd);
    } else {
        my @dir_array = split(/\//, $dir);
        my $tmp_dir = "";
        foreach my $val (@dir_array) {
            next if ($val =~ /^(\s)*$/);

            if (! $tmp_dir) { 
		if ($dir =~ m/^\/.*$/) {
		    $tmp_dir = "/" . $val;
		} else {
		    $tmp_dir = $val;
		}
            } else {
		$tmp_dir .= "/" . $val; 
	    }

            if (! -e $tmp_dir) {
                my $mk_cmd = "mkdir $tmp_dir";
  
		mssg("\# Building '$tmp_dir'\n");
                system($mk_cmd);
            }
        }
    }
}

sub array2tsv
{
    my ($ref) = @_;

    if (@{$ref} == 0) { return(""); }

    my $str = "$$ref[0]";
    for (my $ii = 1; $ii < @{$ref}; $ii++) {
	$str .= "\t";
	if (defined($$ref[$ii])) { $str .= $$ref[$ii]; }
    }
    return($str);
}

sub writeStats
{
    my ($fh, $prefix) = @_;

    foreach my $key (sort keys %stats) {
	## Determine if element is an array
	if (ref($stats{$key}) eq 'ARRAY') {
	    for (my $idx = 0; $idx < @{$stats{$key}}; $idx++) {
		if (defined($prefix)) { print $fh $prefix; }
		print $fh "$key\t$idx\t$stats{$key}[$idx]\n";
	    }
	} else {
	    if (defined($prefix)) { print $fh $prefix; }
	    print $fh "$key\t$stats{$key}\n";
	}
    }
}

sub writeParams
{
    ## Build params file name
    $options{'params-file'}= "$options{'out-dir'}/$options{'output-prefix'}";
    for (my $nploid = 0; $nploid < $options{'n-ploid'}; $nploid++) {
	$options{'params-file'} .= "-" . $options{'n-ploid-names'}[$nploid];
    }
    $options{'params-file'} .= ".transformGenome-params.txt";

    mssg("\# Writing params: $options{'params-file'}\n");
    open(OUT, ">$options{'params-file'}")
	or die "Could not open $options{'out-dir'} for writing: $!";

    foreach my $key (sort keys %options) {
	## Determine if element is an array
	if (ref($options{$key}) eq 'ARRAY') {
	    for (my $idx = 0; $idx < @{$options{$key}}; $idx++) {
		print OUT "$key\t$idx\t$options{$key}[$idx]\n";
	    }
	} else {
	    print OUT "$key\t$options{$key}\n";
	}
    }
    writeStats(\*OUT);
    close(OUT);
}

sub mssg
{
    my ($line) = @_;

    if ($options{'verbose'}) {
	print STDERR "$line";
    }
}

sub warnMssg {
    my ($line) = @_;

    my @fields = split(/\n/, $line);
    foreach my $mssgStr (@fields) {
        print STDERR "[Warning]: $mssgStr\n";
    }
}

sub errMssg {
    my ($mssgStr, $exitStatus) = @_;

    if (! defined($exitStatus)) {
        $exitStatus = 1;
    }
    my @mssgArray = split(/\n/, $mssgStr);

    my $prgm = $0;
    $prgm =~ s/^.*\/([^\/]+)$/$1/;

    foreach my $str (@mssgArray) {
        print STDERR "[ERROR]: {$prgm}: $str\n";
    }

    if ($exitStatus ne "-1") {
        print STDERR "[ERROR]: {$prgm}: Exiting...\n";
        exit($exitStatus);
    }
}

## End of file
