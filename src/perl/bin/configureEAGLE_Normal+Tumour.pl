#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2014 Illumina, Inc.

This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
covered by the "BSD 2-Clause License" (see accompanying LICENSE file)

=head1 NAME

Enhanced Artificial Genome Engine (EAGLE)

=head1 VERSION

EAGLE vB<@EAGLE_VERSION@>

=head1 DESCRIPTION

B<configureEAGLE_Normal+Tumour.pl> - Create and initialize the EAGLE folder.

=head1 SYNOPSIS

  configureEAGLE_Normal+Tumour.pl [<path/to/EAGLE-dir>] \
                    --reference-genome=<path/to/ReferenceDir> \
                    --run-info=<path/to/RunInfo.xml> \
                    [options]

=head1 ARGUMENTS

=over 4

=item B<path/to/EAGLE-dir>    I<(=./EAGLE)>

Output path of the configuration

=back

=head1 OPTIONS

=over 4

=item B<--run-info=<path/to/RunInfo.xml>>, B<-i <path/to/RunInfo.xml>>

Path to the RunInfo.xml to be used as input (mandatory)

=item B<--reference-genome=<path/to/ReferenceDir>>, B<-r <path/to/ReferenceDir>>

Path to the Reference (whole dir)

=item B<--strelka-reference-genome=<path/to/ReferenceDir>>, B<-r <path/to/ReferenceDir>>

Path to the Reference (whole dir)

=item B<--isaac-sorted-reference=<path/to/ReferenceDir>>, B<-r <path/to/ReferenceDir>>

Path to the Reference (xml file)

=item B<--variant-list=<path/to/VariantList>>, B<-v </path/to/VariantList>>

Path to (*.vcf) files that describe the variations (multiple files allowed)

=item B<--sample-genome=<SAMPLE_NAME>>, B<-s <SAMPLE_NAME>>

Override default sample genome dir name with SAMPLE_NAME

=item B<--run-folder=<RUN_FOLDER_NAME>>, B<-r <RUN_FOLDER_NAME>>

Override default RunFolder name with RUN_FOLDER_NAME

=item B<--generate-run-id>

Generate the RunFolder name based on current date and other parameters

=item B<--quality-table=<PATH>>, B<-q <PATH>>    I<(=@EAGLE_FULL_DATADIR@/QualityTables/DefaultQualityTable.qval)>

Path to an alternative quality table (multiple files allowed)

=item B<--template-length-table=<PATH>>, B<-t <PATH>>    I<(=@EAGLE_FULL_DATADIR@/TemplateLengthTables/DefaultTemplateLengthTable.tsv)>

Path to an alternative template length table

=item B<--mismatch-table=<PATH>>, B<-m <PATH>>    I<(=@EAGLE_FULL_DATADIR@/MismatchTables/DefaultMismatchTable.tsv)>

Path to an alternative mismatch table

=item B<--coverage-depth=<NUM>>, B<-d <NUM>>    I<(=30)>

Coverage depth

=item B<--random-seed=<NUM>>

Random seed (should be > 1)

=item B<--workflow=<WKF>>, B<-w <WKF>>    I<(=default)>

Valid workflows are: {default}

=item B<--genome-mutator-options=<EXTRA_OPTIONS>>

Pass EXTRA_OPTIONS to the I<Genome Mutator> tool

=item B<--fragments-allocator-options=<EXTRA_OPTIONS>>

Pass EXTRA_OPTIONS to the I<Fragments Allocator> tool

=item B<--run-folder-generator-options=<EXTRA_OPTIONS>>

Pass EXTRA_OPTIONS to the I<Run Folder Generator> tool

=item B<--sequencer-simulator-options=<EXTRA_OPTIONS>>

Pass EXTRA_OPTIONS to the I<Sequencer Simulator> tool

=item B<--force>

Overwrite existing files (at configuration time)

=item B<--always-force>

Configure Makefiles to support multiple writings to the same targets
    (default behaviour is to require an explicit "make clean")

=item B<--help>, B<-h>, B<-?> 

Brief help message

=item B<--version>, B<-V> 

Version information

=item B<--man>

Full documentation

=back

=head1 NOTE

=head2 For more information, including examples, run:

  @EAGLE_FULL_BINDIR@/configureEAGLE_Normal+Tumour.pl --man

=cut

use lib '@EAGLE_FULL_PERL_LIBDIR@';

use warnings FATAL => 'all';
use strict;
use Cwd qw(abs_path cwd);
use POSIX qw(strftime);
use File::Path qw(mkpath);
use File::Spec;
use File::Spec::Functions qw(rel2abs);
use File::Glob;
use IO::File;
use Carp;

use Pod::Usage;
use Getopt::Long;
use Data::Dumper;

use RunFolder::RunInfo;
use RunFolder::Config;

my $DEBUG            = '@EAGLE_DEBUG_MODE@';
my $EAGLE_DATADIR    = '@EAGLE_FULL_DATADIR@';

=head1 EXAMPLES


=cut

my $programName = (File::Spec->splitpath(abs_path($0)))[2];
my $programPath = (File::Spec->splitpath(abs_path($0)))[1];


# load the default configuration parameters before parsing the command line
my $man             = 0;
my $help            = !(scalar @ARGV);
my $version         = 0;
my $outputDir       = File::Spec->catdir(Cwd::cwd(),'EAGLE');
my $runInfoFile     = undef;
my @referenceGenome = ();
my $strelkaReferenceGenome = undef;
my $isaacSortedReference = undef; # '/illumina/scratch/eagle/References/forISAAC/HumanChr21-6bit/chr21.fa-32mer-6bit-SortedReference.xml';
my @variantList     = ();
my $sampleGenome    = 'sample_genome';
my $runFolder       = undef;
my $generateRunId   = 0;
my $workflow        = 'default';
my $ploidyArgs      = '';
my $force           = 0;
my $alwaysForce     = 0;
my @qualityTable    = ();
my $templateLengthTable = File::Spec->catfile( $EAGLE_DATADIR, 'TemplateLengthTables', 'DefaultTemplateLengthTable.tsv');
my $mismatchTable    = File::Spec->catfile( $EAGLE_DATADIR, 'MismatchTables', 'DefaultMismatchTable.tsv');
my $coverageDepth   = undef;  # FragmentAllocator uses 30 by default, but I want to concentrate default values in just one place
my $randomSeed      = undef;
my $genomeMutatorOptions      = '';
my $fragmentsAllocatorOptions  = '';
my $runFolderGeneratorOptions   = '';
my $sequencerSimulatorOptions = '';

my @sharedVariants = ();
my @tumourOnlyVariants = ();
my $tumourOnlyCnv = '';
my $normalCoverage = 0;
my $tumourOverallCoverage = 0;
my $tumourPurity = 1.0;

my $useQsub = 0;
my $useSge = 0;


my %options = 
(
 'shared-variants'         => \@sharedVariants,
 'tumour-only-variants'    => \@tumourOnlyVariants,
 'tumour-only-cnv'         => \$tumourOnlyCnv,
 'normal-coverage'         => \$normalCoverage,
 'tumour-overall-coverage' => \$tumourOverallCoverage,
 'tumour-purity'           => \$tumourPurity,

    'use-qsub'                => \$useQsub,
    'use-sge'                => \$useSge,
    'run-info'                   => \$runInfoFile,
    'reference-genome'           => \@referenceGenome,
    'strelka-reference-genome'   => \$strelkaReferenceGenome,
    'isaac-sorted-reference'     => \$isaacSortedReference,
    'variant-list'               => \@variantList,
    'sample-genome'              => \$sampleGenome,
    'run-folder'                 => \$runFolder,
    'quality-table'              => \@qualityTable,
    'template-length-table'      => \$templateLengthTable,
    'mismatch-table'             => \$mismatchTable,
    'coverage-depth'             => \$coverageDepth,
    'random-seed'                => \$randomSeed,
    'generate-run-id'            => \$generateRunId,
    'workflow'                   => \$workflow,
    'genome-mutator-options'       => \$genomeMutatorOptions,
    'fragments-allocator-options'  => \$fragmentsAllocatorOptions,
    'run-folder-generator-options' => \$runFolderGeneratorOptions,
    'sequencer-simulator-options'  => \$sequencerSimulatorOptions,
    'force'                      => \$force,
    'always-force'               => \$alwaysForce,
    'help'                       => \$help,
    'version'                    => \$version,
    'man'                        => \$man
);

my $result = GetOptions(\%options,
 'shared-variants=s@',
 'tumour-only-variants=s@',
 'tumour-only-cnv=s',
 'normal-coverage=f',
 'tumour-overall-coverage=f',
 'tumour-purity=f',
 'use-qsub!',
 'use-sge!',
                        'run-info|i=s',
                        'reference-genome|r=s@', 'strelka-reference-genome=s', 'isaac-sorted-reference=s', 'variant-list|v=s@', 'sample-genome|s=s', 'run-folder|r=s',
                        'quality-table|q=s@', 'template-length-table|t=s', 'mismatch-table|m=s', 'coverage-depth|d=i', 'random-seed=i', 'workflow|w=s',
                        'genome-mutator-options=s', 'fragments-allocator-options=s', 'run-folder-generator-options=s', 'sequencer-simulator-options=s',
                        'generate-run-id!', 'force!', 'always-force!',  # 'make!',
                        'help|h|?!', 'version|V!', 'man!')
             or pod2usage(2);

pod2usage(-exitstatus => 0, -verbose => 99, -sections => 'NAME|VERSION|DESCRIPTION|SYNOPSIS|ARGUMENTS|OPTIONS|NOTE') if $help;
pod2usage(-exitstatus => 0, -verbose => 99, -sections => 'NAME|VERSION')                                             if $version;
pod2usage(-exitstatus => 0, -verbose => 2,  -input => $1)  if ($man and $0 =~ /(.*)/);
pod2usage("$0: Too many positional arguments.\n")          if (1 < @ARGV);
pod2usage("$0: reference-genome needs to be specified.\n") if (@referenceGenome == 0);
pod2usage("$0: strelka-reference-genome needs to be specified.\n") if (! defined $strelkaReferenceGenome);
pod2usage("$0: isaac-sorted-reference needs to be specified.\n")   if (! defined $isaacSortedReference);

$outputDir = $ARGV[0]  if @ARGV;
my $fullOutputDir = rel2abs($outputDir);

my @fullSharedVariants = map { rel2abs $_ } @sharedVariants;
my @fullTumourOnlyVariants = map { rel2abs $_ } @tumourOnlyVariants;
my $fullTumourOnlyCnv = rel2abs( $tumourOnlyCnv );
my $fullTemplateLengthTable = rel2abs( $templateLengthTable );


if (! -e $fullOutputDir)
{
    print "Creating directory '$fullOutputDir'\n";
    File::Path::mkpath( $fullOutputDir );
    croak "ERROR: *** failed to create directory '$fullOutputDir' ***"  unless -d $fullOutputDir;
}
else {
#    croak "ERROR: *** already existing directory '$fullOutputDir' ***";
}

if (! -f $fullTemplateLengthTable || ! -r $fullTemplateLengthTable) { print "ERROR: File $fullTemplateLengthTable is not readable\n"; exit -1; }


my $fullMakefile = File::Spec->catfile( $fullOutputDir, "Makefile" );
open my $makefileHandle, '>', "$fullMakefile"          or croak("ERROR: *** failed to open '$fullMakefile' for writing: $! ***\n   ");

#my $referenceGenome = '/illumina/scratch/eagle/References/iGenomes_hg19_with_fai/chr21.fa';

my $commonOptions0 = '--run-folder=RunFolder';
if (defined $runInfoFile) {
  $commonOptions0 .= ' --run-info=' . ${runInfoFile};
}

$commonOptions0 .= ' --template-length-table=' . $fullTemplateLengthTable;

my $commonOptionsPart1 = $commonOptions0;
foreach my $refGenome (@referenceGenome) {
  $commonOptionsPart1 .= ' --reference-genome=' . rel2abs( ${refGenome} );
}
#$commonOptionsPart1 .= ' --genome-mutator-options="--organism-ploidy=2 --ploidy-chromosome=X --ploidy-level=1 --ploidy-chromosome=Y --ploidy-level=1 --ploidy-chromosome=MT --ploidy-level=100';
$commonOptionsPart1 .= ' --genome-mutator-options="--organism-ploidy=2 --ploidy-chromosome=chrX --ploidy-level=1 --ploidy-chromosome=chrY --ploidy-level=1 --ploidy-chromosome=chrM --ploidy-level=100 --no-translocation-error';
my $commonOptionsPart2 = '"';
my $commonOptions = $commonOptionsPart1 . $commonOptionsPart2;

# Add when it goes faster, use: --reference-genome=/illumina/scratch/eagle/RealisticInputs/hs37d5.fa --reference-genome=/illumina/scratch/eagle/RealisticInputs/EagleDecoys.fa --variant-list=/illumina/scratch/eagle/RealisticInputs/dnSnp_set1.vcf

my $tumourCoverageForRequestedPurity = ${tumourOverallCoverage} * ${tumourPurity};
my $normalCoverageForRequestedPurity = ${tumourOverallCoverage} * (1-${tumourPurity});
my $normalDir = "$fullOutputDir/EAGLE_normal";
my $tumourDir = "$fullOutputDir/EAGLE_tumour_purity_" . ${tumourPurity};

#my $EAGLE_LIBEXEC = "/home/ljanin/install/EAGLE-release/libexec/EAGLE";
my $EAGLE_LIBEXEC = ${programPath} . "/../libexec/EAGLE";

my $makeOptions = '-j 8';

my $prefixedSharedVariants = '';
foreach my $sharedVariant (@fullSharedVariants) {
  $prefixedSharedVariants .= ' --variant-list=' . ${sharedVariant};
}

my $prefixedTumourOnlyVariants = '';
foreach my $tumourOnlyVariant (@fullTumourOnlyVariants) {
  $prefixedTumourOnlyVariants .= ' --variant-list=' . ${tumourOnlyVariant};
}

my $QSUB_PREFIX_GENOME_MUTATOR = "";
my $QSUB_PREFIX_ISAAC = "";
my $QSUB_SUFFIX = "";
if (${useQsub}) {
  $QSUB_PREFIX_GENOME_MUTATOR = "qsub -N \"EAGLE-GenomeMutator\" -cwd -sync y -@ /illumina/gridware/SGE_6.1U5/req/automation -b yes \"";
  $QSUB_PREFIX_ISAAC          = "qsub -N \"EAGLE-iSAAC\"         -cwd -sync y -@ /illumina/gridware/SGE_6.1U5/req/automation -b yes \"";
  $QSUB_SUFFIX = "\"";
}
my $SGE_STRING = "";
if (${useSge}) {
  $SGE_STRING="sge";
}

print $makefileHandle "AND = &&\n\n";
print $makefileHandle ".PHONY: all subtraction normal_prep tumour_prep tumour_normal_prep tumout_tumour_prep tumour_normal_sim tumour_tumour_sim tumour_subs tumour_merged tumour normal_with_sample_genome_link normal align_normal align_tumour strelka\n\n";

print $makefileHandle "builds: normal tumour\n";
print $makefileHandle "all: subtraction\n";
print $makefileHandle "subtraction: strelka\n";

print $makefileHandle "\n# Prepare tumour at " . ${tumourPurity}*100 . "% purity\n";
print $makefileHandle "tumour_prep: ${tumourDir}/Makefile\n";
print $makefileHandle "${tumourDir}/Makefile:\n";
print $makefileHandle "\t${programPath}/configureEAGLE.pl ${commonOptions} --coverage-depth=" . ${tumourOverallCoverage}/2 . " ${tumourDir}\n";

print $makefileHandle "\n# Prepare tumour.EAGLE_normalForPurityMix\n";
print $makefileHandle "tumour_normal_prep: ${tumourDir}/EAGLE_normalForPurityMix/Makefile\n";
print $makefileHandle "${tumourDir}/EAGLE_normalForPurityMix/Makefile: ${tumourDir}/Makefile\n";
print $makefileHandle "\t${programPath}/configureEAGLE.pl ${commonOptionsPart1} --prefix=normal_ ${commonOptionsPart2} ${prefixedSharedVariants} --coverage-depth=" . ${normalCoverageForRequestedPurity}/2 . " ${tumourDir}/EAGLE_normalForPurityMix\n";

if ("${tumourOnlyCnv}" eq "") {
  print $makefileHandle "\n# Prepare tumour.EAGLE_tumourForPurityMix\n";
  print $makefileHandle "tumour_tumour_prep: ${tumourDir}/EAGLE_tumourForPurityMix/Makefile\n";
  print $makefileHandle "${tumourDir}/EAGLE_tumourForPurityMix/Makefile: ${tumourDir}/Makefile\n";
  print $makefileHandle "\t${programPath}/configureEAGLE.pl ${commonOptionsPart1} --prefix=tumour_ ${commonOptionsPart2} ${prefixedSharedVariants} ${prefixedTumourOnlyVariants} --coverage-depth=" . ${tumourCoverageForRequestedPurity}/2 . " ${tumourDir}/EAGLE_tumourForPurityMix\n";
}
else {
  print $makefileHandle "\n# Prepare tumour.EAGLE_tumourForPurityMixBeforeCnv\n";
  print $makefileHandle "${programPath}/configureEAGLE.pl ${commonOptions} ${prefixedSharedVariants} ${prefixedTumourOnlyVariants} --coverage-depth=" . ${tumourCoverageForRequestedPurity}/2 . " ${tumourDir}/EAGLE_tumourForPurityMixBeforeCnv\n";
  print $makefileHandle "cd ${tumourDir}/EAGLE_tumourForPurityMixBeforeCnv && \${MAKE} $makeOptions sample\n";

  print $makefileHandle "\n# Prepare tumour.EAGLE_tumourForPurityMix\n";
  print $makefileHandle "${EAGLE_LIBEXEC}/applyCopyNumber.pl --vcf=${tumourDir}/CNV.vcf --copy-number=" . rel2abs(${fullTumourOnlyCnv}) . "\n";
  print $makefileHandle "${programPath}/configureEAGLE.pl ${commonOptions0} --genome-mutator-options=\"--prefix=tumour_\" --reference-genome=${tumourDir}/EAGLE_tumourForPurityMixBeforeCnv/sample_genome --variant-list=${tumourDir}/CNV.vcf --coverage-depth=" . ${tumourCoverageForRequestedPurity}/2 . " ${tumourDir}/EAGLE_tumourForPurityMix\n";
}

print $makefileHandle "\n# Simulate tumour purity sub-datasets\n";
print $makefileHandle "tumour_normal_sim: ${tumourDir}/EAGLE_normalForPurityMix/fragments/fragments.pos\n";
print $makefileHandle "${tumourDir}/EAGLE_normalForPurityMix/fragments/fragments.pos: ${tumourDir}/EAGLE_normalForPurityMix/Makefile\n";
#print $makefileHandle "\t\$(MAKE) -C ${tumourDir}/EAGLE_normalForPurityMix fragments\n";
print $makefileHandle "\tcd ${tumourDir}/EAGLE_normalForPurityMix \\\n";
  print $makefileHandle "\t\$(AND) time ${QSUB_PREFIX_GENOME_MUTATOR}\${MAKE} -C ${tumourDir}/EAGLE_normalForPurityMix fragments${QSUB_SUFFIX}\n";

print $makefileHandle "\n";
print $makefileHandle "tumour_tumour_sim: ${tumourDir}/EAGLE_tumourForPurityMix/fragments/fragments.pos\n";
print $makefileHandle "${tumourDir}/EAGLE_tumourForPurityMix/fragments/fragments.pos: ${tumourDir}/EAGLE_tumourForPurityMix/Makefile\n";
#print $makefileHandle "\t\$(MAKE) -C ${tumourDir}/EAGLE_tumourForPurityMix fragments\n";
print $makefileHandle "\tcd ${tumourDir}/EAGLE_tumourForPurityMix \\\n";
print $makefileHandle "\t\$(AND) time ${QSUB_PREFIX_GENOME_MUTATOR}\${MAKE} -C ${tumourDir}/EAGLE_tumourForPurityMix fragments${QSUB_SUFFIX}\n";

print $makefileHandle "\n# Merge tumour purity sub-datasets to ${tumourDir}\n";
print $makefileHandle "tumour_merged: ${tumourDir}/fragments/fragments.pos\n";
print $makefileHandle "${tumourDir}/fragments/fragments.pos: ${tumourDir}/EAGLE_normalForPurityMix/fragments/fragments.pos ${tumourDir}/EAGLE_tumourForPurityMix/fragments/fragments.pos\n";
print $makefileHandle "\tcd ${tumourDir} \\\n";
print $makefileHandle "\t\$(AND) $EAGLE_LIBEXEC/mergeSampleGenomes.pl -i EAGLE_normalForPurityMix -j EAGLE_tumourForPurityMix \\\n";
print $makefileHandle "\t\$(AND) $EAGLE_LIBEXEC/mergeFragments.pl -i EAGLE_normalForPurityMix -j EAGLE_tumourForPurityMix -a \"`grep CHROMOSOME_ALLELES Makefile | cut -d ' ' -f 3-`\"\n";

print $makefileHandle "\n# Prepare normal\n";
print $makefileHandle "normal_prep: ${normalDir}/Makefile\n";
print $makefileHandle "${normalDir}/Makefile: ${tumourDir}/EAGLE_normalForPurityMix/fragments/fragments.pos\n";
print $makefileHandle "\t${programPath}/configureEAGLE.pl ${commonOptions} ${prefixedSharedVariants} --genome-mutator-options=\"--prefix=normal_\" --coverage-depth=" . ${normalCoverage}/2 . " ${normalDir} \\\n";
print $makefileHandle "\t\$(AND) rm -rf ${normalDir}/reference_genome \\\n";
print $makefileHandle "\t\$(AND) ln -s ${tumourDir}/EAGLE_normalForPurityMix/reference_genome ${normalDir}/reference_genome \\\n";
print $makefileHandle "\t\$(AND) ln -s ${tumourDir}/EAGLE_normalForPurityMix/sample_genome ${normalDir}/sample_genome\n";

print $makefileHandle "\n# Simulate tumour\n";
print $makefileHandle "tumour: ${tumourDir}/RunFolder/RunInfo.xml\n";
print $makefileHandle "${tumourDir}/RunFolder/RunInfo.xml: ${tumourDir}/fragments/fragments.pos\n";
print $makefileHandle "\t\$(MAKE) -C ${tumourDir} ${SGE_STRING}\n";

print $makefileHandle "\n# Simulate normal\n";
print $makefileHandle "normal: ${normalDir}/RunFolder/RunInfo.xml\n";
print $makefileHandle "${normalDir}/RunFolder/RunInfo.xml: ${normalDir}/Makefile\n";
print $makefileHandle "\t\$(MAKE) -C ${normalDir} ${SGE_STRING}\n";

#print $makefileHandle "\n# Align and call variants for normal\n";
print $makefileHandle "\n# Align normal\n";
print $makefileHandle "align_normal: ${normalDir}/Aligned/default/default/sorted.bam\n";
print $makefileHandle "${normalDir}/Aligned/default/default/sorted.bam: ${normalDir}/RunFolder/RunInfo.xml\n";
print $makefileHandle "\tcd ${normalDir} \\\n";
print $makefileHandle "\t\$(AND) time ${QSUB_PREFIX_ISAAC}/illumina/development/iSAAC/latest/bin/isaac-align -r ${isaacSortedReference} -b ${normalDir}/RunFolder/Data/Intensities/BaseCalls/ --stop-at=Bam -m 46${QSUB_SUFFIX}\n";

#print $makefileHandle "\n# Align and call variants for tumour\n";
print $makefileHandle "\n# Align tumour\n";
print $makefileHandle "align_tumour: ${tumourDir}/Aligned/default/default/sorted.bam\n";
print $makefileHandle "${tumourDir}/Aligned/default/default/sorted.bam: ${tumourDir}/RunFolder/RunInfo.xml\n";
print $makefileHandle "\tcd ${tumourDir} && time ${QSUB_PREFIX_ISAAC}/illumina/development/iSAAC/latest/bin/isaac-align -r ${isaacSortedReference} -b RunFolder/Data/Intensities/BaseCalls/ --stop-at=Bam -m 46${QSUB_SUFFIX}\n";

print $makefileHandle "\n# Run Strelka\n";
print $makefileHandle "STRELKA_INSTALL_DIR=/home/csaunders/opt/x86_64-linux/strelka_workflow-0.4.3\n";
print $makefileHandle "strelka: ${fullOutputDir}/StrelkaComparison/task.complete\n";
print $makefileHandle "${fullOutputDir}/StrelkaComparison/task.complete: ${normalDir}/Aligned/default/default/sorted.bam ${tumourDir}/Aligned/default/default/sorted.bam\n";
print $makefileHandle "\t\${STRELKA_INSTALL_DIR}/configureStrelkaWorkflow.pl --normal=${normalDir}/Aligned/default/default/sorted.bam --tumor=${tumourDir}/Aligned/default/default/sorted.bam --ref=${strelkaReferenceGenome} --config=\${STRELKA_INSTALL_DIR}/strelka/etc/strelka_config_eland_default.ini --output-dir=${fullOutputDir}/StrelkaComparison \\\n";
print $makefileHandle "\t\$(AND) \$(MAKE) -C ${fullOutputDir}/StrelkaComparison\n";



print "Successfully configured '$fullOutputDir'. Type 'make' inside this directory to kick-off simulation.\n\n";


__END__

=head1 DIAGNOSTICS

=head2 Exit status

=over 4

=item B<0:> successful completion

=item B<1:> abnormal completion

=item B<2:> fatal error

Z<>

=item B<Errors:> All error messages are prefixed with "ERROR: ".

Z<>

=item B<Warnings:> All warning messages generated by EAGLE are prefixed with "WARNING: ".

Z<>

=back

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=head1 AUTHOR

Lilian Janin

=cut

