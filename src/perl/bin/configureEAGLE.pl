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

B<configureEAGLE.pl> - Create and initialize the EAGLE folder.

=head1 SYNOPSIS

  configureEAGLE.pl [<path/to/EAGLE-dir>] \
                    --reference-genome=<path/to/ReferenceDir> \
                    --variant-list=<path/to/Variants.vcf> \
                    --run-info=<path/to/RunInfo.xml> \
                    [options]

=head1 ARGUMENTS

=over 4

=item B<path/to/EAGLE-dir>    I<(=./EAGLE)>

Output path of the configuration

=back

=head1 OPTIONS

=over 4

=item B<--run-info=<path/to/RunInfo.xml>>, B<-i <path/to/RunInfo.xml>>    I<(=@EAGLE_FULL_DATADIR@/RunInfo/RunInfo_PairedReadsBarcode8x32Tiles.xml)>

Path to the RunInfo.xml to be used as input (mandatory)

=item B<--reference-genome=<path/to/ReferenceDir>>, B<-g <path/to/ReferenceDir>>

Path to the Reference (whole dir)

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

=item B<--qq-table=<PATH>>    I<(=@EAGLE_FULL_DATADIR@/QualityTables/DefaultQQTable.tsv)>

Path to an alternative QQ table

=item B<--template-length-table=<PATH>>, B<-t <PATH>>    I<(=@EAGLE_FULL_DATADIR@/TemplateLengthTables/DefaultTemplateLengthTable.tsv)>

Path to an alternative template length table

=item B<--mismatch-table=<PATH>>, B<-m <PATH>>    I<(=@EAGLE_FULL_DATADIR@/MismatchTables/DefaultMismatchTable.tsv)>

Path to an alternative mismatch table

=item B<--homopolymer-indel-table=<PATH>>    I<(=@EAGLE_FULL_DATADIR@/MismatchTables/DefaultHomopolymerIndelTable.tsv)>

Path to an alternative homopolymer indel table

=item B<--motif-quality-drop-table=<PATH>>    I<(=@EAGLE_FULL_DATADIR@/MotifQualityDropTables/ZeroMotifQualityDropTable.tsv)>

Path to an alternative motif quality drop table

=item B<--gc-coverage-fit-table=<PATH>>

If specified, calculate GC content to affect coverage. Should be a tab-separated file with 2 columns: gc% and coverage multiplier. Intermediate values are interpolated linearly.

=item B<--coverage-depth=<NUM>>    F<(=30)>

Coverage depth (deprecated, replaced by --allele-coverage-depth)

=item B<--allele-coverage-depth=<NUM>>, B<-d <NUM>>    F<(=30)>

Coverage depth for each allele

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

  @EAGLE_FULL_BINDIR@/configureEAGLE.pl --man

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
use File::Basename;
use List::Util qw(max);
use IO::File;
use Carp;

use Pod::Usage;
use Getopt::Long;
use Data::Dumper;

use RunFolder::RunInfo;
use RunFolder::Config;

my $DEBUG            = '@EAGLE_DEBUG_MODE@';
my $EAGLE_DATADIR    = '@EAGLE_FULL_DATADIR@';
#my $EAGLE_ETCDIR     = '@EAGLE_FULL_ETCDIR@';

=head1 EXAMPLES

=head2 Running end-to-end generation of E-coli data with no variants:

 @EAGLE_FULL_BINDIR@/configureEAGLE.pl ./EagleEcoli \
	          --reference-genome=@EAGLE_FULL_DATADIR@/Genomes/E_coli.fa \
	          --variant-list=@EAGLE_FULL_DATADIR@/Variants/None.vcf \
	          --run-info=@EAGLE_FULL_DATADIR@/RunInfo/RunInfo_SingleRead1x1Tiles.xml

 make -C EagleEcoli # use make -j <parallel jobs> to speedup the conversion

=cut


# load the default configuration parameters before parsing the command line
my $originalPath    = Cwd::cwd();
my $originalCmdLine = join " ", $0, @ARGV;
my $man             = 0;
my $help            = !(scalar @ARGV);
my $version         = 0;
my $outputDir       = File::Spec->catdir(Cwd::cwd(),'EAGLE');
my $runInfoFile     = File::Spec->catfile( $EAGLE_DATADIR, 'RunInfo', 'RunInfo_PairedReadsBarcode8x32Tiles.xml');
my @referenceGenome = ();
my @variantList     = ();
my $sampleGenome    = 'sample_genome';
my $runFolder       = undef;
my $generateRunId   = 0;
#my $make           = 0;
my $workflow        = 'default';
my $ploidyArgs      = '';
my $force           = 0;
my $alwaysForce     = 0;
my @qualityTable    = ();
my $qqTable         = File::Spec->catfile( $EAGLE_DATADIR, 'QualityTables', 'DefaultQQTable.tsv');
my $templateLengthTable = File::Spec->catfile( $EAGLE_DATADIR, 'TemplateLengthTables', 'DefaultTemplateLengthTable.tsv');
my $mismatchTable   = File::Spec->catfile( $EAGLE_DATADIR, 'MismatchTables', 'DefaultMismatchTable.tsv');
my $homopolymerIndelTable = File::Spec->catfile( $EAGLE_DATADIR, 'MismatchTables', 'DefaultHomopolymerIndelTable.tsv');
my $motifQualityDropTable = File::Spec->catfile( $EAGLE_DATADIR, 'MotifQualityDropTables', 'ZeroMotifQualityDropTable.tsv');
my $gcCoverageFitTable = undef;
my $coverageDepth   = 30;  # FragmentAllocator uses 30 by default, but I want to concentrate default values in just one place
my $alleleCoverageDepth   = undef;
my $randomSeed      = undef;
my @genomeMutatorOptionsArray = ();
my $genomeMutatorOptions      = '';
my $fragmentsAllocatorOptions = '';
my $runFolderGeneratorOptions = '';
my $sequencerSimulatorOptions = '';
my @errorModelOptions = ();


my %options = 
(
    'run-info'                   => \$runInfoFile,
    'reference-genome'           => \@referenceGenome,
    'variant-list'               => \@variantList,
    'sample-genome'              => \$sampleGenome,
    'run-folder'                 => \$runFolder,
    'quality-table'              => \@qualityTable,
    'qq-table'                   => \$qqTable,
    'template-length-table'      => \$templateLengthTable,
    'mismatch-table'             => \$mismatchTable,
    'homopolymer-indel-table'    => \$homopolymerIndelTable,
    'motif-quality-drop-table'   => \$motifQualityDropTable,
    'gc-coverage-fit-table'      => \$gcCoverageFitTable,
    'error-model-options'        => \@errorModelOptions,
    'coverage-depth'             => \$coverageDepth,
    'allele-coverage-depth'      => \$alleleCoverageDepth,
    'random-seed'                => \$randomSeed,
    'generate-run-id'            => \$generateRunId,
#    'make'                       => \$make,
    'workflow'                   => \$workflow,
    'genome-mutator-options'       => \@genomeMutatorOptionsArray,
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
                        'run-info|i=s',
                        'reference-genome|g=s@', 'variant-list|v=s@', 'sample-genome|s=s', 'run-folder|r=s', 'error-model-options=s@',
                        'quality-table|q=s@', 'qq-table=s', 'template-length-table|t=s', 'mismatch-table|m=s', 'homopolymer-indel-table=s',
                        'motif-quality-drop-table=s', 'gc-coverage-fit-table=s', 'coverage-depth|d=f', 'random-seed=i', 'workflow|w=s',
                        'genome-mutator-options=s@', 'fragments-allocator-options=s', 'run-folder-generator-options=s', 'sequencer-simulator-options=s',
                        'generate-run-id!', 'force!', 'always-force!',  # 'make!',
                        'help|h|?!', 'version|V!', 'man!')
             or pod2usage(2);

pod2usage(-exitstatus => 0, -verbose => 99, -sections => 'NAME|VERSION|DESCRIPTION|SYNOPSIS|ARGUMENTS|OPTIONS|NOTE') if $help;
pod2usage(-exitstatus => 0, -verbose => 99, -sections => 'NAME|VERSION')                                             if $version;
pod2usage(-exitstatus => 0, -verbose => 2,  -input => $1)  if ($man and $0 =~ /(.*)/);
pod2usage("$0: Too many positional arguments.\n")          if (1 < @ARGV);

if (defined $coverageDepth && !defined $alleleCoverageDepth) {
  $alleleCoverageDepth = $coverageDepth;
}
croak "ERROR: *** At least one Reference file or dir is needed. Please use --reference-genome ***\n   "  unless @referenceGenome;
croak "ERROR: *** Path to RunInfo.xml not given. The --run-info switch is not optional! ***\n   "        unless defined $runInfoFile;
croak "ERROR: *** Allele coverage depth must be specified ***\n   "                                             unless defined $alleleCoverageDepth;
print "WARNING: *** No variant list specified ***\n   "                                                  unless @variantList;

$genomeMutatorOptions = join( ' ', @genomeMutatorOptionsArray );

$outputDir = $ARGV[0]  if @ARGV;
my $fullOutputDir = rel2abs($outputDir);

my @fullReferenceGenome = map { rel2abs $_ } @referenceGenome;
my @fullVariantList = map { rel2abs $_ } @variantList;
my $fullQQTable = rel2abs( $qqTable );
my $fullTemplateLengthTable = rel2abs( $templateLengthTable );
my $fullMismatchTable = rel2abs( $mismatchTable );
my $fullHomopolymerIndelTable = rel2abs( $homopolymerIndelTable );
my $fullMotifQualityDropTable = rel2abs( $motifQualityDropTable );
my $fullGcCoverageFitTable = rel2abs( $gcCoverageFitTable ) if defined $gcCoverageFitTable;

# Checking that all input files exist
if (! -f $fullQQTable || ! -r $fullQQTable) { print "ERROR: File $fullQQTable is not readable\n"; exit -1; }
if (! -f $fullTemplateLengthTable || ! -r $fullTemplateLengthTable) { print "ERROR: File $fullTemplateLengthTable is not readable\n"; exit -1; }
if (! -f $fullMismatchTable || ! -r $fullMismatchTable) { print "ERROR: File $fullMismatchTable is not readable\n"; exit -1; }
if (! -f $fullHomopolymerIndelTable || ! -r $fullHomopolymerIndelTable) { print "ERROR: File $fullHomopolymerIndelTable is not readable\n"; exit -1; }
if (! -f $fullMotifQualityDropTable || ! -r $fullMotifQualityDropTable) { print "ERROR: File $fullMotifQualityDropTable is not readable\n"; exit -1; }
if (defined $gcCoverageFitTable && (! -f $fullGcCoverageFitTable || ! -r $fullGcCoverageFitTable)) { print "ERROR: File $fullGcCoverageFitTable is not readable\n"; exit -1; }
foreach my $filename (@fullReferenceGenome) {
  if (! -r $filename) { print "ERROR: File $filename is not readable\n"; exit -1; }
}
my $totalVariantFilesLength = 0;
foreach my $filename (@fullVariantList) {
  if (! -f $filename || ! -r $filename) { print "ERROR: File $filename is not readable\n"; exit -1; }
  $totalVariantFilesLength += -s "$filename";
}


my @fullErrorModelOptions = ();
foreach my $errorModelOption (@errorModelOptions)
{
  my @filenameMatches = ($errorModelOption =~ m/[fF]ile=([^:]*)/g);
  foreach my $rel (@filenameMatches)
  {
    (-e "$rel") || die "Non-existing file: $rel";
    my $abs = rel2abs($rel);
    $errorModelOption =~ s/([fFile]=)$rel/$1$abs/;
  }
  push( @fullErrorModelOptions, $errorModelOption );
}

my $referenceMode;
if (1 == @fullReferenceGenome && -d $fullReferenceGenome[0])
{
    $referenceMode = 'whole-genome';
} else {
    foreach (@fullReferenceGenome)
    {
        croak "ERROR: *** Could not use $_ as reference. Expected a file. ***\n   "   unless -f $_;
    }
    $referenceMode = 'reference-genome';
}

my $runInfo = RunFolder::RunInfo->new();
{
  open my $runInfoHandle, '<', "$runInfoFile"  or croak("ERROR: *** failed to open '$runInfoFile' for reading: $! ***\n   ");
  $runInfo->load($runInfoHandle)               or croak("ERROR: *** could not load '$runInfoFile' in memory ***\n   ");
  close $runInfoHandle                         or croak("ERROR: *** failed to close '$runInfoFile': $! ***\n   ");
}

if (@qualityTable == 0)
{
    my $readNum = 1;
    foreach my $i (sort $runInfo->reads)
    {
        if ((defined $runInfo->read($i)->{IsIndexedRead} && "Y" eq $runInfo->read($i)->{IsIndexedRead})
            || exists $runInfo->read($i)->{Index})
        {
            my $filename = File::Spec->catfile( $EAGLE_DATADIR, "QualityTables", "DefaultQualityTable.barcode.length" . $runInfo->read($i)->{NumCycles} . ".qval");
            (-e ${filename}) or die "ERROR: Cannot find quality table for the requested number of cycles. You may need to manually use the `scaleQualityTable` tool to generate one (${filename})";
            push( @qualityTable, ${filename} );
        }
        else
        {
            my $filename = File::Spec->catfile( $EAGLE_DATADIR, "QualityTables", "DefaultQualityTable.read${readNum}.length" . $runInfo->read($i)->{NumCycles} . ".qval");
            (-e ${filename}) or die "ERROR: Cannot find quality table for the requested number of cycles. You may need to manually use the `scaleQualityTable` tool to generate one (${filename})";
            push( @qualityTable, ${filename} );
            $readNum++;
        }
    }
}
my @fullQualityTable = map { rel2abs $_ } @qualityTable;

my $basesPerCluster = 0;
foreach my $i ($runInfo->reads)
{
    $basesPerCluster += $runInfo->read($i)->{NumCycles}
                        unless (defined $runInfo->read($i)->{IsIndexedRead} && "Y" eq $runInfo->read($i)->{IsIndexedRead})
                            || exists $runInfo->read($i)->{Index};
}
print "... total number of bases: $basesPerCluster\n";

if ($generateRunId)
{
    my $date = strftime "%y%m%d", localtime;
    my @existingRunFolders = <${fullOutputDir}/${date}_*>;
    $runInfo->date($date);
    $runInfo->number(1 + @existingRunFolders);
}
$runFolder = $runInfo->runId()  unless defined $runFolder;

if (! -e $fullOutputDir)
{
    print "Creating directory '$fullOutputDir'\n";
    File::Path::mkpath( $fullOutputDir );
    croak "ERROR: *** failed to create directory '$fullOutputDir' ***"  unless -d $fullOutputDir;
}
else {
    croak "ERROR: *** already existing directory '$fullOutputDir' ***";
}

{
  #my $fullRunInfo = File::Spec->catfile( $fullOutputDir, $RunFolder::RunInfo::RunInfoXml );
  my $fullRunInfo = File::Spec->catfile( $fullOutputDir, "${runFolder}.xml" );
  open my $runInfoHandle, '>', "$fullRunInfo"  or croak("ERROR: *** failed to open '$fullRunInfo' for writing: $! ***\n   ");
  $runInfo->save($runInfoHandle)               or croak("ERROR: *** could not save '$fullRunInfo' in disk ***\n   ");
  close $runInfoHandle                         or croak("ERROR: *** failed to close '$fullRunInfo': $! ***\n   ");
}

my $layout = $runInfo->layout;
croak "ERROR: *** could not retrieve element 'FlowcellLayout' from ${runInfoFile} ***\n   "  unless defined $layout;
croak "ERROR: *** could not retrieve property 'SurfaceCount' from ${runInfoFile} ***\n   "  unless defined $layout->{SurfaceCount};


# detect samtools
my $samtools = `which samtools`;
chomp $samtools;
(-f $samtools) or die "ERROR: samtools needs to be accessible in your PATH";


# creating reference_genome directory
my @bamChromosomes = ();
my @chromosomeAlleles = ();
my $genomeSize = 0;
my $maxChromosomeLength = 0;
my $refGenomeDir = "${fullOutputDir}/reference_genome";
File::Path::mkpath( $refGenomeDir );
open GENOMESIZE_XML, ">${refGenomeDir}/genome_size.xml" or die "Can't open ${refGenomeDir}/genome_size.xml for writing";
print GENOMESIZE_XML "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n<sequenceSizes>\n";
foreach my $genomeEntry (@fullReferenceGenome) {
  my @fastaFiles = ();
  if (-d "$genomeEntry") {
    @fastaFiles = glob "${genomeEntry}/*.fa";
  }
  else {
    push( @fastaFiles, ${genomeEntry} );
  }

  foreach my $fastaFile ( @fastaFiles ) {
    my $baseName = basename( ${fastaFile} );

    # Link fasta file
    #print "Adding ${fastaFile} to reference genome\n";
    system( "ln -s \"${fastaFile}\" \"${refGenomeDir}\"" );

    # Link or generate fasta index file
    if (-f "${fastaFile}.fai") {
      #print "Adding ${fastaFile}.fai to reference genome\n";
      system( "ln -s \"${fastaFile}.fai\" \"${refGenomeDir}\"" );
    }
    else {
      my $cmd = "samtools faidx \"${refGenomeDir}/${baseName}\"";
      print "Generating ${baseName}.fai with: $cmd\n";
      system( $cmd );
    }

    # Add each fasta index entry to genome_size.xml
    my $faiFilename = "${refGenomeDir}/${baseName}.fai";
    open FASTA_INDEX, "<${faiFilename}" or die "Can't open ${faiFilename}";
    while (<FASTA_INDEX>) {
      if (/^([^\t]*)\t([0-9]*)/) {
        my $contigName = $1;
        my $totalBases = $2;
        print "Adding contig ${contigName} of size ${totalBases}\n";
        print GENOMESIZE_XML "    <chromosome fileName=\"$baseName\" contigName=\"$contigName\" totalBases=\"$totalBases\"/>\n";
        if ($contigName =~ /:/) {
          print "Error: Sorry, contig names shouldn't include ':' characters (which appears in \"$contigName\")\n";
          exit 1;
        }
        push @bamChromosomes, "$contigName";

        # stats
        $genomeSize += $totalBases;
        $maxChromosomeLength = max( $maxChromosomeLength, $totalBases);
      }
    }
    close FASTA_INDEX;
  }
}
print GENOMESIZE_XML "</sequenceSizes>\n";
close GENOMESIZE_XML;


my $config = RunFolder::Config->new();
$config->declare( '#ORIGINAL_PATH', $originalPath );
$config->declare( '#ORIGINAL_CMD_LINE', $originalCmdLine );
$config->declare( 'MAKEFILES_DIR',     File::Spec->catdir($EAGLE_DATADIR,"Workflows") );
$config->declare( 'REFERENCE_MODE',    "whole-genome");
$config->declare( 'REFERENCE_GENOME',  "$refGenomeDir");
$config->declare( 'VARIANT_LIST',      "@fullVariantList");
$config->declare( 'BAM_CHROMOSOMES',   "@bamChromosomes");
$config->declare( 'EAGLE_FORCE',       '--force')       if $alwaysForce;
$config->declare( 'QUALITY_TABLE',     "@fullQualityTable");
$config->declare( 'QQ_TABLE',          $fullQQTable);
$config->declare( 'TEMPLATE_LENGTH_TABLE', $fullTemplateLengthTable);
$config->declare( 'MISMATCH_TABLE',    $fullMismatchTable);
$config->declare( 'HOMOPOLYMER_INDEL_TABLE', $fullHomopolymerIndelTable);
$config->declare( 'MOTIF_QUALITY_DROP_TABLE', $fullMotifQualityDropTable);
$config->declare( 'GC_COVERAGE_FIT_TABLE', $fullGcCoverageFitTable) if defined $gcCoverageFitTable;
$config->declare( 'ERROR_MODEL_OPTIONS',"@fullErrorModelOptions");
$config->declare( 'COVERAGE_DEPTH',    $alleleCoverageDepth)  if defined $alleleCoverageDepth;
$config->declare( 'RANDOM_SEED',       $randomSeed)     if defined $randomSeed;
$config->declare( 'BASES_PER_CLUSTER', $basesPerCluster);
$config->declare( 'EAGLE_OUTDIR',      $fullOutputDir);
$config->declare( 'SAMPLE_GENOME',     $sampleGenome);
$config->declare( 'RUN_FOLDER',        $runFolder);
$config->declare( 'GENOME_MUTATOR_OPTIONS',       $genomeMutatorOptions)       if defined $genomeMutatorOptions;
$config->declare( 'FRAGMENTS_ALLOCATOR_OPTIONS',  $fragmentsAllocatorOptions)  if defined $fragmentsAllocatorOptions;
$config->declare( 'RUN_FOLDER_GENERATOR_OPTIONS', $runFolderGeneratorOptions)    if defined $runFolderGeneratorOptions;
$config->declare( 'SEQUENCER_SIMULATOR_OPTIONS',  $sequencerSimulatorOptions)  if defined $sequencerSimulatorOptions;
$config->declare( 'SAMTOOLS',          $samtools);

$config->declare( 'LANES', (join ' ', 1..$layout->{LaneCount}));
if (1 == $layout->{SurfaceCount} || 2 == $layout->{SurfaceCount})
{
    my @tiles = ();
    foreach my $surface (1..$layout->{SurfaceCount})
    {
        foreach my $swath (1..$layout->{SwathCount})
        {
            push @tiles, map{ sprintf "%1d%1d%02d",$surface,$swath,$_ } (1..$layout->{TileCount});
        }
    }
    $config->declare( 'TILES', "@tiles" );
} else {
    croak "ERROR: *** wrong 'SurfaceCount' in ${runInfoFile}. Expected {1,2}, got '$layout->{SurfaceCount}' ***\n   "
}

{
  my $fullMakefile = File::Spec->catfile( $fullOutputDir, "Makefile" );
  open my $makefileHandle, '>', "$fullMakefile.tmp"          or croak("ERROR: *** failed to open '$fullMakefile' for writing: $! ***\n   ");
  $config->generateMakefile($makefileHandle,$workflow)   or croak("ERROR: *** could not save '$fullMakefile' in disk ***\n   ");
  close $makefileHandle                                  or croak("ERROR: *** failed to close '$fullMakefile': $! ***\n   ");

  my $cmd = "make -f $fullMakefile.tmp print-output-contig-names 2> /dev/null | grep -P \"^Chromosome allele: \" | sed 's/^Chromosome allele: //g' | tr '\n' ' '";
  if ($DEBUG) {
    print "Running: $cmd\n";
  }
  my $chromosomeAlleles = `$cmd`;
  print "Alleles: $chromosomeAlleles\n";
  if ($chromosomeAlleles eq "") {
    print "Error: Chromosomes alleles were not detected properly by command $cmd\n";
    exit -1;
  }
  $config->declare( 'CHROMOSOME_ALLELES',  $chromosomeAlleles)  if defined $sequencerSimulatorOptions;

  open $makefileHandle, '>', "$fullMakefile"             or croak("ERROR: *** failed to open '$fullMakefile' for writing: $! ***\n   ");
  $config->generateMakefile($makefileHandle,$workflow)   or croak("ERROR: *** could not save '$fullMakefile' in disk ***\n   ");
  close $makefileHandle                                  or croak("ERROR: *** failed to close '$fullMakefile': $! ***\n   ");
}

# Calculate and report expected RAM usage, to help the user choose an appropriate machine or queue
print "\nExpected RAM usage:\n";
my $genomeMutatorExpectedRamUsage = $genomeSize + $maxChromosomeLength + 20 * $totalVariantFilesLength;
my $totalTileCount = $layout->{TileCount};
if (defined $layout->{SwathCount} && $layout->{SwathCount} > 0) { $totalTileCount *= $layout->{SwathCount}; }
if (defined $layout->{SurfaceCount} && $layout->{SurfaceCount} > 0) { $totalTileCount *= $layout->{SurfaceCount}; }
if (defined $layout->{LaneCount} && $layout->{LaneCount} > 0) { $totalTileCount *= $layout->{LaneCount}; }
my $sequencersimulatorExpectedRamUsage = $genomeSize * $alleleCoverageDepth / $totalTileCount + $maxChromosomeLength;
if ($genomeMutatorExpectedRamUsage == 0) { $genomeMutatorExpectedRamUsage = "unknown"; }
if ($sequencersimulatorExpectedRamUsage == 0) { $sequencersimulatorExpectedRamUsage = "unknown"; }
print " - Genome Mutator (sequential process)    : ~" . int($genomeMutatorExpectedRamUsage/(1024*1024)+1) . " MB\n";
print " - Sequencer Simulator (parallel per tile): ~" . int($sequencersimulatorExpectedRamUsage/(1024*1024)+1) . " MB per process\n";


print "\nSuccessfully configured '$fullOutputDir'. Type 'make' inside this directory to kick-off simulation.\n\n";


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

Mauricio Varea

=cut

