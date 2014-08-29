#!/usr/bin/env perl

=head1 LICENSE

Copyright (c) 2014 Illumina, Inc.

This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
covered by the "BSD 2-Clause License" (see accompanying LICENSE file)

=head1 NAME

RunFolder::RunInfo - Processing of the RunFolder parameters

Library for the extraction of various information from the RunInfo.xml file.

=head1 VERSION

EAGLE vB<@EAGLE_VERSION@>

=head1 DESCRIPTION

This library is intended to work on several versions of the RunFolder.xml.

=head1 AUTHOR

Mauricio Varea

=head1 SYNOPSIS

  use RunFolder::RunInfo;
  my $runInfo = RunFolder::RunInfo->new();
  $runInfo->load($inputHandle);
  ... do something ...
  $runInfo->save($outputHandle);

=cut

package RunFolder::RunInfo;

=head1 DESCRIPTION

=head2 Overview

Provides a number of configuration-related services for EAGLE;

=head2 Dependencies

  strict, warnings, Carp, Exporter, XML::Simple, Scalar::Util, Data::Dumper

=cut

use strict;
use warnings "all";
use Carp;
use Exporter 'import';
use XML::Simple;
use Scalar::Util qw(looks_like_number);
use Data::Dumper;
$Data::Dumper::Terse  = 1;
$Data::Dumper::Indent = 1;

=head2 Exports

  &new, $RunInfoXml

=cut

our @EXPORT_OK = qw(&new $RunInfoXml);
our $RunInfoXml  = "RunInfo.xml";

my $DEBUG        = '@EAGLE_DEBUG_MODE@';

=head2 Methods

=item new(;$)

Constructor

=cut

sub new(;$)
{
    my $proto = shift;
    my $class = ref($proto) || $proto;
    my $self  = {};

    $self->{VERSION}       = undef;
    $self->{FLOWCELL}      = undef;
    $self->{INSTRUMENT}    = undef;
    $self->{DATE}          = undef;
    $self->{NUMBER}        = undef;
    $self->{METHOD}        = undef;
    $self->{READS}         = {};
    $self->{LAYOUT}        = {};
    # $self->{ALIGN_TO_PHIX} = {Lane => []};  ## TODO: Not yet implemented

    bless ($self, $class);
    return $self;
}

=item version($;$)

Either obtains or fills in the XML version.

=cut

sub version($;$)
{
    my $self = shift;
    if (@_) { $self->{VERSION} = shift; }
    return $self->{VERSION};
}

=item method($;$)

Either obtains or fills in the TileNameMethod attribute.

=cut

sub method($;$)
{
    my $self = shift;
    if (@_) { $self->{METHOD} = shift; }
    return $self->{METHOD};
}

=item flowcell($;$)

Either obtains or fills in the flowcell ID.

=cut

sub flowcell($;$)
{
    my $self = shift;
    if (@_) { $self->{FLOWCELL} = shift; }
    return $self->{FLOWCELL};
}

=item instrument($;$)

Either obtains or fills in the instrument name.

=cut

sub instrument($;$)
{
    my $self = shift;
    if (@_) { $self->{INSTRUMENT} = shift; }
    return $self->{INSTRUMENT};
}

=item date($;$)

Either obtains or fills in the date of a particular run.

=cut

sub date($;$)
{
    my $self = shift;
    if (@_)
    {
        $self->{DATE} = shift;
        croak "ERROR: *** ".__PACKAGE__."::date() : Date must be YYMMDD ".$self->_given($self->{DATE})." ***\n   "
              unless defined $self->{DATE} && looks_like_number($self->{DATE});
    }
    return $self->{DATE};
}

=item number($;$)

Either obtains or fills in the run number.

=cut

sub number($;$)
{
    my $self = shift;
    if (@_)
    {
        $self->{NUMBER} = shift;
        croak "ERROR: *** ".__PACKAGE__."::number() : Run number must be numeric ".$self->_given($self->{NUMBER})." ***\n   "
              unless defined $self->{NUMBER} && looks_like_number($self->{NUMBER});
    }
    return $self->{NUMBER};
}

=item runId($;$)

Either obtains or fills in the run ID.

B<Parameters:>

  $stringId - Run ID  "<DATE>_<INSTRUMENT>_<NUMBER>_<FLOWCELL>"

B<Returns:>

    Run ID  "<DATE>_<INSTRUMENT>_<NUMBER>_<FLOWCELL>"
    or  UNDEF if Run ID was not complete.

=cut

sub runId($;$)
{
    my $self = shift;
    if (@_)
    {
        my $stringId = shift;
        my @fields = split '_', $stringId;
        if( 4 == @fields )
        {
            $self->date(       $fields[0] );
            $self->instrument( $fields[1] );
            $self->number(     $fields[2] );
            $self->flowcell(   $fields[3] );
        } else {
            croak "ERROR: *** ".__PACKAGE__."::runId() : Could not parse Run ID: $stringId ***\n   ";
        }
    }
    return undef unless defined $self->{DATE}
                     && defined $self->{INSTRUMENT}
                     && defined $self->{NUMBER}
                     && defined $self->{FLOWCELL};
    return join '_', ($self->{DATE}, $self->{INSTRUMENT}, $self->{NUMBER}, $self->{FLOWCELL});
}

=item reads($)

Obtain the number of reads.

=cut
 
sub reads($)
{
    my $self = shift;
    return sort keys %{$self->{READS}};
}

=item read($$;\%)

Either obtains or fills information about read number N.

B<Parameters:>

  $readNum - Read Number
  $hashRef - {FirstCycle    => (int > 0), 
              LastCycle     => (int > 0),
              NumCycles     => (int > 0),
              IsIndexedRead => (bool, Y/N)}

B<Returns:>

    {FirstCycle    => (int > 0), 
     LastCycle     => (int > 0),
     NumCycles     => (int > 0),
     IsIndexedRead => (bool, Y/N)}
    or  UNDEF if $readNum was not found.

=cut

sub read($$;\%)
{
    my $self = shift;
    my $readNum = shift;
    croak "ERROR: *** ".__PACKAGE__."::read() : Read number must be numeric ".$self->_given($readNum)." ***\n   "
          unless defined $readNum && looks_like_number($readNum);
    if (@_)
    {
        my $hashRef = shift;
        if( 3 <= keys(%{$hashRef})
        &&  defined( $hashRef->{FirstCycle} )  && $hashRef->{FirstCycle} > 0
        &&  defined( $hashRef->{LastCycle} )   && $hashRef->{LastCycle} > 0
        &&  defined( $hashRef->{NumCycles} )   && $hashRef->{NumCycles} > 0 )
        {
            #carp "WARNING: Could not gather whether read is indexed or not\n     "
            #     unless defined $hashRef->{IsIndexedRead} && $hashRef->{IsIndexedRead} =~ /^[YN]$/;
            $self->{READS}->{$readNum} = $hashRef;
        } else {
            croak "ERROR: *** ".__PACKAGE__."::read() : Read info does not follow the specs ***\n   ".Dumper($hashRef);
        }
    }
    return undef unless exists $self->{READS}->{"$readNum"};
    return $self->{READS}->{"$readNum"};
}

=item layout($;\%)

Either obtains or fills in the layout information.

B<Parameters:>

  $hashRef - {LaneCount    => (int > 0), 
              SurfaceCount => (int > 0),
              SwathCount   => (int > 0),
              TileCount    => (int > 0)}

B<Returns:>

    {LaneCount    => (int > 0), 
     SurfaceCount => (int > 0),
     SwathCount   => (int > 0),
     TileCount    => (int > 0)}
    or  UNDEF if layout was empty.

=cut

sub layout($;\%)
{
    my $self = shift;
    if (@_)
    {
        my $hashRef = shift;
        if( 4 == keys(%{$hashRef})
        &&  defined( $hashRef->{LaneCount} )    && $hashRef->{LaneCount} > 0
        &&  defined( $hashRef->{SurfaceCount} ) && $hashRef->{SurfaceCount} > 0
        &&  defined( $hashRef->{SwathCount} )   && $hashRef->{SwathCount} > 0
        &&  defined( $hashRef->{TileCount} )    && $hashRef->{TileCount} > 0)
        {
            $self->{LAYOUT} = $hashRef;
        } else {
            croak "ERROR: *** ".__PACKAGE__."::layout() : Layout info does not follow the specs ***\n   ".Dumper($hashRef);
        }
    }
    return undef unless exists $self->{LAYOUT};
    return $self->{LAYOUT};
}

=item load($$)

Fills class with information from a RunInfo.xml file.

B<Parameters:>

  $handle - Handle to RunInfo.xml

B<Returns:>

    N/A

=cut

sub load($$)
{
    my $self = shift;
    my $handle = shift;
    my $xml = XMLin($handle, ForceArray => 1, KeepRoot => 1) or croak "ERROR: *** couldn't load $RunInfoXml: $! ***\n   ";
    croak "ERROR: *** no top level 'RunInfo' element in $RunInfoXml ***\n   " unless exists $xml->{RunInfo};
    croak "ERROR: *** no 'Run' element in $RunInfoXml. At least 1 required ***\n   " unless exists $xml->{RunInfo}->[0]->{Run};
    carp  "WARNING: More than one 'Run' element given in $RunInfoXml.\n         Will only process the first one!\n"
          if 1 < @{$xml->{RunInfo}->[0]->{Run}};
    $self->version( $xml->{RunInfo}->[0]->{Version} );
    my $run       = $xml->{RunInfo}->[0]->{Run}->[0];
    
    print "Loading $RunInfoXml for '$run->{Id}' ...\n";
    $self->runId($run->{Id});
    # $self->number( $run->{Number} )               if defined $run->{Number};
    # $self->flowcell( $run->{Flowcell}->[0] )      if $run->{Flowcell};
    # $self->instrument( $run->{Instrument}->[0] )  if $run->{Instrument};
    # DON'T overwrite DATE, as some RunInfo.xml have a weird string there
    $self->method( $run->{TileNameMethod} )       if defined $run->{TileNameMethod};
    carp "WARNING: no 'Read' element in <READS /> of $RunInfoXml\n     "  unless @{$run->{Reads}};
    my ($readNum,$prevCycle) = (0,0);
    foreach my $read (@{$run->{Reads}->[0]->{Read}})
    {
        $readNum++;
        $readNum = $read->{Number}  if defined $read->{Number};
        my %readCfg = %{$read};
        delete $readCfg{Number}     if exists $readCfg{Number};

        if (!exists $readCfg{NumCycles} 
        && (exists $readCfg{FirstCycle} && exists $readCfg{LastCycle}) )
        {
            $readCfg{NumCycles} = $readCfg{LastCycle} - $readCfg{FirstCycle} + 1;
        } elsif (exists $readCfg{NumCycles}
             && !(exists $readCfg{FirstCycle} || exists $readCfg{LastCycle}) ) {
                 $readCfg{FirstCycle} = $prevCycle + 1;
                 $readCfg{LastCycle} = $readCfg{FirstCycle} + $readCfg{NumCycles} - 1;
                 $prevCycle = $readCfg{LastCycle};
        }
        croak "ERROR: *** could not parse $RunInfoXml ***\n   "
              unless (exists $readCfg{NumCycles} && exists $readCfg{FirstCycle} && exists $readCfg{LastCycle});
        print "... read #$readNum: [$readCfg{FirstCycle},$readCfg{LastCycle}] ($readCfg{NumCycles} cycles)\n";
        $self->read( 0 + $readNum, \%readCfg );
    }
    carp "WARNING: no 'FlowcellLayout' element in $RunInfoXml\n     "  unless @{$run->{FlowcellLayout}};
    $self->layout( $run->{FlowcellLayout}->[0] );
    return $self;
}

=item save($$)

Dumps class into a RunInfo.xml file.

B<Parameters:>

  $handle - Handle to RunInfo.xml

B<Returns:>

    N/A

=cut

sub save($$)
{
    my $self = shift;
    my $handle = shift;

    my $run = {};
    $run->{Id} = $self->runId;
    $run->{Number} = $self->number;
    $run->{TileNameMethod} = $self->method  if defined $self->method;
    $run->{Flowcell} = [ $self->flowcell ];
    $run->{Instrument} = [ $self->instrument ];
    $run->{Date} = [ $self->date ];
    $run->{FlowcellLayout} = [ $self->layout ];
    foreach my $r ($self->reads)
    {
        my $read = $self->read($r);
        $read->{Number}=$r;
        push @{$run->{Reads}->[0]->{Read}}, $read;
    }
    my $xml = {};
    push @{$xml->{RunInfo}}, {'Version'   => $self->version, 
                              'Run'       => [ $run ],
                              'xmlns:xsi' => 'http://www.w3.org/2001/XMLSchema-instance',
                              'xmlns:xsd' => 'http://www.w3.org/2001/XMLSchema'};
    print $handle XMLout($xml, KeepRoot => 1, NoSort => 1);
    return $self;
}

##
## Internal Methods:
##

sub _given($$)
{
    my $self = shift;
    my $var  = shift;
    return (defined $var ? "(given: $var)" : "(not given)");
}


1;
__END__

=head1 BUGS AND LIMITATIONS

There are no known bugs in this module.

All documented features are fully implemented.

Please report problems to Illumina Technical Support (support@illumina.com)

Patches are welcome.

=cut
