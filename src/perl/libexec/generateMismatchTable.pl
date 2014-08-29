#!/usr/bin/env perl
#
# Copyright (c) 2014 Illumina, Inc.
#
# This file is part of Illumina's Enhanced Artificial Genome Engine (EAGLE),
# covered by the "BSD 2-Clause License" (see accompanying LICENSE file)
#
# \author fnielsen@illumina.com
# 
# Usage: mutationRate.pl -o outputdir -e exportfile
# 
# export file format, example line below: 
## 
# GA208   61      6       1       9447    949     ACAGTG  1       NTTTTTACTGCATTATGTCAGAAATTAAAGATACCAGATCATGTCAGAGAGAGAGCTTGGTTAACTTGGGAGAAA     BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB     c13.fa          48881444        F       GA73    24      136                     174     R       Y
#
# export file format documentation: 
#  http://skaros.chuk.illumina.com:8090/display/docs/export+file+record+fields
#

use strict;
use Getopt::Std;
use File::Basename;
use Compress::Zlib;

my $gzerrno; 

my $DEBUG=1; 

my $output_dir;
my $exportfile;

my %opts;
getopt('oev', \%opts);
$output_dir = $opts{'o'};
$exportfile = $opts{'e'};
$DEBUG = $opts{'v'};

unless (defined($output_dir) && defined($exportfile)){

  print "Usage: $0 -o outputdir -e exportfile\n";
  exit 1;
}

my %muttable; 

sub init_muttable{
  my @bases = ("A", "C", "G", "T", "N", "_");
  foreach my $refbase (@bases){
    foreach my $readbase (@bases){
      my $key = $refbase."->".$readbase;
      $muttable{$key} = 0;
    }
  }
}

sub update_muttable {
  my ($refbase, $readbase) = @_; 
  my $key = $refbase."->".$readbase;
  if (defined $muttable{$key} ){
    $muttable{$key}=$muttable{$key}+1; 
  }
  else { 
    $muttable{$key}=1;
  }
}

sub print_muttable {
  # Precondition: $output_dir is defined and $output_dir exists
  # TODO: ASSERT preconditions

  # open output file
  open(OUTFILE, ">", $output_dir."/mutationTable.tsv");

  # print header
  print OUTFILE "# file created by mutationRate.pl\n#\n";
  print OUTFILE "shorthand\trefbase\treadbase\tcount\n"; 

  # print table
  foreach my $k ( sort (keys(%muttable)) ){ 
    
    printf OUTFILE "%s\t%s\t%4d\n", $k, join("\t",split(/->/, $k)), $muttable{$k};
  }

  # close output file
  close OUTFILE;
}

# trunc_str (string, x)
# 
# truncate a string with x characters from the front of the string
# e.g. trunc_str("ABCD",2) : "CD"
#
sub trunc_str {
  my $str2truncate = shift(@_);
  my $truncLength = shift(@_);

  my $strlength = length($str2truncate);
  my $offset = -($strlength-$truncLength);
  my $returnstr = substr $str2truncate, $offset;

  if ($truncLength == $strlength){
    $returnstr="";
  }

  if ($DEBUG){
    print "'$str2truncate' - $truncLength (strlength:$strlength,offset:$offset) = '$returnstr'\n";
  }

  return $returnstr; 
}

my $gz = gzopen($exportfile, "rb") 
  or die "Cannot open $exportfile: '$gzerrno'\n" ;

print "init muttable...\n" if $DEBUG;
init_muttable();

print "looping over file...\n" if $DEBUG; 

while ($gz->gzreadline($_) > 0) {
  my @cols = split(/\t/); 
  if (scalar @cols < 15){

    print "PARSE ERROR: check the format of your export file, fields must be tab-separated\n"; 
    print "error in parsing cols: ". $_ if $DEBUG;

  } else {

    my $read = $cols[8]; # column 9
    my $matchDesc = $cols[14]; # column 15 (colcount starts at 0)

    my $readlength = length($read);
    my $i = 0; 

    print "\nread: $read, ".length($read)."\n" if $DEBUG;
    print "matchDesc: $matchDesc\n" if $DEBUG; 

    while (length($matchDesc)>0 && $i < $readlength){
    # loop over matchDesc, 
    #   with a pointer to the position in the read
    #   adding substituted bases to hashtable count 
      if ($matchDesc =~ /^(\d+)/ ){
	# move pointer x bases into read
	my $x = $1; 
	$i += $x; 
	my $matchlength = length($x);

	print "matchDesc: $matchDesc, ".length($matchDesc)."\n" if $DEBUG;
	print "x: $x\n" if $DEBUG;
	print "i: $i\n" if $DEBUG;
	print "length x: ".length($x)."\n" if $DEBUG;
	
	# truncate matchDesc with x bases
	$matchDesc = trunc_str($matchDesc,$matchlength);

      } elsif ($matchDesc =~ /^([NAGCT])/){ 
	# match just _one_ substituted base
	my $refbase = $1; 

	# read current pointer base
	my $readbase = substr $read, $i, 1; 

	# save this substitution to table
	update_muttable($refbase, $readbase);

	print "matchDesc: $matchDesc, ".length($matchDesc)."\n" if $DEBUG;
	print "refbase: $refbase\n" if $DEBUG;
	print "i: $i\n" if $DEBUG;
	print "readbase: $readbase\n" if $DEBUG;

	# move pointer 1 base forward
	$i++;

	# truncate matchDesc with 1 base
	$matchDesc = trunc_str($matchDesc,1);

      } elsif ($matchDesc =~ /^(\^(\d+)\$)/){
      	# match an insertion event
      	my $insertlen = $2;
      	my $matchlength = length($1);
      	my $insertedbases = substr $read, $i, $insertlen;

      	# save to substitution table
	for (my $c=0; $c<$insertlen; $c++){	 
	  # loop over bases
	  # update substitution table for each base in insert 
	  my $refbase = "_";
	  my $readbase = substr $insertedbases, $c, 1; 

	  update_muttable($refbase,$readbase);
	}

      	# move pointer
      	$i += $insertlen; 

      	# truncate matchDesc
	print "INS:matchDesc '".$matchDesc."'\nINS:matchlength:".$matchlength."\n" if $DEBUG;
	$matchDesc = trunc_str($matchDesc,$matchlength);
	print "INS:truncated matchDesc:".$matchDesc."\n" if $DEBUG;

      } elsif ($matchDesc =~ /^(\^([NACGT]+)\$)/ ){
      	# match a deletion event 
	my $deletedbases = $2; 
	my $deletionlength = length($2);
	my $matchlength = length($1); 

	# save to substitution table for each base in deletion
	for (my $c=0; $c<$deletionlength; $c++){
	  # loop over bases
	  # update substitution table for each base
	  my $refbase = substr $deletedbases, $c, 1;
	  my $readbase = "_"; 

	  update_muttable($refbase,$readbase);
	}

	# move read pointer nowhere
	$i += 0; 

	# truncate matchDesc 
	$matchDesc = trunc_str($matchDesc,$matchlength);
      } else {
	print "no match found: '".$matchDesc."'\n";
	exit;
      }

    } # end while, ending looping over read

  }
} 

print "\nFinished looping over file\n";

print_muttable_for_eagle();

die "Error reading from $exportfile: '$gzerrno'\n" 
  if (defined $gzerrno && ($gzerrno != Z_STREAM_END));

$gz->gzclose() ;



sub print_muttable_for_eagle {
  # Precondition: $output_dir is defined and $output_dir exists
  # TODO: ASSERT preconditions

  # open output file
  open(OUTFILE, ">", $output_dir."/mismatchTable.tsv");

  # print header
  print OUTFILE "# file created by $0\n";
  print OUTFILE "# x\tx->A\tx->C\tx->G\tx->T\tdel\tx->insertedA\tx->insertedC\tx->insertedG\tx->insertedT\n";

  # print table
  foreach my $i ("A","C","G","T") {
    print OUTFILE "${i}" . "\t";
    foreach my $j ("A","C","G","T") {
      print OUTFILE $muttable{"${i}->${j}"} . "\t";
      delete $muttable{"${i}->${j}"}; # destroy processed entry
    }

    print OUTFILE $muttable{"${i}->_"};
    delete $muttable{"${i}->_"}; # destroy processed entry

    foreach my $j ("A","C","G","T") {
      print OUTFILE "\t" . ($muttable{"_->${j}"} / 4);
    }
    print OUTFILE "\n";
  }

  # destroy processed entries
  foreach my $j ("A","C","G","T") {
    delete $muttable{"_->${j}"};
  }

  # print non-processed entries
  print OUTFILE "# entries not processed:\n";
  foreach my $k ( sort (keys(%muttable)) ){
    printf OUTFILE "# %s\t%s\t%4d\n", $k, join("\t",split(/->/, $k)), $muttable{$k};
  }

  # close output file
  close OUTFILE;
}

