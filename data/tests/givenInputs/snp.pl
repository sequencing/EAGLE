#!/usr/bin/env perl
use strict;
use warnings "all";

print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILT\tINFO\n";
my @chr = (['chr1','A'],['chr2','C'],['chr3','G'],['chr4','T']);
for (my $i = 0; $i < 4; $i++)
{
    foreach my $pos (1..280)
    {
        my $offset =  ($pos - 1) % 4;
        print join("\t",( $chr[$i][0], $pos, ".", $chr[$i][1], $chr[$offset][1], "0", "PASS", "SVTYPE=SNP\n" )) 
           if ($offset != $i);
    }
}
