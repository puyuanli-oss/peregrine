#!/usr/local/bin/perl -w
#use strict;

#celltypev2.pl updated 10/11/2018

use Time::HiRes qw(usleep ualarm gettimeofday tv_interval nanosleep clock_gettime clock_getres clock_nanosleep clock stat);

my $start_time = [Time::HiRes::gettimeofday()];

#ARGV[0] is the ChIA-PET file fresh out of Mango
#ARGV[1] is the cell type and target protein from the original assay i.e. cell-POLR2A
#ARGV[2] is the output file

my ($chr1, $start1, $end1, $chr2, $start2, $end2, $number, $pval, $name, $column);

my $tadcell=$ARGV[1];

open OUT, "> $ARGV[2]" or die $!;

open (FH,'<',$ARGV[0]) or die $!;
#chr1	28843710	28847368	chr1	28973949	28976946	48	1.19176621106095e-05
while(my $line=<FH>){
	chomp $line;
	($chr1, $start1, $end1, $chr2, $start2, $end2, $number, $pval)=(split /\t/, $line);
	$name="$chr1:$start1..$end1-$chr2:$start2..$end2";
	$column="$tadcell\_$pval";
	print OUT "$chr1\t$start1\t$end1\t$name\t$column\n";
	print OUT "$chr2\t$start2\t$end2\t$name\t$column\n";
}

my $diff = Time::HiRes::tv_interval($start_time);

print "\n\n$diff\n"; 