#!/usr/local/bin/perl -w
#use strict;

#orderlinks.pl updated 10/22/2018
#takes files containing links and sorts them by enhancer then gene then tissue

#ARGV[0] is the file e.g. "prelinksDBtad"
#ARGV[1] the output file in nicer order e.g. "linksDBtad"

use Time::HiRes qw(usleep ualarm gettimeofday tv_interval nanosleep clock_gettime clock_getres clock_nanosleep clock stat);

use POSIX qw(ceil floor);

my $start_time = [Time::HiRes::gettimeofday()];

my ($enhID, $gene, $tissue, $assay, %hash, %genes, %tissues, $count, $pval);

open (ENH,'<',$ARGV[0]) or die $!; #linkschia
#64525	HUMAN|HGNC=18786|UniProtKB=Q14320	65	3
while(my $line=<ENH>){
	chomp $line;
	($enhID, $gene, $tissue, $assay)=(split /\t/, $line);
	$hash{$enhID}{$gene}{$tissue}=1;
}

open OUT, "> $ARGV[1]" or die $!;  #linksDBchia

foreach my $enhancer (sort keys %hash){
	foreach my $geen (sort keys %{$hash{$enhancer}}){
		$genes{$geen}=1;
		foreach my $cell (sort keys %{$hash{$enhancer}{$geen}}){
			$tissues{$cell}=1;
			$count++;
			print OUT "$enhancer\t$geen\t$cell\t$assay\n";
		}
	}
}
close OUT;

print "This file contains ";
print scalar keys %hash;
print " different enhancers linked to ";
print scalar keys %genes;
print " genes in ";
print scalar keys %tissues;
print " tissues.";

my $diff = Time::HiRes::tv_interval($start_time);

print "\n\n$diff\n";