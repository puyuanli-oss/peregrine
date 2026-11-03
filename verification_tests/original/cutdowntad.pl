#!/usr/local/bin/perl -w
#use strict;

#cutdowntad.pl updated 10/22/2018
#prints a reduced amount of TAD links that only includes links that already appear in either ChIA-PET, TADintxn, or eQTL data

#ARGV[0] is linksDBchia
#ARGV[1] is linksDBnumeqtl
#ARGV[2] is PSYCHIClinksDB
#ARGV[3] is linkstTADtissues
#ARGV[4] is selecttTAD

use Time::HiRes qw(usleep ualarm gettimeofday tv_interval nanosleep clock_gettime clock_getres clock_nanosleep clock stat);
my $start_time = [Time::HiRes::gettimeofday()];
use POSIX qw(ceil floor);

my ($id, $enhancer, $panthid, $tissue, $assay, $link, $rest, $count1, $count2, %hash);

open OUT, "> $ARGV[4]" or die $!; #selecttTAD

open (ONE,'<',$ARGV[0]) or die $!; #linksDBchia
#10	HUMAN|HGNC=26062|UniProtKB=Q96HA4	61	0.00828927800673878	1
while (my $line=<ONE>){
	chomp $line;
	($enhancer, $panthid)=((split /\t/, $line)[0,1]);
	$link="$enhancer\_$panthid";
	$hash{$link}=$1;
}
close ONE;

open (TWO,'<',$ARGV[1]) or die $!; #linksDBnumeqtl
#enhancer	gene	tissue	number_of_eQTL	assay
#10	HUMAN|HGNC=28208|UniProtKB=Q5SV97	25	2	2
while (my $line=<TWO>){
	chomp $line;
	($enhancer, $panthid)=((split /\t/, $line)[0,1]);
	$link="$enhancer\_$panthid";
	$hash{$link}=1;
}
close TWO;

open (THREE,'<',$ARGV[2]) or die $!; #PSYCHIClinksDB
#EH37E0105481	HUMAN|HGNC=15846|UniProtKB=Q9NP74	49	5.4e-14	4
while (my $line=<THREE>){
	chomp $line;
	($enhancer, $panthid)=((split /\t/, $line)[0,1]);
	$link="$enhancer\_$panthid";
	$hash{$link}=$1;
}
close THREE;

open (FOUR,'<',$ARGV[3]) or die $!; #linkstTADtissues
#1       HUMAN|Gene=LPPR4|UniProtKB=Q7Z2D5       117     3
while (my $line=<FOUR>){
	chomp $line;
	($enhancer, $panthid, $tissue, $assay)=(split /\t/, $line);
	$link="$enhancer\_$panthid";
	$rest="$enhancer\t$panthid\t$tissue\t$assay";
	if ($link ~~ %hash){
		$count1++;
		print OUT "$rest\n"; #should we scratch making a separate line for each tissue and just make it a list?
	}
	else{
		$count2++;
	}
}
close FOUR;

close OUT;
print "Number of TAD links found in TADintxn, eQTL, or ChIA: $count1\nNumber of links found only in TAD data: $count2\n";

my $diff = Time::HiRes::tv_interval($start_time);
print "\n\n$diff\n";