#!/usr/local/bin/perl -w
#use strict;

#chiageneTSS.pl
#prepares Ensembl files to be overlapped with via bedtools (bed formatting, searches for the 600bp before the TSS)

#ARGV[0] is the gene input file e.g. hg19genes
#ARGV[1] is pantherGenelist.txt
#ARGV[2] is the output file containing the bed formatted genes to be overlapped with i.e. TSSgenesbed

use Time::HiRes qw(usleep ualarm gettimeofday tv_interval nanosleep clock_gettime clock_getres clock_nanosleep clock stat);

use POSIX qw(ceil floor);

my $start_time = [Time::HiRes::gettimeofday()];

my ($gene, $tss, $chr, $start, $end, $hash, $panther, $panth, $ensg, $strand, $count1, $count2, %final);

open (FH,'<',$ARGV[0]) or die $!; #or hg19genes
#ENSG00000233440	23708313	23708703	23708313	1	13
#ENSG00000207157	23726825	23726825	23726725	-1	13

while(my $line=<FH>){
	chomp $line;
	($gene, $tss, $end, $start, $strand, $chr)=(split /\t/, $line);
	if ($strand eq "1"){
		$end=$tss;
		$start=($tss-600);
		$chr="chr$chr";
		$key="$chr\t$start\t$end\t$gene";
		$hash{$gene}=$key;
	}
	elsif ($strand eq "-1"){
		$start=$tss;
		$end=($tss+600);
		$chr="chr$chr";
		$key="$chr\t$start\t$end\t$gene";
		$hash{$gene}=$key;
	}
}

open (PAN,'<',$ARGV[1]) or die $!; #or pantherGenelist.txt
#HUMAN|HGNC=24745|UniProtKB=Q8NA69	ENSG00000198723	Uncharacterized protein C19orf45;C19orf45;ortholog	MCG11564 (PTHR34828:SF1)		Homo sapiens
open OUT, "> $ARGV[2]" or die $!; 
while(my $line=<PAN>){
	chomp $line;
	($panth, $ensg)=((split /\t/, $line)[0,1]);
	$panther{$ensg}=$panth;
}
foreach my $ensgene (sort keys %panther){
	if ($ensgene ~~ %hash){
		$line=$hash{$ensgene};
		$final{$line}=1;
	}
}
foreach my $line (sort keys %final){
	print OUT "$line\n";
}
close OUT;

my $diff = Time::HiRes::tv_interval($start_time);

print "\n\n$diff\n";