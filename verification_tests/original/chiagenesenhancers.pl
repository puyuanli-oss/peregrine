#!/usr/local/bin/perl -w
#use strict;

#chiagenesenhancers.pl updated 10/11/2018
#determines which genes' TSS are in ChIA-PET regions with interacting partners containing enhancer overlap

#ARGV[0] is the gene input file e.g. chiaTHUenhancers
#ARGV[1] is chiaTHUgenes
#ARGV[2] is the output file containing the links e.g. linksTHUchia
#ARGV[3] is pantherGenelist.txt
#ARGV[4] is tissuetable

use Time::HiRes qw(usleep ualarm gettimeofday tv_interval nanosleep clock_gettime clock_getres clock_nanosleep clock stat);

use POSIX qw(ceil floor);

my $start_time = [Time::HiRes::gettimeofday()];

my ($chr, $start, $end, $enhID, $one, $two, $three, $chia, $cell, %chia, $region, $thing1, $thing2, %enhancer, $enhancer, $pvalue);

open (ENH,'<',$ARGV[0]) or die $!; #or chiaHCT116POLR2Aenhancers
#chr5    784651  785023  50016   chr5    784065  785788  chr5:692660..694438-chr5:784065..785788 HCT116-POLR2A_0.00426355355394808

#chia pet region->enhancer or gene
#for every gene, search its interacting chia pet regions for enhancers

while(my $line=<ENH>){
	chomp $line;
	($chr, $start, $end, $enhID, $one, $two, $three, $chia, $cell)=(split /\t/, $line);
	($thing1, $thing2)=(split /\-/, $chia);
	($cell, $pvalue)=(split /\_/, $cell);
	$thing1="$thing1\t$cell";
	$thing2="$thing2\t$cell";
	$chia{$thing1}{$thing2}=$pvalue;
	$chia{$thing2}{$thing1}=$pvalue;
	$region="$one:$two..$three\t$cell";
	$enhancer="$chr\t$start\t$end\t$enhID";
	$enhancer{$region}{$enhID}=1;
}

my ($geneID, %genes, $gene, $tissue, %final, $count1, $count2, $panth, %panther, $ensg, $now, $tiss, %tissues, $tissueID);

open (GENE,'<',$ARGV[1]) or die $!; #or chiaTHUgenes
#chr5    693510  694110  ENSG00000171368 chr5    692660  694438  chr5:692660..694438-chr5:784065..785788 HCT116-POLR2A_0.00426355355394808

while(my $line=<GENE>){
	chomp $line;
	($chr, $start, $end, $geneID, $one, $two, $three, $chia, $cell)=(split /\t/, $line);
	($thing1, $thing2)=(split /\-/, $chia);
	#chr1:1014455..1014972
	($cell, $pvalue)=(split /\_/, $cell);
	$thing1="$thing1\t$cell";
	$thing2="$thing2\t$cell";
	$chia{$thing1}{$thing2}=$pvalue;
	$chia{$thing2}{$thing1}=$pvalue;
	$region="$one:$two..$three\t$cell";
	$gene="$chr\t$start\t$end\t$geneID";
	$genes{$geneID}{$region}=1;
}

open OUT, "> $ARGV[2]" or die $!; 
#enhancer gene tissue assay
foreach my $ensg (sort keys %genes){
	foreach my $reg (sort keys %{$genes{$ensg}}){
		if ($reg ~~ %chia){
			foreach my $partner (sort keys %{$chia{$reg}}){
			$count1++;
				if ($partner ~~ %enhancer){
					foreach my $enh (sort keys %{$enhancer{$partner}}){
					$count2++;
					$tissue=((split /\t/, $reg)[1]);
					$cell=((split /\-/, $tissue)[0]);
					$pvalue=$chia{$partner}{$reg};
					print "$pvalue\t";
					$final{$enh}{$ensg}{$cell}=1;
					}
				}
			}
		}
	}
}
open (PAN,'<',$ARGV[3]) or die $!; #or pantherGenelist.txt
#HUMAN|HGNC=24745|UniProtKB=Q8NA69	ENSG00000198723	Uncharacterized protein C19orf45;C19orf45;ortholog	MCG11564 (PTHR34828:SF1)		Homo sapiens
while(my $line=<PAN>){
	chomp $line;
	($panth, $ensg)=((split /\t/, $line)[0,1]);
	$panther{$ensg}=$panth;
}
open (TISS,'<',$ARGV[4]) or die $!; #tissuetable
while (my $line=<TISS>){
#6	Artery Tibial
	chomp $line;
	($tissueID, $tissue)=(split /\t/, $line); 
	$tissue =~ s/\-//g;
	$tissues{$tissue}=$tissueID;
}
print scalar keys %tissues;
print "\n$count1\t$count2";
foreach my $one (sort keys %final){
	foreach my $two (sort keys %{$final{$one}}){
		foreach my $three (sort keys %{$final{$one}{$two}}){
			$now=$panther{$two};
			$tiss=$tissues{$three};
			print OUT "$one\t$now\t$tiss\t1\n";
			#11412	HUMAN|HGNC=17840|UniProtKB=Q9BSY4	55	1
		}
	}
}
close OUT;

my $diff = Time::HiRes::tv_interval($start_time);

print "\n\n$diff\n";