#!/usr/bin/env perl
use strict;
use warnings;

open (INPUT, $ARGV[0]);
while(<INPUT>){
	chomp;
	if($_ =~ /#/){
##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|ALLELE_NUM|DISTANCE|STRAND|FLAGS|PICK|VARIANT_CLASS|MINIMISED|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|RefSeq|GENE_PHENO|SIFT|PolyPhen|DOMAINS|HGVS_OFFSET|GMAF|AFR_MAF|AMR_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE">
		print "$_\n";
		next;
	}
	my @line=split("\t", $_);
	my @info;
	if($line[7] =~ /CSQ=(.*)/){
		@info=split('\|', $1);
	}
	#print "$line[0]\t$line[1]\t$line[2]\t$line[3]\t$line[4]\n";
	my @A_F;
	for (my $i=40; $i<=47; $i++){
		if(defined $info[$i] and $info[$i] =~ /:/){
			my @new_F=split("&", $info[$i]);
			for(my $j=0; $j<=$#new_F; $j++){
				my ($allele, $freq) =split(":", $new_F[$j]);
				if($allele eq $line[4]){
					push @A_F, $freq;
				}
			}
		}
	}
	@A_F = sort { $b <=> $a } @A_F;
	if (defined $A_F[0]){
		if ($A_F[0] <=0.01){
			print "$_\n";
		}
	}
	else{
		print "$_\n";
	}
}
