#!/usr/bin/perl
use strict;
use warnings;


open(IN, "$ARGV[0]") or die "cannot open file $ARGV[0]:$!\n";
my $comment =<IN>;
my $title = <IN>;
my @title = split(/\t/,$title);
#print "$comment"."$title";
my ($Chromosome, $Start_Position, $Tumor_Sample_Barcode, $Reference_Allele, $Tumor_Seq_Allele1, $Tumor_Seq_Allele2, $gene, $vtype);
for(0..$#title){
	$gene 			= $_ if $title[$_] =~ "Hugo_Symbol";
	$Chromosome 		= $_ if $title[$_] =~ "Chromosome";
	$Start_Position 	= $_ if $title[$_] =~ "Start_Position";
	$Tumor_Sample_Barcode 	= $_ if $title[$_] eq "Tumor_Sample_Barcode";
	$Reference_Allele	= $_ if $title[$_] eq "Reference_Allele";
	$Tumor_Seq_Allele1 	= $_ if $title[$_] eq "Tumor_Seq_Allele1";
	$Tumor_Seq_Allele2 	= $_ if $title[$_] eq "Tumor_Seq_Allele2";
	$vtype			= $_ if $title[$_] eq "Variant_Type";
}
#print "$Chromosome, $Start_Position, $Tumor_Sample_Barcode, $Reference_Allele, $Tumor_Seq_Allele1, $Tumor_Seq_Allele2\n";
#=cut
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	next unless defined $line[$Tumor_Seq_Allele1];
	if($line[$vtype] =~/DEL/ or $line[$vtype] =~/INS/ or $line[$vtype] =~/TNP/ or $line[$vtype] =~/DNP/){
		next;
	}
	my $key;
	if($line[$Start_Position] =~ /^[0-9]+$/){
		my $end = $line[$Start_Position] + length($line[$Tumor_Seq_Allele1]) - 1;
		if($line[$Tumor_Seq_Allele1] eq $line[$Reference_Allele]){
			$key="$line[$Chromosome]\t$line[$Start_Position]\t$end\t$line[$Reference_Allele]\t$line[$Tumor_Seq_Allele2]\t$line[$gene]";
		}else{
			$key="$line[$Chromosome]\t$line[$Start_Position]\t$end\t$line[$Reference_Allele]\t$line[$Tumor_Seq_Allele1]\t$line[$gene]";
		}
		print "$key\n";
	}
}

close IN;
