#!/usr/bin/perl
use strict;
use warnings;


open(IN, "$ARGV[0]") or die "cannot open file $ARGV[0]:$!\n";
my $comment =<IN>;
my $title = <IN>;
my @title = split(/\t/,$title);
chomp $title;
print "$comment"."$title\ti_TumorVAF\n";
#print "$title\n";
my ($Gene, $Chromosome, $Start_Position, $Tumor_Sample_Barcode, $Reference_Allele, $Tumor_Seq_Allele1, $Tumor_Seq_Allele2, $Variant_Classification, $Filter, $HGVS ,$Total, $Ref, $Alt);
for(0..$#title){
	$Gene = $_ if $title[$_] =~ "Hugo_Symbol";
	$Chromosome = $_ if $title[$_] =~ "Chromosome";
	$HGVS = $_ if $title[$_] =~ "HGVSp_Short";
	$Start_Position = $_ if $title[$_] =~ "Start_Position";
	$Tumor_Sample_Barcode = $_ if $title[$_] eq "Tumor_Sample_Barcode";
	$Reference_Allele = $_ if $title[$_] eq "Reference_Allele";
	$Tumor_Seq_Allele1 = $_ if $title[$_] eq "Tumor_Seq_Allele1";
	$Tumor_Seq_Allele2 = $_ if $title[$_] eq "Tumor_Seq_Allele2";
	$Variant_Classification = $_ if $title[$_] eq "Variant_Classification";
	$Filter = $_ if $title[$_] eq "FILTER";
	$Total=$_ if $title[$_] eq "t_depth";
	$Ref=$_ if $title[$_] eq "t_ref_count";
	$Alt=$_ if $title[$_] eq "t_alt_count";

}
#print "Chr\tStart\tEnd\tRef\tAlt\tSample\tGene\tVariant_Classification\tHGVSp\tTotalDepth\tRefDepth\tAltDepth\tVAF\n";
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	next unless defined $line[$Tumor_Seq_Allele1];
	my $key;
	if($line[$Start_Position] =~ /^[0-9]+$/ and $line[$Filter] !~ /common_variant/ and ($line[$Variant_Classification] =~ /Missense/ or $line[$Variant_Classification] =~ /Nonsense_Mutation/ or $line[$Variant_Classification] =~ /Nonstop_Mutation/ or $line[$Variant_Classification] =~ /Frame/ or $line[$Variant_Classification] =~ /Splice_Site/) and $line[$Total]=~ /^[0-9]+$/ and $line[$Total] >0){
		my $vaf=0;
		if ($line[$Alt] =~/^[0-9]+$/ and $line[$Total]=~ /^[0-9]+$/){
			$vaf=sprintf("%.2f", $line[$Alt]/$line[$Total]);
		}
		print "$_\t$vaf\n";
	}
}

close IN;
