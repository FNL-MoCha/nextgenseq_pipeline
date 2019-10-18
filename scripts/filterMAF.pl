#!/usr/bin/perl
use strict;
use warnings;
$| = 1;

open(IN, "$ARGV[0]") or die "cannot open file $ARGV[0]:$!\n";
my $comment =<IN>;
my $title  = <IN>;
chomp $title;
my @title = split(/\t/,$title);
print "$comment"."$title\ti_TumorVAF\n";
my ($class, $Filter, $Total, $Ref, $Alt);
for(0..$#title){
	$class  =$_ if $title[$_] eq "Variant_Classification";
	$Filter =$_ if $title[$_] eq "gnomAD_AF";
	$Total  =$_ if $title[$_] eq "t_depth";
	$Ref    =$_ if $title[$_] eq "t_ref_count";
	$Alt    =$_ if $title[$_] eq "t_alt_count";

}
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	if(!defined $line[$Filter]){
		$line[$Filter] =0;
	}
	if($line[$Filter] <=0.01){
		if($line[$class] =~ /Missense/ or $line[$class] =~ /Nonsense_Mutation/ or $line[$class] =~ /Nonstop_Mutation/ or $line[$class] =~ /Frame/ or $line[$class] =~ /Splice_Site/){
			my $vaf=0;
			if ($line[$Alt] =~/^[0-9]+$/ and $line[$Total]=~ /^[0-9]+$/ and $line[$Total] >0 and $line[$Alt] >0){
				$vaf=sprintf("%.4f",$line[$Alt]/$line[$Total]);
			}
			print "$_\t$vaf\n";
		}
	}
}
close IN;
