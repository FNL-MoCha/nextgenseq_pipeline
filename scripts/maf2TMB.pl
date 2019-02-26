#!/usr/bin/perl
use strict;
use warnings;
my $count=0;
open(IN, "$ARGV[0]") or die "cannot open file $ARGV[0]:$!\n";
my $comment =<IN>;
my $title = <IN>;
chomp $title;
my @title = split(/\t/,$title);
my ($Variant_Classification, $Variant_Type, $VAF, $Total, $Ref, $Alt);
for(0..$#title){
	$Variant_Classification =$_ if $title[$_] eq "Variant_Classification";
	$Variant_Type           =$_ if $title[$_] eq "Variant_Type";
	$Total                  =$_ if $title[$_] eq "t_depth";
	$Ref                    =$_ if $title[$_] eq "t_ref_count";
	$Alt                    =$_ if $title[$_] eq "t_alt_count";
	$VAF                    =$_ if $title[$_] eq "i_TumorVAF";

}
while(<IN>){
	chomp;
	my @line = split("\t", $_);
	if($line[$VAF] >=0.05 and $line[$Total] >=25 and $line[$Alt] >=3){
		if($line[$Variant_Classification] =~ /[Missense|Nonsense_Mutation|Nonstop_Mutation|Frame]/){
			$count = $count+1;
		}
	}
}
close IN;
my $tmb=sprintf("%.2f",$count/$ARGV[1]);
print "$tmb\n";
