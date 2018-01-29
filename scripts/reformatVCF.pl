#!/usr/local/bin/perl
use warnings;
use strict;
use File::Basename;
use List::Util qw(first);
use List::Util qw(sum);
use 5.010;
local $SIG{__WARN__} = sub {
	my $message =shift;
	die $message;
};
#####################################################
# Version: v0.1
# Author: Rajesh Patidar rajbtpatidar@gmail.com
# How to run:
# 		$0 <vcf FILE Name> >File.txt
#
#####################################################
# Take Temporary location from command line
my $input=$ARGV[0];
my $output=$ARGV[1];
open (VAR, $input)or die "Error: cannot read variant file $input: $!\n";
if($output =~/annovar/){
	print "Chr\tStart\tEnd\tRef\tAlt\tQual\tGenotype\tTotalDepth\tRefDepth\tAltDepth\tVAF\n";
}
while (<VAR>) {
	chomp;
	if (m/^#/){
		if ($output=~ /vcf/) {
			if ($. ==1){
				print "$_\n";
				print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
				print "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n";
				print "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">\n";
				print "##FORMAT=<ID=AF,Number=1,Type=Integer,Description=\"Allelic frequency for the alt alleles in the order listed\">\n";
			}
			else{
				if($_ !~ /^##FORMAT/){
					print "$_\n";
				}
			}
			
		}
	}
	else{
		my @field=split(/\t/,$_);
		my ($chr, $position, $rsID, $ref, $alt, $qual, $filter, $info, $format, $sample) = @field;
		my @arr = split(":", $sample);
		my $vaf="";
		my @alleles =split(",", $alt);	
		if($format =~ /^GT:AD:DP:GQ:PL$/){ # Haplotype Caller
			my @AD=split(",",$arr[1]);
			for (my $idx=1; $idx <=@alleles; ++$idx){
				if($AD[$idx] ==0){$vaf ="$vaf,0";}else{
					my $allele_vaf=sprintf("%.2f", $AD[$idx]/$arr[2]);
					$vaf ="$vaf,$allele_vaf";
				}
			}
			$vaf=~ s/^,//s;
				if($output=~ /vcf/){
					print "$chr\t$position\t$rsID\t$ref\t$alt\t$qual\t.\t$info\tGT:DP:AD:AF\t$arr[0]:$arr[2]:$arr[1]:$vaf\n";
				}
				else{
					printAnnovarInput($chr, $position, $ref, $alt, $qual, $arr[0], $arr[2], $arr[1], $vaf);
				}
		}
		elsif($format =~ /^GT:GL:GOF:GQ:NR:NV$/){ # Platypus
			if($alt =~/,/){
				my @AD=split(",",$arr[5]);
				my @DP=split(",",$arr[4]);
				for (my $idx=0; $idx <=$#alleles; ++$idx){
					my $vaf=sprintf("%.2f", $AD[$idx]/$DP[$idx]);
					my $refCount=$DP[$idx] - $AD[$idx];
					$refCount= "$refCount,$AD[$idx]";
					if($output=~ /vcf/){
						print "$chr\t$position\t$rsID\t$ref\t$alleles[$idx]\t$qual\t.\t$info\tGT:DP:AD:AF\t0/1:$DP[$idx]:$refCount:$vaf\n";
					}
					else{
						printAnnovarInput($chr, $position, $ref, $alt, $qual, $arr[0], $arr[2], $arr[1], $vaf);
					}
				}
			}
			else{
				my $refCount=$arr[4] - $arr[5];
				$vaf=sprintf("%.2f", $arr[5]/$arr[4]);
				$arr[5]= "$refCount,$arr[5]";
				if ($output=~ /vcf/){
					print "$chr\t$position\t$rsID\t$ref\t$alt\t$qual\t.\t$info\tGT:DP:AD:AF\t$arr[0]:$arr[4]:$arr[5]:$vaf\n";
				}
				else{
					printAnnovarInput($chr, $position, $ref, $alt, $qual, $arr[0], $arr[2], $arr[1], $vaf);
						
				}
			}
		}
		else{
			print STDERR "Sorry format $format is not recognized\n";
			die;
		}
	}
}
close VAR;
sub printAnnovarInput{
	my ($chr, $position, $ref, $alt, $qual, $genotype, $refCount, $altCount, $vaf) =(@_);
	my @AD=split(",",$altCount);
	if(length($ref)==1 && length($alt)==1) { # Output single Allelic SNPS
		print "$chr\t$position\t$position\t$ref\t$alt\t$qual\t`$genotype\t$refCount\t$AD[0]\t$AD[1]\t$vaf\n";
	}
	elsif (length($ref) > 1 || length($alt) > 1){ # Everything Else
		my @alleles=split(",", $alt);
		my @vafs=split(",", $vaf);
		for(my $idx=0; $idx <=$#alleles; ++$idx){
			if(length($ref) > length($alleles[$idx])){ # deletion or block substitution
				my $head = substr($ref, 0, length ($alleles[$idx]));
				if ($head eq $alleles[$idx]) {
					my $ref_allele1 = substr ($ref, length ($alleles[$idx]));
					print $chr,"\t",$position+length($head),"\t",$position+length($ref)-1,"\t","$ref_allele1\t-\t$qual\t`$genotype\t$refCount\t$AD[0]\t",$AD[$idx + 1],"\t",$vafs[$idx],"\n";
				} else {
					print $chr,"\t",$position,"\t",$position+length($ref)-1,"\t",$ref,"\t",$alleles[$idx],"\t$qual\t`$genotype\t$refCount\t$AD[0]\t",$AD[$idx + 1],"\t",$vafs[$idx],"\n";
				}
			}
			elsif(length($alleles[$idx]) >= length ($ref)){
				my $head = substr ($alleles[$idx], 0, length ($ref));
				if ($head eq $ref) {
					my $mut = substr ($alleles[$idx], length ($ref));
					print $chr,"\t",$position+length($ref)-1,"\t",$position+length($ref)-1,"\t","-\t",$mut,"\t$qual\t`$genotype\t$refCount\t$AD[0]\t",$AD[$idx + 1],"\t",$vafs[$idx],"\n";
				} else {
					print $chr,"\t",$position,"\t",$position+length($ref)-1, "\t", $ref, "\t", $alleles[$idx], "\t$qual\t`$genotype\t$refCount\t$AD[0]\t",$AD[$idx + 1],"\t",$vafs[$idx],"\n";
				}
			}
			else{
				print "Unexpected Record:1 $chr, $position, $ref, $alt, $qual, $genotype, $refCount, $altCount, $vaf\n";
				die;
			}
		}
	}
	else{
		print "Unexpected Record:2 $chr, $position, $ref, $alt, $qual, $genotype, $refCount, $altCount, $vaf\n";
		die;
	}
}
