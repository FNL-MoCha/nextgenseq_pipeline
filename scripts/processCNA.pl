#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {my $message =shift; die $message;};
open(SAM, $ARGV[0]);
open(HUGO, $ARGV[1]);
my %LIST;
while(<SAM>){
	chomp;
	my @t =split("\t", $_);
	my @exon=split(",", $t[3]);
	foreach my $ex(@exon){
		my ($gene, $e)=split("___", $ex);
		$LIST{"$gene"} =$t[4];
	}
}
#close SAM;
while(<HUGO>){
	chomp;
	if($_ =~ /^#version 2.4/ or $_ =~ /^Hugo_Symbol/){
		print "$ARGV[2]\n"; 
		next;
	}
	my @a = split ("\t", $_);
	if (exists $LIST{$a[0]}){
		#my $out=(2**$LIST{$a[0]})*2;
		my $out=sprintf("%.0f",((2**$LIST{$a[0]})*2));
		if($out ==0 ){
			$out = -2;	
		}
		elsif($out ==1){
			$out= -1;
		}
		elsif($out==2){
			$out=0;
		}
		elsif($out >2 and $out <5){
			$out=1;
		}
		else{
			$out=2;
		}
		print "$out\n";
	}
	else{
		print "0\n";
	}
}
close HUGO;
