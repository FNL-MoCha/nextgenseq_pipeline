#!/usr/bin/perl
use strict;
use warnings;
use 5.010;
local $SIG{__WARN__} = sub {my $message =shift; die $message;};
open(SAM, $ARGV[0]);
while(<SAM>){
	chomp;
	if($_ =~ /^Hugo_Symbol/){
		print "Hugo_Symbol\tEntrez_Gene_Id\t$ARGV[1]\n";
		next;
	}
	my ($hugo, $entrez, @sample)= split("\t", $_);
	@sample=sort{$a<=>$b}@sample;
	if($sample[-1] ==0){
	
	}
	elsif($sample[0] == $sample[-1]){
		print "$hugo\t$entrez\t$sample[0]\n";
	}
}
close SAM;
