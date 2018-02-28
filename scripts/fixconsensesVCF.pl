#!/usr/local/bin/perl -sw

open(FH, $ARGV[0]);
while(<FH>){
	chomp;
	if($_ =~/##/){
		print "$_\n";
		next;
	}
	elsif($_ =~ /^#CHROM/){
		print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$ARGV[1]\n";
		next;
	}
	else{
		print "$_\n";
	}
}
close FH;
