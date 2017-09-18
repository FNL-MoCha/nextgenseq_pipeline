#!/usr/local/bin/perl -ws
local $SIG{__WARN__} = sub {
        my $message =shift;
        die $message;
};
# Remove block substitutions longer ten 3 bp and failed calls
open(FH, $ARGV[0]);
while(<FH>){
	chomp;
	if ($_ =~ /^#/){
		print "$_\n";
		next;
	}
	my @line = split("\t", $_);
	if ($line[6] =~ /PASS/){
		if (length($line[3]) eq length($line[4]) and length($line[3]) <=3){
			print "$_\n";
		}
		elsif(length($line[3]) ne length($line[4])){	
			print "$_\n";	
		}
	}
}
