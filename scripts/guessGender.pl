#!/usr/local/bin/perl -s
no if $] >= 5.018, warnings => "experimental::smartmatch";
open(FH, $ARGV[0]);
my $het=0;
my $homo=0;
while(<FH>){
	chomp;
	my @a=split("\t", $_);
	if ($a[0] =~ /chrX/ and $a[2] >0.2){
		if ($a[1] ~~ [60001..2699520] or $a[1] ~~ [154931044..155260560]){
		}
		else{
			if ($a[1] >0.8){
				$homo++;
			}
			else{
				$het++;
			}
		}
	}
}
my $result = $het/($het+$homo);

if ($result <0.3){
	print "Male\n";
}
elsif($result >0.7){
	print "Female";
}
else{
	print "$het\t$homo\t$result\n";
}
