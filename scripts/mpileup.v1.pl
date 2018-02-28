#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use 5.010;
local $SIG{__WARN__} = sub {
        my $message =shift;
        die $message;
};
my $FILE = $ARGV[0]; # Vcf like file Name
my $BAM  = $ARGV[1];
unless (open(FH, "$FILE")){
	print STDERR "Can not find the file $FILE\n";
}
my @a = split(/[.]/, basename($BAM));
while(<FH>){
	chomp;
	my $line = $_;
	my @d = split("\t", $_);
	if($_ =~ /^Chr/ or $_ =~ /^#/){
		print "$a[0]\n";
		next;
	}
	my $rnaseq;		
	$rnaseq = `samtools mpileup -d1000000000000 -r $d[0]:$d[1]-$d[1] "$BAM" 2>/dev/null |cut -f 3-5`;
	chomp $rnaseq;
	if($rnaseq =~ /\d/){ # Have coverage 
		my @s = split("\t", $rnaseq);
		if ($s[1] <1){
			print "0\n";
		}
		else{
			$s[2] =uc($s[2]);
			my $bases = join '', sort, split //, $s[2];
			my $Ref = () = ($bases =~ /$d[3]/g); # 4th Column is the Ref Column
				my $Alt = () = ($bases =~ /$d[4]/g); # 5th Column is the Alt Column 
				my $totalRef = $Ref + $Alt;
			if($totalRef >=1){ 
				if($Alt >=1){
					my $vaf = ($Alt/$totalRef);
					print "$vaf\n";
				}
				else{
					print "-1\n";
				}
			}
			else{
				print "-1\n";
			}
		}
	}
	else{ # no output from mpileup no COVERAGE
		print "-1\n";
	}
}
close FH;
