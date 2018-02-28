#!/usr/bin/perl -sw

#mutation_id     ref_counts      var_counts      normal_cn       minor_cn        major_cn        variant_case    variant_freq    genotype
#NA12156:BB:chr1:70820008        2538    91      2       0       2       NA12156 0.034613921643210345    BB
#NA18507:BB:chr22:32587251       544     1276    2       0       2       NA18507 0.701098901098901       BB
#NA19240:BB:chr5:135228165       2470    20      2       0       2       NA19240 0.008032128514056224    BB
print "mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn\tmajor_cn\tvariant_freq\n";
my $gender=`cat $ARGV[2]`;
chomp $gender;
my $normalCN=2;
open(FH, $ARGV[0]); # Somatic file
while(<FH>){
	chomp;
	my @line=split("\t", $_);
	if($_=~ /^Chr/){
		next;
	}
	my ($minor, $major)= getCN($line[0], $line[1]);
	if ($line[0] =~ /chrX/ and $gender =~/Male/){
		$normalCN=1;
	}
	if($line[0]=~ /chrY/ and $gender =~/Male/){
		$normalCN=1;
	}
	if($line[0] =~ /chrY/ and $gender=~/Female/){
		next;
	}
	print "$line[6]:$line[0]:$line[1]\t$line[89]\t$line[91]\t$normalCN\t$minor\t$major\t$line[92]\n";
}
close FH;

sub getCN{
	my ($chr, $pos) =@_;
	my $minor=0;
	my $major=2;
	open(FH1, $ARGV[1]);
	while(<FH1>){
		chomp;
		my @a =split("\t", $_);
		if ($chr eq $a[0] and $pos >=$a[1] and $pos <=$a[2]){
			if ($a[11] <$a[12]){
				$minor=$a[11];
				$major=$a[12];
			}
			else{
				$minor=$a[12];
				$major=$a[11];
			}
		}
	}
	return($minor, $major);
	close FH1;
}
