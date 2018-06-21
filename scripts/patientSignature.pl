#!/usr/local/bin/perl -sw


############################
my $patient=`echo $ARGV[0] |cut -f1 -d '\/'`;
chomp $patient;

my $list=join(" ", @ARGV);

getAncestry();


#############################
sub getFile{
	my ($type) = @_;
	foreach my $file(@ARGV){
		if($file =~ /$type/i){
			return ($file);
		}
	}
}
sub getAncestry{
	my $data;
	print "Sample\tMean\tSignature.1\tSignature.10\tSignature.11\tSignature.12\tSignature.13\tSignature.14\tSignature.15\tSignature.16\tSignature.17\tSignature.18\tSignature.19\tSignature.2\tSignature.20\tSignature.21\tSignature.22\tSignature.23\tSignature.24\tSignature.25\tSignature.26\tSignature.27\tSignature.28\tSignature.29\tSignature.3\tSignature.30\tSignature.4\tSignature.5\tSignature.6\tSignature.7\tSignature.8\tSignature.9\n";
	foreach $f(@ARGV){
		open(F, $f);
		while(my $row= <F>){
			chomp $row;
			if($row =~ /^Sample/){next;}
			my @line = split(" ", $row);
			push @$data, [@$_] for ([split("\t", $row)]);
		}
		close F;
	}
	push @$data, ["$patient\tMean", map { my $col = $_; avg( map { $_->[$col] }  @$data) }  1..30];	
	for my $i (@$data){
		if(@$i[0] =~ /^$patient\tMean/){
			print join ("\t", @$i)."\n";
		}
	}	
}
sub avg {
	my (@Arr) = @_;
	my $mean =0;
	my $totalelements = @Arr;
	foreach my $ele(@Arr){
		$mean = $mean + $ele;
	}
	if ($mean <=0){
		return ($mean);
	}
	else{
		$mean = sprintf("%.3f", $mean/$totalelements);
		return ($mean);
	}
}
