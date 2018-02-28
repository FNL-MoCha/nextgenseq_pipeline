#!/usr/local/bin/perl -sw


############################
my $patient=`echo $ARGV[0] |cut -f1 -d '\/'`;
chomp $patient;

my $list=join(" ", @ARGV);

if($list =~ /germline/i){
	my $file=getFile('germline');
	getAncestry('germline', $file);
	exit;
}
elsif($list =~ /ORIGINATOR/i){
	my $file=getFile('originator');
	getAncestry('originator', $file);
	exit;
}
else{
	getAncestry('consensus');
	exit;
}


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
	my ($source, @files) =@_;
	if($source =~/germline/ or $source =~/originator/){
		open(FH, $files[0]);
		while(<FH>){
			chomp;
			if($_ =~ /^#/){
				print "$_\tSource\n";
			}
			else{
				my @line = split(" ", $_);
				print "$patient\t".join("\t", @line[1..9])."\t$source\n";
			}
		}
		close FH;
	}
	else{
		print "#Sample ID\tPopulation label\tNumber of SNPs\tPredicted PC1\tPredicted PC2\tPredicted PC3\t% YRI ancestry\t% CEU ancestry\t% East Asian ancestry\t% Native American Ancestry\tSource\n";
		my @DATA;
		my $data;
		foreach $f(@ARGV){
			open(F, $f);
			while(my $row= <F>){
				chomp $row;
				if($row =~ /^#/){next;}
				my @line = split(" ", $row);
				push @$data, [@$_] for ([split(/ /, $row)]);
				#print $line[3];
			}
			close F;
		}
		push @$data, ["$patient\tPOP", map { my $col = $_; avg( map { $_->[$col] }  @$data) }  2..9];	
		for my $i (@$data){
			if(@$i[0] =~ /^$patient\tPOP/){
				print join ("\t", @$i)."\t$source"."\n";
			}
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
