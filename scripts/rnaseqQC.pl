#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;
no if $] >= 5.017011, warnings => 'experimental::smartmatch';
use 5.010;
#local $SIG{__WARN__} = sub {my $message =shift; die $message;};
my $Patient  =$ARGV[0];
my $Library  =$ARGV[1];
my $Diagnosis=$ARGV[2];
my $pct_human=$ARGV[3];
chomp $pct_human;

#PF_BASES     PF_ALIGNED_BASES  RIBOSOMAL_BASES  CODING_BASES  UTR_BASES  INTRONIC_BASES  INTERGENIC_BASES  IGNORED_READS  CORRECT_STRAND_READS  INCORRECT_STRAND_READS  PCT_RIBOSOMAL_BASES  PCT_CODING_BASES  PCT_UTR_BASES  PCT_INTRONIC_BASES  PCT_INTERGENIC_BASES  PCT_MRNA_BASES  PCT_USABLE_BASES  PCT_CORRECT_STRAND_READS  MEDIAN_CV_COVERAGE  MEDIAN_5PRIME_BIAS  MEDIAN_3PRIME_BIAS  MEDIAN_5PRIME_TO_3PRIME_BIAS  SAMPLE  LIBRARY  READ_GROUP
#10472336775  10471765851       35636850         19039801      2419240    1370225         10413301463       0              0                     0                       0.003403             0.001818          0.000231       0.000131            0.994417              0.002049        0.002049          0   
#PF_BASES        PF_ALIGNED_BASES        RIBOSOMAL_BASES CODING_BASES    UTR_BASES       INTRONIC_BASES  INTERGENIC_BASES        IGNORED_READS   CORRECT_STRAND_READS    INCORRECT_STRAND_READS  NUM_R1_TRANSCRIPT_STRAND_READS  NUM_R2_TRANSCRIPT_STRAND_READS  NUM_UNEXPLAINED_READS   PCT_R1_TRANSCRIPT_STRAND_READS  PCT_R2_TRANSCRIPT_STRAND_READS  PCT_RIBOSOMAL_BASES     PCT_CODING_BASES        PCT_UTR_BASES   PCT_INTRONIC_BASES      PCT_INTERGENIC_BASES    PCT_MRNA_BASES  PCT_USABLE_BASES        PCT_CORRECT_STRAND_READS        MEDIAN_CV_COVERAGE      MEDIAN_5PRIME_BIAS      MEDIAN_3PRIME_BIAS      MEDIAN_5PRIME_TO_3PRIME_BIAS    SAMPLE  LIBRARY READ_GROUP
#9521173114      9520958285      17485747        7474294608      1767837707      213259218       48081287        0       0       0       163086  17587773        19107015        0.009187        0.990813        0.001837        0.785036        0.185679        0.022399        0.00505 0.970715        0.970693        0       1.218971        0.172508        0.232287        0.138931


print "#Patient\tLibrary\tDiagnosis\tTOTOAL_READS\tPCT_HUMAN\tGENES_ABOVE_0_TPM\tALIGNED_READS\tPCT_ALIGNED_READS\tPCT_ALIGNED_Q20_BASES\tPCT_RIBOSOMAL_BASES\tPCT_CODING_BASES\tPCT_UTR_BASES\tPCT_INTRONIC_BASES\tPCT_INTERGENIC_BASES\tPCT_MRNA_BASES\tPCT_USABLE_BASES\n";
my $totalReads =`grep "^Total Sequences" $ARGV[4]|cut -f2`; # fastqc 
chomp $totalReads;
$totalReads=$totalReads*2;

my $aligned_Reads = `grep "^PAIR" $ARGV[5] |cut -f6`;#RNASEQ.AlignmentSummaryMetrics.txt
chomp $aligned_Reads;

my $pct_mapped=sprintf ("%.2f", ($aligned_Reads/$totalReads)*100);

my $hq_bases=`grep "^PAIR" $ARGV[5] |cut -f10`;#RNASEQ.AlignmentSummaryMetrics.txt
chomp $hq_bases;

my $hq_bases_20=`grep "^PAIR" $ARGV[5] |cut -f11`;#RNASEQ.AlignmentSummaryMetrics.txt
chomp $hq_bases_20;
my $pct_hq_20=sprintf ("%.2f", (($hq_bases_20/$hq_bases)*100));
print "$Patient\t$Library\t$Diagnosis\t$totalReads\t$pct_human\t$ARGV[7]\t$aligned_Reads\t$pct_mapped\t$pct_hq_20";

open(FH, $ARGV[6]); # RNASEQ.RnaSeqMetrics.txt
my ($PCT_RIBOSOMAL_BASES,$PCT_CODING_BASES, $PCT_UTR_BASES, $PCT_INTRONIC_BASES, $PCT_INTERGENIC_BASES, $PCT_MRNA_BASES, $PCT_USABLE_BASES);
while(<FH>){
	chomp;
	next if $. <7 or $. >8;
	#next if $_=~ /^$/ or $_ =~ /^#/ or $. >9;
	my @a = split("\t", $_);
	if($_ =~/PF_BASES/ ){
		for(0..$#a){
			$PCT_RIBOSOMAL_BASES	= $_ if $a[$_] =~ "^PCT_RIBOSOMAL_BASES\$";
			$PCT_CODING_BASES 	= $_ if $a[$_] =~ "PCT_CODING_BASES";
			$PCT_UTR_BASES 		= $_ if $a[$_] =~ "PCT_UTR_BASES";
			$PCT_INTRONIC_BASES	= $_ if $a[$_] =~ "PCT_INTRONIC_BASES";
			$PCT_INTERGENIC_BASES	= $_ if $a[$_] =~ "PCT_INTERGENIC_BASES";
			$PCT_MRNA_BASES		= $_ if $a[$_] =~ "PCT_MRNA_BASES";
			$PCT_USABLE_BASES	= $_ if $a[$_] =~ "PCT_USABLE_BASES";
		}
		next;
	}
	foreach	my $idx($PCT_RIBOSOMAL_BASES, $PCT_CODING_BASES, $PCT_UTR_BASES, $PCT_INTRONIC_BASES , $PCT_INTERGENIC_BASES, $PCT_MRNA_BASES, $PCT_USABLE_BASES){
		$a[$idx] = $a[$idx]*100 if $a[$idx] =~ /\d+/;
		print "\t$a[$idx]";
	}
	print "\n";
}
close FH;
