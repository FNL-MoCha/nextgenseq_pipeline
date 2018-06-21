"""Find adjacent SNPs/MNPs in a VCF, and look in the BAM file to see if these
should be joined into a larger variant.

This script can only handle up to 2 internal compound ALT alleles for every
stretch of adjacent SNPs. This script is not designed to test all possible
combinations of internal joined SNPs for SNP stretches longer than 3.
This script also cannot handle multiallelic variants or indels.

"ALERT" will be added to VCF INFO columns for cases that are too complex.
 
 
 
########################################################################## 
 
 module load samtools
 python joinAdjacentSNPs.py -v 128128/20170910/128128~338-R~L43~WES/calls/128128~338-R~L43~WES.HC_DNASeq.raw.vcf -o test.vcf 1 128128/20170910/128128~338-R~L43~WES/128128~338-R~L43~WES.bwa.final.bam /data/MoCha/patidarr/ref/ucsc.hg19.2bit
 
 
 
"""
from __future__ import division
import copy
import operator
import subprocess
import sys


# We're splitting multiallelics now so should not find any multiallelic VCF lines anymore.
# If the adjacent variants contain some multiallelics, should add alert.
# If the adjacent variants are only multiallelics, should not add an alert.
def checkMultiallelic(variants):
    seenPos = []
    countSplit = 0
    countComma = 0
    for var in variants:
        if not [var[0], var[1]] in seenPos:
            seenPos.append([var[0], var[1]])
        if 'SPLITMULTIALLELIC' in var[7]:
            countSplit += 1
        elif ',' in var[4]:
            countComma += 1
    if countSplit == len(variants) and len(seenPos) == 1:
        # This set of variants came from a single original multiallelic variant. don't add ALERT to this set.
        multiallelic = 'All'
    elif countSplit > 0:
        # Not everything in this set came from a single original multiallelic variant.
        multiallelic = 'Some'
    elif countComma > 0:
        multiallelic = 'Some'
    else:
        multiallelic = 'None'
    return multiallelic


def checkIndels(variants):
    indelsPresent = False
    for var in variants:
        for alt in var[4].split(','):
            if len(var[3]) != len(alt):
                indelsPresent = True
                break
    return indelsPresent


def getReference(chrom, start, end, ref2bit):
    # For twoBit2Fa, coords need to be 0-based, left closed, right open.
    fixStart = int(start)-1
    fixEnd = int(end)-1
    coords = "{}:{}-{}".format(chrom, fixStart, fixEnd)
    step_get_ref = subprocess.Popen(['twoBitToFa', '-noMask', ref2bit + ":" + coords, 'stdout'],
                                    stdout=subprocess.PIPE)
    step_remove_header = subprocess.Popen(['grep', '-v', '>'],
                                          stdin=step_get_ref.stdout,
                                          stdout=subprocess.PIPE)
    step_get_ref.stdout.close()
    # Stuff from twoBitToFa should only ever be a few bases, since we are
    # looking for things that are "adjacent" so there shouldn't be multiple
    # lines in output, but set this up so that multiple lines won't break
    # things.
    ref = "".join(step_remove_header.communicate()[0].rstrip('\n\r').splitlines())
    print >> sys.stderr, "Intervening genomic sequence: %s" % ref
    return ref


def getJoinedAlleles(variants, ref2bit):
    # Initialize joined alleles
    refJoined = variants[0][3]
    altJoined = variants[0][4]
    for i in range(1, len(variants)):
        # Is there space between this and the previous variant?
        # If there is, add it to the joined alleles.
        currChrom = variants[i][0]
        currPos = int(variants[i][1])
        prevPos = int(variants[i-1][1])
        prevLength = len(variants[i-1][3])
        prevEnd = prevPos + prevLength - 1
        if currPos - prevEnd > 1:
            intervening_genomic_seq = getReference(currChrom, prevEnd+1,
                                                   currPos, ref2bit)
            refJoined += intervening_genomic_seq
            altJoined += intervening_genomic_seq
        refJoined += variants[i][3]
        altJoined += variants[i][4]
    return refJoined, altJoined


def correctVcfLine(variants, newRef, newAlt):
    if len(variants) > 1:
        correctedLine = copy.copy(variants[0])
        correctedLine[3] = newRef
        correctedLine[4] = newAlt
        correctedLine[7] = correctedLine[7].replace('SNP', 'MNP').replace('snp', 'mnp') + ';JOINED'
    else:
        correctedLine = variants[0]
    return correctedLine


def alertSet(variants):
    alertVars = []
    for var in variants:
        alertLine = copy.copy(var)
        alertLine[7] += ";ALERT"
        alertVars.append("\t".join(map(str, alertLine)))
    return alertVars


def examineAdjacent(adjacentVar, inBam, ref2bit):
    # Minimum number of reads to support an allele
    min_support = 5
    print >> sys.stderr, "Examining adjacent variants:\n%s" % "\n".join(map(str, adjacentVar))
    correctedVars = []
    if len(adjacentVar) > 3:
        print >> sys.stderr, "WARNING: This is a stretch of %s adjacent variants. This script is not designed to test all possible combinations of internal joined SNPs for stretches longer than 3." % len(adjacentVar)
    if checkIndels(adjacentVar) is True:
        print >> sys.stderr, "WARNING: Variant stretch includes indels."
        correctedVars = alertSet(adjacentVar)
    if checkMultiallelic(adjacentVar) == 'Some':
        print >> sys.stderr, "WARNING: Variant stretch includes multiallelic variants."
        correctedVars = alertSet(adjacentVar)
    elif checkMultiallelic(adjacentVar) == 'All':
        print >> sys.stderr, "WARNING: This variant set comes from a single original multiallelic variant and does not represent true adjacent SNPs/MNPs. Leaving variants as is."
        # correctedVars might have already been set if there are indels, so need to reset.
        correctedVars = []
        for var in adjacentVar:
            correctedVars.append("\t".join(map(str, var)))
    if len(correctedVars) == 0:
        refJoined, altJoined = getJoinedAlleles(adjacentVar, ref2bit)
        print >> sys.stderr, "Joined REF allele is %s" % refJoined
        print >> sys.stderr, "Joined ALT allele is %s" % altJoined
        # coordinates of joined region
        chrom = adjacentVar[0][0]
        start = int(adjacentVar[0][1])
        end = int(adjacentVar[-1][1]) + len(adjacentVar[-1][4]) - 1
        coords = chrom + ":" + str(start) + "-" + str(end)

        # get counts of alleles from BAM file.
        observed_alleles = {}
        print >> sys.stderr, "Samtools view -f 2 -F 1792 %s" % coords
        # don't count duplicates.
        step_samtools = subprocess.Popen(['samtools', 'view', '-f', '2', '-F', '1792', inBam, coords], stdout=subprocess.PIPE)
        step_cut = subprocess.Popen(['cut', '-f', '4,6,10'], stdin=step_samtools.stdout, stdout=subprocess.PIPE)
        step_samtools.stdout.close()
        bamData=step_cut.communicate()[0].rstrip('\n\r').split('\n')
        if len(bamData[0]) > 0:
            print >> sys.stderr, "Total reads at this position: %s" % len(bamData)
            fullMatchCount = 0
            spanFullLocusCount = 0
            for read in iter(bamData):
                pos, cigar, seq = read.split("\t")
                fullMatchCigar = str(len(seq)) + "M"
                if cigar == fullMatchCigar:
                    fullMatchCount += 1
                    # VCF is 1-based but python is 0-based
                    left = start - int(pos)
                    right = end - int(pos) + 1
                    allele = seq[left:right]
                    if len(allele) == len(altJoined):
                        spanFullLocusCount += 1
                        if allele in observed_alleles:
                            observed_alleles[allele] += 1
                        else:
                            observed_alleles[allele] = 1
            print >> sys.stderr, "Reads with full match CIGAR (%s): %s" % (fullMatchCigar, fullMatchCount)
            print >> sys.stderr, "Reads with full match CIGAR that span the full locus: %s (%.3f%% of all reads)" % (spanFullLocusCount, 100.0*spanFullLocusCount/len(bamData))
        # Sort alleles by abundance. This produces a list of duples. For every duple, the first element is the allele
        # and the second is the count.
        observed_alleles_sorted = sorted(observed_alleles.items(), key=operator.itemgetter(1), reverse=True)

        # Require at least 40% usable reads
        if (len(observed_alleles.keys()) >= 1
            and observed_alleles_sorted[0][1] >= min_support
            and 1.0*spanFullLocusCount/len(bamData) >= 0.4):
            # For now, script handles at most 2 major alleles
            major_alleles = [observed_alleles_sorted[0][0]]
            # major_allele_read_counts = int(observed_alleles_sorted[0][1])
            if (len(observed_alleles_sorted) > 1
                and (observed_alleles_sorted[1][1] >= min_support
                     or observed_alleles_sorted[1][1] >= 0.6*observed_alleles_sorted[0][1])):
                major_alleles.append(observed_alleles_sorted[1][0])
                # major_allele_read_counts += int(observed_alleles_sorted[1][1])
            # Convert each major allele into a binary array showing wheter it
            # matches REF or ALT at each variant position
            abBinary = []
            for abAllele in major_alleles:
                convBin = []
                for k in range(len(refJoined)):
                    if refJoined[k:k+1] == altJoined[k:k+1]:
                        # This is some intervening genomic sequence between the adjacent SNPs.
                        continue
                    elif abAllele[k:k+1] == refJoined[k:k+1]:
                        convBin.append(0)
                    else:
                        convBin.append(1)
                abBinary.append(convBin)
            print >> sys.stderr, "Major alleles: %s" % major_alleles
            print >> sys.stderr, "Major alleles binary arrays: %s" % abBinary

            # How many reads aren't from the major allele(s)/the reference allele?
            # Only look at variant positions, not intervening genomic sequence
            nonref_minor_allele_read_counts = 0
            for y in range(len(major_alleles), len(observed_alleles_sorted)):
                currAllele = observed_alleles_sorted[y][0]
                nonRef = 0
                for k in range(len(refJoined)):
                    if refJoined[k:k+1] == altJoined[k:k+1]:
                        # This is some intervening genomic sequence between the adjacent SNPs.
                        continue
                    elif refJoined[k:k+1] == currAllele[k:k+1]:
                        continue
                    else:
                        nonRef+=1
                if nonRef > 0:
                    nonref_minor_allele_read_counts += int(observed_alleles_sorted[y][1])
            nonref_minor_allele_percent = 100*float(nonref_minor_allele_read_counts)/float(spanFullLocusCount)

            for obs in observed_alleles_sorted:
                print >> sys.stderr, "Found allele %s, count %s" % (obs[0], obs[1])
            print >> sys.stderr, "Reads from non-reference minor alleles are %.3f%% of total" % nonref_minor_allele_percent

            if nonref_minor_allele_percent <= 5:
                if ((len(major_alleles) == 2 and altJoined in major_alleles and refJoined in major_alleles)
                    or (len(major_alleles) == 1 and altJoined in major_alleles)
                    or (len(major_alleles) == 2 and [0]*len(adjacentVar) in abBinary and [1]*len(adjacentVar) in abBinary)
                    or (len(major_alleles) == 1 and [1]*len(adjacentVar) in abBinary)):
                    # Join the variants
                    print >> sys.stderr, "Joining variants:\n%s" % "\n".join(map(str, adjacentVar))
                    joinedVars = correctVcfLine(adjacentVar, refJoined, altJoined)
                    correctedVars.append("\t".join(map(str, joinedVars)))
                else:
                    if len(adjacentVar) == 2 or len(major_alleles) == 1:
                        print >> sys.stderr, "Data does not support joining variants. Leaving them separate"
                        for origVar in adjacentVar:
                            correctedVars.append("\t".join(map(str, origVar)))
                    else:
                        # Convert major allele sequence to 0,1 for ref/alt.
                        print >> sys.stderr, "Full combined variant %s not found. Looking at internal subsets" % altJoined

                        # Check that at every position, one allele is REF and the other is ALT.
                        # Sum of abBinary[0][q] and abBinary[1][q] should be 1 for all q.
                        if all(abbin0 + abbin1 == 1
                               for abbin0, abbin1 in zip(abBinary[0], abBinary[1])):
                            # Sort orig SNPs into two sets to join.
                            first_set_to_join = []
                            second_set_to_join = []
                            # This puts the allele with the ALT SNP at first position first in list.
                            # Use this to sort orig SNPs into alleles.
                            abBinary.sort(reverse=True)
                            alleleMap = abBinary[0]
                            for z, amz in enumerate(alleleMap):
                                if amz == 1:
                                    first_set_to_join.append(adjacentVar[z])
                                else:
                                    second_set_to_join.append(adjacentVar[z])
                            print >> sys.stderr, "Joining variants."
                            print >> sys.stderr, "Set1:\n%s" % "\n".join(map(str, first_set_to_join))
                            refJoinedSet1, altJoinedSet1 = getJoinedAlleles(first_set_to_join, ref2bit)
                            joinedVarsSet1 = correctVcfLine(first_set_to_join, refJoinedSet1, altJoinedSet1)
                            print >> sys.stderr, "Set2:\n%s" % "\n".join(map(str, second_set_to_join))
                            refJoinedSet2, altJoinedSet2 = getJoinedAlleles(second_set_to_join, ref2bit)
                            joinedVarsSet2 = correctVcfLine(second_set_to_join, refJoinedSet2, altJoinedSet2)
                            correctedVars.append("\t".join(map(str, joinedVarsSet1)))
                            correctedVars.append("\t".join(map(str, joinedVarsSet2)))
                        else:
                            # this case might happen if one or more (but not all) SNPs is homozygous ALT, in which case that SNP would segregate with all other variants.
                            print >> sys.stderr, "WARNING: Two major alleles are not mutually exclusive. This may be a multiallelic locus."
                            correctedVars = alertSet(adjacentVar)
            else:
                if len(adjacentVar) == 2:
                    print >> sys.stderr, "Data does not support joining variants. Leaving them separate"
                    for origVar in adjacentVar:
                        correctedVars.append("\t".join(map(str, origVar)))
                else:
                    print >> sys.stderr, "WARNING: There may be more than two internal subsets of variants."
                    correctedVars = alertSet(adjacentVar)
        else:
            print >> sys.stderr, "WARNING: Not enough reads to process variants. Poor mapping region."
            # consider: don't print ALERT here?  If there aren't enough reads this is just a bad mapping region, and looking in the BAM file won't help.
            correctedVars = alertSet(adjacentVar)
    return correctedVars


def peek_line(f):
    currplace = f.tell()
    nextline = f.readline()
    f.seek(currplace)
    return nextline


def main(distance, inFile, inBam, ref2bit, outFile):
    adjacentVar = []
    #for line in inFile:
    while True:
        line = inFile.readline()
        line = line.rstrip('\n\r')
        if line.startswith("##INFO=<ID=VT,"):
            # Rewrite the VT header tag
            print >> outFile, '##INFO=<ID=VT,Number=1,Type=String,Description="Variant type, can be SNP, MNP, INS or DEL">'
        elif line.startswith("##"):
            # Copy the rest of the VCF header
            print >> outFile, line
        elif line.startswith("#"):
            # Add custom tags to the end of the VCF header, just before #CHROM...
            print >> outFile, '##joinAdjacentSNPs_command=%s' % " ".join(map(str, sys.argv))
            print >> outFile, '##INFO=<ID=JOINED,Number=0,Type=Flag,Description="Variant is the result of joining adjacent SNPs/MNPs. Please use caution when interpreting sample-specific values.">'
            print >> outFile, '##INFO=<ID=ALERT,Number=0,Type=Flag,Description="Variant is adjacent to other variants, but case was too complex for joinAdjacentSNPs to handle. Please review.">'
            print >> outFile, line
        else:
            fields = line.split('\t')
            chrom = fields[0]
            pos = int(fields[1])
            end = pos + len(fields[3]) - 1

            # Check the next line to see if the variant is adjacent.
            nextLine = peek_line(inFile)
            if nextLine:
                nextFields = nextLine.rstrip('\n\r').split('\t')
                nextChrom = nextFields[0]
                nextPos = int(nextFields[1])
                if chrom == nextChrom and end >= nextPos - distance:
                    # This is adjacent to the next variant
                    if len(adjacentVar) == 0:
                        print >> sys.stderr, "----- Found adjacent variants -----"
                        adjacentVar.append(fields)
                    adjacentVar.append(nextFields)
                else:
                    if len(adjacentVar) == 0:
                        print >> outFile, line
                    else:
                        # Deal with adjacent variants
                        correctedVars = examineAdjacent(adjacentVar, inBam, ref2bit)
                        print >> outFile, "\n".join(map(str, correctedVars))
                        # Reset
                        adjacentVar = []
            else:
		break
                # EOF
                if len(adjacentVar) == 0:
                    print >> outFile, line
                else:
                    # Deal with adjacent variants
                    correctedVars = examineAdjacent(adjacentVar, inBam, ref2bit)
                    print >> outFile, "\n".join(map(str, correctedVars))


if __name__ == '__main__':
    import argparse
    AP = argparse.ArgumentParser(description=__doc__)
    AP.add_argument('distance', type=int,
                    help="""Distance in bp between SNVs/MNVs to be considered
                    'adjacent'. In other words, POS2-POS1 must be <= [distance]
                    for the SNPs to be adjacent.""")
    AP.add_argument('bam',
                    help="""BAM file of mapped reads from which the VCF calls
                    were made.""")
    AP.add_argument('ref2bit',
                    help="Reference genome in UCSC 2bit format.")
    AP.add_argument('-v', '--vcf',
                    type=argparse.FileType('r'), default=sys.stdin,
                    help="Variants in VCF format, sorted by genomic position.")
    AP.add_argument('-o', '--output',
                    type=argparse.FileType('w'), default=sys.stdout,
                    help="Output VCF filename. [Default: standard output]")
    args = AP.parse_args()
    main(args.distance, args.vcf, args.bam, args.ref2bit, args.output)
