#!/bin/sh

patient=$1
rna=`../patidarr/ngs_pipeline/scripts/addAnnotations2vcf.pl /data/MoCha/patidarr/ref/aMOI.bed $patient/20170910/$patient/db/${patient}.rnaseq |cut -f 1-6 |sort |uniq -c |sed -e 's/^[ \t]*//' |sed -e 's/ /\t/' |sort -r`

if [ -f $patient/20170910/$patient/db/${patient}.variants ]
then
	dna=`../patidarr/ngs_pipeline/scripts/addAnnotations2vcf.pl /data/MoCha/patidarr/ref/aMOI.bed $patient/20170910/$patient/db/${patient}.variants |grep -v fibroblast|cut -f 1-6 |sort |uniq -c|sed -e 's/^[ \t]*//' |sed -e 's/ /\t/'|sort -r `
else
	dna=`../patidarr/ngs_pipeline/scripts/addAnnotations2vcf.pl /data/MoCha/patidarr/ref/aMOI.bed $patient/20170910/$patient/db/${patient}.germline |grep -v germline |grep -v fibroblast|cut -f 1-6 |sort |uniq -c|sed -e 's/^[ \t]*//' |sed -e 's/ /\t/'|sort -r `
fi

WES=`/bin/ls -1d $patient/20170910/*WES |grep -v germline |grep -v fibroblast|wc -l`
RNASEQ=`/bin/ls -1d $patient/20170910/*RNASEQ |wc -l`
echo -e "\n\n\e[1;31mPatient ID:\t$patient\n \e[0m"
echo -e "\e[1;31m$WES\e[0m\tTotal number of samples with Exome Seq"
echo -e "\n$dna\n\n"
echo -e "\e[1;31m$RNASEQ\e[0m\tTotal number of samples with RNASeq"
echo -e "\n$rna\n\n"

if [ -f $patient/20170910/$patient/db/${patient}.variants ]
then
	echo "../patidarr/ngs_pipeline/scripts/addAnnotations2vcf.pl /data/MoCha/patidarr/ref/aMOI.bed $patient/20170910/$patient/db/${patient}.variants |cut -f 1-7,13,14"
else
	echo "../patidarr/ngs_pipeline/scripts/addAnnotations2vcf.pl /data/MoCha/patidarr/ref/aMOI.bed $patient/20170910/$patient/db/${patient}.germline |cut -f 1-7,13,14"
fi
echo "../patidarr/ngs_pipeline/scripts/addAnnotations2vcf.pl /data/MoCha/patidarr/ref/aMOI.bed $patient/20170910/$patient/db/${patient}.rnaseq |cut -f 1-7,13,14"

echo "${patient}/20170910/variantHeatmap/${patient}.matrix"
