RNASEQ_BAM =[]
RNASEQ_FUSION1=[]
RNASEQ_FUSION =[]
EXPRESSION=[]
RNA_CALLS =[]
SUB2RNA = {}
SUB_RNASEQ=[]
SUB_FUSION={}
ActionableFiles1 =[]
for subject,samples in config['RNASeq'].items():
	SUB_RNASEQ.append(subject)
	for sample in samples:
		SUB2RNA[sample]=subject
for subject  in config['RNASeq'].keys():
	DBFiles         +=[subject+"/"+subject+"/db/"+subject+".rnaseq"]
	ActionableFiles +=[subject+ACT_DIR+subject+".fusion.actionable.txt"]
	ActionableFiles +=[subject+ACT_DIR+subject+".rnaseq.actionable.txt"]
	for sample in config['RNASeq'][subject]:
		RNASEQ_BAM += [subject+"/"+sample+"/"+sample+".star.final.bam"]
		RNASEQ_BAM += [subject+"/"+sample+"/"+sample+".tophat.final.bam"]
		ALL_FASTQC += [subject+"/"+sample+"/qc/fastqc/"+sample+"_R2_fastqc.html"]
		RNASEQ_FUSION += [subject+"/"+sample+"/fusion/tophat-fusion.txt"]
		RNASEQ_FUSION += [subject+"/"+sample+"/fusion/fusion-catcher.txt"]
		RNASEQ_FUSION += [subject+"/"+sample+"/fusion/defuse.filtered.txt"]
		ALL_QC      += [subject+"/"+sample+"/qc/"+sample+".star.flagstat.txt"]
		ALL_QC      += [subject+"/"+sample+"/qc/"+sample+".star.hotspot.depth"]
		ALL_QC      += [subject+"/"+sample+"/qc/"+sample+".star.gt"]
		add_to_SUBJECT_ANNO(subject, "rnaseq", [subject+"/"+sample+"/calls/"+sample+".hapCaller.annotated.txt"])
		EXPRESSION += [subject+"/"+sample+"/exonExp_UCSC/"+sample+".exonExpression.UCSC.txt"]
		EXPRESSION += [subject+"/"+sample+"/exonExp_ENS/"+sample+".exonExpression.ENS.txt"]
		for gtf in config['GTF']:
			EXPRESSION += [subject+"/"+sample+"/cufflinks_"+gtf+"/genes.fpkm_tracking_log2"]
	RNA_CALLS  += ["{subject}/{sample}/calls/{sample}.hapCaller.raw.vcf".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	SUB_FUSION[subject] = ["{subject}/{sample}/fusion/{sample}.actionable.fusion.txt".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	if subject in SUBJECT_VCFS:
		SUBJECT_VCFS[subject] += ["{subject}/{sample}/calls/{sample}.hapCaller.snpEff.txt".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	else:
		SUBJECT_VCFS[subject] = []
		SUBJECT_VCFS[subject] += ["{subject}/{sample}/calls/{sample}.hapCaller.snpEff.txt".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	if subject in SUB_HOT:
		SUB_HOT[subject] += ["{subject}/{sample}/qc/{sample}.star.hotspot.depth".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_LOH[subject] += ["{subject}/{sample}/qc/{sample}.star.loh".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_COV[subject] += ["{subject}/{sample}/qc/{sample}.star.coverage.txt".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_GT[subject]  += ["{subject}/{sample}/qc/{sample}.star.gt".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.star.final.bam".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.star.final.bam.tdf".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.tophat.final.bam".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.tophat.final.bam.tdf".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
	else:
		SUB_HOT[subject] = []
		SUB_LOH[subject] = []
		SUB_COV[subject] = []
		SUB_GT[subject]  = []
		SUB_IGV[subject] = []
		SUB_HOT[subject] += ["{subject}/{sample}/qc/{sample}.star.hotspot.depth".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_LOH[subject] += ["{subject}/{sample}/qc/{sample}.star.loh".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_COV[subject] += ["{subject}/{sample}/qc/{sample}.star.coverage.txt".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_GT[subject]  += ["{subject}/{sample}/qc/{sample}.star.gt".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.star.final.bam".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.star.final.bam.tdf".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.tophat.final.bam".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
		SUB_IGV[subject] += ["{subject}/{sample}/{sample}.tophat.final.bam.tdf".format(subject=SUB2RNA[s], sample=s) for s in config['RNASeq'][subject]]
############
#       RNASeq All
############
rule RNASeq:
	input:
		RNASEQ_BAM,
		EXPRESSION,
		RNA_CALLS,
		SUBJECT_VCFS,
		expand("{subject}"+ACT_DIR+"{subject}.fusion.actionable.txt", subject=config['RNASeq']),
		expand("{subject}/qc/{subject}.hotspot_coverage.png", subject=config['RNASeq']),
		expand("{subject}/qc/{subject}.coveragePlot.png", subject=config['RNASeq']),
		expand("{subject}/qc/{subject}.circos.png", subject=config['RNASeq']),
	output: temp("rnaseqDone")
	params:
		rulename  = "RNASeq_final",
		batch     = config[config['host']]["job_default"]

	shell: """
	#######################
	echo "This is the end of RNASeq" >{output}
	#######################
	"""
############
#	Tophat
############
rule tophat:
	input:
		R=lambda wildcards: FQ[wildcards.sample],
#		R1=DATA_DIR + "/{sample}/{sample}_R1.fastq.gz",
#		R2=DATA_DIR + "/{sample}/{sample}_R2.fastq.gz",
	output:
		"{base}/{sample}/tophat_{sample}/accepted_hits.bam",
		"{base}/{sample}/tophat_{sample}/accepted_hits.bam.bai"
	version: config["tophat"]
	params:
		rulename  = "tophat",
		batch     = config[config['host']]["job_tophat"],
		ref=config['Bowtie2Index']
	shell: """
	#######################
	module load tophat/{version}
	module load samtools
	tophat -p ${{THREADS}} -o ${{LOCAL}} --keep-fasta-order --rg-id {wildcards.sample} --no-coverage-search --rg-sample {wildcards.sample} --rg-library {wildcards.sample} --rg-platform ILLUMINA --fusion-search --fusion-min-dist 100000 --mate-inner-dist 84 --mate-std-dev 74 {params.ref} {input.R[0]} {input.R[1]}
	cp -rf ${{LOCAL}}/* {wildcards.base}/{wildcards.sample}/tophat_{wildcards.sample}/
	samtools index {wildcards.base}/{wildcards.sample}/tophat_{wildcards.sample}/accepted_hits.bam
	#######################
	"""
############
#       Link Tophat bam file
############
rule symlink_tophatBam:
	input:
		bam="{base}/{sample}/tophat_{sample}/accepted_hits.bam",
		bai="{base}/{sample}/tophat_{sample}/accepted_hits.bam.bai",
	output:
		bam="{base}/{sample}/{sample}.tophat.final.bam",
		bai="{base}/{sample}/{sample}.tophat.final.bam.bai"
	params:
		rulename  = "tophat_link",
		batch     = config[config['host']]["job_default"]
	shell: """
	#######################
	cd {wildcards.base}/{wildcards.sample}/
	ln -sf tophat_{wildcards.sample}/accepted_hits.bam {wildcards.sample}.tophat.final.bam
	ln -sf tophat_{wildcards.sample}/accepted_hits.bam.bai {wildcards.sample}.tophat.final.bam.bai
	#######################
	"""
############
#       Tophat-fusion
############
rule tophat_fusion:
	input:
		"{base}/{sample}/tophat_{sample}/accepted_hits.bam",
		"{base}/{sample}/tophat_{sample}/accepted_hits.bam.bai"
	output:
		"{base}/{sample}/tophatfusion_out/result.txt",
		"{base}/{sample}/fusion/tophat-fusion.txt"
	version: config["tophat"]
	params:
		rulename = "tp",
		batch    =config[config['host']]['job_tophatPost'],
		ref      =config['BowtieIndex'],
		bowtie   =config['bowtie'],
		tp_ref   =config['tophat_post_ref']
	shell: """
	#######################
	module load tophat/{version}
	module load bowtie/{params.bowtie}
	module load blast
	cd {wildcards.base}/{wildcards.sample}/
	rm -f blast ensGene.txt ensGtp.txt mcl refGene.txt
	ln -s {params.tp_ref}/* .
	tophat-fusion-post -p ${{THREADS}} --num-fusion-pairs 1 {params.ref}
	rm blast ensGene.txt ensGtp.txt mcl refGene.txt
	sed -i  '1s/^/Sample\\tGene_left\\tChr_left\\tCoordinate_left\\tGene_right\\tChr_right\\tCoordinate_right\\t#SpanningReads\\t#SpanningMatePairs\\t#SpanningMateEndOfPair\\tScore\\n/' tophatfusion_out/result.txt
	ln -sf ../tophatfusion_out/result.html fusion/tophat-fusion.html
	ln -sf ../tophatfusion_out/result.txt  fusion/tophat-fusion.txt
	#######################
	"""
############
#       Cufflinks
############
rule cufflinks:
	input:
		bam="{base}/{sample}/tophat_{sample}/accepted_hits.bam",
		bai="{base}/{sample}/tophat_{sample}/accepted_hits.bam.bai",
		convertor =NGS_PIPELINE + "/scripts/transformlog2_FPKM.py",
		ref=lambda wildcards: config['GTF'][wildcards.gtf]
	output:
		"{base}/{sample}/cufflinks_{gtf}/genes.fpkm_tracking_log2",
		"{base}/{sample}/cufflinks_{gtf}/isoforms.fpkm_tracking_log2"
	version: config['cufflinks']
	params:
		rulename   = "cuff",
		batch      =config[config['host']]['job_cufflinks']
	shell: """
	#######################
	module load cufflinks/{version}
	cufflinks -p ${{THREADS}} -G {input.ref} --max-bundle-frags 8000000000000 --max-bundle-length 10000000 -o {wildcards.base}/{wildcards.sample}/cufflinks_{wildcards.gtf} {input.bam}
	/usr/bin/python {input.convertor} genes {wildcards.base}/{wildcards.sample}/cufflinks_{wildcards.gtf}/genes.fpkm_tracking {wildcards.base}/{wildcards.sample}/cufflinks_{wildcards.gtf}/genes.fpkm_tracking_log2
	python {input.convertor} isoforms {wildcards.base}/{wildcards.sample}/cufflinks_{wildcards.gtf}/isoforms.fpkm_tracking {wildcards.base}/{wildcards.sample}/cufflinks_{wildcards.gtf}/isoforms.fpkm_tracking_log2
	#######################
	"""
############
#       Exon Expression
############
rule exon_exp:
	input:
		bam="{base}/{sample}/tophat_{sample}/accepted_hits.bam",
		bai="{base}/{sample}/tophat_{sample}/accepted_hits.bam.bai",
		convertor =NGS_PIPELINE + "/scripts/exon_exp.sh",
		ucsc=config["exon_Bed_UCSC"],
		ens=config["exon_Bed_ENS"],
	output:
		ucsc="{base}/{sample}/exonExp_UCSC/{sample}.exonExpression.UCSC.txt",
		ens="{base}/{sample}/exonExp_ENS/{sample}.exonExpression.ENS.txt"
	version: config['samtools']
	params:
		rulename   = "exonExp",
		batch      =config[config['host']]['job_exonExp']
	shell: """
	#######################
	module load samtools/{version}
	totalReads=`samtools flagstat {input.bam} |head -1 | sed 's/\s/\\t/g' | cut -f1`

        split -d -l 12000 {input.ucsc} ${{LOCAL}}/ucsc
        for file in ${{LOCAL}}/ucsc*
        do
		sh {input.convertor} ${{totalReads}} ${{file}} {input.bam} ${{file}}.out &
        done
        wait;
        cat ${{LOCAL}}/ucsc*.out >{output.ucsc}
	rm -rf ${{LOCAL}}/*

	split -d -l 15000 {input.ens} ${{LOCAL}}/ens
	for file in ${{LOCAL}}/ens*
	do
		sh {input.convertor} ${{totalReads}} ${{file}} {input.bam} ${{file}}.out &
	done
	wait
	cat ${{LOCAL}}/ens*.out >{output.ens}
	#######################
        """
############
#       Fusioncatcher
############
rule fusioncatcher:
	input:
		R=lambda wildcards: FQ[wildcards.sample]
#		R1=DATA_DIR + "/{sample}/{sample}_R1.fastq.gz",
#		R2=DATA_DIR + "/{sample}/{sample}_R2.fastq.gz",
	output:
		"{base}/{sample}/fusion/fusion-catcher.txt"
	version: config['fusioncatcher']
	params:
		rulename = "fc",
		batch    = config[config['host']]['job_fusioncatch']
	shell: """
	#######################
	module load fusioncatcher/{version}
	fusioncatcher -p ${{THREADS}} -i {input.R[0]},{input.R[1]} -o ${{LOCAL}}/
	cp ${{LOCAL}}/final-list_candidate-fusion-genes.GRCh37.txt {wildcards.base}/{wildcards.sample}/fusion/fusion-catcher.txt
	#######################
	"""
############
#       deFuse
############
rule deFuse:
	input:
		R=lambda wildcards: FQ[wildcards.sample],
#		R1=DATA_DIR + "/{sample}/{sample}_R1.fastq.gz",
#		R2=DATA_DIR + "/{sample}/{sample}_R2.fastq.gz",
		config=config["defuse_config"],
	output:
		"{base}/{sample}/fusion/defuse.raw.txt",
		"{base}/{sample}/fusion/defuse.filtered.txt",
		"{base}/{sample}/fusion/defuse.Reads/defuse.done"
	version: config["defuse"]
	params:
		rulename = "deFuse",
		batch    = config[config['host']]["job_deFuse"]
	shell: """
	#######################
	module load defuse/{version}
	module load R
	export TMPDIR="{wildcards.base}/{wildcards.sample}/fusion/temp/R"
	export TMP="{wildcards.base}/{wildcards.sample}/fusion/temp/R"
	export TEMP="{wildcards.base}/{wildcards.sample}/fusion/temp/R"

	gunzip -c {input.R[0]} >{wildcards.base}/{wildcards.sample}/fusion/{wildcards.sample}_R1.fastq &
	gunzip -c {input.R[1]} >{wildcards.base}/{wildcards.sample}/fusion/{wildcards.sample}_R2.fastq &
	wait
	defuse.pl -c {input.config} \
		-1 {wildcards.base}/{wildcards.sample}/fusion/{wildcards.sample}_R1.fastq\
		-2 {wildcards.base}/{wildcards.sample}/fusion/{wildcards.sample}_R2.fastq\
		-p ${{THREADS}} \
		-n {wildcards.sample}\
		-o {wildcards.base}/{wildcards.sample}/fusion/temp\
		-s direct
	cp {wildcards.base}/{wildcards.sample}/fusion/temp/results.filtered.tsv  {wildcards.base}/{wildcards.sample}/fusion/defuse.filtered.txt
	cp {wildcards.base}/{wildcards.sample}/fusion/temp/results.tsv           {wildcards.base}/{wildcards.sample}/fusion/defuse.raw.txt

	for ID in `cut -f1 {wildcards.base}/{wildcards.sample}/fusion/defuse.filtered.txt|grep -v cluster_id`;
	do
		echo "get_reads.pl -c {input.config} -o {wildcards.base}/{wildcards.sample}/fusion/temp/ -i ${{ID}} >{wildcards.base}/{wildcards.sample}/fusion/defuse.Reads/${{ID}}.txt"
	done >{wildcards.base}/{wildcards.sample}/fusion/cmd.swarm
	cat {wildcards.base}/{wildcards.sample}/fusion/cmd.swarm | parallel -j ${{THREADS}} --no-notice
	touch {wildcards.base}/{wildcards.sample}/fusion/defuse.Reads/defuse.done
	rm -rf {wildcards.base}/{wildcards.sample}/fusion/{wildcards.sample}_R1.fastq {wildcards.base}/{wildcards.sample}/fusion/{wildcards.sample}_R2.fastq
	rm -rf {wildcards.base}/{wildcards.sample}/fusion/temp {wildcards.base}/{wildcards.sample}/fusion/cmd.swarm
	#######################
	"""
############
#       STAR
############
rule STAR:
	input:
		R=lambda wildcards: FQ[wildcards.sample],
#		R1=DATA_DIR + "/{sample}/{sample}_R1.fastq.gz",
#		R2=DATA_DIR + "/{sample}/{sample}_R2.fastq.gz",
		ref=config["reference"],
	output:
		temp("{base}/{sample}/{sample}.star.bam"),
		temp("{base}/{sample}/{sample}.star.bam.bai")
	version: config["STAR"]
	params:
		rulename  = "STAR",
		batch     = config[config['host']]['job_STAR'],
		star_ref  = config['STAR_ref'],
		awk       = NGS_PIPELINE + "/scripts/SJDB.awk",
		home      = WORK_DIR,
		picard    = config['picard']
	shell: """
	#######################
	module load STAR/{version}
	cd ${{LOCAL}}/
	# run 1st pass
	STAR --outTmpDir STEP1 \
		--genomeDir {params.star_ref} \
		--readFilesIn {input.R[0]} {input.R[1]} \
		--readFilesCommand zcat\
		--outFileNamePrefix {wildcards.sample} \
		--runThreadN ${{THREADS}} \
		--outFilterMismatchNmax 2
	echo "Finished Step 1"

	# make splice junctions database file out of SJ.out.tab, filter out non-canonical junctions
	mkdir GenomeForPass2
	awk -f {params.awk} {wildcards.sample}SJ.out.tab > GenomeForPass2/{wildcards.sample}.out.tab.Pass1.sjdb
	echo "Finished Step 2"

	# generate genome with junctions from the 1st pass
	STAR --outTmpDir STEP1\
		--genomeDir GenomeForPass2\
		--runMode genomeGenerate\
		--genomeSAindexNbases 8\
		--genomeFastaFiles {input.ref}\
		--sjdbFileChrStartEnd GenomeForPass2/{wildcards.sample}.out.tab.Pass1.sjdb\
		--sjdbOverhang 100\
		--runThreadN ${{THREADS}}
	echo "Finished Step 3"

	# run 2nd pass with the new genome
	STAR --outTmpDir STEP1\
		--genomeDir GenomeForPass2\
		--runThreadN ${{THREADS}}\
		--outSAMattributes All\
		--readFilesIn {input.R[0]} {input.R[1]}\
		--readFilesCommand zcat\
		--genomeLoad NoSharedMemory\
		--outFileNamePrefix {wildcards.sample}_pass2
	echo "Finished Step 4"

	module load picard/{params.picard}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $PICARD_JARPATH/AddOrReplaceReadGroups.jar\
	VALIDATION_STRINGENCY=SILENT\
	INPUT={wildcards.sample}_pass2Aligned.out.sam\
	OUTPUT={params.home}/{wildcards.base}/{wildcards.sample}/{wildcards.sample}.star.bam\
	SORT_ORDER=coordinate RGLB={wildcards.sample} RGPU={wildcards.sample} RGPL=ILLUMINA RGSM={wildcards.sample} RGCN=khanlab

	echo "Finished Step 5"
	module load samtools
	samtools index {params.home}/{wildcards.base}/{wildcards.sample}/{wildcards.sample}.star.bam
	#######################
	"""
############
#       GATK_RNASeq
############
############
#       GATK Best Practices
############
rule GATK_RNASeq:
	input:  bam="{base}/{sample}/{sample}.star.dd.bam",
		bai="{base}/{sample}/{sample}.star.dd.bam.bai",
		ref=config["reference"],
		phase1=config["1000G_phase1"],
		mills=config["Mills_and_1000G"]
	output:
		bam="{base}/{sample}/{sample}.star.final.bam",
		index="{base}/{sample}/{sample}.star.final.bam.bai",
	version: config["GATK"]
	params:
		rulename  = "gatk",
		batch     = config[config['host']]["job_gatk"]
	shell: """
	#######################
	module load GATK/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T SplitNCigarReads -R {input.ref} -I {input.bam} -o ${{LOCAL}}/{wildcards.sample}.trim.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T RealignerTargetCreator -R {input.ref} -known {input.phase1} -known {input.mills} -I ${{LOCAL}}/{wildcards.sample}.trim.bam -o ${{LOCAL}}/{wildcards.sample}.realignment.intervals

	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T IndelRealigner -R {input.ref} -known {input.phase1} -known {input.mills} -I ${{LOCAL}}/{wildcards.sample}.trim.bam --targetIntervals ${{LOCAL}}/{wildcards.sample}.realignment.intervals -o ${{LOCAL}}/{wildcards.sample}.lr.bam

	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T BaseRecalibrator -R {input.ref} -knownSites {input.phase1} -knownSites {input.mills} -I ${{LOCAL}}/{wildcards.sample}.lr.bam -o ${{LOCAL}}/{wildcards.sample}.recalibration.matrix.txt

	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T PrintReads -R {input.ref} -I ${{LOCAL}}/{wildcards.sample}.lr.bam -o {output.bam} -BQSR ${{LOCAL}}/{wildcards.sample}.recalibration.matrix.txt

	mv {wildcards.base}/{wildcards.sample}/{wildcards.sample}.star.final.bai {output.index}
	######################
	"""
############
# RNASeq Hapcaller
############
rule HapCall_RNASeq:
	input:
		bam="{base}/{sample}/{sample}.star.final.bam",
		bai="{base}/{sample}/{sample}.star.final.bam.bai",
		ref=config["reference"],
		dbsnp=config["dbsnp"]
	output:
		vcf="{base}/{sample}/calls/{sample}.hapCaller.raw.vcf"
	version: config["GATK"]
	params:
		rulename = "HC",
		batch    = config[config['host']]["job_gatk"]
	shell: """
	#######################
	module load GATK/{version}
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T HaplotypeCaller -R {input.ref} -I {input.bam} -o ${{LOCAL}}/{wildcards.sample}.vcf --dbsnp {input.dbsnp} -dontUseSoftClippedBases -stand_call_conf 30 -stand_emit_conf 30
	java -Xmx${{MEM}}g -Djava.io.tmpdir=${{LOCAL}} -jar $GATK_JAR -T VariantFiltration -R {input.ref} -V ${{LOCAL}}/{wildcards.sample}.vcf -window 35 -cluster 3 --filterExpression "FS > 30.0 || QD < 2" -filterName "RNASeqFilters_FS_QD" -o {output.vcf}
	#######################
	"""
############
# Coverage
############
rule coverage:
	input:
		bam="{subject}/{sample}/{sample}.star.final.bam",
		bai="{subject}/{sample}/{sample}.star.final.bam.bai",
		interval=config["refSeq_bed"]
	output:
		"{subject}/{sample}/qc/{sample}.star.coverage.txt"
	version: config["bedtools"]
	params:
		rulename = "coverage",
		batch    = config[config['host']]["job_bedtools"]
	shell: """
	#######################
	module load bedtools/{version}
	bedtools coverage -abam {input.bam} -b {input.interval} -hist |grep "^all" > {output}
	#######################
	"""
############
# Filter fusion for every library
############
rule sub_fusion:
	input:
		tophat="{subject}/{sample}/fusion/tophat-fusion.txt",
		fc="{subject}/{sample}/fusion/fusion-catcher.txt",
		defuse="{subject}/{sample}/fusion/defuse.filtered.txt",
		convertor = NGS_PIPELINE + "/scripts/" + config['Actionable_fusion'],
		mitelman  = config["Mitelman"],
		omim	  = config["omim_fusion"],
		TCGA      = config["TCGA_fusion"]
	output:
		"{subject}/{sample}/fusion/{sample}.actionable.fusion.txt"
	params:
		rulename = "Sub_F",
		batch    = config[config['host']]["job_default"]
	shell: """
	#######################
	mkdir -p {wildcards.subject}/Actionable
	perl {input.convertor} {input.mitelman} {input.omim} {input.TCGA} {wildcards.sample} {input.defuse} {input.tophat} {input.fc} {wildcards.subject}/{ACT_DIR}/ |sort >{output}
	#######################
	"""
############
# Combine filtered fusions to actionable.
############
rule Actionable_fusion:
	input:
		fusion=lambda wildcards: SUB_FUSION[wildcards.subject]
	output:
		"{subject}/{ACT_DIR}/{subject}.fusion.actionable.txt"
	params:
		rulename = "Actionable_F",
		batch    = config[config['host']]["job_default"]
	shell: """
	#######################
	cat {input.fusion} |sort |uniq >{output}
	#######################
	"""