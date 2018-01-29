#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("gridExtra"))

option_list <- list(
	make_option("--input",  help="bcftools covrage compare metrix, required"),
	make_option("--patient",  help="patient ID"),
	make_option("--output", help="output pdf file name.")
	)

opt <- parse_args(OptionParser(option_list=option_list))

if( ! is.element('input', names(opt)) ) stop("Maggie result is required. ")
if( ! is.element('patient', names(opt)) ) stop("Patient name is required. ")
if( ! is.element('output', names(opt)) ) stop("Output file name is required. ")
patient=opt$patient

a <-read.table(opt$input, header=T,sep="\t")
a$Score <-round(as.numeric(a$Score/a$Sites),digits=5)
a <- a[order(a$Score),]
#a$Score[a$Score==0]<-NA

a$Color[grep(".*WES.*",a$Sample,perl=TRUE,value = FALSE)]="green"
a$Color[grep(".*RNASEQ.*",a$Sample,perl=TRUE,value = FALSE)]="darkred"
a$Color[a$Sites <=15000]="orange"
a$Color[grep(paste(patient,".*WES.*",sep=""),a$Sample,perl=TRUE,value = FALSE)]="blue"
a$Color[grep(paste(patient,".*RNASEQ.*",sep=""),a$Sample,perl=TRUE,value = FALSE)]="red"

a$shape<-19
a$shape[grep(patient,a$Sample,perl=TRUE,value = FALSE)]=8

a$size<-0.60
a$size[grep(patient,a$Sample,perl=TRUE,value = FALSE)]=1.2

d<-nrow(a[a$Sites<=15000,])
pdf(opt$output, height=10, width=10, points=20)
plot(a$Score,
	log="x", 
	col=a$Color,
	pch=a$shape,
	cex=a$size,
	las=1,
	ylab="Discordance",
	xlab="Ordered Samples", 
	main=paste(patient,"\nMAGGIE Result",sep="\n")
)

legend("bottomright",
	legend = c("WES", "RNASEQ", "Patient WES", "Patient RNASEQ",paste(d,"Samples <15K sites",sep=" ")),
	text.col = c("green","darkred","blue","red","orange"),
	bty = "n"
)
plot.new()
c<-a[c("Sample","Score")]
grid.table(head(c,n=20),rows = NULL)

dev.off()
