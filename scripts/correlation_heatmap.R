#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(require(ConsensusClusterPlus))
suppressPackageStartupMessages(require(gplots))
suppressPackageStartupMessages(require(RColorBrewer))
suppressPackageStartupMessages(require(mixtools))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(reshape2))

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

# run assemble_rnaseq_gene.py to generate the rnaseq gene expression matrix 
data <- read.csv(opt$input,check.name=FALSE,row.names=1)

cormat <- cor(data)
melted_cormat <- melt(cormat,na.rm=FALSE)
pdf(opt$output)

#ggplot(melted_cormat,aes(Var2,Var1,fill=value)) + labs(title = patient) + theme(plot.title = element_text(hjust = 0.5)) + geom_tile(color="white") + scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=0.8,limit=c(0.5,1), na.value = "black", space="Lab",name="Spearman\nCorrelation") +
ggplot(melted_cormat,aes(Var2,Var1,fill=value)) + labs(title = patient) + theme(plot.title = element_text(hjust = 0.5)) + geom_tile(color="white") + scale_fill_gradient2(low="blue",high="red",mid="white",midpoint=0.6,limit=c(0.2,1), na.value = "black", space="Lab",name="Spearman\nCorrelation") +
  theme_minimal() + theme(axis.text.x=element_text(angle=45,vjust=1,size=6,hjust=1)) + theme(axis.text.y=element_text(hjust=0,size=6)) + 
  coord_fixed()
dev.off()
