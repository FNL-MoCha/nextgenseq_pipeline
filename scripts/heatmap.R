#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("RColorBrewer"))
suppressPackageStartupMessages(library("pheatmap"))
suppressPackageStartupMessages(library("gplots"))
option_list <- list(
        make_option(c("-s", "--subject"),help="Title of the plot(subject name),required"),
        make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
        help="to output some information about a job.  [default: %default]")
    )
opt <- parse_args(OptionParser(option_list=option_list))

subject<-opt$subject

input <- read.table(paste("/data/MoCha/processedDATA/",subject,"/20170910/variantHeatmap/",subject,".matrix",sep=""), header = T,sep="\t",check.names=FALSE)
pdf(paste("/data/MoCha/processedDATA/",subject,"/20170910/",subject,".variantsHeatMap.pdf", sep=""), width=30, height=15, onefile=FALSE)
par(mar=c(2,2,5,0))

rownames(input) <- paste(input$Gene,input$Start,input$Alt,sep="__")
mat <- as.matrix(input[7:ncol(input)])

pairs.breaks <- seq(0, 1, by=0.2)
brk <- length(pairs.breaks)-1
mycol <- colorpanel(n=brk,low="blue",high="red")

mat <-t(mat)
pheatmap(mat,
	color = c("#FFFFFF","#0000FF", "#4000BF", "#800080", "#BF0040", "#FF0000"),
	border_color = NA,
	breaks = c(-1,0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
	cluster_rows = T,
	cluster_cols = T,
	treeheight_row = 0, 
	treeheight_col = 0,
	fontsize = 20, 
	fontsize_row = 10, 
	fontsize_col = 2,	
	main=paste(subject,"\n",sep=""),
	legend = TRUE
)
dev.off()
