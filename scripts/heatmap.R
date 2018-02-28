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
	color = mycol,
	border_color = NA,
	breaks = pairs.breaks,
	cluster_rows = T,
	cluster_cols = T,
	fontsize = 20, 
	fontsize_row = 10, 
	fontsize_col = 2,	
	main=paste(subject,"\n",sep=""),
	legend = TRUE
)
dev.off()
