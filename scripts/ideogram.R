library(stringr)
library(RColorBrewer)

########################################
ideogram=function(segmentedData, color=NULL, title=NULL){
  if(is.null(title)){
    title= "B Allele Frequency Plot"
  }
  chr.lens = c(249250621, 243199373, 198022430, 191154276, 
               180915260, 171115067, 159138663, 146364022, 
               141213431, 135534747, 135006516, 133851895, 
               115169878, 107349540, 102531392, 90354753, 
               81195210,  78077248,  59128983,  63025520, 
               48129895,  51304566,  155270560, 59373566)

  chr.mid = c( 124625310, 370850307, 591461209, 786049562,
               972084330, 1148099494,1313226359,1465977701,
               1609766428,1748140517,1883411148,2017840354,
               2142351240,2253610949,2358551415,2454994488,
               2540769469,2620405698,2689008814,2750086065,
               2805663773,2855381003,2958668566,3065990629)
segmentedData=data.table::setDT(segmentedData)
segmentedData[,Position := as.numeric(as.character(Position))]
segmentedData$Chr = gsub(pattern = 'chr', replacement = '', x = segmentedData$Chr, fixed = TRUE)
segmentedData$Chr = gsub(pattern = 'X',   replacement = '23', x = segmentedData$Chr, fixed = TRUE)
segmentedData$Chr = gsub(pattern = 'Y',   replacement = '24', x = segmentedData$Chr, fixed = TRUE)
segmentedData$Chr = factor(x = segmentedData$Chr, levels = 1:24, labels = 1:24)
segmentedData = segmentedData[order(Chr, Position, decreasing = FALSE)]
seg.spl = split(segmentedData, segmentedData$Chr)
seg.spl.transformed = seg.spl[[1]]
if(nrow(seg.spl.transformed) > 0){
  seg.spl.transformed$Position_updated = seg.spl.transformed$Position
}
chr.lens.sumsum = cumsum(chr.lens)

for(i in 2:length(seg.spl)){
  x.seg = seg.spl[[i]]
  if(nrow(x.seg) > 0){
    x.seg$Position = x.seg$Position + chr.lens.sumsum[i-1]
  }
  seg.spl.transformed = rbind(seg.spl.transformed, x.seg, fill = TRUE)
}

plot(NA, NA,
     axes =FALSE,
     xlab =NA,
     ylab =NA,
     xaxs ="i",
     ylim =c(0,1),
     xlim =c(0, max(chr.lens.sumsum)),
     main =title,
     cex  =2
     )
abline(h=0.5, col="lightgrey")
points(x = seg.spl.transformed$Position, y = seg.spl.transformed$LOH, col = color, pch = 16, cex = 0.5)
segments(x0 = chr.lens.sumsum, y0 = 0, x1 = chr.lens.sumsum, y1 = 1, col = 'black', lwd = 1)
axis(side = 2, at = c(0,.5,1), las = 2)
axis(side = 1, at = chr.mid, labels = c(1:22, "X", "Y"), tick = FALSE, line = -0.7, cex.axis = 1.5, col=c('blue', 'red'))
box()
}
########################################

args<-commandArgs(TRUE)
DIR = str_trim(args[1])
FILE=str_trim(args[2])
SAM=str_trim(args[3])

files <- list.files(path = DIR, pattern=".loh$")
for (i in 1:length(files)){
	if (length(files) >1){
		if (length(grep("_P.*.bwa.loh",files[i]))>0){
			files <- files[-grep("_P.*.bwa.loh",files,perl=TRUE,value = FALSE)]
		}
	}
}


labs <- paste("", gsub("Sample_|\\~WES|\\~RNASEQ|\\.bwa|\\.star|\\.loh", "", files, perl=TRUE), sep="")



cols <-c('#26294a','#01545a','#bd544f','#017351',
	'#03c383','#b8bd4f','#aad962','#fbbf45',
	'#bd8b4f','#ef6a32','#ed0346','#d76e60',
	'#a12a5e','#710162','#26294a','#01545a',
	'#bd544f','#017351','#03c383','#b8bd4f',
	'#aad962','#fbbf45','#bd8b4f','#ef6a32',
	'#ed0346','#d76e60','#a12a5e','#710162'
       )

options(stringsAsFactors = FALSE)
if (length(files) >4){
	plot_height=length(files)*80
}else{
	plot_height=length(files)*150
}
png(FILE ,width = 1200, height = plot_height)
#png(FILE ,width = 1200, height = plot_height, res=100, points=12, type=c("cairo"))
par(mfrow=c(length(files),1), mar=c(1,4,2,1))

for (i in 1:length(files)){
	LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T)
	ideogram(LOH.data, col=cols[i], title=labs[i])
}
dev.off()
