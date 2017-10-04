suppressPackageStartupMessages(library(OmicCircos))
library(stringr)
library(RColorBrewer)

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


labs <- paste("", gsub("Sample_|\\.bwa|\\.star|\\.loh", "", files, perl=TRUE), sep="")



data("UCSC.hg19.chr", package="OmicCircos");
data("UCSC.hg19", package="OmicCircos");
ref     <- UCSC.hg19.chr;
ref[,1] <- gsub("chr", "", ref[,1]);
chr     <- unique(ref[,1]);
chr.l   <- c();
for (ch in chr){
  dat.i <- which(ref[,1]==ch);
  m     <- max(ref[dat.i,3])/2;
  chr.l <- rbind(chr.l, c(ch, m, ch));
}
angle <- function(x, po.min, po.max, angle.start, angle.end){ 
  ag <- (x - po.min)/(po.max-po.min) * (angle.end - angle.start) + angle.start;
  ag;
}

angles    <- angle(as.numeric(chr.l[,2]), as.numeric(UCSC.hg19[,6]), 
                   as.numeric(UCSC.hg19[,7]), as.numeric(UCSC.hg19[,2]),
                   as.numeric(UCSC.hg19[,3]));
r         <- 418;
coord.x   <- 400 + cos((angles)/180*pi) * r;
coord.y   <- 400 - sin((angles)/180*pi) * r;


cols <-c('#26294a','#01545a','#bd544f','#017351',
	'#03c383','#b8bd4f','#aad962','#fbbf45',
	'#bd8b4f','#ef6a32','#ed0346','#d76e60',
	'#a12a5e','#710162','#26294a','#01545a',
	'#bd544f','#017351','#03c383','#b8bd4f',
	'#aad962','#fbbf45','#bd8b4f','#ef6a32',
	'#ed0346','#d76e60','#a12a5e','#710162'
       )

options(stringsAsFactors = FALSE)
set.seed(1234)
png(FILE ,width = 1000, height = 1000, res=100, points=12, type=c("cairo"))
par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="", main="")

circos(R=400, cir="hg19", type="chr",    mapping=UCSC.hg19.chr ,print.chr.lab=FALSE, W=10, lwd=5, cex=1.5)
text(coord.x, coord.y, chr.l[,1], cex=1.4, col=c("red", "blue"));

r=350
if (length(files) >6){
	#for (i in 1:length(files)){
	for (i in 1:6){
        	LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T)
       		circos(cir="hg19", R=r, W=50, type="s", mapping=LOH.data, col.v=3, col=cols[i], B=FALSE, cex=0.0001, lwd=1)
		r=r-45
	}
} else{
	for (i in 1:length(files)){
                LOH.data   <-read.table(paste(DIR,files[i] ,sep = ""), sep="\t", quote="", head=T)
                circos(cir="hg19", R=r, W=50, type="s", mapping=LOH.data, col.v=3, col=cols[i], B=FALSE, cex=0.0001, lwd=1)
                r=r-45
        }

}
dev.off()
newname <-gsub("circos.png", "legend.circos.png", FILE)
png(newname ,width = 1000, height = 1000, res=100, points=12, type=c("cairo"))
par(mar=c(2, 2, 2, 2))
plot(c(1,800), c(1,800), type="n", axes=FALSE, xlab="", ylab="")
if (length(files) >6){
	legend("center", legend=labs[0:6], col=cols, lty=1, lwd=4, cex=2)
}else{
	legend("center", legend=labs, col=cols, lty=1, lwd=4, cex=2)
}
dev.off()

