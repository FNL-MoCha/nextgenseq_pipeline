#EV Plot with Percent Ancestry Overlay
data=read.table("ancestry.txt", as.is=T, header=F)
names(data) <- c("ID", "Case", "SNPs", "EV1", "EV2", "EUR", "AFR", "ASN")

plot(data$EV1, data$EV2, pch=20, col="gray", xlab="EV1", ylab="EV2")
text(data$EV1, data$EV2,labels=round(data$EUR,2)*100, cex=0.4, offset=0.1, pos=3)
text(data$EV1, data$EV2,labels=round(data$AFR,2)*100, cex=0.4, offset=0.1, pos=2)
text(data$EV1, data$EV2,labels=round(data$ASN,2)*100, cex=0.4, offset=0.1, pos=1)

# Triangle Plot
data$total=data$EUR+data$AFR+data$ASN         # Need to account
data$European=data$EUR/data$total             # for slight rounding
data$African=data$AFR/data$total              # in the ancestry
data$Asian=data$ASN/data$total                # estimation file for
data_p=data[c("European","Asian","African")]  # triax.plot to work

library(plotrix)
triax.plot(data_p, pch=20, cc.axes=T, show.grid=T)
