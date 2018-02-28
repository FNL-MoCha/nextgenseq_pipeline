cols=c("dodgerblue","green","pink","red","purple","darkgreen","gold", "coral2","grey")
input<-read.table("../signature2", sep="\t", header=T, check.names =FALSE)
rownames(input) <-input$Signature
input<-input[2:ncol(input)]
colnames(input) <-paste(colnames(input)," (",colSums(input),")",sep="")
data.perc <- apply(input, 2, function(x){x/sum(x)*100})
op <- par(mar=c(6,14,2,10))
xx <-barplot(data.perc, 
        col=cols,
        space=0.01,
        border=NA,
        las=1,
        horiz = T,
        #xaxt='n',
        xlab="Percent of cases with \ndominant Signature",
        legend = rownames(input),
        args.legend = list(x = "topright",
                           bty = "o",
                           box.lwd = 0,
                           inset=c(-0.6, 0.3))
        )
#Signature	Digestive-Gastrointestinal	Respiratory-Thoracic	Head-Neck	Musculoskeletal	Genitourinary	Skin	Gynecologic	NeurologicEndocrine-Neuroendocrine	Unknown-Primary	Breast	Hematologic-Blood
#Age	36	4	7	4	4	1	1	2	0	0	0	0
#APOBEC	0	1	1	0	4	0	0	0	0	0	0	0
#BRCA1/2	10	13	13	20	13	5	6	2	3	1	1	1
#MMR	0	0	0	0	0	0	1	0	0	0	0	0
#MSI	7	1	1	0	1	0	4	0	0	0	0	0
#Smoking	0	11	1	0	0	0	0	0	0	0	0	0
#TMZ	0	1	0	0	0	0	0	0	0	0	0	0
#UV	0	0	2	0	0	14	0	0	0	0	0	0
#Other	4	0	5	2	3	0	0	0	1	1	0	0