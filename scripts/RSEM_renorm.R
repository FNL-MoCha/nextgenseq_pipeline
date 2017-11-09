#!/usr/bin/env Rscript
# Convert expected counts to transcripts per million (TPM).
# 
#    Lior Pachter. Models for transcript quantification from RNA-Seq.
#    arXiv:1104.3889v2 
#    
#    Wagner, et al. Measurement of mRNA abundance using RNA-seq data:
#    RPKM measure is inconsistent among samples. Theory Biosci. 24 July 2012.
#    doi:10.1007/s12064-012-0162-3

# use the results from RSEM with expected counts, efflength, to get TPM in a sample
#/data/MoCha/patidarr/ngs_pipeline/scripts/RSEM_renorm.R SA0350~F605~M540M786~RNASEQ.isoforms.results SA0350~F605~M540M786~RNASEQ.genes.results SA0350~F605~M540M786~RNASEQ.isoforms.results.out SA0350~F605~M540M786~RNASEQ.genes.results.out

args <- commandArgs(T)
if (length(args) != 4) {
  stop("two file needs to be the input argument", call. = FALSE)
} 


###first argument: isoform.txt, second argument: gene.txt
Iso <- read.table(args[1], header = T, sep = "\t", as.is = T, check.names = F)
Gene <- read.table(args[2], header = T, sep = "\t", as.is = T, check.names = F)
message(sprintf("Reading counts file: %s, %s\n", args[1], args[2]),appendLF = FALSE)

#remove SNORs
Iso.noSNOR <- subset(Iso, !grepl("^SNOR", Iso$gene_id))
Gene.noSNOR <- subset(Gene, !grepl("^SNOR", Gene$gene_id))

#function for re-normalization
Cal_exp <- function(dat){
  rate <- dat$expected_count/dat$effective_length
  #assign 0 rate & TPM to genes with effLen <=1 !!!!
  rate[dat$effective_length <= 1] <- 0 
  sum_rate <- sum(rate)
  tpm <- rate/sum_rate*1e6
  fpkm <- rate/sum(dat$expected_count)*1e6*1e3
  return(list(TPM = tpm, FPKM = fpkm))
}

#process the data
Iso.tpm <- round(Cal_exp(Iso.noSNOR)[[1]], 2)
Iso.fpkm <- round(Cal_exp(Iso.noSNOR)[[2]], 2)
Iso.noSNOR <- data.frame(Iso.noSNOR[,1:5], TPM = Iso.tpm, FPKM = Iso.fpkm, IsoPct = Iso.noSNOR$IsoPct)

#get the real expected counts from isoform table
library("data.table")
Iso.noSNOR.mod <- Iso.noSNOR
Iso.noSNOR.mod$expected_count[Iso.noSNOR.mod$effective_length <= 1] <- 0
DT <- data.table(Iso.noSNOR.mod, key = "gene_id")
iso2gen <- as.data.frame(DT[, list("expected_count" = sum(expected_count)), by = "gene_id"])

#implement the expected counts into gene table
Gene.noSNOR$expected_count <- iso2gen$expected_count
Gene.tpm <- round(Cal_exp(Gene.noSNOR)[[1]], 2)
Gene.fpkm <- round(Cal_exp(Gene.noSNOR)[[2]], 2)
Gene.noSNOR <- data.frame(Gene.noSNOR[,1:5], TPM = Gene.tpm, FPKM = Gene.fpkm)


#write new file (assume args[1] name: *.genes.results.old)
write.table(Iso.noSNOR, args[3], sep = "\t", quote = F, row.names = F)
write.table(Gene.noSNOR, args[4], sep = "\t", quote = F, row.names = F)
