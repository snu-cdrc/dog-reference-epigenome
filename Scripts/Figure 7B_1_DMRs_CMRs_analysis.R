# This are necessary libraries to use
library(DESeq)
library(GenomicRanges)
library(GenomeInfoDb)
library(data.table)
library(goldmine)

# Load required code (methylaction source code)
source("SourceCode_Rev_line310.R")

# Set working directory and load the sample list
samp <- readSampleInfo("BeagelRefMBD.csv")
print(samp)

# Set basic parameters
chrs <- paste0("chr",c(1:38,"X"))
fragsize <- 200
winsize <- 100
ncore <- 1
freq <- 1 #2/3

# Get reads from BAM files
print(Sys.time())
#2021-07-12 13:02:40 KST
reads <- getReads(samp=samp, chrs=chrs, fragsize=fragsize, ncore=ncore);print(Sys.time())
#2021-07-12 13:30:47 KST
save(reads, file = paste0(gsub("-","",Sys.Date()),"_1_reads_BAMfiles_100bp.Rdata"))
#assign("OldReads", get(load("../20200624_MethylAction_New/20200624_1_reads_BAMfiles_100bp.Rdata")))

# Get counts from above reads file with proper winsize (default: 50 bp) 
print(Sys.time())
#2021-07-12 13:45:42 KST
counts <- getCounts(samp=samp, reads=reads, chrs=chrs, winsize=winsize, ncore=ncore);print(Sys.time())
#2021-07-12 13:48:10 KST
save(counts, file = paste0(gsub("-","",Sys.Date()),"_2_counts_BAMfiles_100bp.Rdata"))

# MethylAction (DMR analysis without bootstrapping or permutations)
print(Sys.time())
#2021-07-12 13:48:58 KST
ma <- methylaction(samp=samp, counts=counts, reads=reads, poifdr=0.1, stageone.p=0.01, anodev.p=0.01, post.p=0.05, freq=freq, minsize=100, joindist=0, nperms=0, perm.boot=F, ncore=ncore);print(Sys.time())
#2021-07-13 05:03:14 KST
save(ma, file = paste0(gsub("-","",Sys.Date()),"_3_ma_100bp.Rdata"))

# Summaries for MethylAction
maPat <- as.data.frame(table(ma$dmr$pattern,ma$dmr$frequent))
maSum <- maSummary(ma)
maFDR <- as.data.frame(ma$fdr)

# Export outputs
write.csv(maPat, row.names = F, file=paste0(gsub("-","",Sys.Date()),"_3_maPat_tsDMR.csv"))
write.csv(maSum, row.names = F, file=paste0(gsub("-","",Sys.Date()),"_3_maSum_tsDMR.csv"))
write.csv(maFDR, row.names = F, file=paste0(gsub("-","",Sys.Date()),"_3_maFDR_tsDMR.csv"))
write.csv(makeDT(ma$dmr), row.names=FALSE, file=paste0(gsub("-","",Sys.Date()),"_3_dmrs_tsDMR.csv"))

# Export reads and dmrs to BED format
maBed(ma,file=paste0(gsub("-","",Sys.Date()),"_3_dmrs_tsDMR.bed"))

# Draw heatmap with frequent DMRs
library(gplots)
pdf(paste0(gsub("-","",Sys.Date()),"_3_ma_tsDMR_heatmap.pdf"))
maHeatmap(ma)
dev.off()

# Draw karyogram with frequent DMRs
library(ggplot2)
pdf(paste0(gsub("-","",Sys.Date()),"_3_ma_tsDMR_karyogam.pdf", 8,6))
maKaryogram(ma=ma, reads=reads)
dev.off()

# Export normalized counts
norm <- as.data.frame(ma[["data"]][["windows"]][["signal.norm"]])
norm$dmrid <- paste0(norm$seqnames,"_",norm$start)
colnames(norm)[1] <- "chr"
save(norm, file = paste0(gsub("-","",Sys.Date()),"_4_ma_100bp_normcounts.Rdata"))
