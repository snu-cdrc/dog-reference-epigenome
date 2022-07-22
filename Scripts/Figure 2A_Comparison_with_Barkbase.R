library(DESeq2)
library(tximport)
library(ggplot2)
library(ggrepel)

library(BiocParallel)
register(MulticoreParam(4))
## ---------------------------------------
path_rsem = "output_rsem"

in_sample_info = "sample_info_simple.txt"
suffix_rsem = ".str.rs.Ev102.genes.results"

out_name= "DTE_PCA-vst.RSEM.pdf"
## ---------------------------------------
setwd(path_rsem)
# read sample info file
sampleData <- read.table(in_sample_info, sep = "\t", header = TRUE) 
rownames(sampleData) = sampleData$label
# save rsem output files path
files <- file.path(path_rsem, paste0(sampleData$label, suffix_rsem)) 
# gather rsem files to merged table
names(files) <- sampleData$label
files
txi <- tximport(files, type = "rsem", txIn = FALSE, txOut = FALSE) 
summary(txi)
# fix "Error: all(lengths > 0) is not TRUE" error
txi$length[txi$length == 0] = 0.01
# Filter out lowly expressed genes
#filter = apply(txi$counts, 1, function(row) length(row[row>5])>=2)
filter = apply(txi$counts, 1, function(row) sum(row) >= 10)
txi$length = txi$length[filter,]
txi$counts = txi$counts[filter,]
rm(filter)

# read in 
ddsTxi <- DESeqDataSetFromTximport(txi, colData = sampleData, ~ type)
# set the reference
ddsTxi = DESeq(ddsTxi, parallel = T) # DESeq
ddsTxi

## Visualization -------------------------
vsd = vst(ddsTxi)
vsd
# vst
pcaData = plotPCA(vsd, intgroup=c("type"),
                  returnData=TRUE)
percentVar = round(100 * attr(pcaData, "percentVar"))

pdf(out_name)
#png(out_name, width=7, height=7, units = "in", res = 300)
ggplot(pcaData, aes(PC1, PC2, colour = type)) +
  geom_point(size = 2) + theme_bw() +
  geom_text_repel(aes(label = type), nudge_x = -1, nudge_y = 0.2, size = 3) +
  ggtitle("Principal Component Analysis (PCA)", subtitle = "vst transformation") +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))
dev.off()

