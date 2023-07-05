library(TissueEnrich)
library(SummarizedExperiment)

rm(list=ls())
setwd("TissueEnrich")

# Raw data
expressionData<-read.table('Dog_11tissue_RNA_FPKM-UQ.txt',header=TRUE,row.names=1,sep='\t')
se = SummarizedExperiment(assays = SimpleList(as.matrix(expressionData)),
                          rowData = row.names(expressionData),
                          colData = colnames(expressionData))

output = teGeneRetrieval(se, foldChangeThreshold=4, maxNumberOfTissues=3)

head(assay(output))
write.csv(assay(output), 'Dog_11tissue_RNA_FPKM-UQ_tiEnrich_4fold_group2-3.txt', row.names=FALSE, quote=FALSE)
