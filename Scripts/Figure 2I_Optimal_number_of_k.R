library(factoextra)

rm(list = ls())
#Change my working Dir
path = "05_Merge_LMM_out_220127_kmean"
inFile_name = "LMM_2-3.txt" # Tissue or Species

outFile_name = "HM_LMM"
Title = "HM_LMM"
#------------------------------------------------------------------------
# Read Data
setwd(path)
inFile = read.table(inFile_name, sep="\t" ,header=T, row.names = 1)
#inFile = inFile[,-c(1:4)]
inFile_rowCol = t(inFile)

## Functions
##--------------------------------------------------------------------------
outPlot <- function(outName){
  outFile_png = paste0(outName, ".png")
  outFile_pdf = paste0(outName, ".pdf")
  dev.print(png, file = outFile_png,
    width=800,
    height=600)
  dev.print(pdf, file = outFile_pdf)
}

## Optimal # of k ----------------------------------------------------------
# Optimal number of clusters for k-means
fviz_nbclust(inFile_rowCol, kmeans, method = "gap_stat",
             diss = dist(inFile_rowCol, method= "euclidean"),
             k.max = 25,
             nboot = 1000)
outName = paste0(outFile_name, "_optiman_k_n1000")
outPlot(outName)



