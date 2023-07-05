library(factoextra)

rm(list = ls())
#Change my working Dir
path = "Merge_LMM_out"
inFile_name = "HM_LMM.txt" # Tissue or Species

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

## PCA analysis -------------------------------------------------------------
res.pca <- prcomp(inFile_rowCol, scale = TRUE)
# Visualize eigenvalues/variances
fviz_screeplot(res.pca, addlabels = TRUE, main = Title)#, ylim = c(0, 50))

outName = paste0(outFile_name, "_var_percent")
outPlot(outName)
##--------------------------------------------------------------------------
fviz_pca_ind(res.pca,
             col.ind = "cos2", # Color by the quality of representation
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)
outName = paste0(outFile_name, "_PCA")
outPlot(outName)



## Dendrogram --------------------------------------------------------------
# Compute hierarchical clustering and cut into 4 clusters
res.dist <- factoextra::get_dist(inFile_rowCol, method = "euclidean")
res.hclust <- stats::hclust(res.dist, method = "average")

# Visualize
fviz_dend(res.hclust, rect = TRUE, cex = 0.7,
          k = 3,
          k_colors = c("#2E9FDF", "#E7B800", "#FC4E07", "#00AFBB", "darkorchid3"))

outName = paste0(outFile_name, "_Dend_kmean_euc_ave")
outPlot(outName)