library(readr)
library(stringr)
library(EnrichedHeatmap)
library(circlize)

##-----##
## CGI ##                                
##-----##
CGI <- read_csv("CGIregions_CanFam3_unmasked.csv")
CGI <- CGI[which(CGI$name == "island"),]
colnames(CGI) <- c("Chr", "Start", "End", "width", "strand", "name")
cpg <- CGI[which(!str_detect(CGI$Chr, "ChrUn")),]
cpg$Chr <- as.character(cpg$Chr)
cpggr <- GRanges(cpg$Chr, ranges = IRanges(cpg$Start, cpg$End))

##------------##
## ChIP / MBD ##
##------------##
tissList = c("CL","CO","CR","KI","LI","LU","MG","OV","PA","SP","ST")
assign('k27ac',get(load(paste0("k27ac_",tissue,".bedgraph_GRanges.Rdata"))))
assign('k27me3',get(load(paste0("k27me3_",tissue,".bedgraph_GRanges.Rdata"))))
assign('k4me1',get(load(paste0("k4me1_",tissue,".bedgraph_GRanges.Rdata"))))
assign('k4me3',get(load(paste0("k4me3_",tissue,".bedgraph_GRanges.Rdata"))))
assign('k9me3',get(load(paste0("k9me3_",tissue,".bedgraph_GRanges.Rdata"))))
assign('methyl',get(load(paste0("MBD_",tissue,".bedgraph_GRanges.Rdata"))))

##-----------------##
## generate matrix ##
##-----------------##
mat_k27ac = normalizeToMatrix(k27ac, cpggr,
                             extend = 4000, mean_mode = "coverage", smooth = TRUE, w = 50,
                             target_ratio = 0.1, empty_value = 0)
mat_k27me3 = normalizeToMatrix(k27me3, cpggr,
                              extend = 4000, mean_mode = "coverage", smooth = TRUE, w = 50,
                              target_ratio = 0.1, empty_value = 0)
mat_k4me1 = normalizeToMatrix(k4me1, cpggr,
                             extend = 4000, mean_mode = "coverage", smooth = TRUE, w = 50,
                             target_ratio = 0.1, empty_value = 0)
mat_k4me3 = normalizeToMatrix(k4me3, cpggr,
                             extend = 4000, mean_mode = "coverage", smooth = TRUE, w = 50,
                             target_ratio = 0.1, empty_value = 0)
mat_k9me3 = normalizeToMatrix(k9me3, cpggr,
                             extend = 4000, mean_mode = "coverage", smooth = TRUE, w = 50,
                             target_ratio = 0.1, empty_value = 0)
mat_mbd = normalizeToMatrix(methyl, cpggr,
                           extend = 4000, mean_mode = "coverage", smooth = TRUE, w = 50,
                           target_ratio = 0.1, empty_value = 0)
save(mat_k27ac,mat_k27me3,mat_k4me1,mat_k4me3,mat_k9me3,mat_mbd,
    file = paste0(gsub("-","",Sys.Date()),"_1_EH(mat)_hmarks&methyl_",tissue,"_inCGI.Rdata"))

for (i in 1:length(tissList)) {
  tissue <- tissList[i]
  load(file = paste0("20211004_1_EH(mat)_hmarks&methyl_",tissue,"_inCGI.Rdata"))
  
  ##------------##
  ## set colors ##
  ##------------##
  k27ac_col = colorRamp2(c(0,quantile(mat_k27ac, 0.25),
                           quantile(mat_k27ac, 0.95, na.rm = TRUE)),
                         c("white", "white", "purple"))
  k4me1_col = colorRamp2(c(0,quantile(mat_k4me1, 0.25), 
                           quantile(mat_k4me1, 0.95, na.rm = TRUE)),
                         c("white", "white", "purple"))
  k4me3_col = colorRamp2(c(0,quantile(mat_k4me3, 0.25), 
                           quantile(mat_k4me3, 0.95, na.rm = TRUE)),
                         c("white", "white", "purple"))
  k27me3_col = colorRamp2(c(0,quantile(mat_k27me3, 0.25), 
                            quantile(mat_k27me3, 0.95, na.rm = TRUE)),
                          c("white", "white", "darkgreen"))
  k9me3_col = colorRamp2(c(0,quantile(mat_k9me3, 0.25), 
                           quantile(mat_k9me3, 0.95, na.rm = TRUE)),
                         c("white", "white", "darkgreen"))
  mbd_col = colorRamp2(c(0,quantile(mat_mbd, 0.25), 
                         quantile(mat_mbd, 0.95, na.rm = TRUE)),
                       c("white", "white", "indianred"))
  
  ##------------------##
  ## enriched heatmap ##
  ##------------------##
  ht_list = EnrichedHeatmap(mat_mbd, use_raster = TRUE, col = mbd_col, name="Methylation",
                            width = unit(25, "mm"),
                            column_title = paste0("[",tissue,"] Methyl")) +
    EnrichedHeatmap(mat_k27me3, use_raster = TRUE, col = k27me3_col, name="H3K27me3",
                    width = unit(25, "mm"), column_title = "H3K27me3") +
    EnrichedHeatmap(mat_k9me3, use_raster = TRUE, col = k9me3_col, name="H3K9me3",
                    width = unit(25, "mm"), column_title = "H3K9me3") +
    EnrichedHeatmap(mat_k27ac, use_raster = TRUE, col = k27ac_col, name="H3K27ac",
                    width = unit(25, "mm"), column_title = "H3K27ac") +
    EnrichedHeatmap(mat_k4me1, use_raster = TRUE, col = k4me1_col, name="H3K4me1",
                    width = unit(25, "mm"), column_title = "H3K4me1") +
    EnrichedHeatmap(mat_k4me3, use_raster = TRUE, col = k4me3_col, name="H3K4me3",
                    width = unit(25, "mm"), column_title = "H3Kme3")
  
  pdf(paste0(gsub("-","",Sys.Date()),"_1_EH_methyl&hmarks_",tissue,"_inCGI_col(0,Q0.25,Q0.95).pdf"), width = 8, height = 10)
  draw(ht_list)
  dev.off()
}
