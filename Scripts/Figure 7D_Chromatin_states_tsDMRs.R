rm(list = ls())
library(readr)
library(stringr)
library(GenomicRanges)
library(ComplexHeatmap)
library(EnrichedHeatmap)
library(circlize)

##-------##
## tsDMR ##
##-------##
tsDMR <- read_csv("20210719_1_tsDMR20377_anno.csv")
tsDMRgr <- GRanges(seqnames = tsDMR$chr,
                   ranges = IRanges(tsDMR$start, tsDMR$end),
                   strand = tsDMR$Strand)
names(tsDMRgr) <- tsDMR$dmrid
head(tsDMRgr)

##----------##
## ChromHMM ##
##----------##
## Chromatin states _ simplify
map = c("1"  = "Promoter", "2"  = "Promoter",
        "3"  = "Promoter", "4"  = "Promoter",
        "5"  = "Enhancer", "6"  = "Enhancer",
        "7"  = "Enhancer", "8"  = "TssBivalent",
        "9"  = "Heterochromatin", "10" = "Heterochromatin",
        "11" = "Heterochromatin", "12" = "Heterochromatin",
        "13" = "Quiescent")

## Chromatin states _ per tissue
dir <- c("./chromHMM/")
dir_file <- fs::dir_ls(dir, regexp = "\\.bed$")
tissList <- c("CL","CO","CR","KI","LI","LU","MG","OV","PA","SP","ST")

for (i in 1:length(tissList)) { 
  ch <- read.delim(dir_file[i], sep = "\t", header = FALSE)[1:4]
  gr <- GRanges(ch$V1,
                ranges = IRanges(ch$V2, ch$V3),
                states = ch$V4)
  gr$states_simplified = map[gr$states]
  assign(paste0("ch_",tissList[i]), gr)
}

##--------------------##
## color (per states) ##
##--------------------##
states_col = c("Promoter"        = "#f0512f",
               "Enhancer"        = "#4287c6",
               "Heterochromatin" = "#7e5aa5",
               "TssBivalent"     = "yellow") # "Quiescent"       = "gray"

states_name = names(states_col)
n_states = length(states_col)

##-----------------##
## generate matrix ##
##-----------------##
ch_merged <- c(ch_CL, ch_CO, ch_CR, ch_KI, ch_LI,
               ch_LU, ch_MG, ch_OV, ch_PA, ch_SP, ch_ST) ##**
head(ch_merged) 
ch_merged <- ch_merged[which(ch_merged$states_simplified != "Quiescent"),] ##**
mat_merged = normalizeToMatrix(ch_merged, tsDMRgr, value_column = "states_simplified",
                               extend = 5000, mean_mode = "w0", smooth = TRUE, w = 50,
                               empty_value = 0) ##**

##-----------------------##
## draw enriched heatmap ##
##-----------------------##
pdf(paste0(gsub("-","",Sys.Date()), "_1_EH_Total_tsDMR_chromHMM(merged).pdf"))
print(EnrichedHeatmap(mat_merged, name="states", col=states_col, use_raster = TRUE,
                      width = unit(25, "mm"), height = unit(80, "mm"), 
                      column_title = "tsDMR"))
dev.off()
