#===================================================#
## 1. Filtering DMRs with tissue-specific patterns ##
#===================================================#
## all DMR list
setwd("/media/vetbio/Seagate Backup Plus Drive/CCRCD-Beagle/Tissue-specific DMRs/20210711_MethylAction_tsMBD")
library(readr)
tsDMR <- read_csv("20210713_3_dmrs_tsDMR.csv")                                               
colnames(tsDMR)
#
tsDMR$dmrid <- paste0(tsDMR$chr, "_", tsDMR$start)
data <- tsDMR[c("dmrid","chr","start","end","anodev.padj","pattern", "frequent")]

## normalized counts
library(dplyr)
load("D:/CCRCD-Beagle/Tissue-specific DMRs/20210711_MethylAction_tsMBD/20210713_4_ma_100bp_normcounts.Rdata")
colnames(norm)
data <- left_join(data, norm[c("dmrid","CL_679","CL_jj8","CO_679","CO_jj8","CR_679","CR_jj8","KI_679","KI_jj8","LI_679","LI_jj8","LU_679","LU_jj8","MG_244","MG_679","OV_244","OV_679","PA_679","PA_jj8","SP_679","SP_jj8","ST_679","ST_jj8")])

## tsDMR with tissue-specific pattern
library(readxl)
  tsPat <- read_xlsx("20210713_3_dmrs_tsDMR_trimmed_edit.xlsx", sheet = "Tissue-specific")
colnames(tsPat)
colnames(tsPat)[1] <- "pattern"

## save
data <- data[which(data$pattern %in% tsPat$pattern),]
data <- left_join(data, tsPat[c("pattern","tsDMR")])
data <- data[c("dmrid","chr","start","end","anodev.padj","pattern","tsDMR","frequent",
               "CL_679","CL_jj8","CO_679","CO_jj8","CR_679","CR_jj8","KI_679","KI_jj8",
               "LI_679","LI_jj8","LU_679","LU_jj8","MG_244","MG_679","OV_244","OV_679",
               "PA_679","PA_jj8","SP_679","SP_jj8","ST_679","ST_jj8")];print(nrow(data))
data_freq <- data[which(data$frequent == TRUE),];print(nrow(data_freq))

##
write_csv(data, file = "20210715_1_tsDMR_norm.csv")
write_csv(data_freq, file = "20210715_1_tsDMR_norm_freq.csv")
