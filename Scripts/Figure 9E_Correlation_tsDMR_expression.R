library(readr)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(ggformula)
library(ggpmisc)
library(ggsci)
library(gridExtra)

# tsDMR list
tsDMR <- read_csv("20210927_1_tsDMR20377_anno5.csv")
tsDMR <- tsDMR[-grep("CRCL", tsDMR$tsDMR),] ##**except 'CRCL1'
colnames(tsDMR)
tsDMR$CpG_annotation <- ifelse(tsDMR$CpG_annotation == "-", "Non-CGI", tsDMR$CpG_annotation)
tsDMR$CpG_annotation <- ifelse(tsDMR$CpG_annotation == "CpG Island", "CGI", tsDMR$CpG_annotation)
tsDMR$CpG_annotation <- factor(tsDMR$CpG_annotation,
                               levels = c("CGI", "Shore", "Shelf", "Non-CGI"))

# target tissue
unique(tsDMR$tsDMR)
tisList <- c("CL","CO","CR","KI","LI","LU","MG","OV","PA","SP","ST")

for (i in 1:length(tisList)) {
  
  tissue <- tisList[i]
  print(tissue)

  ##-----------##
  # methylation #
  ##-----------##
  dat1 <- tsDMR[grep(tissue, tsDMR$tsDMR),] ##**
  dat1["meanM"] <- rowMeans(dat1[grep(paste0("^",tissue), names(dat1))]) ##**
  colnames(dat1)
  dat1["meanM_others"] <- rowMeans(dat1[-c(1:8, ##**
                                           grep(paste0("^",tissue), names(dat1)), ##**
                                           31:ncol(dat1))]) ##**
  # except intergenic
  dat2 <- dat1[which(dat1$Annotation2 != "Intergenic"),]
  dat2 <- dat2[c("Ensembl","Annotation2","CpG_annotation",
                 "meanM", "meanM_others")] ##**
  
  # aggregate
  colnames(dat2)
  dat3 <- aggregate(cbind(meanM, meanM_others) ~ ., ##**
                    data=dat2, sum)  ##**
  
  # calculate log2FC
  dat3["log2FC_M"] <- log2((dat3$meanM+1)/(dat3$meanM_others+1)) ##** counts+1
  
  ##---------------##
  # gene expression #
  ##---------------##
  exp <- read.delim("Dog_11tissue_RNA_FPKM_raw.txt")
  colnames(exp)
  colnames(exp)[1] <- "Ensembl"
  exp["meanE"] <- exp[grep(paste0("^",tissue), names(exp))] ##**
  exp["meanE_others"] <- rowMeans(exp[-c(1,grep(paste0("^",tissue), names(exp)),ncol(exp))]) ##**
  exp["log2FC_E"] <- log2((exp$meanE+0.1)/(exp$meanE_others+0.1)) ##** fpkm+0.01
  colnames(exp)
  exp2 <- exp[c("Ensembl", "meanE", "meanE_others","log2FC_E")] ##**
  
  ##-----------------------------##
  # methylation + gene expression #
  ##-----------------------------##
  colnames(dat3);colnames(exp2)
  assign(tissue, left_join(dat3, exp2, by = "Ensembl"))
  
}

dat <- rbind(CL,CO,CR,KI,LI,LU,MG,OV,PA,SP,ST) ##**
dat$Annotation2 <- factor(dat$Annotation2,
                          levels = c("promoter-TSS", "intron", "exon", "TTS")) ##**

##--------##
# plotting #
##--------##
source("20190417_Stat_Smooth_Func_with_Pvalue(2).R")
corPlot <- ggplot(dat, aes(x=log2FC_M, y=log2FC_E)) +
  #geom_text(aes(label=Ensembl), vjust=-1.1, colour="grey35") +
  #theme_wsj() +
  #theme_calc(base_size=20) +
  theme_bw() +
  scale_color_npg() + 
  
  ggtitle(paste0("tsDMGs (",nrow(dat),")")) + ##**
  
  theme(plot.title = element_text(hjust = 0.5, size = 20)) +
  xlab("Methylation (log2(F.C))") +
  ylab("Gene Expression (log2(F.C))") +
  #geom_point(aes(colour=Annotation2, shape=CpG_annotation), size=1.5) +
  geom_point(aes(colour=Annotation2), size=1.5) +
  geom_smooth(method = "lm", colour = "red") +
  stat_smooth_func_with_pval(geom="text", method="lm",
                             hjust=0.25, vjust=1,
                             parse=TRUE, size=2.5)

png(paste0(gsub("-","",Sys.Date()),"_1_corPlot_total-tsDMGs_log2FC.png"),
    width = 3000, height = 2400, res = 360)
print(corPlot+facet_grid(Annotation2~CpG_annotation, scales = "free"))
dev.off()

pdf(paste0(gsub("-","",Sys.Date()),"_1_corPlot_total-tsDMGs_log2FC.pdf"),
    width = 8, height = 6)
print(corPlot+facet_grid(Annotation2~CpG_annotation, scales = "free"))
dev.off()

