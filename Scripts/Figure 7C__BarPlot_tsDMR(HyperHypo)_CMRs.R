library(readr)
library(stringr)
library(ggplot2)
library(cowplot)

##-------##
## tsDMR ##
##-------##
tsDMR_anno <- read_csv("20210719_1_tsDMR20377_anno.csv")
tsDMR_anno$FC <- ifelse(str_sub(tsDMR_anno$tsDMR,3,3)=="1","hyper","hypo")

met_col <- c("indianred3","lightblue4") # Hyper/Hypo

plotBar <- function(dat, X, lp, title){
  dat2 <- as.data.frame(with(dat, eval(parse(text=paste0("table(",X," = ",X,", FC = FC)")))))
  aes2 <- eval(parse(text=paste0("aes(",X,", y=Freq, fill=FC)")))
  ggp <- ggplot(dat2, aes2) + 
    geom_col(position = "dodge", col = "black") + 
    geom_text(aes(label=Freq), position = position_dodge(0.9), vjust=-1, size=3) + 
    scale_fill_manual(values = met_col) + 
    labs(list(x="",y="Count", title=title)) + 
    #coord_cartesian(xlim=c(0.5,length(levels(dat[,X]))+0.5), ylim=c(0,max(table(dat[,c(X, "group")]))*1.1), expand = F) + 
    theme_classic() + 
    theme(legend.title=element_blank(), 
          legend.position = "top", 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_text(size=rel(1.2)))
  return(ggp)
}

png(paste0(gsub("-","",Sys.Date()),"_2_BarPlot_GR_tsDMR.png"), width=600, height=500, res=100)
print(plot_grid(
  plotBar(tsDMR_anno, X="Annotation2", lp=c(0.2,0.8), title="tissue-specific DMRs")))
dev.off()

png(paste0(gsub("-","",Sys.Date()),"_2_BarPlot_Cpg&Rp_tsDMR.png"), width=600, height=600, res=120)
print(plot_grid(
  plotBar(tsDMR_anno, X="Repeat", lp=c(0.8,0.8), title="Repeats"),
  plotBar(tsDMR_anno, X="CpG_annotation", lp=c(0.8,0.8), title="CpG region")))
dev.off()

##-----------##
## common MR ##
##-----------##
assign("common_anno",get(load("20210721_1_CommonMeth_anno.Rdata")))

plotBar2 <- function(dat, X, lp, title){
  dat2 <- as.data.frame(with(dat, eval(parse(text=paste0("table(",X," = ",X,")")))))
  aes2 <- eval(parse(text=paste0("aes(",X,", y=Freq)")))
  ggp <- ggplot(dat2, aes2) + 
    geom_col(position = "dodge", col = 'black', fill = "mediumorchid4") + 
    geom_text(aes(label=Freq), position = position_dodge(0.9), vjust=-1, size=3) + 
    scale_fill_manual(values = "mediumorchid4") + 
    labs(list(x="",y="Count", title=title)) + 
    #coord_cartesian(xlim=c(0.5,length(levels(dat[,X]))+0.5), ylim=c(0,max(table(dat[,c(X, "group")]))*1.1), expand = F) + 
    theme_classic() + 
    theme(legend.title=element_blank(), 
          legend.position = "top", 
          axis.ticks.x = element_blank(), 
          axis.text.x = element_text(size=rel(1.2)))
  return(ggp)
}

png(paste0(gsub("-","",Sys.Date()),"_2_BarPlot_GR_commonDMR.png"), width=600, height=500, res=100)
print(plot_grid(
  plotBar2(common_anno, X="Annotation2", lp=c(0.2,0.8), title="common methylated regions")))
dev.off()

png(paste0(gsub("-","",Sys.Date()),"_2_BarPlot_Cpg&Rp_commonDMR.png"), width=600, height=600, res=120)
print(plot_grid(
  plotBar2(common_anno, X="Repeat", lp=c(0.8,0.8), title="Repeats"),
  plotBar2(common_anno, X="CpG_annotation", lp=c(0.8,0.8), title="CpG region")))
dev.off()


##-------------##
# Export to pdf #
##-------------##
pdf(paste0(gsub("-","",Sys.Date()),"_2_BarPlot_GR_common&tsDMR.pdf"), width=12, height=12)
print(plot_grid(
  nrow = 2, ncol = 2,
  plotBar2(common_anno, X="Annotation2", lp=c(0.2,0.8), title="common methylated regions"),
  plotBar(tsDMR_anno, X="Annotation2", lp=c(0.2,0.8), title="tissue-specific DMRs"),
  
  plot_grid(
    ncol = 2,
    plotBar2(common_anno, X="Repeat", lp=c(0.8,0.8), title="Repeats"),
    plotBar(tsDMR_anno, X="Repeat", lp=c(0.8,0.8), title="Repeats")),
  
  plot_grid(
    ncol = 2,
    plotBar2(common_anno, X="CpG_annotation", lp=c(0.8,0.8), title="CpG region"),
    plotBar(tsDMR_anno, X="CpG_annotation", lp=c(0.8,0.8), title="CpG region"))))

dev.off()
