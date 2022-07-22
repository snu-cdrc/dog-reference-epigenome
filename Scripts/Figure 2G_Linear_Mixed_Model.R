library(reshape2)
library(stringr)
library(lme4)
library(ggplot2)



## Dog vs Human

rm(list = ls())
### Input file --------------------------------------------
path = "04_LMM_anlaysis"
in_file_name = "210511_LMM_input.txt"
out_plot_name = "210511_LMM_out_D_H.png"
out_plot_name_pdf = "210511_LMM_out_D_H.pdf"
out_table_name = "210511_LMM_out_D_H.txt"

plot_title = "D_H species"
plot_x_axis = "Fraction of variance across tissues"
plot_y_axis = "Fraction of variance across species"

col = c("seagreen", "goldenrod2", "gray91") # Species / Tissue / None

## Read Data ---------------------------------------------
setwd(path)
in_file = read.table(in_file_name, sep="\t" , header=T, stringsAsFactors = F)

## LMM processing  ---------------------------------------
in_3col <- melt(in_file, id="Dog") # format convert (row x col table -> row_id, col_id, value format)

label <- as.data.frame(str_split(in_3col$variable, "_", simplify = T)) # Human_Colon -> Human, Colon (to split label)
colnames(label) <- c("Species", "Tissue")
in_3col <- cbind(in_3col, label)

## Subsampling
in_3col_sub <- subset(in_3col, Species %in% c("D", "H"))

list_lmm <- list() # make empty list
for(i in 1:length(unique(in_file$Dog))){ # the number of genes
  name <- unique(in_file$Dog)[i]
  mod <- lmer(value~(1|Species) + (1|Tissue), subset(in_3col_sub, Dog==name))
  vardat <- as.data.frame(summary(mod)$varcor)
  vardat$vcov2 <- vardat$vcov/sum(vardat$vcov)
  list_lmm[[i]] <- data.frame(Name=name,
                    Tissue=vardat$vcov2[vardat$grp=="Tissue"],
                    Species=vardat$vcov2[vardat$grp=="Species"])
}

## annotation -------------------------------------------
list_lmm_T <- do.call(rbind, list_lmm) # format convert (change row -> col)
list_lmm_T$group <- "N"
list_lmm_T$group[list_lmm_T$Tissue >= quantile(list_lmm_T$Tissue, 0.75) & list_lmm_T$Tissue > list_lmm_T$Species] <- "Tissue"
list_lmm_T$group[list_lmm_T$Species >= quantile(list_lmm_T$Species, 0.75)] <- "Species"
list_lmm_T$group <- factor(list_lmm_T$group, levels=c("Species", "Tissue", "N"))

## Make plot and table ----------------------------------
pdf(out_plot_name_pdf)
print(
ggplot(list_lmm_T, aes(x=Tissue, y=Species, color=group)) + geom_point() + 
  scale_color_manual(values=col, labels=c("High across species", "High across tissues", "None")) + 
  guides(color=guide_legend("")) + 
  labs(list(x = plot_x_axis, y = plot_y_axis, title = plot_title)) + 
  theme_bw() + 
  theme(legend.position = c(0.7, 0.8), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
)
dev.off()

write.table(list_lmm_T, out_table_name, sep="\t", row.names = F, quote = F)



## Dog vs Mouse

rm(list = ls())
### Input file --------------------------------------------
in_file_name = "210511_LMM_input.txt"
out_plot_name = "210511_LMM_out_D_M.png"
out_plot_name_pdf = "210511_LMM_out_D_M.pdf"
out_table_name = "210511_LMM_out_D_M.txt"

plot_title = "D_M species"
plot_x_axis = "Fraction of variance across tissues"
plot_y_axis = "Fraction of variance across species"

col = c("seagreen", "goldenrod2", "gray91") # Species / Tissue / None

## Read Data ---------------------------------------------
setwd(path)
in_file = read.table(in_file_name, sep="\t" , header=T, stringsAsFactors = F)

## LMM processing  ---------------------------------------
in_3col <- melt(in_file, id="Dog") # format convert (row x col table -> row_id, col_id, value format)

label <- as.data.frame(str_split(in_3col$variable, "_", simplify = T)) # Human_Colon -> Human, Colon (to split label)
colnames(label) <- c("Species", "Tissue")
in_3col <- cbind(in_3col, label)

## Subsampling
in_3col_sub <- subset(in_3col, Species %in% c("D", "M"))

list_lmm <- list() # make empty list
for(i in 1:length(unique(in_file$Dog))){ # the number of genes
  name <- unique(in_file$Dog)[i]
  mod <- lmer(value~(1|Species) + (1|Tissue), subset(in_3col_sub, Dog==name))
  vardat <- as.data.frame(summary(mod)$varcor)
  vardat$vcov2 <- vardat$vcov/sum(vardat$vcov)
  list_lmm[[i]] <- data.frame(Name=name,
                              Tissue=vardat$vcov2[vardat$grp=="Tissue"],
                              Species=vardat$vcov2[vardat$grp=="Species"])
}

## annotation -------------------------------------------
list_lmm_T <- do.call(rbind, list_lmm) # format convert (change row -> col)
list_lmm_T$group <- "N"
list_lmm_T$group[list_lmm_T$Tissue >= quantile(list_lmm_T$Tissue, 0.75) & list_lmm_T$Tissue > list_lmm_T$Species] <- "Tissue"
list_lmm_T$group[list_lmm_T$Species >= quantile(list_lmm_T$Species, 0.75)] <- "Species"
list_lmm_T$group <- factor(list_lmm_T$group, levels=c("Species", "Tissue", "N"))

## Make plot and table ----------------------------------
pdf(out_plot_name_pdf)
print(
  ggplot(list_lmm_T, aes(x=Tissue, y=Species, color=group)) + geom_point() + 
    scale_color_manual(values=col, labels=c("High across species", "High across tissues", "None")) + 
    guides(color=guide_legend("")) + 
    labs(list(x = plot_x_axis, y = plot_y_axis, title = plot_title)) + 
    theme_bw() + 
    theme(legend.position = c(0.7, 0.8), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
)
dev.off()

write.table(list_lmm_T, out_table_name, sep="\t", row.names = F, quote = F)



## Human vs Mouse

rm(list = ls())
### Input file --------------------------------------------
in_file_name = "210511_LMM_input.txt"
out_plot_name = "210511_LMM_out_H_M.png"
out_plot_name_pdf = "210511_LMM_out_H_M.pdf"
out_table_name = "210511_LMM_out_H_M.txt"

plot_title = "H_M species"
plot_x_axis = "Fraction of variance across tissues"
plot_y_axis = "Fraction of variance across species"

col = c("seagreen", "goldenrod2", "gray91") # Species / Tissue / None

## Read Data ---------------------------------------------
setwd(path)
in_file = read.table(in_file_name, sep="\t" , header=T, stringsAsFactors = F)

## LMM processing  ---------------------------------------
in_3col <- melt(in_file, id="Dog") # format convert (row x col table -> row_id, col_id, value format)

label <- as.data.frame(str_split(in_3col$variable, "_", simplify = T)) # Human_Colon -> Human, Colon (to split label)
colnames(label) <- c("Species", "Tissue")
in_3col <- cbind(in_3col, label)

## Subsampling
in_3col_sub <- subset(in_3col, Species %in% c("H", "M"))

list_lmm <- list() # make empty list
for(i in 1:length(unique(in_file$Dog))){ # the number of genes
  name <- unique(in_file$Dog)[i]
  mod <- lmer(value~(1|Species) + (1|Tissue), subset(in_3col_sub, Dog==name))
  vardat <- as.data.frame(summary(mod)$varcor)
  vardat$vcov2 <- vardat$vcov/sum(vardat$vcov)
  list_lmm[[i]] <- data.frame(Name=name,
                              Tissue=vardat$vcov2[vardat$grp=="Tissue"],
                              Species=vardat$vcov2[vardat$grp=="Species"])
}

## annotation -------------------------------------------
list_lmm_T <- do.call(rbind, list_lmm) # format convert (change row -> col)
list_lmm_T$group <- "N"
list_lmm_T$group[list_lmm_T$Tissue >= quantile(list_lmm_T$Tissue, 0.75) & list_lmm_T$Tissue > list_lmm_T$Species] <- "Tissue"
list_lmm_T$group[list_lmm_T$Species >= quantile(list_lmm_T$Species, 0.75)] <- "Species"
list_lmm_T$group <- factor(list_lmm_T$group, levels=c("Species", "Tissue", "N"))

## Make plot and table ----------------------------------
pdf(out_plot_name_pdf)
print(
  ggplot(list_lmm_T, aes(x=Tissue, y=Species, color=group)) + geom_point() + 
    scale_color_manual(values=col, labels=c("High across species", "High across tissues", "None")) + 
    guides(color=guide_legend("")) + 
    labs(list(x = plot_x_axis, y = plot_y_axis, title = plot_title)) + 
    theme_bw() + 
    theme(legend.position = c(0.7, 0.8), panel.grid = element_blank(), plot.title = element_text(hjust = 0.5))
)
dev.off()

write.table(list_lmm_T, out_table_name, sep="\t", row.names = F, quote = F)
