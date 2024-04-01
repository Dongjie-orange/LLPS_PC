rm(list = ls())
gc()
library(ComplexHeatmap)
library(maftools)
library(tibble)
library(ggsci)


# SNP ------------------------------

load("../A01_Data_Processing/PAAD/TCGA/tcga_tmb_snp_cnv.rds")
load("../A01_Data_Processing/PAAD/TCGA/mutsig_for_heatmap.rds")
load("train_tcga.rds")
clin_tcga <- clin_tcga[colnames(TCGA_SNP), ]
identical(rownames(clin_tcga), TCGA_TMB$ID)

dat_top_anno <- data.frame(TMB = TCGA_TMB$TMB, Group = clin_tcga$LLPS_subtype,
                           row.names = rownames(clin_tcga))
dat_top_anno <- dat_top_anno[order(dat_top_anno$Group), ]

dat <- TCGA_SNP
dat <- dat[, rownames(dat_top_anno)]

col <- c(Missense_Mutation = "#a6bddb", Splice_Site = "#beaed4",
         Nonsense_Mutation = "#fdc086", Frame_Shift_Del = "#a6d854", Frame_Shift_Ins = "#e78ac3",
         Multi_Hit = "#fc8d62", In_Frame_Del = "#0BA087", In_Frame_Ins = "#E69E85",
         Gain = "#E64A35", Loss = "#4EBBD6")
alter_fun <- list(background = alter_graphic("rect", fill = "#e9e9e9"),
                  Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
                  Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
                  Nonsense_Mutation = alter_graphic("rect", fill = col["Nonsense_Mutation"],
                  ), Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
                  Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
                  Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"]),
                  In_Frame_Del = alter_graphic("rect", fill = col["In_Frame_Del"]),
                  In_Frame_Ins = alter_graphic("rect", fill = col["In_Frame_Ins"]),
                  Gain = alter_graphic("rect", fill = col["Gain"]), Loss = alter_graphic("rect",
                                                                                         fill = col["Loss"]))

#top annotation
top_anno <- HeatmapAnnotation(TMB = anno_barplot(dat_top_anno$TMB, border = F,
                                                 bar_width = 0.6, axis = T, gp = gpar(fill = "#3070B0", col = NA),
                                                 axis_param = list(gp = gpar(fontsize = 15))), Group = dat_top_anno$Group,
                              col = list(Group = c(LS1 = "#709ECD", LS2 = "#4F9595", LS3 = "#E5928E")),
                              annotation_height = unit(c(3.8, 0.8), "cm"), show_legend = F, show_annotation_name = T,
                              gap = unit(1, "mm"), annotation_label = gt_render(c("TMB", "LLPS subtype"),
                                                                                gp = gpar(fontsize = 20, fontface = "bold")))

lgd_list <- list(Legend(labels = c("LS1", "LS2", "LS3"), title = "LLPS subtype",
                        legend_gp = gpar(fill = c("#709ECD", "#4F9595", "#E5928E")),
                        labels_gp = gpar(fontsize = 16), grid_height = unit(0.5, "cm"),
                        grid_width = unit(0.5, "cm"), title_gp = gpar(fontsize = 16,
                                                                      fontface = "bold")), Legend(labels = c("Missense_Mutation",
                                                                                                             "Multi_Hit", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                                                                                                             "Splice_Site", "In_Frame_Del", "In_Frame_Ins"), title = "Alterations",
                                                                                                  legend_gp = gpar(fill = c("#a6bddb", "#beaed4", "#fdc086", "#a6d854",
                                                                                                                            "#e78ac3", "#fc8d62", "#0BA087", "#E69E85")), labels_gp = gpar(fontsize = 16),
                                                                                                  grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"),
                                                                                                  title_gp = gpar(fontsize = 16, fontface = "bold")))

dat_right_anno <- dat
dat_right_anno[dat_right_anno == ""] <- 0
dat_right_anno[!dat_right_anno == 0] <- 1
dimnames <- list(rownames(dat_right_anno), colnames(dat_right_anno))
dat_right_anno <- matrix(as.numeric(as.matrix(dat_right_anno)),
                         nrow = nrow(dat_right_anno), dimnames = dimnames)
table(dat_top_anno$Group)

df <- data.frame(row.names = rownames(dat), LS1 = apply(dat_right_anno[, 1:83],
                                                        1, sum)/apply(dat_right_anno, 1, sum), LS2 = apply(dat_right_anno[, 84:92],
                                                                                                           1, sum)/apply(dat_right_anno, 1, sum), LS3 = apply(dat_right_anno[, 93:143],
                                                                                                                                                              1, sum)/apply(dat_right_anno, 1, sum))

right_anno <- rowAnnotation(anno = anno_barplot(df, border = F, bar_width = 0.6,
                                                axis = T, gp = gpar(fill = c("#709ECD", "#4F9595", "#E5928E"), col = NA),
                                                axis_param = list(gp = gpar(fontsize = 18))), show_legend = F, show_annotation_name = F,
                            width = unit(4, "cm"))

# plot
pdf("SNP_oncoplot.pdf", width = 18, height = 18)
ht_list <- oncoPrint(dat, alter_fun = alter_fun, alter_fun_is_vectorized = FALSE,
                     width = unit(6 * 4, "cm"), height = unit(4 * 4, "cm"), pct_gp = gpar(fontsize = 16),
                     col = col, row_title = NULL, column_title = NULL, row_gap = unit(3, "mm"),
                     column_gap = unit(2, "mm"), row_order = 1:nrow(dat), column_order = NULL,
                     row_split = NULL, column_split = dat_top_anno$Group, top_annotation = top_anno,
                     right_annotation = right_anno, show_heatmap_legend = F)
draw(ht_list, heatmap_legend_list = lgd_list)
dev.off()


# CNV ------------------------------
rm(list = ls())
gc()

load("../A01_Data_Processing/PAAD/TCGA/tcga_tmb_snp_cnv.rds")
load("../A01_Data_Processing/PAAD/TCGA/mutsig_for_heatmap.rds")
load("train_tcga.rds")
clin_tcga <- clin_tcga[colnames(TCGA_SNP), ]
identical(rownames(clin_tcga), TCGA_TMB$ID)

dat_top_anno <- data.frame(TMB = TCGA_TMB$TMB, Group = clin_tcga$LLPS_subtype,
                           row.names = rownames(clin_tcga))
dat_top_anno <- dat_top_anno[order(dat_top_anno$Group), ]

dat <- TCGA_CNV
dat <- dat[, rownames(dat_top_anno)]
dat[dat == 0] <- ""
dat[dat == 1] <- "Gain"
dat[dat == -1] <- "Loss"
col <- c(Gain = "#E64A35", Loss = "#4EBBD6")
alter_fun <- list(background = alter_graphic("rect", fill = "#e9e9e9"),
                  Gain = alter_graphic("rect", fill = col["Gain"]), Loss = alter_graphic("rect",
                                                                                         fill = col["Loss"]))

#top annotation
top_anno <- HeatmapAnnotation(TMB = anno_barplot(dat_top_anno$TMB, border = F,
                                                 bar_width = 0.6, axis = T, gp = gpar(fill = "#3070B0", col = NA),
                                                 axis_param = list(gp = gpar(fontsize = 15))), Group = dat_top_anno$Group,
                              col = list(Group = c(LS1 = "#709ECD", LS2 = "#4F9595", LS3 = "#E5928E")),
                              annotation_height = unit(c(3.8, 0.8), "cm"), show_legend = F, show_annotation_name = T,
                              gap = unit(1, "mm"), annotation_label = gt_render(c("TMB", "LLPS subtype"),
                                                                                gp = gpar(fontsize = 20, fontface = "bold")))

lgd_list <- list(Legend(labels = c("LS1", "LS2", "LS3"), title = "LLPS subtype",
                        legend_gp = gpar(fill = c("#709ECD", "#4F9595", "#E5928E")),
                        labels_gp = gpar(fontsize = 16), grid_height = unit(0.5, "cm"),
                        grid_width = unit(0.5, "cm"), title_gp = gpar(fontsize = 16,
                                                                      fontface = "bold")), Legend(labels = c("Gain", "Loss"), title = "CNV (arm-level)",
                                                                                                  legend_gp = gpar(fill = c("#E64A35", "#4EBBD6")), labels_gp = gpar(fontsize = 16),
                                                                                                  grid_height = unit(0.5, "cm"), grid_width = unit(0.5, "cm"),
                                                                                                  title_gp = gpar(fontsize = 16, fontface = "bold")))

dat_right_anno <- dat
dat_right_anno[dat_right_anno == ""] <- 0
dat_right_anno[!dat_right_anno == 0] <- 1
dimnames <- list(rownames(dat_right_anno), colnames(dat_right_anno))
dat_right_anno <- matrix(as.numeric(as.matrix(dat_right_anno)),
                         nrow = nrow(dat_right_anno), dimnames = dimnames)
table(dat_top_anno$Group)

df <- data.frame(row.names = rownames(dat), LS1 = apply(dat_right_anno[, 1:83],
                                                        1, sum)/apply(dat_right_anno, 1, sum), LS2 = apply(dat_right_anno[, 84:92],
                                                                                                           1, sum)/apply(dat_right_anno, 1, sum), LS3 = apply(dat_right_anno[, 93:143],
                                                                                                                                                              1, sum)/apply(dat_right_anno, 1, sum))

right_anno <- rowAnnotation(anno = anno_barplot(df, border = F, bar_width = 0.6,
                                                axis = T, gp = gpar(fill = c("#709ECD", "#4F9595", "#E5928E"), col = NA),
                                                axis_param = list(gp = gpar(fontsize = 18))), show_legend = F, show_annotation_name = F,
                            width = unit(4, "cm"))

# plot
pdf("CNV_oncoplot.pdf", width = 20, height = 18)
ht_list <- oncoPrint(dat, alter_fun = alter_fun, alter_fun_is_vectorized = FALSE,
                     width = unit(6 * 4, "cm"), height = unit(9 * 4, "cm"), pct_gp = gpar(fontsize = 16),
                     col = col, row_title = NULL, column_title = NULL, row_gap = unit(3, "mm"),
                     column_gap = unit(2, "mm"), row_order = 1:nrow(dat), column_order = NULL,
                     row_split = NULL, column_split = dat_top_anno$Group, top_annotation = top_anno,
                     right_annotation = right_anno, show_heatmap_legend = F)
draw(ht_list, heatmap_legend_list = lgd_list)
dev.off()



# CNV ------------------------------
rm(list = ls())
gc()
library(ggplot2)
library(ggpubr)
library(ggsignif)

load("../Case12 ES_PC/dat_tmb_amp_del_arm_focal.rds")
load("train_tcga.rds")

clin_tcga <- clin_tcga[rownames(dat), ]
dat$LLPS_subtype <- clin_tcga$LLPS_subtype

i <- "TMB"

for (i in colnames(dat)[c(8:11)])
{
  
  p <- ggplot(data = dat) + geom_boxplot(mapping = aes(x = LLPS_subtype,
                                                       y = dat[, i], colour = LLPS_subtype), alpha = 0.5,
                                         size = 1, width = 0.5, outlier.shape = 1, outlier.size = 1) +
    geom_jitter(mapping = aes(x = LLPS_subtype, y = dat[,
                                                        i], colour = LLPS_subtype), alpha = 0.3, size = 3) +
    scale_color_manual(limits = c("LS1", "LS2", "LS3"),
                       values = c("#709ECD", "#4F9595", "#E5928E")) +
    scale_y_continuous(labels = scales::label_comma(accuracy = 0.1)) +
    geom_signif(mapping = aes(x = LLPS_subtype, y = dat[,
                                                        i]), comparisons = list(c("LS1", "LS2"), c("LS1",
                                                                                                   "LS3"), c("LS2", "LS3")), map_signif_level = T,
                tip_length = 0, y_position = c(max(dat[, i]) +
                                                 0.2, max(dat[, i]) + 0.4, max(dat[, i]) + 0.6),
                size = 0.6, textsize = 4, test = "wilcox.test") +
    theme_classic() + labs(title = i, x = "", y = "CNV burden (log2)") +
    theme(plot.title = element_text(size = 15, colour = "black",
                                    hjust = 0.5), axis.title = element_text(size = 15,
                                                                            color = "black", vjust = 1.9, hjust = 0.5, angle = 90),
          legend.position = "none", axis.text = element_text(size = 13,
                                                             color = "black", vjust = 0.5, hjust = 0.5,
                                                             angle = 0), panel.grid = element_blank())
  
  pdf(file = paste0(i, "_LLPS_subtype.pdf"), width = 3.5,
      height = 5.5, onefile = T)
  print(p)
  dev.off()
}












