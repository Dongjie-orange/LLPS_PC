rm(list = ls())
gc()
jco <- c("#709ECD", "#4F9595", "#E5928E", "#8F67B7")

## immune fraction ---------------------------------------
rm(list = ls())
gc()
library(GSEABase)
library(GSVA)
library(tibble)
library(tidyr)

load("train_tcga.rds")
load("merge_gene.rds")
load("../Case12 ES_PC/cellMarker_ssGSEA.Rdata")

ssgseaScore <- gsva(as.matrix(exp_tcga), cellMarker, method = "ssgsea", kcdf = "Gaussian",
                    abs.ranking = TRUE)
ssgseaScore <- t(scale(t(ssgseaScore)))
ssgseaOut <- rbind(id = colnames(ssgseaScore), ssgseaScore)
ssgseaOut <- ssgseaOut[-1, ]
a <- ssgseaOut %>%
  t() %>%
  as.data.frame()
a$LLPS_subtype <- clin_tcga$LLPS_subtype
a <- a %>%
  rownames_to_column("sample")
b <- gather(a, key = ssGSEA, value = Expression, -c(LLPS_subtype, sample))
b$Expression <- as.numeric(b$Expression)

pdf("immune_abundance_boxplot_LLPS_subtype.pdf", width = 12, height = 4.5, onefile = F)
ggboxplot(b, x = "ssGSEA", y = "Expression", fill = "LLPS_subtype",
          palette = jco[1:3], outlier.shape = NA, size = 0.3) + theme_bw() +
  stat_compare_means(aes(group = LLPS_subtype), method = "wilcox.test",
                     label = "p.signif", symnum.args = list(cutpoints = c(0,
                                                                          0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*",
                                                                                                             ""))) + theme(panel.grid.major = element_line(colour = NA),
                                                                                                                           panel.background = element_rect(fill = "transparent", colour = NA),
                                                                                                                           plot.background = element_rect(fill = "transparent", colour = NA),
                                                                                                                           panel.grid.minor = element_blank(), text = element_text(size = 11),
                                                                                                                           axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
dev.off()




## cancer immune cycle ---------------------------------------

load("../Case12 ES_PC/Cancer_Immunity_Cycle_marker_ssGSEA.Rdata")

ssgseaScore <- gsva(as.matrix(exp_tcga), cellMarker, method = "ssgsea", kcdf = "Gaussian",
                    abs.ranking = TRUE)

ssgseaScore <- t(scale(t(ssgseaScore)))
ssgseaOut <- rbind(id = colnames(ssgseaScore), ssgseaScore)
ssgseaOut <- ssgseaOut[-1, ]
a <- ssgseaOut %>%
  t() %>%
  as.data.frame()
a$LLPS_subtype <- clin_tcga$LLPS_subtype
a <- a %>%
  rownames_to_column("sample")
b <- gather(a, key = ssGSEA, value = Expression, -c(LLPS_subtype, sample))
b$Expression <- as.numeric(b$Expression)
b$ssGSEA <- factor(b$ssGSEA, levels = unique(b$ssGSEA))

col <- rev(c("#E64A35", "#4EBBD6", "#0BA087", rep("#3E5388", 17), "#9FCFC1", "#E69E85",
             "#8691B1"))

pdf("Cancer_Immunity_Cycle_boxplot_LLPS_subtype.pdf", width = 12, height = 6, onefile = F)
ggboxplot(b, x = "ssGSEA", y = "Expression", fill = "LLPS_subtype",
          xlab = "Cancer Immunity Cycle", ylab = "Expression",
          palette = jco[1:3], outlier.shape = NA, size = 0.3) +
  theme_bw() + stat_compare_means(aes(group = LLPS_subtype),
                                  method = "wilcox.test", label = "p.signif", symnum.args = list(cutpoints = c(0,
                                                                                                               0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*",
                                                                                                                                                  ""))) + theme(panel.grid.major = element_line(colour = NA),
                                                                                                                                                                panel.background = element_rect(fill = "transparent",
                                                                                                                                                                                                colour = NA), plot.background = element_rect(fill = "transparent",
                                                                                                                                                                                                                                             colour = NA), panel.grid.minor = element_blank(),
                                                                                                                                                                text = element_text(size = 15), legend.position = "right",
                                                                                                                                                                axis.text.x = element_text(angle = 45, hjust = 1, colour = col))
dev.off()


## ICB expression ---------------------------------------

icb_gene <- c("CD274", "CD276", "CTLA4", "HHLA2", "ICOS", "ICOSLG", "PDCD1", "PDCD1LG2",
              "TMIGD2", "VTCN1", "BTLA", "CD27", "CD40", "CD40LG", "CD70", "TNFRSF18", "TNFRSF4",
              "TNFRSF9", "TNFRSF14", "ENTPD1", "HAVCR2", "IDO1", "LAG3", "NCR3", "NT5E", "SIGLEC15")
icb_gene[icb_gene %in% rownames(exp_tcga)]
icb_exp <- as.data.frame(exp_tcga[icb_gene, ])

d <- icb_exp %>%
  t() %>%
  as.data.frame()
d$LLPS_subtype <- clin_tcga$LLPS_subtype
d <- d %>%
  rownames_to_column("sample")

e <- gather(d, key = ICB, value = Expression, -c(LLPS_subtype, sample))
e$Expression <- as.numeric(e$Expression)

pdf("ICB_exp_boxplot_LLPS_subtype.pdf", width = 12, height = 4.5, onefile = F)
ggboxplot(e, x = "ICB", y = "Expression", fill = "LLPS_subtype",
          palette = jco[1:3], outlier.shape = NA, size = 0.3) + theme_bw() +
  stat_compare_means(aes(group = LLPS_subtype), method = "wilcox.test",
                     label = "p.signif", symnum.args = list(cutpoints = c(0,
                                                                          0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*",
                                                                                                             ""))) + theme(panel.grid.major = element_line(colour = NA),
                                                                                                                           panel.background = element_rect(fill = "transparent", colour = NA),
                                                                                                                           plot.background = element_rect(fill = "transparent", colour = NA),
                                                                                                                           panel.grid.minor = element_blank(), text = element_text(size = 11),
                                                                                                                           axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
dev.off()


## Estimate ---------------------------------------
library(stringr)
library(estimate)
library(ggsignif)
library(ggplot2)
library(see)
library(plyr)

write.table(exp_tcga, file = "exp_tcga.txt", sep = "\t", quote = F, row.names = T)
filterCommonGenes(input.f = "exp_tcga.txt", output.f = "gene.gct", id = "GeneSymbol")
estimateScore("gene.gct", "estimat_score.gct", platform = "affymetr")
ESTIMATE_score <- read.table("estimat_score.gct", skip = 2, header = TRUE, row.names = 1,
                             check.names = F)
ESTIMATE_score <- as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)]))
rownames(ESTIMATE_score) <- gsub(".", "-", rownames(ESTIMATE_score), fixed = T)
identical(rownames(ESTIMATE_score), rownames(clin_tcga))
ESTIMATE_score$LLPS_subtype <- factor(clin_tcga$LLPS_subtype,
                                      levels = unique(clin_tcga$LLPS_subtype))

colnames(ESTIMATE_score)
#[1] 'StromalScore'  'ImmuneScore'   'ESTIMATEScore' 'TumorPurity'   'LLPS_subtype' 
pdf("ESTIMATEScore_boxplot_LLPS_subtype.pdf", width = 4.5, height = 4.5, onefile = F)
ggplot(data = ESTIMATE_score, aes(x = LLPS_subtype, y = ESTIMATEScore,
                                  color = LLPS_subtype)) + geom_jitter(alpha = 0.5, size = 2,
                                                                       position = position_jitterdodge(jitter.width = 0.35, jitter.height = 0,
                                                                                                       dodge.width = 0.8)) + geom_boxplot(alpha = 0.2, width = 0.45,
                                                                                                                                          position = position_dodge(width = 0.8), size = 0.75, outlier.colour = NA) +
  geom_violin(alpha = 0.2, width = 0.9, position = position_dodge(width = 0.8),
              size = 0.75) + scale_color_manual(values = jco[1:3]) + theme_classic() +
  theme(legend.position = "none") + theme(text = element_text(size = 15)) +
  ylab("ESTIMATE Score") + ggsignif::geom_signif(comparisons = list(c("LS1",
                                                                      "LS2"), c("LS1", "LS3"), c("LS2", "LS3")), y_position = c(4000,
                                                                                                                                4800, 5600), tip_length = 0, size = 0.5, test = "wilcox.test",
                                                 color = "black", map_signif_level = c(`***` = 0.001, `**` = 0.01,
                                                                                       `*` = 0.05))
dev.off()


## TIDE ---------------------------------------

tide <- read.csv("../Case12 ES_PC/TIDE_result_TCGA.csv")
tide <- tide %>%
  column_to_rownames("Patient")
tide <- tide[rownames(clin_tcga), ]
tide$LLPS_subtype <- factor(clin_tcga$LLPS_subtype,
                            levels = unique(clin_tcga$LLPS_subtype))

pdf("TIDEScore_boxplot_LLPS_subtype.pdf", width = 4.5, height = 4.5, onefile = F)
ggplot(data = tide, aes(x = LLPS_subtype, y = TIDE, color = LLPS_subtype)) +
  geom_jitter(alpha = 0.5, size = 2, position = position_jitterdodge(jitter.width = 0.35,
                                                                     jitter.height = 0, dodge.width = 0.8)) + geom_boxplot(alpha = 0.2, width = 0.45,
                                                                                                                           position = position_dodge(width = 0.8), size = 0.75, outlier.colour = NA) +
  geom_violin(alpha = 0.2, width = 0.9, position = position_dodge(width = 0.8),
              size = 0.75) + scale_color_manual(values = jco[1:3]) + theme_classic() +
  theme(legend.position = "none") + theme(text = element_text(size = 15)) +
  ylab("TIDE") + ggsignif::geom_signif(comparisons = list(c("LS1", "LS2"),
                                                          c("LS1", "LS3"), c("LS2", "LS3")), y_position = c(2, 2.4, 2.8), tip_length = 0,
                                       size = 0.5, test = "t.test", color = "black", map_signif_level = c(`***` = 0.001,
                                                                                                          `**` = 0.01, `*` = 0.05))
dev.off()



## MSI ---------------------------------------

library(cBioPortalData)
cbio <- cBioPortal()
studies <- getStudies(cbio)
head(studies$studyId)
id <- "ucec_tcga_pan_can_atlas_2018"
clinical <- clinicalData(cbio, id)
table(clinical$patientId %in% substr(rownames(clin_tcga), 1, 12))



# IPS (RJ/TCGA)-----------------------------------------------------------------------
#RJ

ips <- read.csv("../Case12 ES_PC/IPS.txt", header = TRUE, sep = "\t", dec = ".",
                check.names = FALSE)
colnames(ips)[1] <- "ID"
ips <- ips %>%
  column_to_rownames("ID")
ips <- ips[rownames(clin_tcga), ]
identical(rownames(clin_tcga), rownames(ips))
ips$LLPS_subtype <- factor(clin_tcga$LLPS_subtype,
                           levels = unique(clin_tcga$LLPS_subtype))

pdf("IPS_boxplot_LLPS_subtype.pdf", width = 4.5, height = 4.5, onefile = F)
ggplot(data = ips, aes(x = LLPS_subtype, y = IPS, color = LLPS_subtype)) +
  geom_jitter(alpha = 0.5, size = 2, position = position_jitterdodge(jitter.width = 0.35,
                                                                     jitter.height = 0, dodge.width = 0.8)) + geom_boxplot(alpha = 0.2,
                                                                                                                           width = 0.45, position = position_dodge(width = 0.8), size = 0.75,
                                                                                                                           outlier.colour = NA) + geom_violin(alpha = 0.2, width = 0.9,
                                                                                                                                                              position = position_dodge(width = 0.8), size = 0.75) +
  scale_color_manual(values = jco[1:3]) + theme_classic() +
  theme(legend.position = "none") + theme(text = element_text(size = 15)) +
  ylab("IPS") + ggsignif::geom_signif(comparisons = list(c("LS1",
                                                           "LS2"), c("LS1", "LS3"), c("LS2", "LS3")), y_position = c(10,
                                                                                                                     10.4, 10.8), tip_length = 0, size = 0.5, test = "t.test",
                                      color = "black", map_signif_level = c(`***` = 0.001, `**` = 0.01,
                                                                            `*` = 0.05))
dev.off()

## submap ---------------------------------------
rm(list = ls())
gc()
library(ComplexHeatmap)
load("train_tcga.rds")

library(readxl)
library(tibble)
library(tidyverse)
library(cmapR)

df <- read_xls("../A02_Templet  File/1.0_TCGA_master/Submap/NIHMS927815-supplement-Tables_S1-S18.xls",
               sheet = 9)
df[df == "NA"] <- NA
df <- df %>%
  column_to_rownames("sample")
marker_gene <- colnames(df)[5:ncol(df)]

df_c <- df %>%
  drop_na("aCTLA4_response")
exp_c <- cbind(Gene = df_c[, 1], df_c[, 5:ncol(df_c)])
pd_c <- df_c$aCTLA4_response
df_p <- df %>%
  drop_na("aPD1_response")
exp_p <- cbind(Gene = df_p[, 1], df_p[, 5:ncol(df_p)])
pd_p <- df_p$aPD1_response
samegene <- intersect(rownames(exp_tcga), marker_gene)

exp_tcga <- exp_tcga[samegene, ]
identical(rownames(clin_tcga), colnames(exp_tcga))

exp_c <- as.data.frame(t(exp_c[, samegene]))
exp_p <- as.data.frame(t(exp_p[, samegene]))

#Expression dataset-
rj <- data.frame(NAME = row.names(exp_tcga), Description = row.names(exp_tcga))
exp_tcga <- cbind(rj, exp_tcga)
row.names(exp_tcga) <- NULL
write.table(exp_tcga, file = "exp_tcga.txt", sep = "\t", row.names = F, col.names = T,
            quote = F)

group <- ifelse(clin_tcga$LLPS_subtype == "LS1", 1, ifelse(clin_tcga$LLPS_subtype ==
                                                             "LS2", 2, 3))
group <- paste(group, collapse = " ")
group <- c(paste(c(176, 3, 1), collapse = " "), "# 1 2 3", group)
write.table(group, file = "group_tcga.cls", col.names = F, row.names = F, quote = F,
            sep = "\t")

#plot
df <- read_xls("submap_result_TCGA.xls", sheet = 1)

submap_plot <- function(df, name)
{
  df <- df %>%
    column_to_rownames("group")
  col_fun <- colorRampPalette(c("#FFFFE3", "#4FA1D0"))(5)
  id <- factor(c(rep("Nominal p value", 3), rep("Bonferroni corrected",
                                                3)), levels = c("Nominal p value", "Bonferroni corrected"))
  left <- rowAnnotation(df = data.frame(pvalue = id),
                        show_annotation_name = T, col = list(pvalue = c(`Nominal p value` = "#E5928E",
                                                                        `Bonferroni corrected` = "#4F9595")),
                        annotation_legend_param = list(labels_gp = gpar(fontsize = 12),
                                                       title_gp = gpar(fontsize = 12, fontface = "bold"),
                                                       ncol = 1), annotation_name_gp = gpar(fontsize = 12,
                                                                                            fontface = "bold"), simple_anno_size = unit(4,
                                                                                                                                        "mm"))
  p <- Heatmap(as.matrix(df), col = col_fun, cell_fun = function(j,
                                                                 i, x, y, width, height, fill)
  {
    if (df[i, j] < 0.05)
    {
      grid.text(paste0("p=", sprintf("%.4f",
                                     df[i, j])), x, y, gp = gpar(fontsize = 12))
    }
  }, name = "value", rect_gp = gpar(col = "white",
                                    lwd = 2), row_split = id, row_title = NULL,
  row_gap = unit(2, "mm"), width = ncol(df) *
    unit(20, "mm"), height = nrow(df) * unit(20,
                                             "mm"), left_annotation = left, row_names_side = "right",
  cluster_rows = F, cluster_columns = F, border = F,
  heatmap_legend_param = list(labels_gp = gpar(fontsize = 12),
                              at = seq(0, 1, 0.2), title_gp = gpar(fontsize = 12,
                                                                   fontface = "bold")))
  pdf(file = paste0(name, "_submap_plot.pdf"), width = 6,
      height = 6, onefile = T)
  print(p)
  dev.off()
}
submap_plot(df, "TCGA")








