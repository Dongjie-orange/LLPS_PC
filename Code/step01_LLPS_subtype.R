rm(list = ls())
gc()
library(data.table)
library(readxl)
library(NMF)
library(tibble)
library(stringr)
library(ggVennDiagram)
library(ggplot2)
library(dplyr)

## identifcation of LLPS genes ----------

# llps gene
llps_gene <- read_xlsx("相分离基因列表.xlsx", 1, col_names = F)
llps_gene <- llps_gene$...1
# cancer-driven gene
cancer_gene <- read_xlsx("driver_gene.xlsx", 1, col_names = T)
cancer_gene <- cancer_gene$Symbol
# deg
load("../Case11 immune_MS_bulk_sc_seq/DEG_T_vs_N.Rdata")
deg <- rownames(DEG[DEG$adj.P.Val < 0.05, ])
merge_gene <- Reduce(intersect, list(llps_gene, deg))
save(llps_gene, merge_gene, file = "merge_gene.rds")
# veen plot
x1 <- list(llps_gene, deg)
name <- c("LLPS gene", "DEG")
ggVennDiagram(x1, category.names = name, label = "count", label_alpha = 0,
              set_size = 5, edge_lty = "solid", edge_size = 1) + scale_fill_gradient(low = "white",
                                                                                     high = "white", guide = "none") + scale_color_manual(values = c("#709ECD",
                                                                                                                                                     "#E5928E"))


## classification via NMF ----------

# NMF
rm(list = ls())
gc()
load("merge_gene.rds")
load("train_tcga.rds")
exp_tcga <- exp_tcga[apply(exp_tcga, 1, function(x) sum(x == 0) < 0.25 *
                             ncol(exp_tcga)), ]
rt <- exp_tcga[rownames(exp_tcga) %in% llps_gene, ]

ranks <- 2:5
estim.coad <- nmf(rt, ranks, nrun = 50, method = "brunet")
pdf("NMF_rank_survey.pdf", width = 6.5, height = 5.5, onefile = F)
plot(estim.coad)
dev.off()

seed <- 20230330
nmf.rank <- nmf(rt, rank = 3, nrun = 50, seed = seed, method = "brunet")
group <- predict(nmf.rank)
table(group)
LLPS_subtype <- as.data.frame(group)
LLPS_subtype$group <- ifelse(group == 1, "LS1", ifelse(group == 2, "LS2",
                                                       "LS3"))
table(LLPS_subtype$group)
# LS1 LS2 LS3 90 21 65
save(LLPS_subtype, file = "LLPS_subtype.rds")

jco <- c("#709ECD", "#4F9595", "#E5928E", "#8F67B7")
pdf("Consensus_map_of_NMF.pdf", width = 6.5, height = 6.5, onefile = F)
consensusmap(nmf.rank, labRow = NA, labCol = NA, annCol = data.frame(cluster = group),
             annColors = list(cluster = c(`1` = jco[1], `2` = jco[2], `3` = jco[3])))
dev.off()

## PCA plot ----------
library(Rtsne)

tsne_out <- Rtsne(t(rt), perplexity = 20)
pdat <- data.frame(tsne_out$Y, factor(LLPS_subtype$group))
colnames(pdat) <- c("Y1", "Y2", "group")

pdf("tSNE_plot.pdf", width = 5, height = 5, onefile = F)
ggplot(data = pdat, aes(x = Y1, y = Y2, color = LLPS_subtype$group)) +
  geom_point(size = 2) + labs(x = "tSNE_1", y = "tSNE_2", color = "LLPS_subtype") +
  guides(fill = "none") + theme_bw() + scale_fill_manual(values = jco[1:3]) +
  scale_colour_manual(values = jco[1:3]) + stat_ellipse(aes(color = LLPS_subtype$group,
                                                            fill = LLPS_subtype$group), geom = "polygon", alpha = 0.1, linetype = 1) +
  theme(plot.title = element_text(hjust = 0.5, size = 15), axis.text = element_text(size = 11),
        axis.title = element_text(size = 13), legend.text = element_text(size = 11),
        plot.margin = unit(c(0.4, 0.4, 0.4, 0.4), "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), legend.position = c(0.85, 0.12))
dev.off()


## KM plot ----------
identical(rownames(LLPS_subtype), rownames(clin_tcga))
clin_tcga$LLPS_subtype <- LLPS_subtype$group
save(exp_tcga, clin_tcga, file = "train_tcga_paad.rds")
library(survival)
library(survminer)
sfit <- survfit(Surv(OS.time, OS) ~ LLPS_subtype, data = clin_tcga)

pdf("KM_plot_LLPS_subtype.pdf", width = 6.5, height = 6.5, onefile = F)
ggsurvplot(sfit, pval = T, palette = jco[1:3], risk.table = T)
dev.off()


## Gene heatmap ----------
library(pheatmap)
clin_tcga <- clin_tcga[order(clin_tcga$LLPS_subtype), ]

dat <- sweep(rt, 1, apply(rt, 1, median, na.rm = T))
dat[dat < -2] <- -2
dat <- dat[, rownames(clin_tcga)]

bk <- c(seq(-2, -0.1, by = 0.01), seq(0, 2, by = 0.01))
annotation_col <- data.frame(row.names = rownames(clin_tcga), LLPS_subtype = clin_tcga$LLPS_subtype)
ann_colors <- list(LLPS_subtype = c(LS1 = jco[1], LS2 = jco[2], LS3 = jco[3]))

pdf("LLPS_gene_heatmap.pdf", width = 6.5, height = 5, onefile = F)
pheatmap(dat, scale = "row", show_rownames = F, show_colnames = F, annotation_colors = ann_colors,
         annotation_col = annotation_col, cluster_row = T, cluster_col = F,
         color = c(colorRampPalette(colors = c("#354D7B", "white"))(length(bk)/2),
                   colorRampPalette(color = c("white", "#ED7D31"))(length(bk)/2)),
         legend_breaks = seq(-2, 2, 1), breaks = bk)
dev.off()



