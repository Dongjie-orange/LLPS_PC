rm(list = ls())
gc()


## clinical --------------
library(plyr)
library(ggsci)
library(ggplot2)

load("train_tcga.rds")

x <- clin_tcga

for (i in colnames(x)[c(4:15)])
{
  
  x <- x %>%
    drop_na(i)
  ggdata <- table(x$LLPS_subtype, x[, which(colnames(x) == i)]) %>%
    as.data.frame()
  ggdata2 <- ddply(ggdata, "Var1", transform, percent_Freq = Freq/sum(Freq) *
                     100)
  ggdata2 <- ddply(ggdata2, "Var1", transform, lable = cumsum(percent_Freq) -
                     0.5 * percent_Freq)
  ggdata2$ll <- paste0(round(ggdata2$percent_Freq/100, 2) * 100,
                       "%")
  jco <- c("#EFDDBA", "#77BC82", "#9FCFC1", "#C6D7A3")
  
  p <- ggplot(ggdata2, aes(Var1, percent_Freq, fill = Var2)) +
    geom_bar(stat = "identity", width = 0.8) + xlab(NULL) +
    ylab("Fraction (%)") + geom_text(aes(label = ll), size = 3.8,
                                     position = position_stack(vjust = 0.5), color = "black") +
    scale_fill_manual(values = jco[1:2]) + theme_classic() +
    scale_x_discrete(expand = c(0.3, 0.2)) + scale_y_continuous(expand = c(0.01,
                                                                           0)) + theme(legend.position = "right", plot.title = element_text(hjust = 0.5,
                                                                                                                                            face = "plain"), axis.title = element_text(size = 13),
                                                                                       axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12)) +
    labs(fill = i)
  
  print(p)
  ggsave(filename = paste0(i, "_bar.pdf"), height = 5.5, width = 4)
  dev.off()
  
}

## pathway difference --------------
rm(list = ls())
gc()
library(GSEABase)
library(GSVA)
library(GSVAdata)
data(c2BroadSets)
library(tibble)
library(tidyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(limma)

jco <- c("#709ECD", "#4F9595", "#E5928E", "#8F67B7")
load("train_tcga.rds")

hallmarkSet <- getGmt("h.all.v2022.1.Hs.symbols.gmt")

entrez_id <- mapIds(x = org.Hs.eg.db, keys = rownames(exp_tcga), keytype = "SYMBOL",
                    column = "ENTREZID")
entrez_id <- na.omit(entrez_id)
entrez_id <- as.data.frame(entrez_id)
exp2 <- exp_tcga[rownames(entrez_id), ]
rownames(exp2) <- entrez_id$entrez_id

ssgseaScore <- gsva(as.matrix(exp_tcga), hallmarkSet, kcdf = "Gaussian", min.sz = 10,
                    max.sz = 500, verbose = TRUE)


ssgseaOut <- rbind(id = colnames(ssgseaScore), ssgseaScore)
ssgseaOut <- ssgseaOut[-1, ]
source("../A02_Templet  File/0.0_R_Encyclopedia/transfer_to_matrix.R")
ssgseaOut <- numeric_Matrix(ssgseaOut)
identical(colnames(ssgseaOut), rownames(clin_tcga))

group <- clin_tcga$LLPS_subtype
design <- model.matrix(~0 + factor(group))
colnames(design) <- levels(factor(group))
rownames(design) <- rownames(clin_tcga)

cont.wt <- makeContrasts("LS1-LS2", "LS1-LS3", "LS2-LS1", "LS2-LS3", "LS3-LS1", "LS3-LS2",
                         levels = design)

fit <- lmFit(ssgseaOut, design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
tT <- topTableF(fit2, adjust = "BH", sort.by = "F", n = Inf)
tT_diff <- tT[tT$P.Value < 0.05, ]
tT_diff <- tT[order(tT$P.Value), ]

s1_dif <- rownames(tT_diff[tT_diff$LS1.LS2 > 0 & tT_diff$LS1.LS3 > 0, ])
s2_dif <- rownames(tT_diff[tT_diff$LS2.LS1 > 0 & tT_diff$LS2.LS3 > 0, ])
s3_dif <- rownames(tT_diff[tT_diff$LS3.LS1 > 0 & tT_diff$LS3.LS2 > 0, ])

index <- c(s1_dif, s2_dif, s3_dif)
ssgseaOut_sig <- ssgseaOut[index, ]
rownames(ssgseaOut_sig) <- str_split_fixed(rownames(ssgseaOut_sig), "HALLMARK_", 2)[, 2]

top_anno <- HeatmapAnnotation(df = data.frame(LLPS_subtype = clin_tcga$LLPS_subtype),
                              show_annotation_name = F, show_legend = F, col = list(LLPS_subtype = c(LS1 = jco[1],
                                                                                                     LS2 = jco[2], LS3 = jco[3])), simple_anno_size = unit(3, "mm"))
left_anno <- rowAnnotation(df = data.frame(Pathway = c(rep("LS1", 25),
                                                       rep("LS2", 4), rep("LS3", 21))), show_annotation_name = F, show_legend = F,
                           col = list(Pathway = c(LS1 = jco[1], LS2 = jco[2], LS3 = jco[3])),
                           simple_anno_size = unit(3, "mm"))

col_fun <- colorRamp2(c(-1, 0, 1), c("#354D7B", "white", "#ED7D31"))
rowsplit <- factor(c(rep("LS1", 25), rep("LS2", 4), rep("LS3", 21)), levels = c("LS1",
                                                                                "LS2", "LS3"))

dat <- as.matrix(t(scale(t(ssgseaOut_sig))))
dat[dat > 1] <- 1
dat[dat < -1] <- -1

pdf("Heatmap_pathway.pdf", width = 10, height = 8, onefile = F)
Heatmap(dat, name = "Expression", col = col_fun, row_split = rowsplit,
        column_split = clin_tcga$LLPS_subtype, row_title = NULL,
        column_title = c("LS1", "LS2", "LS3"), column_gap = unit(1.5,
                                                                 "mm"), row_gap = unit(1.5, "mm"), row_names_gp = grid::gpar(fontsize = 10),
        top_annotation = top_anno, left_annotation = left_anno, row_names_side = "right",
        cluster_rows = F, row_order = NULL, column_order = NULL,
        show_column_names = FALSE, border = T, border_gp = gpar(col = "black",
                                                                lty = 1), heatmap_legend_param = list(labels_gp = gpar(fontsize = 12),
                                                                                                      title_gp = gpar(fontsize = 12, fontface = "bold")))
dev.off()




# vocalno plot --------------
rm(list = ls())
gc()
library(DESeq2)
library(edgeR)
library(ImageGP)
library(ggplot2)
library(ggpubr)
library(egg)
library(ggrepel)

load("train_tcga.rds")
table(clin_tcga$LLPS_subtype)

ls1 <- rownames(clin_tcga)[1:90]
ls2 <- rownames(clin_tcga)[91:111]
ls3 <- rownames(clin_tcga)[112:176]

# DESeq2
gset <- exp_tcga[rowMeans(exp_tcga) > 0, ]
group_list <- c(rep("B", 90), rep("B", 21), rep("A", 65))
group_list <- factor(group_list, levels = c("A", "B"))
data <- apply(gset, 2, as.integer)
row.names(data) <- row.names(gset)
condition <- group_list
coldata <- data.frame(row.names = colnames(data), condition)
dds <- DESeqDataSetFromMatrix(countData = data, colData = coldata, design = ~condition)
dds$condition <- relevel(dds$condition, ref = "B")
dds <- DESeq(dds)
nrDEG <- as.data.frame(results(dds))
nrDEG_LS3 <- nrDEG
save(nrDEG_LS3, file = "nrDEG_LS3.rds")


# plot
baseA <- counts(dds, normalized = T)[, dds$condition == "A"]
if (is.vector(baseA))
{
  baseMeanA <- as.data.frame(baseA)
} else
{
  baseMeanA <- as.data.frame(rowMeans(baseA))
}
colnames(baseMeanA) <- "A"
head(baseMeanA)

baseB <- counts(dds, normalized = T)[, dds$condition == "B"]
if (is.vector(baseB))
{
  baseMeanB <- as.data.frame(baseB)
} else
{
  baseMeanB <- as.data.frame(rowMeans(baseB))
}
colnames(baseMeanB) <- "B"
head(baseMeanB)

nrDEG$A <- log2(baseMeanA$A + 1)
nrDEG$B <- log2(baseMeanB$B + 1)
nrDEG$level <- ifelse(nrDEG$padj < 0.05, ifelse(nrDEG$log2FoldChange >= 1, "Up",
                                                ifelse(nrDEG$log2FoldChange <= -1, "Down", "NoSig")), "NoSig")
head(nrDEG)


nrDEG$label <- ""
nrDEG <- nrDEG[order(nrDEG$padj), ]
up.genes <- head(rownames(nrDEG)[which(nrDEG$level == "Up")], 5)
down.genes <- head(rownames(nrDEG)[which(nrDEG$level == "Down")], 5)
top10genes <- c(as.character(up.genes), as.character(down.genes))
nrDEG$label[match(top10genes, rownames(nrDEG))] <- top10genes

# plot
pdf(paste0("LS3 vs LS1+LS2", "_volcano.pdf"), width = 5, height = 5, onefile = F)
sp_scatterplot(nrDEG, xvariable = "A", yvariable = "B", color_variable = "level",
               title = "LS3 vs LS1+LS2", color_variable_order = c("NoSig", "Up", "Down"),
               manual_color_vector = c("#AFB6B6", "#E5928E", "#709ECD")) + coord_fixed(1) +
  labs(x = "LS3", y = "LS1+LS2") + theme_bw() + theme(panel.grid = element_blank(),
                                                      axis.text = element_text(size = 11), axis.title = element_text(size = 13),
                                                      legend.title = element_blank(), legend.text = element_text(size = 11),
                                                      legend.position = c(0.995, 0.012), legend.justification = c(1, 0)) +
  geom_text_repel(data = nrDEG, aes(label = label), color = "black", size = 4,
                  fontface = "italic", arrow = arrow(ends = "first", length = unit(0.01,
                                                                                   "npc")), box.padding = 0.2, point.padding = 0.3, segment.color = "black",
                  segment.size = 0.3, force = 1, max.iter = 3000, max.overlaps = Inf)
dev.off()


# vocalno plot --------------
rm(list = ls())
gc()
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)

load("nrDEG_LS1.rds")
load("nrDEG_LS2.rds")
load("nrDEG_LS3.rds")

ls1_gene <- rownames(nrDEG_LS1[nrDEG_LS1$log2FoldChange > 0 & nrDEG_LS1$padj < 0.05, ])
ls2_gene <- rownames(nrDEG_LS2[nrDEG_LS2$log2FoldChange > 0 & nrDEG_LS2$padj < 0.05, ])
ls3_gene <- rownames(nrDEG_LS3[nrDEG_LS3$log2FoldChange > 0 & nrDEG_LS3$padj < 0.05, ])

gene.df <- bitr(ls3_gene, fromType = "SYMBOL", toType = c("ENTREZID", "SYMBOL"),
                OrgDb = org.Hs.eg.db)
gene <- gene.df$ENTREZID
head(gene.df)

# GO
ego <- enrichPathway(gene = gene, organism = "human", pvalueCutoff = 0.05)
go <- ego@result
go_diff <- go[go$pvalue < 0.05, ]
go_diff <- go_diff[order(go_diff$pvalue), ]
go_diff$Description

bardata <- go_diff[c(1, 2, 4, 5, 9, 14, 25, 33, 36, 46), ]

bardata <- bardata[order(bardata$Count, decreasing = T), ]
bardata$Description <- factor(bardata$Description, levels = rev(bardata$Description))

pdf("Reactome pathway in LS2.pdf", width = 6.5, height = 4, onefile = F)
ggplot(bardata, aes(y = Description, x = Count)) + geom_bar(stat = "identity",
                                                            aes(fill = qvalue), width = 0.8) + scale_fill_gradientn(colours = c("#4D6381",
                                                                                                                                "#BAC2CC")) + theme_bw() + theme(panel.grid = element_blank(),
                                                                                                                                                                 axis.text = element_text(size = 11), axis.title = element_text(size = 13)) +
  labs(title = "Reactome pathway in LS2")
dev.off()







