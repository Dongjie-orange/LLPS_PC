rm(list = ls())
gc()
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(ggalluvial)


############################  clinical  ########################################
rm(list = ls())
gc()
load("train_tcga.rds")
clin_tcga <- clin_tcga[order(clin_tcga$riskscore), ]

# distribution plot
dat_plot1 <- clin_tcga
dat_plot1$id2 <- paste("p", 1:176, sep = "")
dat_plot1$id2 <- factor(dat_plot1$id2, levels = paste("p", 1:176, sep = ""))
p1 <- ggplot(data = dat_plot1, aes(x = id2, y = riskscore,
                                   fill = id2)) + geom_bar(stat = "identity", fill = ifelse(dat_plot1$id2 %in%
                                                                                              paste("p", 1:88, sep = ""), "#709ECD", "#E5928E"),
                                                           width = 2, position = position_dodge(width = 0.9)) +
  theme(panel.background = element_rect(fill = "transparent",
                                        color = "gray"), axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) + ylim(0, 1) +
  theme(panel.grid = element_blank()) + guides(fill = FALSE) +
  xlab(NULL)
pdf("riskscore_disribution.pdf", width = 7.5, height = 2.5, onefile = F)
print(p1)
dev.off()

#clinical plot
dat_plot2 <- clin_tcga[, c(1, 4:16)]
dat_plot2$OS <- ifelse(dat_plot2$OS == 0, "Alive", "Dead")
dat_plot2$patient <- paste("P", seq(1:176), sep = "")
dat_plot2 <- melt(dat_plot2, id = "patient")
dat_plot2$patient <- factor(dat_plot2$patient, levels = paste("P", seq(1:176), sep = ""))
dat_plot2$variable <- factor(dat_plot2$variable, levels = rev(unique(dat_plot2$variable)))

unique(dat_plot2$variable)
#Levels: Age Alcohol Diabetes Grade M N T Relapse Radiation Gender Site Stage LLPS_subtype

cols <- c(`NA` = "lightgrey", Dead = "#CD3022", Alive = "#78B489", `<=65` = "#EE8784",
          `>65` = "#8BB2F9", Yes = "#F9A96C", No = "#F1F892", `G1+G2` = "#C998C1",
          `G3+G4` = "#456093", M0 = "#78B489", M1 = "#CD3022", N0 = "#D80EC1", N1 = "#A1D5EE",
          `T1+T2` = "#DCD55B", `T3+T4` = "#EE1299", Progression = "#F09E38", Stable = "#67C9FA",
          Response = "#2A64F6", Male = "#39807E", Female = "#D86D2B", `Body/Tail` = "#351D7D",
          Head = "#DCD55B", `I+II` = "#56B93F", `III+IV` = "#BD382F", LS1 = "#709ECD",
          LS2 = "#4F9595", LS3 = "#E5928E")

p2 <- dat_plot2 %>%
  ggplot(aes(x = patient, y = variable)) + geom_tile(aes(fill = value)) +
  scale_x_discrete("") + scale_y_discrete("") + scale_fill_manual(values = cols) +
  theme(axis.text.x.bottom = element_text(size = 0),
        axis.text.y.left = element_text(size = 12), axis.ticks = element_blank(),
        legend.position = "bottom")
pdf("clinical_disribution.pdf", width = 8.5, height = 6, onefile = F)
print(p2)
dev.off()

# sankey plot

sanky <- clin_tcga[, c(18, 16, 1)]
sanky$OS <- ifelse(sanky$OS == 0, "Alive", "Dead")
corLodes <- to_lodes_form(sanky, axes = 1:ncol(sanky), id = "Cohort")
colors <- c("#709ECD", "#E5928E", "#709ECD", "#4F9595", "#E5928E", "#709ECD", "#E5928E")
mycol <- colors

ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort, fill = stratum,
                     label = stratum)) + scale_x_discrete(expand = c(0, 0)) + geom_flow(width = 1/10,
                                                                                        aes.flow = "forward") + geom_stratum(alpha = 0.8, width = 1/10, color = NA) +
  scale_fill_manual(values = mycol) + geom_text(stat = "stratum", size = 4,
                                                color = "black", angle = 90) + xlab("") + ylab("") + theme_bw() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.y = element_blank()) + theme(panel.grid = element_blank()) +
  theme(panel.border = element_blank()) + ggtitle("") + guides(fill = FALSE)




############################  pathway correlation  ########################################

rm(list = ls())
gc()
library(GSEABase)
library(GSVA)
library(GSVAdata)
library(psych)
library(corrplot)
load("train_tcga.rds")

# GSVA
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
ssgseaOut <- as.data.frame(t(ssgseaOut))

mantel <- corr.test(clin_tcga$riskscore, ssgseaOut, method = "spearman")
res <- as.data.frame(cbind(t(as.data.frame(mantel$r)),
                           t(as.data.frame(mantel$p.adj)))) %>%
  rownames_to_column("env")
colnames(res)[2:3] <- c("r", "p")
res <- res[order(res$r, decreasing = T), ]
res <- res[c(1, 3:6, 9:17, 19), ]

ssgseaOut <- ssgseaOut[, res$env]
colnames(ssgseaOut) <- str_split_fixed(colnames(ssgseaOut), "HALLMARK_", 2)[, 2]
ssgseaOut$Riskscore <- clin_tcga$riskscore


res1 <- cor.mtest(ssgseaOut, conf.level = 0.95)
pdf(file = "hallmark_corrplot.pdf", width = 10, height = 10, onefile = F)
corrplot(cor(ssgseaOut, method = "spearman"), type = "upper", tl.col = "black",
         p.mat = res1$p, insig = "label_sig", sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.8,
         pch.col = "black", col = colorRampPalette(c("white", "#4FA1D0"))(10))
dev.off()


############################  immune correlation  ########################################
rm(list = ls())
gc()
library(GSEABase)
library(GSVA)
library(GSVAdata)
library(psych)
library(corrplot)
load("train_tcga.rds")
load("../Case12 ES_PC/cellMarker_ssGSEA.Rdata")


ssgseaScore <- gsva(as.matrix(exp_tcga), cellMarker, kcdf = "Gaussian", min.sz = 10,
                    max.sz = 500, verbose = TRUE)

ssgseaOut <- rbind(id = colnames(ssgseaScore), ssgseaScore)
ssgseaOut <- ssgseaOut[-1, ]
source("../A02_Templet  File/0.0_R_Encyclopedia/transfer_to_matrix.R")
ssgseaOut <- numeric_Matrix(ssgseaOut)
identical(colnames(ssgseaOut), rownames(clin_tcga))
ssgseaOut <- as.data.frame(t(ssgseaOut))
ssgseaOut$Riskscore <- as.numeric(clin_tcga$riskscore)


res1 <- cor.mtest(ssgseaOut, conf.level = 0.95)
pdf(file = "immune_corrplot.pdf", width = 15, height = 15, onefile = F)
corrplot(cor(ssgseaOut, method = "spearman"), type = "upper", tl.col = "black",
         p.mat = res1$p, insig = "label_sig", sig.level = c(0.001, 0.01, 0.05), pch.cex = 0.8,
         pch.col = "black", col = colorRampPalette(c("white", "#4FA1D0"))(10))
dev.off()



############################  ESTIMATE  ########################################

rm(list = ls())
gc()
library(stringr)
library(estimate)
library(ggsignif)
library(ggplot2)
library(see)
library(plyr)

load("train_tcga.rds")

ESTIMATE_score <- read.table("estimat_score.gct", skip = 2, header = TRUE, row.names = 1,
                             check.names = F)
ESTIMATE_score <- as.data.frame(t(ESTIMATE_score[, 2:ncol(ESTIMATE_score)]))
rownames(ESTIMATE_score) <- gsub(".", "-", rownames(ESTIMATE_score), fixed = T)
identical(rownames(ESTIMATE_score), rownames(clin_tcga))
ESTIMATE_score$Group <- clin_tcga$group

colnames(ESTIMATE_score)
# [1] 'StromalScore'  'ImmuneScore'   'ESTIMATEScore' 'TumorPurity'   'Riskscore'  



p <- ggstatsplot::ggbetweenstats(data = ESTIMATE_score, x = Group,
                                 y = TumorPurity, type = "np", mean.ci = TRUE, pairwise.comparisons = TRUE,
                                 pairwise.display = "all", ylab = "TumorPurity", messages = F, ggstatsplot.layer = F,
                                 ggtheme = ggplot2::theme_classic(), package = "ggsci", palette = "uniform_startrek",
                                 results.subtitle = F, subtitle = "p=1.45e-03")
pdf(file = paste0("TumorPurity", "risk_group.pdf"), width = 3.5, height = 5.5,
    onefile = T)
print(p)
dev.off()

