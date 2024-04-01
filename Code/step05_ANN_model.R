rm(list = ls())
gc()
library(DESeq2)
library(edgeR)
library(ImageGP)
library(ggplot2)
library(ggpubr)
library(egg)
library(ggrepel)
library(survival)
library(glmnet)
library(neuralnet)
library(tibble)
library(tidyverse)
library(sampling)
library(pROC)
library(NeuralNetTools)
library(nnet)
library(caret)
library(RColorBrewer)
library(mlbench)
library(survival)
library(survminer)

############################  DEG  ########################################
load("train_tcga.rds")
gset <- exp_tcga[rowMeans(exp_tcga) > 0, ]
group <- clin_tcga$LLPS_subtype
design <- model.matrix(~0 + factor(group))
colnames(design) <- levels(factor(group))
rownames(design) <- rownames(clin_tcga)
cont.wt <- makeContrasts("LS1-LS2", "LS1-LS3", "LS2-LS3", levels = design)

fit <- lmFit(gset, design)
fit2 <- contrasts.fit(fit, cont.wt)
fit2 <- eBayes(fit2)
tT <- topTableF(fit2, adjust = "BH", sort.by = "F", n = Inf)

s1_dif <- rownames(tT[abs(tT$LS1.LS2) > 2 & tT$adj.P.Val < 0.01, ])
s2_dif <- rownames(tT[abs(tT$LS1.LS3) > 2 & tT$adj.P.Val < 0.01, ])
s3_dif <- rownames(tT[abs(tT$LS2.LS3) > 2 & tT$adj.P.Val < 0.01, ])

deg <- unique(c(s1_dif, s2_dif, s3_dif))


############################  univariate Cox & bootstrap  ########################################

data <- cbind(clin_tcga[, 1:2], as.data.frame(t(exp_tcga[deg, ])))
coxR <- data.frame()
coxf <- function(x)
{
  fmla1 <- as.formula(Surv(OS.time, OS) ~ data[, x])
  mycox <- coxph(fmla1, data = data)
}

for (a in colnames(data[, 3:ncol(data)]))
{
  mycox <- coxf(a)
  coxResult <- summary(mycox)
  coxR <- rbind(coxR, cbind(id = a, HR = coxResult$coefficients[, "exp(coef)"],
                            P = coxResult$coefficients[, "Pr(>|z|)"]))
}

sigGene <- coxR[as.numeric(as.character(coxR$P)) < 0.01, ]
data1 <- cbind(data[, 1:2], data[, as.character(sigGene$id)])

### bootstrapping method 
patients <- rownames(data1)
outTab <- data.frame()
for (gene in colnames(data1[, 3:ncol(data1)]))
{
  Mboot <- replicate(1000, expr = {
    indices <- sample(patients, size = nrow(data1) * 0.7, replace = TRUE)
    data <- data1[indices, ]
    fmla1 <- as.formula(Surv(data[, "OS.time"], data[, "OS"]) ~ data[, gene])
    mycox <- coxph(fmla1, data = data)
    coxResult <- summary(mycox)
    P <- coxResult$coefficients[, "Pr(>|z|)"]
  })
  times <- length(Mboot[which(Mboot < 0.05)])
  outTab <- rbind(outTab, cbind(gene = gene, times = times))
}

prog_gene <- outTab[as.numeric(as.character(outTab$times)) > 900, ]$gene
save(prog_gene, file = "bootstrap_prog_gene.rds")


############################  lasso  ########################################

set.seed(20230402)
load("bootstrap_prog_gene.rds")
est_dd2 <- cbind(data[, 1:2], data[, prog_gene])
fit <- glmnet(as.matrix(est_dd2[, 3:ncol(est_dd2)]), as.matrix(Surv(est_dd2$OS.time,
                                                                    est_dd2$OS)), family = "cox")

pdf("lasso_lambda.pdf", width = 4.5, height = 5, onefile = F)
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

cv.fit <- cv.glmnet(as.matrix(est_dd2[, 3:ncol(est_dd2)]), as.matrix(Surv(est_dd2$OS.time,
                                                                          est_dd2$OS)), family = "cox", nfolds = 10)

pdf("lasso_cvfit.pdf", width = 4.5, height = 5, onefile = F)
plot(cv.fit)
abline(v = log(c(cv.fit$lambda.min, cv.fit$lambda.1se)), lty = "dashed")
dev.off()

coef <- coef(fit, s = cv.fit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene <- row.names(coef)[index]  #8
save(lassoGene, file = "lassoGene.rds")


############################  ANN  ########################################

rm(list = ls())
gc()
load("train_tcga.rds")
load("lassoGene.rds")

exp_lasso <- as.data.frame(t(exp_tcga[lassoGene, ]))
temp <- cbind(exp_lasso, Survival = clin_tcga$OS)
data <- temp %>%
  mutate_at("Survival", as.factor)

##data loading and cleaning----
#random sampling
set.seed(1221)
select <- sample(1:nrow(data), nrow(data) * 0.7)
train <- data[select, ]
test <- data[-select, ]

# standardized data
train[, 1:8] <- scale(train[, 1:8])
test[, 1:8] <- scale(test[, 1:8])
data[, 1:8] <- scale(data[, 1:8])
##Implementation of BP neural network using packet-----

mynnet <- nnet(Survival ~ ., linout = F, size = 2, decay = 0.04, maxit = 200,
               data = train)
#save(mynnet,file = 'mynnet.rds')

#linout judge whether the output is linear
# size is Number of nodes 
#decay is attenuation rate
#decay and learning rate are similar
#Maximum number of iterations
pdf("neural interpretation diagram.pdf", width = 11, height = 7.5, onefile = F)
plotnet(mynnet, circle_col = "#4FA1D0", bord_col = "white", pos_col = "#4F9595",
        neg_col = colorRampPalette(c("#4F9595", "white"))(10)[5], circle_cex = 6)
dev.off()

pdf("variable importance with garson.pdf", width = 6, height = 5, onefile = F)
garson(mynnet)
dev.off()

olden(mynnet)
lekprofile(mynnet, group_vals = 5)
lekprofile(mynnet, group_vals = 5, group_show = TRUE)

##result evaluation----

##model predict
out <- predict(mynnet, test)
out[out < 0.5] <- 0
out[out >= 0.5] <- 1
##calculation accuracy
rate <- sum(out == test$Survival)/length(test$Survival)
rate

####Predict on the training set and test set respectively
##The construction of ROC curve function ROC () is convenient.
##Note that this function will add columns to the data frame called
ROC <- function(model, train, test, objcolname, ifplot = TRUE)
{
  library(ROCR, quietly = T)
  train$p <- predict(model, train)
  test$p <- predict(model, test)
  
  predTr <- prediction(train$p, train[, objcolname])
  perfTr <- performance(predTr, "tpr", "fpr")
  
  predTe <- prediction(test$p, test[, objcolname])
  perfTe <- performance(predTe, "tpr", "fpr")
  
  tr_auc <- round(as.numeric(performance(predTr, "auc")@y.values), 3)
  te_auc <- round(as.numeric(performance(predTe, "auc")@y.values), 3)
  
  if (ifplot == T)
  {
    plot(perfTr, col = "#4F9595", main = "ROC", lwd = 2)
    plot(perfTe, col = "#E5928E", lty = 1, add = TRUE, lwd = 2)
    abline(0, 1, lty = 1, col = "#AFB6B6")
    
    tr_str <- paste("Train-AUC:", tr_auc, sep = "")
    legend(0.6, 0.35, tr_str, border = NA, box.col = NA, fill = "#4F9595", col = NA)
    te_str <- paste("Test-AUC:", te_auc, sep = "")
    legend(0.6, 0.15, te_str, border = NA, box.col = NA, fill = "#E5928E", col = NA)
  }
  auc <- data.frame(tr_auc, te_auc)
  return(auc)
}
ROC(model = mynnet, train = train, test = test, objcolname = "Survival", ifplot = T)
# 6 6.5



##Adjusting parameters----
##The predictive variables of the input data must be binary. And all variables only include model input and output variables.
##When adjusting parameters, if there are too many variables, the size should not be too large.
##Build the parameter adjustment function network().
network <- function(formula, data, size, adjust, decay = 0, maxit = 200, scale = TRUE,
                    samplerate = 0.7, seed = 1, linout = FALSE, ifplot = TRUE)
{
  
  library(nnet)
  ##The specification output variable is 0,1
  yvar <- colnames(data) == (all.vars(formula)[1])
  levels(data[, yvar]) <- c(0, 1)
  ##Establish training set and test set by sampling
  set.seed(seed)
  select <- sample(1:nrow(data), nrow(data) * samplerate)
  train <- data[select, ]
  test <- data[-select, ]
  ##Standardize according to given judgment
  if (scale == T)
  {
    xvar <- colnames(data) != (all.vars(formula)[1])
    train[, xvar] <- scale(train[, xvar])
    test[, xvar] <- scale(test[, xvar])
  }
  ##Recycle NNET training parameters
  obj <- eval(parse(text = adjust))
  auc <- data.frame()
  for (i in obj)
  {
    if (adjust == "size")
    {
      mynnet <- nnet(formula, size = i, linout = linout, decay = decay, maxit = maxit,
                     trace = FALSE, data = train)
    } else if (adjust == "decay")
    {
      mynnet <- nnet(formula, size = size, linout = linout, decay = i, maxit = maxit,
                     trace = FALSE, data = train)
    }
    ##Call the previous roc() to get the AUC value of the corresponding parameterֵ
    objcolname <- all.vars(formula)[1]
    auc0 <- ROC(model = mynnet, train = train, test = test, objcolname = objcolname,
                ifplot = F)
    ##Output data frames corresponding to different values of specified parameters
    out <- data.frame(i, auc0)
    auc <- rbind(auc, out)
  }
  
  names(auc) <- c(adjust, "Train_auc", "Test_auc")
  # if(ifplot==T){
  #   library(plotrix)
  #   # twoord.plot(auc1[,1] , auc1$Train_auc , auc1[,1] , auc1$Test_auc , lcol=4 , rcol=2 , xlab=adjust , 
  #   #             ylab='Train_auc' , rylab='Test_auc' , type=c('l','b'),lab=c(15,5,10))
  # }
  return(auc)
}
auc <- network(Survival ~ ., data = data, size = 1:16, adjust = "size", decay = 1e-04,
               maxit = 200, scale = T)
auc$size[which(auc$Test_auc == max(auc$Test_auc))]  #2

auc <- network(Survival ~ ., data = data, size = 14, adjust = "decay", decay = c(0,
                                                                                 seq(1e-04, 0.01, 3e-04)), maxit = 200)
auc$decay[which(auc$Test_auc == max(auc$Test_auc))]  #0.004

plot(auc$decay, auc$Train_auc, ylim = c(0.5, 1))
lines(auc$decay, auc$Train_auc)
points(auc$decay, auc$Test_auc)
lines(auc$decay, auc$Test_auc)




#risk score----

KM <- function(data, clin, name)
{
  
  risk <- as.data.frame(predict(mynnet, data))
  risk$group <- ifelse(risk$V1 > median(risk$V1), "High_risk", "Low_risk")
  clin_data <- clin[rownames(data), ]
  if (identical(rownames(risk), rownames(clin_data)) == T)
  {
    clin_data$group <- risk$group
  }
  sfit <- survfit(Surv(OS.time, OS) ~ group, data = clin_data)
  pp <- ggsurvplot(sfit, pval = T, palette = c("#E5928E", "#709ECD"), risk.table = T)
  pdf(paste0(name, "_KM_plot_ANN.pdf"), width = 6.5, height = 6.5, onefile = F)
  print(pp)
  dev.off()
  
}
KM(train, clin_tcga, "train")
KM(test, clin_tcga, "test")



# load('../A01_Data_Processing/PAAD/all_15_set.rds')
# #[1] 'AU_Array'    'TCGA'        'AU_Seq'      'CA_Seq'      'E_MTAB_6134' 'GSE62452'    'GSE28735'    'GSE78229'   
# #[9] 'GSE79668'    'GSE85916'    'CELL'        'GSE21501'    'GSE57495'    'GSE71729'    'RJ'  


load("/Users/djchen/Desktop/R data/A01_Data_Processing/PAAD/ICGC-AU-Array/CDJ_PDAC_ICGC_AU_Array_mRNA_T_data.Rdata")
load("/Users/djchen/Desktop/R data/A01_Data_Processing/PAAD/ICGC-AU-Array/CDJ_PDAC_ICGC_AU_Array_Clinical_T_data.Rdata")
traning_set <- merge(clin_array_icgc_au[, c("os_time", "os_status")],
                     as.data.frame(t(exp_array_icgc_au)), by = "row.names", all = T)
colnames(traning_set)[1:3] <- c("ID", "OS.time", "OS")


# KM_2 <- function(data,name){
#   data <- data %>% column_to_rownames('ID')
#   risk <- as.data.frame(predict(mynnet, scale(data[,lassoGene]) ))
#   risk$group <- ifelse(risk$V1 > median(risk$V1),'High_risk','Low_risk')
#   if (identical(rownames(risk),rownames(data))== T)
# {
  #     data$group <- risk$group
  #   }
  #   sfit <- survfit(Surv(OS.time, OS) ~ group,
  #                   data = data)
  #   pp <- ggsurvplot(sfit,
  #                    data = data,
  #                    pval = T,
  #                    palette = c('#E5928E','#709ECD'),
  #                    risk.table = T)
  #   pdf(paste0(name,'_KM_plot_ANN.pdf'),width = 6.5,height = 6.5,onefile = F)
  #   print(pp)
  #   dev.off()
  # }
  # KM_2(traning_set,'ICGC_AU')
  
  
  i <- 1
  KM_2(all_15_set[[i]], names(all_15_set)[i])
  
  
  # TCGA riskscore ----
  
  risk <- as.data.frame(predict(mynnet, data))
  risk$group <- ifelse(risk$V1 > median(risk$V1), "High_risk", "Low_risk")
  risk <- risk[rownames(clin_tcga), ]
  clin_tcga$riskscore <- risk$V1
  clin_tcga$group <- risk$group
  save(exp_tcga, clin_tcga, file = "train_tcga.rds")
  
  
  # forest plot ----
  
  exp <- all_15_set[[2]]
  rt <- cbind(exp[, 2:3], exp[, lassoGene])
  colnames(rt)[1:2] <- c("futime", "fustat")
  #单因素独立预后分析
  outTab <- data.frame()
  for (i in colnames(rt[, 3:ncol(rt)]))
  {
    cox <- coxph(Surv(futime, fustat) ~ rt[, i], data = rt)
    coxSummary <- summary(cox)
    coxP <- coxSummary$coefficients[, "Pr(>|z|)"]
    outTab <- rbind(outTab, cbind(id = i, HR = coxSummary$conf.int[, "exp(coef)"],
                                  HR.95L = coxSummary$conf.int[, "lower .95"], HR.95H = coxSummary$conf.int[,
                                                                                                            "upper .95"], pvalue = coxSummary$coefficients[, "Pr(>|z|)"]))
  }
  write.table(outTab, file = "uniCox_results.txt", sep = "\t", row.names = F, quote = F)
  bioForest <- function(coxFile = null, forestCol = null, forestFile = null)
  {
    rt <- read.table(coxFile, header = T, sep = "\t", row.names = 1, check.names = F)
    gene <- rownames(rt)
    hr <- sprintf("%.3f", rt$HR)
    hrLow <- sprintf("%.3f", rt$HR.95L)
    hrHigh <- sprintf("%.3f", rt$HR.95H)
    Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
    pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
    pdf(file = forestFile, width = 6.5, height = 5)
    n <- nrow(rt)
    nRow <- n + 1
    ylim <- c(1, nRow)
    layout(matrix(c(1, 2), nc = 2), width = c(3, 2.2))
    xlim <- c(0, 3)
    par(mar = c(4, 2.5, 2, 1))
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, xlab = "", ylab = "")
    text.cex <- 0.8
    text(0, n:1, gene, adj = 0.1, cex = text.cex)
    text(0.5, n + 1, "IRS", cex = text.cex, font = 2, adj = 1)
    text(1.5 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
    text(1.5 - 0.5 * 0.2, n + 1, "pvalue", cex = text.cex, font = 2, adj = 1)
    text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
    text(3, n + 1, "Hazard ratio", cex = text.cex, font = 2, adj = 1, )
    par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
    xlim <- c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
    plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, ylab = "", xaxs = "i",
         xlab = "Hazard ratio")
    arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3,
           length = 0.05, col = "#4DBBD6", lwd = 2.5)
    abline(v = 1, col = "black", lty = 2, lwd = 2)
    boxcolor <- ifelse(as.numeric(hr) > 1, forestCol, forestCol)
    points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.3)
    axis(1)
    dev.off()
  }
  bioForest(coxFile = "uniCox_results.txt", forestCol = "#E64B35",
            forestFile = "IRS_Forest.pdf")
  
  
  
  
  # riskscore mutation ########################################
  rm(list = ls())
  gc()
  library(tidyverse)
  library(gmodels)
  library(ggpubr)
  library(ggthemes)
  library(grid)
  library(gridExtra)
  library(ggplot2)
  library(circlize)
  library(ComplexHeatmap)
  library(survival)
  library(survminer)
  library(plyr)
  library(ggsci)
  library(ggsignif)
  library(ComplexHeatmap)
  load("train_tcga.rds")
  load("../A01_Data_Processing/PAAD/TCGA/tcga_tmb_snp_cnv.rds")
  
  
  pd_select <- clin_tcga
  mut <- TCGA_SNP
  merge_id <- intersect(rownames(pd_select), colnames(mut))
  pd_select <- pd_select[merge_id, ]
  pd_select <- pd_select[order(pd_select$group), ]
  mut <- mut[, rownames(pd_select)]
  table(pd_select$group)
  
  if (T)
  {
    pd_select_1 <- pd_select[1:82, ]
    mut_1 <- mut[, 1:82]
    pd_select <- pd_select_1
    mut <- mut_1
  }
  
  if (F)
  {
    pd_select_2 <- pd_select[83:nrow(pd_select), ]
    mut_2 <- mut[, 83:nrow(pd_select)]
    pd_select <- pd_select_2
    mut <- mut_2
  }
  
  #网格颜色
  col <- c(Missense_Mutation = "#F1A6A4", Splice_Site = "#BF9E9A",
           Nonsense_Mutation = "#D385BB", Frame_Shift_Del = "#2BA9B5", Frame_Shift_Ins = "#EEBD8E",
           Multi_Hit = "#A1D599", In_Frame_Del = "#C7B5D1", In_Frame_Ins = "#B4CADF")
  alter_fun <- list(background = alter_graphic("rect", fill = "#e9e9e9"),
                    Missense_Mutation = alter_graphic("rect", fill = col["Missense_Mutation"]),
                    Splice_Site = alter_graphic("rect", fill = col["Splice_Site"]),
                    Nonsense_Mutation = alter_graphic("rect", fill = col["Nonsense_Mutation"]),
                    Frame_Shift_Del = alter_graphic("rect", fill = col["Frame_Shift_Del"]),
                    Frame_Shift_Ins = alter_graphic("rect", fill = col["Frame_Shift_Ins"]),
                    Multi_Hit = alter_graphic("rect", fill = col["Multi_Hit"]),
                    In_Frame_Del = alter_graphic("rect", fill = col["In_Frame_Del"]),
                    In_Frame_Ins = alter_graphic("rect", fill = col["In_Frame_Ins"]))
  
  top_anno <- HeatmapAnnotation(Group = pd_select$group,
                                col = list(Group = c(High_risk = "#E5928E", Low_risk = "#709ECD")),
                                annotation_height = unit(c(3.8, 0.8), "cm"), show_legend = F,
                                show_annotation_name = T, gap = unit(1, "mm"))
  
  lgd_list <- list(Legend(labels = c("High_risk", "Low_risk"),
                          title = "", legend_gp = gpar(fill = c("#E5928E", "#709ECD")),
                          labels_gp = gpar(fontsize = 13), grid_height = unit(0.5,
                                                                              "cm"), grid_width = unit(0.5, "cm"), title_gp = gpar(fontsize = 13,
                                                                                                                                   fontface = "bold")), Legend(labels = c("Missense_Mutation",
                                                                                                                                                                          "Multi_Hit", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
                                                                                                                                                                          "Splice_Site", "In_Frame_Del", "In_Frame_Ins"), title = "Alterations",
                                                                                                                                                               legend_gp = gpar(fill = c("#F1A6A4", "#BF9E9A", "#D385BB",
                                                                                                                                                                                         "#2BA9B5", "#EEBD8E", "#A1D599", "#C7B5D1", "#B4CADF")),
                                                                                                                                                               labels_gp = gpar(fontsize = 13), grid_height = unit(0.5,
                                                                                                                                                                                                                   "cm"), grid_width = unit(0.5, "cm"), title_gp = gpar(fontsize = 13,
                                                                                                                                                                                                                                                                        fontface = "bold")))
  
  #right annotation
  #data_processing
  dat_right_anno <- mut
  dat_right_anno[dat_right_anno == ""] <- 0
  dat_right_anno[!dat_right_anno == 0] <- 1
  dimnames <- list(rownames(dat_right_anno), colnames(dat_right_anno))
  dat_right_anno <- matrix(as.numeric(as.matrix(dat_right_anno)),
                           nrow = nrow(dat_right_anno), dimnames = dimnames)
  
  # plot
  if (T)
  {
    pdf("SNP_oncoplot_riskscore.pdf", width = 12, height = 8)
    ht_list_1 <- oncoPrint(mut, alter_fun = alter_fun, alter_fun_is_vectorized = FALSE,
                           width = unit(5 * 3, "cm"), height = unit(9 * 1, "cm"), pct_gp = gpar(fontsize = 13),
                           col = col, row_title = NULL, column_title = NULL, row_gap = unit(3,
                                                                                            "mm"), column_gap = unit(2, "mm"), row_order = 1:nrow(mut), column_order = NULL,
                           row_split = NULL, top_annotation = top_anno, show_heatmap_legend = F)
    draw(ht_list_1)
    ht_list_2 <- oncoPrint(mut, alter_fun = alter_fun, alter_fun_is_vectorized = FALSE,
                           width = unit((61 * 5 * 3)/82, "cm"), height = unit(9 * 1, "cm"),
                           pct_gp = gpar(fontsize = 13), col = col, row_title = NULL, column_title = NULL,
                           row_gap = unit(3, "mm"), column_gap = unit(2, "mm"), row_order = 1:nrow(mut),
                           column_order = NULL, row_split = NULL, top_annotation = top_anno,
                           show_heatmap_legend = F)
    ht <- ht_list_1 + ht_list_2
    draw(ht, heatmap_legend_list = lgd_list)
    dev.off()
  }
  
  