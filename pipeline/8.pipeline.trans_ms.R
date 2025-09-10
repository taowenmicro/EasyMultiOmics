# 转录与代谢---
rm(list=ls())
# library(EasyMultiOmics)
ps02 = EasyMultiOmics::ps.ms
ps01 = EasyMultiOmics::ps.trans
ps.tem = ps.ms %>% filter_OTU_ps(500)

# 创建转录组-代谢组分析目录
trans_metab_path =  "./result/transcriptome_metabolite/"
dir.create(trans_metab_path, recursive = TRUE)

# 创建基因筛选与排序分析目录
gene_selection_path <- file.path(trans_metab_path, "gene_selection")
dir.create(gene_selection_path, recursive = TRUE)

# 创建基因筛选与排序分析总表
gene_selection_wb <- createWorkbook()

#0 基因筛选-------
res = loadingPCA.trans(ps =ps.trans)
# p = res[[1]]
# p
dat = res[[2]]
dat$id = row.names(dat)
id = dat %>% arrange(desc(PCone)) %>% head(15) %>% .$id
id

ftab = ps.trans %>% scale_micro() %>%
  subset_taxa(
    row.names(tax_table(ps.trans)) %in% c(id)
  ) %>% vegan_otu() %>%
  as.data.frame()

head(ftab)

# 保存基因筛选结果
addWorksheet(gene_selection_wb, "gene_selection_PCA")
writeData(gene_selection_wb, "gene_selection_PCA", dat, rowNames = TRUE)
addWorksheet(gene_selection_wb, "selected_genes")
writeData(gene_selection_wb, "selected_genes", ftab, rowNames = TRUE)

#1 RDA_CCA:对代谢物影响的基因共排序-----
result = RDA_CCA(ps = ps.tem,
                 env = ftab,
                 # path = RDApath,
                 chose.env = FALSE
)

p1 = result[[1]]
p1
# 提取作图数据
dataplot = result[[2]]
dataplot
# 提取带有标记的图片
p2 = result[[3]]
p2
aov = result[[4]]
aov

# 保存RDA/CCA结果
ggsave(file.path(gene_selection_path, "RDA_CCA_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(gene_selection_path, "RDA_CCA_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(gene_selection_path, "RDA_CCA_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(gene_selection_path, "RDA_CCA_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(gene_selection_wb, "RDA_CCA_plotdata")
writeData(gene_selection_wb, "RDA_CCA_plotdata", dataplot, rowNames = TRUE)
addWorksheet(gene_selection_wb, "RDA_CCA_anova")
writeData(gene_selection_wb, "RDA_CCA_anova", aov, rowNames = TRUE)

#2 RDA_CCA_explain_percent: 对代谢物的解释比例-----
result = RDA_CCA_explain_percent(ps = ps.tem,
                                 env.dat = ftab)
sam = result[[1]]
sam
all = result[[2]]
all

# 保存解释比例结果
addWorksheet(gene_selection_wb, "RDA_CCA_explain_sample")
writeData(gene_selection_wb, "RDA_CCA_explain_sample", sam, rowNames = TRUE)
addWorksheet(gene_selection_wb, "RDA_CCA_explain_all")
writeData(gene_selection_wb, "RDA_CCA_explain_all", all, rowNames = TRUE)

#3 rdacca.hp.micro:    层次分割对代谢物影响的基因#---------
library(rdacca.hp)
library(vegan)
library(ggClusterNet)
library(tidyverse)
res = rdacca.hp.micro(OTU =ps.tem %>% filter_OTU_ps(500) %>%vegan_otu() %>% as.data.frame(),
                      env = ftab[,1:5],
                      cca = FALSE
)
p3_rda = res[[1]]
p3_rda
dat1 = res[[2]]
dat1
p3_db_rda  = res[[3]]
p3_db_rda
dat2 = res[[4]]
dat2

# 保存层次分割结果
ggsave(file.path(gene_selection_path, "rdacca_hp_plot1.png"), plot = p3_rda, width = 10, height = 8, dpi = 300)
ggsave(file.path(gene_selection_path, "rdacca_hp_plot1.pdf"), plot = p3_rda, width = 10, height = 8)
ggsave(file.path(gene_selection_path, "rdacca_hp_plot2.png"), plot = p3_db_rda, width = 10, height = 8, dpi = 300)
ggsave(file.path(gene_selection_path, "rdacca_hp_plot2.pdf"), plot = p3_db_rda, width = 10, height = 8)
addWorksheet(gene_selection_wb, "rdacca_hp_data1")
writeData(gene_selection_wb, "rdacca_hp_data1", dat1, rowNames = TRUE)
addWorksheet(gene_selection_wb, "rdacca_hp_data2")
writeData(gene_selection_wb, "rdacca_hp_data2", dat2, rowNames = TRUE)

#4 rdacca.hp.micro.p: 层次分割对代谢物影响的基因#---------
res <-rdacca.hp.micro.p(
  OTU = ps.tem%>% filter_OTU_ps(200) %>%vegan_otu() %>% as.data.frame(),
  env = ftab[,1:5],
  cca = FALSE,
  dbRDA = FALSE
)
p3_rda = res[[1]]
p3_rda
dat1 = res[[2]]
dat1

# 保存显著层次分割结果
ggsave(file.path(gene_selection_path, "rdacca_hp_sig.png"), plot = p3_rda, width = 10, height = 8, dpi = 300)
ggsave(file.path(gene_selection_path, "rdacca_hp_sig.pdf"), plot = p3_rda, width = 10, height = 8)
addWorksheet(gene_selection_wb, "rdacca_hp_sig_data")
writeData(gene_selection_wb, "rdacca_hp_sig_data", dat1, rowNames = TRUE)

# 保存基因筛选与排序分析总表
saveWorkbook(gene_selection_wb, file.path(gene_selection_path, "transcriptome_metabolome_gene_selection.xlsx"), overwrite = TRUE)

# 相关性分析-------
# 创建相关性分析目录
trans_metab_cor_path <- file.path(trans_metab_path, "correlation")
dir.create(trans_metab_cor_path, recursive = TRUE)

# 创建相关性分析总表
trans_metab_cor_wb <- createWorkbook()

#5 MetalTast: mantal 代谢物与基因的相关性计算----
otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps.tem)))
head(otu1)

tabOTU1 = list(bac = otu1 )

rep = MetalTast (env.dat = ftab, tabOTU = tabOTU1,distance = "bray",method = "metal")
repR = rep[c(-seq(from=1,to=dim(rep)[2],by=2)[-1])]
repP = rep[seq(from=1,to=dim(rep)[2],by=2)]
head( repR)
head( repP)
mat = cbind(repR,repP)
head(mat)

# 保存MetalTast结果
addWorksheet(trans_metab_cor_wb, "metalTast_results")
writeData(trans_metab_cor_wb, "metalTast_results", mat, rowNames = TRUE)

#6 cor_env_ggcorplot_top:高丰度代谢物和基因的关系探索----
result = cor_env_ggcorplot_top(ps =  ps.tem,
                               env1 = ftab,
                               label =  TRUE,
                               col_cluster = TRUE,
                               row_cluster = TRUE,
                               method = "spearman",
                               r.threshold= 0,
                               p.threshold= 1  )

p1 <- result[[1]]
p1
p2 <- result[[2]]
p2
topotu =  result[[3]]
topotu
cordata = result[[4]]
cordata

# 保存高丰度相关性结果
ggsave(file.path(trans_metab_cor_path, "cor_top_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_cor_path, "cor_top_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(trans_metab_cor_path, "cor_top_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_cor_path, "cor_top_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(trans_metab_cor_wb, "cor_top_otu")
writeData(trans_metab_cor_wb, "cor_top_otu", topotu, rowNames = TRUE)
addWorksheet(trans_metab_cor_wb, "cor_top_data")
writeData(trans_metab_cor_wb, "cor_top_data", cordata, rowNames = TRUE)

#7 cor_env_ggcorplot_top:丰度差异大的代谢物和基因关系探索----
result = cor_env_ggcorplot_rcv(ps =  ps.tem,
                               env1 = ftab,
                               label =  TRUE,
                               col_cluster = TRUE,
                               row_cluster = TRUE,
                               method = "spearman",
                               r.threshold= 0,
                               p.threshold= 1  )

p1 <- result[[1]]
p1
p2 <- result[[2]]
p2
rcvotu =  result[[3]]
rcvotu
cordata = result[[4]]
cordata

# 保存差异相关性结果
ggsave(file.path(trans_metab_cor_path, "cor_rcv_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_cor_path, "cor_rcv_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(trans_metab_cor_path, "cor_rcv_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_cor_path, "cor_rcv_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(trans_metab_cor_wb, "cor_rcv_otu")
writeData(trans_metab_cor_wb, "cor_rcv_otu", rcvotu, rowNames = TRUE)
addWorksheet(trans_metab_cor_wb, "cor_rcv_data")
writeData(trans_metab_cor_wb, "cor_rcv_data", cordata, rowNames = TRUE)

#8 MatCorPlot:代谢物和基因Science组合图表------
otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps.tem)))
head(otu1)
tabOTU1 = list(bac = otu1 )
detach("package:mia", unload = TRUE)
# tabOTU1 = list(b = otu)
# tabOTU1 = list(f = otu2）

p <-  MatCorPlot(env.dat = ftab,
                 tabOTU = tabOTU1,
                 diag = FALSE,
                 range = 0.2,
                 numpoint = 21,
                 sig = FALSE,
                 siglabel = FALSE,
                 shownum = FALSE,
                 curvature = 0,
                 numsymbol = NULL,
                 lacx = "right",
                 lacy = "bottom",
                 p.thur = 1,
                 onlysig = TRUE)

p

p0 <- MatCorPlot2(env.dat = ftab,
                  tabOTU = tabOTU1,
                  corva = -0.05,
                  diag = TRUE,
                  range = 0.05,
                  numpoint = 21,
                  numpoint2 = 22,
                  sig = FALSE,
                  siglabel = FALSE,
                  shownum = FALSE,
                  curvature = 0,
                  numsymbol = NULL,
                  lacx = "right",
                  lacy = "bottom",
                  p.thur = 0.5,
                  onlysig = FALSE
)
p0

# 保存MatCorPlot结果
ggsave(file.path(trans_metab_cor_path, "MatCorPlot.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(trans_metab_cor_path, "MatCorPlot.pdf"), plot = p, width = 12, height = 10)
ggsave(file.path(trans_metab_cor_path, "MatCorPlot2.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(trans_metab_cor_path, "MatCorPlot2.pdf"), plot = p0, width = 12, height = 10)

# 保存相关性分析总表
saveWorkbook(trans_metab_cor_wb, file.path(trans_metab_cor_path, "transcriptome_metabolome_correlation.xlsx"), overwrite = TRUE)

# 多组学整合分析-------
# 创建多组学整合分析目录
trans_metab_omics_path <- file.path(trans_metab_path, "multi_omics")
dir.create(trans_metab_omics_path, recursive = TRUE)

# 创建多组学整合分析总表
trans_metab_omics_wb <- createWorkbook()

#9 heatmap.line.omics---------
res=heatmap.line.omics(ps01=ps.ms,
                       ps02=ps.trans,
                       lab.1 = "ms",
                       lab.2 = "trans")
grid.draw(res[[1]])
res[[2]]
res[[3]]
res[[4]]

# 保存热图线性组学结果
ggsave(file.path(trans_metab_omics_path, "heatmap_line_omics.png"), plot = res[[1]], width = 12, height = 10, dpi = 300)
ggsave(file.path(trans_metab_omics_path, "heatmap_line_omics.pdf"), plot = res[[1]], width = 12, height = 10)
addWorksheet(trans_metab_omics_wb, "heatmap_line_data1")
writeData(trans_metab_omics_wb, "heatmap_line_data1", res[[2]], rowNames = TRUE)
addWorksheet(trans_metab_omics_wb, "heatmap_line_data2")
writeData(trans_metab_omics_wb, "heatmap_line_data2", res[[3]], rowNames = TRUE)
addWorksheet(trans_metab_omics_wb, "heatmap_line_data3")
writeData(trans_metab_omics_wb, "heatmap_line_data3", res[[4]], rowNames = TRUE)

#10 volcano.line.omics: 火山互连 ----
library(ggnewscale)
ps1 = ps.trans %>% filter_OTU_ps(500)
ps2 = ps.ms %>% filter_OTU_ps(500)
res = volcano.line.omics(ps1 = ps1,ps2 = ps2,
                         lab.1 = "tran",
                         lab.2 = "ms",
                         group = "Group",
                         r.threshold = 0.8,
                         p.threshold= 0.05,
                         col.point = NULL,
                         col.line = NULL,
                         method.cor = "pearson",
                         n = 50
)

res[[1]]$WT_OE
res[[1]]$WT_KO
res[[1]]$OE_KO

# 保存volcano.line.omics结果
if (!is.null(res[[1]]$WT_OE)) {
  ggsave(file.path(trans_metab_omics_path, "volcano_line_WT_OE.png"), plot = res[[1]]$WT_OE, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_metab_omics_path, "volcano_line_WT_OE.pdf"), plot = res[[1]]$WT_OE, width = 10, height = 8)
}
if (!is.null(res[[1]]$WT_KO)) {
  ggsave(file.path(trans_metab_omics_path, "volcano_line_WT_KO.png"), plot = res[[1]]$WT_KO, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_metab_omics_path, "volcano_line_WT_KO.pdf"), plot = res[[1]]$WT_KO, width = 10, height = 8)
}
if (!is.null(res[[1]]$OE_KO)) {
  ggsave(file.path(trans_metab_omics_path, "volcano_line_OE_KO.png"), plot = res[[1]]$OE_KO, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_metab_omics_path, "volcano_line_OE_KO.pdf"), plot = res[[1]]$OE_KO, width = 10, height = 8)
}

res= volcano.line2.omics(ps01=ps.ms,
                         ps02=ps.trans,
                         group1= "ms",
                         group2= "trans",
                         r.threshold = 0.6,
                         p.threshold = 0.1,
                         method = "spearman",
                         top = 500)

p=res[[1]]$KO_OE
p
dat=res[[2]]$KO_OE

# 保存volcano.line2.omics结果
if (!is.null(p)) {
  ggsave(file.path(trans_metab_omics_path, "volcano_line2_KO_OE.png"), plot = p, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_metab_omics_path, "volcano_line2_KO_OE.pdf"), plot = p, width = 10, height = 8)
}
if (!is.null(dat)) {
  addWorksheet(trans_metab_omics_wb, "volcano_line2_KO_OE_data")
  writeData(trans_metab_omics_wb, "volcano_line2_KO_OE_data", dat, rowNames = TRUE)
}

#11 quare.line.omics: 九象限火山图 ----
results <- quare.line.omics(ps.16s =ps.ms,
                            ps.ms =  ps.trans,
                            group1 = "ms",
                            group2 = "trans",
                            r.threshold = 0.7,
                            p.threshold = 0.05,
                            method = "spearman",
                            top = 500)
results$KO_WT
results$OE_WT$data

# 保存九象限火山图结果
if (!is.null(results$KO_WT)) {
  ggsave(file.path(trans_metab_omics_path, "quare_line_KO_WT.png"), plot = results$KO_WT, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_metab_omics_path, "quare_line_KO_WT.pdf"), plot = results$KO_WT, width = 10, height = 8)
}
if (!is.null(results$OE_WT$data)) {
  addWorksheet(trans_metab_omics_wb, "quare_line_OE_WT_data")
  writeData(trans_metab_omics_wb, "quare_line_OE_WT_data", results$OE_WT$data, rowNames = TRUE)
}

#12 随机森林寻找对代谢物群落影响较大的基因--------
result <- ms_micro.rf.omics(ps  = ps.tem,env = ftab )
p <- result[[2]]
p
data = result[[1]]
head(data)

# 保存随机森林结果
ggsave(file.path(trans_metab_omics_path, "rf_importance.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_omics_path, "rf_importance.pdf"), plot = p, width = 10, height = 8)
addWorksheet(trans_metab_omics_wb, "rf_importance_data")
writeData(trans_metab_omics_wb, "rf_importance_data", data, rowNames = TRUE)

# 保存多组学整合分析总表
saveWorkbook(trans_metab_omics_wb, file.path(trans_metab_omics_path, "transcriptome_metabolome_multi_omics.xlsx"), overwrite = TRUE)

# 机器学习----
# 创建机器学习分析目录
trans_metab_ml_path <- file.path(trans_metab_path, "machine_learning")
dir.create(trans_metab_ml_path, recursive = TRUE)

# 创建机器学习分析总表
trans_metab_ml_wb <- createWorkbook()

library(randomForest)
library(caret)
library(ROCR)
library(e1071)
library(rpart)
library(mia)
library(TreeSummarizedExperiment)
library("class")
library(kknn)
library(ipred)
library(glmnet)
library(nnet)
library(MLeval)

#13 合并代谢物与基因-----
ps03 = merge.ps(ps1 = ps.ms,
                ps2 = ps.trans,
                N1 = 100,
                N2 = 100,
                scale = TRUE,
                onlygroup = TRUE,#不进行列合并，只用于区分不同域
                dat1.lab = "ms",
                dat2.lab = "trans")

map= sample_data(ps03)
head(map)
i=1
id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
ps.cs = ps03 %>% subset_samples.wt("Group" ,id.g[,i])

# 保存合并数据信息
addWorksheet(trans_metab_ml_wb, "merged_data_info")
writeData(trans_metab_ml_wb, "merged_data_info", as.data.frame(map), rowNames = TRUE)

#14 loadingPCA.omics:载荷矩阵筛选特征代谢物与基因------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(trans_metab_ml_path, "loadingPCA_omics.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_ml_path, "loadingPCA_omics.pdf"), plot = p34.1, width = 10, height = 8)
addWorksheet(trans_metab_ml_wb, "loadingPCA_results")
writeData(trans_metab_ml_wb, "loadingPCA_results", dat, rowNames = TRUE)

#15 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps = ps03 ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存交叉验证结果
ggsave(file.path(trans_metab_ml_path, "rfcv_omics.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_ml_path, "rfcv_omics.pdf"), plot = prfcv, width = 10, height = 8)
addWorksheet(trans_metab_ml_wb, "rfcv_results")
writeData(trans_metab_ml_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#16 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = ps.cs %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
dat1 =  res[[2]]
dat1
dat2 =  res[[3]]
dat2

# 保存ROC结果
ggsave(file.path(trans_metab_ml_path, "ROC_omics_plot1.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_ml_path, "ROC_omics_plot1.pdf"), plot = p33.1, width = 10, height = 8)
addWorksheet(trans_metab_ml_wb, "ROC_AUC")
writeData(trans_metab_ml_wb, "ROC_AUC", dat1, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "ROC_results")
writeData(trans_metab_ml_wb, "ROC_results", dat2, rowNames = TRUE)

#17 randomforest.omics: 随机森林筛选特征基因与代谢物----
res = randomforest.omics( ps = ps03 %>% filter_OTU_ps(100),
                          group  = "Group",
                          optimal = 50)

p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
colnames(dat)
p42.4 =res[[4]]
p42.4

# 保存随机森林结果
ggsave(file.path(trans_metab_ml_path, "randomforest_omics_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_ml_path, "randomforest_omics_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(trans_metab_ml_path, "randomforest_omics_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_ml_path, "randomforest_omics_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(trans_metab_ml_path, "randomforest_omics_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_ml_path, "randomforest_omics_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(trans_metab_ml_wb, "randomforest_results")
writeData(trans_metab_ml_wb, "randomforest_results", dat, rowNames = TRUE)

#18 svm: 筛选特征基因与代谢物------
res <- svm.omics(ps = ps03 %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(trans_metab_ml_wb, "svm_AUC")
writeData(trans_metab_ml_wb, "svm_AUC", AUC, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "svm_importance")
writeData(trans_metab_ml_wb, "svm_importance", importance, rowNames = TRUE)

#19 GLM: 筛选特征基因与代谢物---------
res <- glm.omics(ps = ps.cs %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(trans_metab_ml_wb, "glm_AUC")
writeData(trans_metab_ml_wb, "glm_AUC", AUC, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "glm_importance")
writeData(trans_metab_ml_wb, "glm_importance", importance, rowNames = TRUE)

#20 xgboost: 筛选特征基因与代谢物-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance = importance$importance

# 保存XGBoost结果
addWorksheet(trans_metab_ml_wb, "xgboost_accuracy")
writeData(trans_metab_ml_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "xgboost_importance")
writeData(trans_metab_ml_wb, "xgboost_importance", importance, rowNames = TRUE)

#21 kNN: 筛选特征基因与代谢物--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存kNN结果
addWorksheet(trans_metab_ml_wb, "knn_accuracy")
writeData(trans_metab_ml_wb, "knn_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "knn_importance")
writeData(trans_metab_ml_wb, "knn_importance", importance, rowNames = TRUE)

#22 decisiontree: 筛选特征基因与代谢物------
res =decisiontree.omics(ps=ps03, top =1000,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(trans_metab_ml_wb, "decisiontree_accuracy")
writeData(trans_metab_ml_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "decisiontree_importance")
writeData(trans_metab_ml_wb, "decisiontree_importance", importance, rowNames = TRUE)

#23 bagging: 筛选特征基因与代谢物-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(trans_metab_ml_wb, "bagging_accuracy")
writeData(trans_metab_ml_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "bagging_importance")
writeData(trans_metab_ml_wb, "bagging_importance", importance, rowNames = TRUE)

#24 lasso: 筛选特征基因与代谢物-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(trans_metab_ml_wb, "lasso_accuracy")
writeData(trans_metab_ml_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "lasso_importance")
writeData(trans_metab_ml_wb, "lasso_importance", importance, rowNames = TRUE)

#25 naivebayes: 筛选特征基因与代谢物------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(trans_metab_ml_wb, "naivebayes_accuracy")
writeData(trans_metab_ml_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "naivebayes_importance")
writeData(trans_metab_ml_wb, "naivebayes_importance", importance, rowNames = TRUE)

#26 nnet神经网络: 筛选特征基因与代谢物------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(trans_metab_ml_wb, "nnet_accuracy")
writeData(trans_metab_ml_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_metab_ml_wb, "nnet_importance")
writeData(trans_metab_ml_wb, "nnet_importance", importance, rowNames = TRUE)

# 保存机器学习总表
saveWorkbook(trans_metab_ml_wb, file.path(trans_metab_ml_path, "transcriptome_metabolome_machine_learning.xlsx"), overwrite = TRUE)

# 网络分析----
# 创建网络分析目录
trans_metab_network_path <- file.path(trans_metab_path, "network")
dir.create(trans_metab_network_path, recursive = TRUE)

# 创建网络分析总表
trans_metab_network_wb <- createWorkbook()

#27 代谢物与基因跨域网络-----
library(sna)
library(igraph)
detach("package:mia", unload = TRUE)
otu1 =  otu_table(ps03)
row.names(otu1) = gsub("-","_", row.names(otu1))
otu_table(ps03)  = otu1
res =ggClusterNet:: corBionetwork.st(
  ps.st= ps03,# phyloseq对象
  N= 500,
  g1 = "Group",# 分组1
  g2 = NULL,# 分组2
  g3 = NULL,# 分组3
  ord.g1 = NULL, # 排序顺序
  ord.g2 = NULL, # 排序顺序
  ord.g3 = NULL, # 排序顺序
  order = NULL, # 出图每行代表的变量
  fill = "filed",
  size = "igraph.degree",
  method = "spearman",
  clu_method = "cluster_fast_greedy",
  select_layout = TRUE,
  layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 5,
  scale = TRUE,
  env = NULL,
  bio = TRUE,
  minsize = 4,
  maxsize = 14)

dat = res[[2]]
cor = dat$cortab
p = res[[1]]
p

# 保存跨域网络结果
ggsave(file.path(trans_metab_network_path, "corBionetwork.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(trans_metab_network_path, "corBionetwork.pdf"), plot = p, width = 12, height = 10)
addWorksheet(trans_metab_network_wb, "corBionetwork_data")
writeData(trans_metab_network_wb, "corBionetwork_data", dat$cortab, rowNames = TRUE)

#28 hub.network.omics-------
res <- hub.network.omics(cor = cor, top = 10)
p = res [[1]]
dat = res [[2]]

# 保存Hub网络结果
ggsave(file.path(trans_metab_network_path, "hub_network.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_metab_network_path, "hub_network.pdf"), plot = p, width = 10, height = 8)
addWorksheet(trans_metab_network_wb, "hub_network_data")
writeData(trans_metab_network_wb, "hub_network_data", dat, rowNames = TRUE)

#29 module.compare.net.pip -----
dat = module.compare.net.pip(
  ps = NULL,
  corg = cor,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 3)

res = dat[[1]]

# 保存模块比较结果
addWorksheet(trans_metab_network_wb, "module_compare_results")
writeData(trans_metab_network_wb, "module_compare_results", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(trans_metab_network_wb, file.path(trans_metab_network_path, "transcriptome_metabolome_network_analysis.xlsx"), overwrite = TRUE)

#30 pathway_show.omics 无----

#31 module_show.omics 无----

#32 reaction.show.omics 无-----

#33 reaction_network.omics 无----

#34 pathway.gene.compound.omics 无-----
