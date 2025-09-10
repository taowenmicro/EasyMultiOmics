# 多组学联合----
rm(list=ls())
library(vegan)
library(phyloseq)
library(dplyr)
library(ggClusterNet)
library(openxlsx)
# 导入数据-
ps.trans=ps.trans
ps.ms=ps.ms
ps.micro=ps.micro
ps.meta=ps.kegg
col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")
package.amp()

#0 合并多组学ps对象#-------
ps.all = merge.ps(ps1 = ps.trans,
                  ps2 = ps.ms,
                  N1 = 0,
                  N2 = 0,
                  scale = TRUE,
                  onlygroup = TRUE,#不进行列合并，只用于区分不同域
                  dat1.lab = "trans",
                  dat2.lab = "ms")
ps.all

ps.all2 = merge.ps(ps1 = ps.all,
                   ps2 = ps.micro,
                   N1 = 0,
                   N2 = 0,
                   scale = TRUE,
                   onlygroup = TRUE,#不进行列合并，只用于区分不同域
                   dat1.lab = "",
                   dat2.lab = "micro")

tax = ps.all2 %>% vegan_tax() %>% as.data.frame()
tax$filed %>% unique()

ps03 = merge.ps(ps1 = ps.all2,
                ps2 = ps.kegg,
                N1 = 0,
                N2 = 0,
                scale = TRUE,
                onlygroup = TRUE,#不进行列合并，只用于区分不同域
                dat1.lab = "",
                dat2.lab = "meta")

tax = ps03 %>% vegan_tax() %>% as.data.frame()
tax$filed %>% unique()

# otu文件的变化
otu = ps03 %>% vegan_otu() %>% t() %>% as.data.frame()
head(otu)

#转录与宏基因分组map文件要完全一致
map = ps03 %>% sample_data()
head(map)


# 创建结果保存根目录
multi_omics_path <- file.path("./result/multi_omics_analysis")
dir.create(multi_omics_path, recursive = TRUE, showWarnings = FALSE)

#1 ordinate.omics: 排序分析#----------
# 创建排序分析目录
ordinate_path <- file.path(multi_omics_path, "ordinate")
dir.create(ordinate_path, recursive = TRUE, showWarnings = FALSE)

# 创建排序分析总表
ordinate_wb <- createWorkbook()

result = ordinate.omics(ps = ps03, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                        pair = F)
p3_1 = result[[1]]
p3_1 = p3_1 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +
  theme_classic()
# 保存排序图
ggsave(file.path(ordinate_path, "ordinate_PCoA.png"), plot = p3_1, width = 8, height = 6, dpi = 300)
ggsave(file.path(ordinate_path, "ordinate_PCoA.pdf"), plot = p3_1, width = 8, height = 6)


#带标签图形出图
p3_2 = result[[3]]
p3_2 = p3_2 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +
  theme_classic()
# 保存带标签排序图
ggsave(file.path(ordinate_path, "ordinate_PCoA_with_labels.png"), plot = p3_2, width = 8, height = 6, dpi = 300)
ggsave(file.path(ordinate_path, "ordinate_PCoA_with_labels.pdf"), plot = p3_2, width = 8, height = 6)

#---------排序-精修图
plotdata =result[[2]]
head(plotdata)
# 求均值
cent <- aggregate(cbind(x,y) ~Group, data = plotdata, FUN = mean)
cent
# 合并到样本坐标数据中
segs <- merge(plotdata, setNames(cent, c('Group','oNMDS1','oNMDS2')),
              by = 'Group', sort = FALSE)

library(ggsci)
p3_3 = p3_1 +geom_segment(data = segs,
                          mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow")
p3_3 = p3_3 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +
  theme_classic()
# 保存精修排序图
ggsave(file.path(ordinate_path, "ordinate_PCoA_refined.png"), plot = p3_3, width = 8, height = 6, dpi = 300)
ggsave(file.path(ordinate_path, "ordinate_PCoA_refined.pdf"), plot = p3_3, width = 8, height = 6)


addWorksheet(ordinate_wb, "PCoA_data")
writeData(ordinate_wb, "PCoA_data", plotdata, rowNames = TRUE)

#2 ordinateTest.omics:群落水平差异检测#-------
dat1 = ordinateTest.omics(ps = ps03, Micromet = "adonis", dist = "bray")
dat1
addWorksheet(ordinate_wb, "ordinateTest_adonis")
writeData(ordinate_wb, "ordinateTest_adonis", dat1, rowNames = TRUE)

#3 pairordinateTest.omics:两两分组群落水平差异检测#-------
dat2 = pairordinateTest.omics(ps = ps03, Micromet = "MRPP", dist = "bray")
dat2
addWorksheet(ordinate_wb, "pairordinateTest_MRPP")
writeData(ordinate_wb, "pairordinateTest_MRPP", dat2, rowNames = TRUE)

saveWorkbook(ordinate_wb, file.path(ordinate_path, "ordinate_analysis_results.xlsx"), overwrite = TRUE)


#4 cluster.omics:样品聚类#-----
# 创建聚类分析目录
cluster_path <- file.path(multi_omics_path, "cluster")
dir.create(cluster_path, recursive = TRUE, showWarnings = FALSE)

# 创建聚类分析总表
cluster_wb <- createWorkbook()

res = cluster.omics (ps= ps03,hcluter_method = "complete",
                     dist = "bray",
                     cuttree = 3,
                     row_cluster = TRUE,
                     col_cluster =  TRUE)
p4 = res[[1]]
# 保存聚类热图
ggsave(file.path(cluster_path, "cluster_heatmap.png"), plot = p4, width = 10, height = 8, dpi = 300)
ggsave(file.path(cluster_path, "cluster_heatmap.pdf"), plot = p4, width = 10, height = 8)

p4_1 = res[[2]]
# 保存聚类热图（带树）
ggsave(file.path(cluster_path, "cluster_heatmap_with_tree.png"), plot = p4_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cluster_path, "cluster_heatmap_with_tree.pdf"), plot = p4_1, width = 10, height = 8)

p4_2 = res[[3]]
# 保存聚类图（带树）
ggsave(file.path(cluster_path, "cluster_tree.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cluster_path, "cluster_tree.pdf"), plot = p4_2, width = 10, height = 8)

dat = res[[4]]
addWorksheet(cluster_wb, "cluster_data")
writeData(cluster_wb, "cluster_data", dat, rowNames = TRUE)
saveWorkbook(cluster_wb, file.path(cluster_path, "cluster_analysis_results.xlsx"), overwrite = TRUE)


# 机器学习----
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

# 创建机器学习目录
ml_path <- file.path(multi_omics_path, "machine_learning")
dir.create(ml_path, recursive = TRUE, showWarnings = FALSE)

# 创建机器学习总表
ml_wb <- createWorkbook()


map= sample_data(ps03)
head(map)
i=1
id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
ps.cs = ps03 %>% subset_samples.wt("Group" ,id.g[,i])

#5 loadingPCA.omics------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
# 保存 loadingPCA 图
ggsave(file.path(ml_path, "loading_PCA.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "loading_PCA.pdf"), plot = p34.1, width = 10, height = 8)
dat = res[[2]]
addWorksheet(ml_wb, "loadingPCA_data")
writeData(ml_wb, "loadingPCA_data", dat, rowNames = TRUE)


#6 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps = ps03 ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
# 保存交叉验证图
ggsave(file.path(ml_path, "rfcv_cross_validation.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "rfcv_cross_validation.pdf"), plot = prfcv, width = 10, height = 8)
# result[[2]]# plotdata
rfcvtable = result[[3]]
addWorksheet(ml_wb, "rfcv_data")
writeData(ml_wb, "rfcv_data", rfcvtable, rowNames = TRUE)


#7 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = ps.cs %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
# 保存ROC曲线图
ggsave(file.path(ml_path, "ROC_curve.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "ROC_curve.pdf"), plot = p33.1, width = 10, height = 8)

dat1 =  res[[2]]
addWorksheet(ml_wb, "ROC_AUC")
writeData(ml_wb, "ROC_AUC", dat, rowNames = TRUE)

dat2 =  res[[3]]
addWorksheet(ml_wb, "ROC_data")
writeData(ml_wb, "ROC_data", dat, rowNames = TRUE)

#8 randomforest.omics: 随机森林筛选多组学标志物----
res = randomforest.omics( ps = ps03 %>% filter_OTU_ps(100),
                          group  = "Group",
                          optimal = 50)
p42.1 = res[[1]]
# 保存随机森林重要性图1
ggsave(file.path(ml_path, "randomforest_importance1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "randomforest_importance1.pdf"), plot = p42.1, width = 10, height = 8)
p42.2 =res[[2]]
# 保存随机森林重要性图2
ggsave(file.path(ml_path, "randomforest_importance2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "randomforest_importance2.pdf"), plot = p42.2, width = 10, height = 8)
dat =res[[3]]
addWorksheet(ml_wb, "randomforest_data")
writeData(ml_wb, "randomforest_data", dat, rowNames = TRUE)
colnames(dat)
p42.4 =res[[4]]
# 保存随机森林ROC图
ggsave(file.path(ml_path, "randomforest_ROC.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "randomforest_ROC.pdf"), plot = p42.4, width = 10, height = 8)

#9 svm: 筛选多组学标志物------
res <- svm.omics(ps = ps03 %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
addWorksheet(ml_wb, "svm_AUC")
writeData(ml_wb, "svm_AUC", AUC, rowNames = TRUE)
importance = res[[2]]
addWorksheet(ml_wb, "svm_importance")
writeData(ml_wb, "svm_importance", importance, rowNames = TRUE)

#10 GLM: 筛选多组学标志物---------
res <- glm.omics(ps = ps.cs %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
addWorksheet(ml_wb, "glm_AUC")
writeData(ml_wb, "glm_AUC", AUC, rowNames = TRUE)
importance = res[[2]]
addWorksheet(ml_wb, "glm_importance")
writeData(ml_wb, "glm_importance", importance, rowNames = TRUE)

#11 xgboost: 筛选多组学标志物-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
addWorksheet(ml_wb, "xgboost_accuracy")
writeData(ml_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
importance = res[[2]]
importance = importance$importance
addWorksheet(ml_wb, "xgboost_importance")
writeData(ml_wb, "xgboost_importance", importance, rowNames = TRUE)

#12 kNN: 筛选多组学标志物--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
addWorksheet(ml_wb, "knn_accuracy")
writeData(ml_wb, "knn_accuracy", accuracy, rowNames = TRUE)
importance = res[[2]]
importance = importance$importance
addWorksheet(ml_wb, "knn_importance")
writeData(ml_wb, "knn_importance", importance, rowNames = TRUE)

#13 decisiontree: 筛选多组学标志物------
res =decisiontree.omics(ps=ps03, top =1000,  k = 5)
accuracy = res[[1]]
addWorksheet(ml_wb, "decisiontree_accuracy")
writeData(ml_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
importance = res[[2]]
addWorksheet(ml_wb, "decisiontree_importance")
writeData(ml_wb, "decisiontree_importance", importance, rowNames = TRUE)

#14 bagging: 筛选多组学标志物-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
addWorksheet(ml_wb, "bagging_accuracy")
writeData(ml_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
importance = res[[2]]
addWorksheet(ml_wb, "bagging_importance")
writeData(ml_wb, "bagging_importance", importance, rowNames = TRUE)

#15 lasso: 筛选多组学标志物-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
addWorksheet(ml_wb, "lasso_accuracy")
writeData(ml_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
importance = res[[2]]
addWorksheet(ml_wb, "lasso_importance")
writeData(ml_wb, "lasso_importance", importance, rowNames = TRUE)

#16 naivebayes: 筛选多组学标志物------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
addWorksheet(ml_wb, "naivebayes_accuracy")
writeData(ml_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
importance = res[[2]]
addWorksheet(ml_wb, "naivebayes_importance")
writeData(ml_wb, "naivebayes_importance", importance, rowNames = TRUE)

#17 nnet神经网络: 筛选多组学标志物------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
addWorksheet(ml_wb, "nnet_accuracy")
writeData(ml_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
importance = res[[2]]
addWorksheet(ml_wb, "nnet_importance")
writeData(ml_wb, "nnet_importance", importance, rowNames = TRUE)

saveWorkbook(ml_wb, file.path(ml_path, "machine_learning_results.xlsx"), overwrite = TRUE)


# lavaan 结构方程模型----
# 创建结构方程模型目录
sem_path <- file.path(multi_omics_path, "structural_equation_models")
dir.create(sem_path, recursive = TRUE, showWarnings = FALSE)

# 创建结构方程模型总表
sem_wb <- createWorkbook()

res <- lavaan.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
                    ps.meta=ps.kegg,filter=0.05)
addWorksheet(sem_wb, "lavaan_summary")
writeData(sem_wb, "lavaan_summary", res, rowNames = TRUE)

# piecewiseSEM结构方程模型----
library(piecewiseSEM)
res <- piecewiseSEM.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
                          ps.meta=ps.kegg,filter=0.05)
r_squared= res[[1]]
aic = res[[2]]
addWorksheet(sem_wb, "piecewiseSEM_r_squared")
writeData(sem_wb, "piecewiseSEM_r_squared", r_squared, rowNames = TRUE)
addWorksheet(sem_wb, "piecewiseSEM_aic")
writeData(sem_wb, "piecewiseSEM_aic", aic, rowNames = TRUE)

# plspm结构方程模型----
res <- plspm.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
                   ps.meta=ps.kegg,filter=0.05)
dat = res [[1]]
dat = dat$crossloadings
cor =  res [[2]]
effect =  res [[3]]
addWorksheet(sem_wb, "plspm_dat")
writeData(sem_wb, "plspm_dat", dat, rowNames = TRUE)
addWorksheet(sem_wb, "plspm_cor")
writeData(sem_wb, "plspm_cor", cor, rowNames = TRUE)
addWorksheet(sem_wb, "plspm_effect")
writeData(sem_wb, "plspm_effect", effect, rowNames = TRUE)

saveWorkbook(sem_wb, file.path(sem_path, "structural_equation_models_results.xlsx"), overwrite = TRUE)

