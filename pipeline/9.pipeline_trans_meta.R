# 转录与宏基因组联合------
rm(list=ls())
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(EasyStat)
library(openxlsx)

#转录与宏基因联合------
ps02 = EasyMultiOmics::ps.trans
ps01 = EasyMultiOmics::ps.kegg
sample_data(ps01)
sample_data(ps02)

#0 转录与宏基因数据合并----
# 创建转录与宏基因组联合分析目录
trans_meta_path <- "./result/transcriptome_metagenome/"
dir.create(trans_meta_path, recursive = TRUE)

ps03 = merge.ps(ps1 = ps01,
                ps2 = ps02,
                N1 = 100,
                N2 = 100,
                scale = TRUE,
                onlygroup = TRUE,#不进行列合并，只用于区分不同域
                dat1.lab = "meta",
                dat2.lab = "trans")

ps03
# 注意tax文件所产生的变化，OTU的名字加上了前缀，增加了filed列，定义细菌还是真菌数据
tax= ps03%>% vegan_tax() %>% as.data.frame()
head(tax)

# otu文件的变化
otu = ps03 %>% vegan_otu() %>% t() %>% as.data.frame()
head(otu)

#转录与宏基因分组map文件要完全一致
map = ps03 %>% sample_data()
head(map)

id = sample_data(ps03)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
sample_data(ps03)
ps.cs = pst = ps03 %>% subset_samples(Group %in% group)

#1 ordinate.omics: 排序分析#----------
# 创建Beta多样性分析目录
trans_meta_beta_path <- file.path(trans_meta_path, "beta_diversity")
dir.create(trans_meta_beta_path, recursive = TRUE)

result = ordinate.omics(ps = ps03, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                        pair = F)
p3_1 = result[[1]]
p3_1

#带标签图形出图
p3_2 = result[[3]]
p3_2

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
p3_3

# 保存Beta多样性结果
ggsave(file.path(trans_meta_beta_path, "PCoA_plot.png"), plot = p3_1, width = 8, height = 6, dpi = 300)
ggsave(file.path(trans_meta_beta_path, "PCoA_plot.pdf"), plot = p3_1, width = 8, height = 6)
ggsave(file.path(trans_meta_beta_path, "PCoA_with_labels.png"), plot = p3_2, width = 8, height = 6, dpi = 300)
ggsave(file.path(trans_meta_beta_path, "PCoA_with_labels.pdf"), plot = p3_2, width = 8, height = 6)
ggsave(file.path(trans_meta_beta_path, "PCoA_refined.png"), plot = p3_3, width = 8, height = 6, dpi = 300)
ggsave(file.path(trans_meta_beta_path, "PCoA_refined.pdf"), plot = p3_3, width = 8, height = 6)

# 创建Beta多样性分析总表
trans_meta_beta_wb <- createWorkbook()
addWorksheet(trans_meta_beta_wb, "PCoA_coordinates")
writeData(trans_meta_beta_wb, "PCoA_coordinates", plotdata, rowNames = TRUE)
addWorksheet(trans_meta_beta_wb, "centroids")
writeData(trans_meta_beta_wb, "centroids", cent, rowNames = TRUE)

#2 ordinateTest.omics:群落水平差异检测#-------
dat1 = ordinateTest.omics(ps = ps03, Micromet = "adonis", dist = "bray")
dat1

# 保存总体差异检测结果
addWorksheet(trans_meta_beta_wb, "overall_test")
writeData(trans_meta_beta_wb, "overall_test", dat1, rowNames = TRUE)

#3 pairordinateTest.omics:两两分组群落水平差异检测#-------
dat2 = pairordinateTest.omics(ps = ps03, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较结果
addWorksheet(trans_meta_beta_wb, "pairwise_test")
writeData(trans_meta_beta_wb, "pairwise_test", dat2, rowNames = TRUE)

#4 mantal.omics ：mantal相似性检测#------
dat = mantal.omics(ps01= ps01,ps02= ps02)
dat

# 保存mantel检验结果
addWorksheet(trans_meta_beta_wb, "mantel_test")
writeData(trans_meta_beta_wb, "mantel_test", dat, rowNames = TRUE)

#5 cluster.omics:样品聚类-----
res = cluster.omics (ps= ps03,hcluter_method = "complete",
                     dist = "bray",
                     cuttree = 3,
                     row_cluster = TRUE,
                     col_cluster =  TRUE)

p4 = res[[1]]
p4
p4_1 = res[[2]]
p4_1
p4_2 = res[[3]]
p4_2
dat = res[[4]]
dat

# 保存聚类分析结果
ggsave(file.path(trans_meta_beta_path, "cluster_plot1.png"), plot = p4, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_beta_path, "cluster_plot1.pdf"), plot = p4, width = 10, height = 8)
ggsave(file.path(trans_meta_beta_path, "cluster_plot2.png"), plot = p4_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_beta_path, "cluster_plot2.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(trans_meta_beta_path, "cluster_plot3.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_beta_path, "cluster_plot3.pdf"), plot = p4_2, width = 10, height = 8)
addWorksheet(trans_meta_beta_wb, "cluster_data")
writeData(trans_meta_beta_wb, "cluster_data", dat, rowNames = TRUE)
saveWorkbook(trans_meta_beta_wb, file.path(trans_meta_beta_path, "beta_diversity.xlsx"), overwrite = TRUE)

#6 heatmap.line.omics---------
# 创建相关性分析目录
trans_meta_correlation_path <- file.path(trans_meta_path, "correlation_analysis")
dir.create(trans_meta_correlation_path, recursive = TRUE)

res=heatmap.line.omics(ps01=ps.kegg,
                       ps02=ps.trans,
                       lab.1 = "meta",
                       lab.2 = "trans")
grid.draw(res[[1]])
res[[2]]
res[[3]]
res[[4]]

# 保存热图连线结果
png(file.path(trans_meta_correlation_path, "heatmap_line.png"), width = 12, height = 10, units = "in", res = 300)
grid.draw(res[[1]])
dev.off()

pdf(file.path(trans_meta_correlation_path, "heatmap_line.pdf"), width = 12, height = 10)
grid.draw(res[[1]])
dev.off()

# 创建相关性分析总表
trans_meta_correlation_wb <- createWorkbook()
addWorksheet(trans_meta_correlation_wb, "heatmap_line_data1")
writeData(trans_meta_correlation_wb, "heatmap_line_data1", res[[2]], rowNames = TRUE)
addWorksheet(trans_meta_correlation_wb, "heatmap_line_data2")
writeData(trans_meta_correlation_wb, "heatmap_line_data2", res[[3]], rowNames = TRUE)
addWorksheet(trans_meta_correlation_wb, "heatmap_line_data3")
writeData(trans_meta_correlation_wb, "heatmap_line_data3", res[[4]], rowNames = TRUE)

#7 volcano.line.omics: 火山互连 ----
library(ggnewscale)
ps1 = ps.trans %>% filter_OTU_ps(500)
ps2 = ps.ms %>% filter_OTU_ps(500)
res = volcano.line.omics(ps1 = ps1,
                         ps2 = ps2,
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

# 保存火山互连结果
ggsave(file.path(trans_meta_correlation_path, "volcano_line_WT_OE.png"), plot = res[[1]]$WT_OE, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_correlation_path, "volcano_line_WT_OE.pdf"), plot = res[[1]]$WT_OE, width = 10, height = 8)
ggsave(file.path(trans_meta_correlation_path, "volcano_line_WT_KO.png"), plot = res[[1]]$WT_KO, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_correlation_path, "volcano_line_WT_KO.pdf"), plot = res[[1]]$WT_KO, width = 10, height = 8)
ggsave(file.path(trans_meta_correlation_path, "volcano_line_OE_KO.png"), plot = res[[1]]$OE_KO, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_correlation_path, "volcano_line_OE_KO.pdf"), plot = res[[1]]$OE_KO, width = 10, height = 8)

res= volcano.line2.omics(ps01=ps.kegg,
                         ps02=ps.trans,
                         group1= "meta",
                         group2= "trans",
                         r.threshold = 0.6,
                         p.threshold = 0.1,
                         method = "spearman",
                         top = 500)

p=res[[1]]$KO_OE
p
dat=res[[2]]$KO_OE

# 保存火山互连结果2
ggsave(file.path(trans_meta_correlation_path, "volcano_line2_KO_OE.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_correlation_path, "volcano_line2_KO_OE.pdf"), plot = p, width = 10, height = 8)
addWorksheet(trans_meta_correlation_wb, "volcano_line2_KO_OE")
writeData(trans_meta_correlation_wb, "volcano_line2_KO_OE", dat, rowNames = TRUE)

#8 quare.line.omics: 九象限火山图 ----
results <- quare.line.omics(ps.16s =ps.kegg,
                            ps.ms =  ps.trans,
                            group1 = "kegg",
                            group2 = "trans",
                            r.threshold = 0.7,
                            p.threshold = 0.05,
                            method = "spearman",
                            top = 500)
results$KO_WT
results$OE_WT$data

# 保存九象限图结果
ggsave(file.path(trans_meta_correlation_path, "nine_quadrant_KO_WT.png"), plot = results$KO_WT, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_correlation_path, "nine_quadrant_KO_WT.pdf"), plot = results$KO_WT, width = 10, height = 8)
addWorksheet(trans_meta_correlation_wb, "nine_quadrant_OE_WT")
writeData(trans_meta_correlation_wb, "nine_quadrant_OE_WT", results$OE_WT$data, rowNames = TRUE)
saveWorkbook(trans_meta_correlation_wb, file.path(trans_meta_correlation_path, "correlation_analysis.xlsx"), overwrite = TRUE)

# 机器学习----
# 创建机器学习目录
trans_meta_ml_path <- file.path(trans_meta_path, "machine_learning")
dir.create(trans_meta_ml_path, recursive = TRUE)

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

#9 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps =pst ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存交叉验证结果
ggsave(file.path(trans_meta_ml_path, "rfcv_plot.png"), plot = prfcv, width = 8, height = 6, dpi = 300)
ggsave(file.path(trans_meta_ml_path, "rfcv_plot.pdf"), plot = prfcv, width = 8, height = 6)

# 创建机器学习分析总表
trans_meta_ml_wb <- createWorkbook()
addWorksheet(trans_meta_ml_wb, "rfcv_results")
writeData(trans_meta_ml_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#10 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = pst %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

# 保存ROC曲线结果
ggsave(file.path(trans_meta_ml_path, "ROC_plot1.png"), plot = p33.1, width = 8, height = 6, dpi = 300)
ggsave(file.path(trans_meta_ml_path, "ROC_plot1.pdf"), plot = p33.1, width = 8, height = 6)
ggsave(file.path(trans_meta_ml_path, "ROC_plot2.png"), plot = p33.2, width = 8, height = 6, dpi = 300)
ggsave(file.path(trans_meta_ml_path, "ROC_plot2.pdf"), plot = p33.2, width = 8, height = 6)
addWorksheet(trans_meta_ml_wb, "ROC_data")
writeData(trans_meta_ml_wb, "ROC_data", dat, rowNames = TRUE)

#11 loadingPCA.omics:载荷矩阵筛选特征功能------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(trans_meta_ml_path, "PCA_loading_plot.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_ml_path, "PCA_loading_plot.pdf"), plot = p34.1, width = 10, height = 8)
addWorksheet(trans_meta_ml_wb, "PCA_loading")
writeData(trans_meta_ml_wb, "PCA_loading", dat, rowNames = TRUE)

#12 randomforest.omics: 随机森林筛选特征功能----
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
ggsave(file.path(trans_meta_ml_path, "randomforest_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_ml_path, "randomforest_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(trans_meta_ml_path, "randomforest_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_ml_path, "randomforest_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(trans_meta_ml_path, "randomforest_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_ml_path, "randomforest_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(trans_meta_ml_wb, "randomforest_features")
writeData(trans_meta_ml_wb, "randomforest_features", dat, rowNames = TRUE)

#13 svm: 筛选转录与宏基因组特征基因------
res <- svm.omics(ps = ps.cs %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(trans_meta_ml_wb, "SVM_AUC")
writeData(trans_meta_ml_wb, "SVM_AUC", AUC, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "SVM_importance")
writeData(trans_meta_ml_wb, "SVM_importance", importance, rowNames = TRUE)

#14 GLM: 筛选转录与宏基因组特征基因---------
res <- glm.omics(ps = ps.cs %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(trans_meta_ml_wb, "GLM_AUC")
writeData(trans_meta_ml_wb, "GLM_AUC", AUC, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "GLM_importance")
writeData(trans_meta_ml_wb, "GLM_importance", importance, rowNames = TRUE)

#15 xgboost: 筛选转录与宏基因组特征基因-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(trans_meta_ml_wb, "XGBoost_accuracy")
writeData(trans_meta_ml_wb, "XGBoost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "XGBoost_importance")
writeData(trans_meta_ml_wb, "XGBoost_importance", importance, rowNames = TRUE)

#16 kNN: 筛选转录与宏基因组特征基因--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存kNN结果
addWorksheet(trans_meta_ml_wb, "kNN_accuracy")
writeData(trans_meta_ml_wb, "kNN_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "kNN_importance")
writeData(trans_meta_ml_wb, "kNN_importance", importance, rowNames = TRUE)

#17 decisiontree: 筛选转录与宏基因组特征基因------
res =decisiontree.omics(ps=ps03, top =1000,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(trans_meta_ml_wb, "DecisionTree_accuracy")
writeData(trans_meta_ml_wb, "DecisionTree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "DecisionTree_importance")
writeData(trans_meta_ml_wb, "DecisionTree_importance", importance, rowNames = TRUE)

#18 bagging: 筛选转录与宏基因组特征基因-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(trans_meta_ml_wb, "Bagging_accuracy")
writeData(trans_meta_ml_wb, "Bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "Bagging_importance")
writeData(trans_meta_ml_wb, "Bagging_importance", importance, rowNames = TRUE)

#19 lasso: 筛选转录与宏基因组特征基因-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(trans_meta_ml_wb, "Lasso_accuracy")
writeData(trans_meta_ml_wb, "Lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "Lasso_importance")
writeData(trans_meta_ml_wb, "Lasso_importance", importance, rowNames = TRUE)

#20 naivebayes: 筛选转录与宏基因组特征基因------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(trans_meta_ml_wb, "NaiveBayes_accuracy")
writeData(trans_meta_ml_wb, "NaiveBayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "NaiveBayes_importance")
writeData(trans_meta_ml_wb, "NaiveBayes_importance", importance, rowNames = TRUE)

#21 nnet神经网络: 筛选转录与宏基因组特征基因------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(trans_meta_ml_wb, "NeuralNet_accuracy")
writeData(trans_meta_ml_wb, "NeuralNet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_meta_ml_wb, "NeuralNet_importance")
writeData(trans_meta_ml_wb, "NeuralNet_importance", importance, rowNames = TRUE)

# 保存机器学习分析总表
saveWorkbook(trans_meta_ml_wb, file.path(trans_meta_ml_path, "machine_learning.xlsx"), overwrite = TRUE)

#22 宏基因与转录跨域网络-----
# 创建网络分析目录
trans_meta_network_path <- file.path(trans_meta_path, "network_analysis")
dir.create(trans_meta_network_path, recursive = TRUE)

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

# 保存跨域网络图
ggsave(file.path(trans_meta_network_path, "cross_domain_network.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(trans_meta_network_path, "cross_domain_network.pdf"), plot = p, width = 12, height = 10)

# 创建网络分析总表
trans_meta_network_wb <- createWorkbook()
addWorksheet(trans_meta_network_wb, "network_correlation_data")
writeData(trans_meta_network_wb, "network_correlation_data", dat$cortab, rowNames = TRUE)

#23 hub.network.omics-------
res <- hub.network.omics(cor = cor, top = 10)
p = res [[1]]
dat = res [[2]]

# 保存Hub网络结果
ggsave(file.path(trans_meta_network_path, "hub_network.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_meta_network_path, "hub_network.pdf"), plot = p, width = 10, height = 8)
addWorksheet(trans_meta_network_wb, "hub_network_data")
writeData(trans_meta_network_wb, "hub_network_data", dat, rowNames = TRUE)

#24 module.compare.net.pip -----
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
addWorksheet(trans_meta_network_wb, "module_compare_data")
writeData(trans_meta_network_wb, "module_compare_data", res, rowNames = TRUE)

# 保存总表
saveWorkbook(trans_meta_network_wb, file.path(trans_meta_network_wb, "trans_meta_network.xlsx"), overwrite = TRUE)


