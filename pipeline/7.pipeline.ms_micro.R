# 微生物与代谢-------
rm(list=ls())
library(openxlsx)

ps.micro =  EasyMultiOmics::ps.micro
ps.ms = EasyMultiOmics::ps.ms
ps.tem = ps.micro %>% filter_OTU_ps(500)

#-主题--颜色等
package.amp()
res = theme_my(ps.micro)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

#0 代谢物筛选-------
# 创建微生物与代谢物联合分析目录
microbe_metabolite_path <- "./result/microbe_metabolite/"
dir.create(microbe_metabolite_path, recursive = TRUE)

res = loadingPCA.ms(ps = ps.ms)
dat = res[[2]]
dat$id = row.names(dat)
id = dat %>% arrange(desc(PCone)) %>% head(15) %>% .$id
id

ftab = ps.ms %>% scale_micro() %>%
  subset_taxa(
    row.names(tax_table(ps.ms)) %in% c(id)
  ) %>% vegan_otu() %>%
  as.data.frame()

head(ftab)

#1 RDA_CCA:对微生物影响的代谢物共排序-----
# 创建约束分析目录
microbe_metabolite_constraint_path <- file.path(microbe_metabolite_path, "constraint_analysis")
dir.create(microbe_metabolite_constraint_path, recursive = TRUE)

result = RDA_CCA(ps = ps.tem,
                 env = ftab,
                 # path = RDApath,
                 chose.env = FALSE
)

p1 = result[[1]]
p1_final = p1+theme_classic()
# 提取作图数据
dataplot = result[[2]]
dataplot
# 提取带有标记的图片
p2 = result[[3]]
p2
aov = result[[4]]
aov

# 保存RDA/CCA结果
ggsave(file.path(microbe_metabolite_constraint_path, "RDA_CCA_plot1.png"), plot = p1_final, width = 8, height = 6, dpi = 300)
ggsave(file.path(microbe_metabolite_constraint_path, "RDA_CCA_plot1.pdf"), plot = p1_final, width = 8, height = 6)
ggsave(file.path(microbe_metabolite_constraint_path, "RDA_CCA_plot2.png"), plot = p2, width = 8, height = 6, dpi = 300)
ggsave(file.path(microbe_metabolite_constraint_path, "RDA_CCA_plot2.pdf"), plot = p2, width = 8, height = 6)

# 创建约束分析总表
microbe_metabolite_constraint_wb <- createWorkbook()
addWorksheet(microbe_metabolite_constraint_wb, "RDA_CCA_data")
writeData(microbe_metabolite_constraint_wb, "RDA_CCA_data", dataplot, rowNames = TRUE)
addWorksheet(microbe_metabolite_constraint_wb, "RDA_CCA_anova")
writeData(microbe_metabolite_constraint_wb, "RDA_CCA_anova", aov, rowNames = TRUE)
addWorksheet(microbe_metabolite_constraint_wb, "selected_metabolites")
writeData(microbe_metabolite_constraint_wb, "selected_metabolites", ftab, rowNames = TRUE)

#2 RDA_CCA_explain_percent: 对微生物群落的解释比例-----
result = RDA_CCA_explain_percent(ps = ps.tem,
                                 env.dat = ftab)
sam = result[[1]]
sam
all = result[[2]]
all

# 保存解释比例结果
addWorksheet(microbe_metabolite_constraint_wb, "explain_sample")
writeData(microbe_metabolite_constraint_wb, "explain_sample", sam, rowNames = TRUE)
addWorksheet(microbe_metabolite_constraint_wb, "explain_all")
writeData(microbe_metabolite_constraint_wb, "explain_all", all, rowNames = TRUE)

#3 rdacca.hp.micro:   rdacca.hp.micro.p  层次分割对微生物影响的代谢物#---------
library(rdacca.hp)
library(vegan)

res = rdacca.hp.micro(OTU = ps.tem %>% filter_OTU_ps(500) %>%vegan_otu() %>% as.data.frame(),
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
ggsave(file.path(microbe_metabolite_constraint_path, "hierarchical_partition1.png"), plot = p3_rda, width = 8, height = 6, dpi = 300)
ggsave(file.path(microbe_metabolite_constraint_path, "hierarchical_partition1.pdf"), plot = p3_rda, width = 8, height = 6)
ggsave(file.path(microbe_metabolite_constraint_path, "hierarchical_partition2.png"), plot = p3_db_rda, width = 8, height = 6, dpi = 300)
ggsave(file.path(microbe_metabolite_constraint_path, "hierarchical_partition2.pdf"), plot = p3_db_rda, width = 8, height = 6)
addWorksheet(microbe_metabolite_constraint_wb, "hierarchical_data1")
writeData(microbe_metabolite_constraint_wb, "hierarchical_data1", dat1, rowNames = TRUE)
addWorksheet(microbe_metabolite_constraint_wb, "hierarchical_data2")
writeData(microbe_metabolite_constraint_wb, "hierarchical_data2", dat2, rowNames = TRUE)

#4 rdacca.hp.micro.p: 层次分割对微生物影响的代谢物#---------
res <-rdacca.hp.micro.p(
  OTU = ps.tem %>% filter_OTU_ps(200) %>%vegan_otu() %>% as.data.frame(),
  env = ftab[,1:5],
  cca = FALSE,
  dbRDA = FALSE
)

p3_rda = res[[1]]
p3_rda
dat1 = res[[2]]
dat1

# 保存带p值的层次分割结果
ggsave(file.path(microbe_metabolite_constraint_path, "hierarchical_partition_p.png"), plot = p3_rda, width = 8, height = 6, dpi = 300)
ggsave(file.path(microbe_metabolite_constraint_path, "hierarchical_partition_p.pdf"), plot = p3_rda, width = 8, height = 6)
addWorksheet(microbe_metabolite_constraint_wb, "hierarchical_p_data")
writeData(microbe_metabolite_constraint_wb, "hierarchical_p_data", dat1, rowNames = TRUE)
saveWorkbook(microbe_metabolite_constraint_wb, file.path(microbe_metabolite_constraint_path, "constraint_analysis.xlsx"), overwrite = TRUE)

# 相关性分析-------
# 创建相关性分析目录
microbe_metabolite_correlation_path <- file.path(microbe_metabolite_path, "correlation_analysis")
dir.create(microbe_metabolite_correlation_path, recursive = TRUE)

#7 MetalTast:  单一代谢物与微生物群落的相关性计算----
otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps.tem)))
head(otu1)
tabOTU1 = list(bac = otu1 )
rep = MetalTast (env.dat = ftab, tabOTU = tabOTU1,distance = "bray",method = "metal")
head(rep)

# 保存MetalTast结果
microbe_metabolite_correlation_wb <- createWorkbook()
addWorksheet(microbe_metabolite_correlation_wb, "MetalTast_results")
writeData(microbe_metabolite_correlation_wb, "MetalTast_results", rep, rowNames = TRUE)

#8 cor_env_ggcorplot_top:高丰度微生物和代谢物的关系探索----
result = cor_env_ggcorplot_top(ps =  ps.tem,
                               env1 = ftab,
                               label =  TRUE,
                               col_cluster = TRUE,
                               row_cluster = TRUE,
                               method = "spearman",
                               r.threshold= 0,
                               p.threshold= 0  )

p1 <- result[[1]]
p1
p2 <- result[[2]]
p2
topotu =  result[[3]]
topotu
cordata = result[[4]]
cordata

# 保存高丰度相关性分析结果
ggsave(file.path(microbe_metabolite_correlation_path, "cor_top_heatmap.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "cor_top_heatmap.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(microbe_metabolite_correlation_path, "cor_top_bubble.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "cor_top_bubble.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(microbe_metabolite_correlation_wb, "cor_top_otu")
writeData(microbe_metabolite_correlation_wb, "cor_top_otu", topotu, rowNames = TRUE)
addWorksheet(microbe_metabolite_correlation_wb, "cor_top_data")
writeData(microbe_metabolite_correlation_wb, "cor_top_data", cordata, rowNames = TRUE)

#9 cor_env_ggcorplot_top:丰度差异大的微生物和代谢物关系探索----
result = cor_env_ggcorplot_rcv(ps =  ps.tem,
                               env1 = ftab,
                               label = TRUE,
                               col_cluster = TRUE,
                               row_cluster = TRUE,
                               method = "spearman",
                               r.threshold= 0,
                               p.threshold= 0  )

p1 <- result[[1]]
p1
p2 <- result[[2]]
p2
rcvotu =  result[[3]]
rcvotu
cordata = result[[4]]
cordata

# 保存差异大相关性分析结果
ggsave(file.path(microbe_metabolite_correlation_path, "cor_rcv_heatmap.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "cor_rcv_heatmap.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(microbe_metabolite_correlation_path, "cor_rcv_bubble.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "cor_rcv_bubble.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(microbe_metabolite_correlation_wb, "cor_rcv_otu")
writeData(microbe_metabolite_correlation_wb, "cor_rcv_otu", rcvotu, rowNames = TRUE)
addWorksheet(microbe_metabolite_correlation_wb, "cor_rcv_data")
writeData(microbe_metabolite_correlation_wb, "cor_rcv_data", cordata, rowNames = TRUE)

#10 MatCorPlot:微生物群落和代谢Science组合图表------
otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps.tem)))
head(otu1)

tabOTU1 = list(bac = otu1 )
p <-  MatCorPlot(env.dat = ftab,
                 tabOTU = tabOTU1,
                 diag =  FALSE,
                 range = 0.2,
                 numpoint = 21,
                 sig =  FALSE,
                 siglabel = FALSE,
                 shownum =  FALSE,
                 curvature = 0,
                 numsymbol = NULL,
                 lacx = "right",
                 lacy = "bottom",
                 p.thur = 0.05,
                 onlysig = TRUE)

p

# 保存MatCorPlot结果
ggsave(file.path(microbe_metabolite_correlation_path, "MatCorPlot.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "MatCorPlot.pdf"), plot = p, width = 12, height = 10)

#11 heatmap.line.omics---------
res=heatmap.line.omics(ps01=ps.tem,
                       ps02=ps.ms,
                       lab.1 = "micro",
                       lab.2 = "ms")
grid.draw(res[[1]])
res[[2]]
res[[3]]
res[[4]]

# 保存热图连线结果
png(file.path(microbe_metabolite_correlation_path, "heatmap_line.png"), width = 12, height = 10, units = "in", res = 300)
grid.draw(res[[1]])
dev.off()

pdf(file.path(microbe_metabolite_correlation_path, "heatmap_line.pdf"), width = 12, height = 10)
grid.draw(res[[1]])
dev.off()

addWorksheet(microbe_metabolite_correlation_wb, "heatmap_line_data1")
writeData(microbe_metabolite_correlation_wb, "heatmap_line_data1", res[[2]], rowNames = TRUE)
addWorksheet(microbe_metabolite_correlation_wb, "heatmap_line_data2")
writeData(microbe_metabolite_correlation_wb, "heatmap_line_data2", res[[3]], rowNames = TRUE)
addWorksheet(microbe_metabolite_correlation_wb, "heatmap_line_data3")
writeData(microbe_metabolite_correlation_wb, "heatmap_line_data3", res[[4]], rowNames = TRUE)

#12 volcano.line2.omics: 火山互连 ----
res= volcano.line2.omics(ps01=ps.tem,
                         ps02 = ps.ms,
                         group1= "micro",
                         group2= "ms",
                         r.threshold = 0.6,
                         p.threshold = 0.1,
                         method = "spearman",
                         top = 500)

p = res[[1]]$KO_OE
p
dat = res[[2]]$KO_OE

# 保存火山互连结果
ggsave(file.path(microbe_metabolite_correlation_path, "volcano_line_KO_OE.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "volcano_line_KO_OE.pdf"), plot = p, width = 10, height = 8)
addWorksheet(microbe_metabolite_correlation_wb, "volcano_line_KO_OE")
writeData(microbe_metabolite_correlation_wb, "volcano_line_KO_OE", dat, rowNames = TRUE)

#13 quare.line.omics: 九象限火山图 ----
results <- quare.line.omics(ps.16s = ps.tem,
                            ps.ms =  ps.ms,
                            group1 = "microbe",
                            group2 = "metabolite",
                            r.threshold = 0.7,
                            p.threshold = 0.05,
                            method = "spearman",
                            top = 500)

p = results$KO_WT
p
dat = results$OE_WT$data

# 保存九象限图结果
ggsave(file.path(microbe_metabolite_correlation_path, "nine_quadrant_KO_WT.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "nine_quadrant_KO_WT.pdf"), plot = p, width = 10, height = 8)
addWorksheet(microbe_metabolite_correlation_wb, "nine_quadrant_OE_WT")
writeData(microbe_metabolite_correlation_wb, "nine_quadrant_OE_WT", dat, rowNames = TRUE)

#14 随机森林寻找对微生物群落影响较大的代谢物--------
result <- ms_micro.rf.omics(ps  = ps.tem,env = ftab )
p <- result[[2]]
p
data = result[[1]]
head(data)

# 保存随机森林筛选结果
ggsave(file.path(microbe_metabolite_correlation_path, "rf_importance.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_correlation_path, "rf_importance.pdf"), plot = p, width = 10, height = 8)
addWorksheet(microbe_metabolite_correlation_wb, "rf_importance_data")
writeData(microbe_metabolite_correlation_wb, "rf_importance_data", data, rowNames = TRUE)
saveWorkbook(microbe_metabolite_correlation_wb, file.path(microbe_metabolite_correlation_path, "correlation_analysis.xlsx"), overwrite = TRUE)

# 机器学习----
# 创建机器学习目录
microbe_metabolite_ml_path <- file.path(microbe_metabolite_path, "machine_learning")
dir.create(microbe_metabolite_ml_path, recursive = TRUE)

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

#15 合并微生物与代谢物-----
ps0 = merge.ps(ps1 = ps.micro,
               ps2 = ps.ms,
               N1 = 100,
               N2 = 100,
               scale = TRUE,
               onlygroup = TRUE,#不进行列合并，只用于区分不同域
               dat1.lab = "16S",
               dat2.lab = "compounds")

map= sample_data(ps0)
head(map)
i=1
id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
ps.cs = ps0 %>% subset_samples.wt("Group" ,id.g[,i])
ps03 = ps.cs

#16 loadingPCA.omics:载荷矩阵筛选特征微生物与代谢物------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(microbe_metabolite_ml_path, "PCA_loading_plot.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_ml_path, "PCA_loading_plot.pdf"), plot = p34.1, width = 10, height = 8)

# 创建机器学习分析总表
microbe_metabolite_ml_wb <- createWorkbook()
addWorksheet(microbe_metabolite_ml_wb, "PCA_loading")
writeData(microbe_metabolite_ml_wb, "PCA_loading", dat, rowNames = TRUE)

#17 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps = ps03 ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存交叉验证结果
ggsave(file.path(microbe_metabolite_ml_path, "rfcv_plot.png"), plot = prfcv, width = 8, height = 6, dpi = 300)
ggsave(file.path(microbe_metabolite_ml_path, "rfcv_plot.pdf"), plot = prfcv, width = 8, height = 6)
addWorksheet(microbe_metabolite_ml_wb, "rfcv_results")
writeData(microbe_metabolite_ml_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#18 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = ps03 %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
dat =  res[[3]]
dat

# 保存ROC曲线结果
ggsave(file.path(microbe_metabolite_ml_path, "ROC_plot.png"), plot = p33.1, width = 8, height = 6, dpi = 300)
ggsave(file.path(microbe_metabolite_ml_path, "ROC_plot.pdf"), plot = p33.1, width = 8, height = 6)
addWorksheet(microbe_metabolite_ml_wb, "ROC_data")
writeData(microbe_metabolite_ml_wb, "ROC_data", dat, rowNames = TRUE)

#19 randomforest.omics: 随机森林筛选特征微生物与代谢物----
res = randomforest.omics( ps = ps03 ,
                          group  = "Group",
                          optimal = 50)
p42.1 = res[[1]]
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林结果
ggsave(file.path(microbe_metabolite_ml_path, "randomforest_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_ml_path, "randomforest_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(microbe_metabolite_ml_path, "randomforest_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_ml_path, "randomforest_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(microbe_metabolite_ml_path, "randomforest_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_ml_path, "randomforest_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(microbe_metabolite_ml_wb, "randomforest_features")
writeData(microbe_metabolite_ml_wb, "randomforest_features", dat, rowNames = TRUE)

#20 svm: 筛选特征微生物与代谢物------
res <- svm.omics(ps = ps03 %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(microbe_metabolite_ml_wb, "SVM_AUC")
writeData(microbe_metabolite_ml_wb, "SVM_AUC", AUC, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "SVM_importance")
writeData(microbe_metabolite_ml_wb, "SVM_importance", importance, rowNames = TRUE)

#21 GLM: 筛选特征微生物与代谢物---------
res <- glm.omics(ps = ps03 %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(microbe_metabolite_ml_wb, "GLM_AUC")
writeData(microbe_metabolite_ml_wb, "GLM_AUC", AUC, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "GLM_importance")
writeData(microbe_metabolite_ml_wb, "GLM_importance", importance, rowNames = TRUE)

#22 xgboost: 筛选特征微生物与代谢物-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(microbe_metabolite_ml_wb, "XGBoost_accuracy")
writeData(microbe_metabolite_ml_wb, "XGBoost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "XGBoost_importance")
writeData(microbe_metabolite_ml_wb, "XGBoost_importance", importance, rowNames = TRUE)

#23 kNN: 筛选特征微生物与代谢物--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存kNN结果
addWorksheet(microbe_metabolite_ml_wb, "kNN_accuracy")
writeData(microbe_metabolite_ml_wb, "kNN_accuracy", accuracy, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "kNN_importance")
writeData(microbe_metabolite_ml_wb, "kNN_importance", importance, rowNames = TRUE)

#24 decisiontree: 筛选特征微生物与代谢物------
res =decisiontree.omics(ps=ps03, top =1000,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(microbe_metabolite_ml_wb, "DecisionTree_accuracy")
writeData(microbe_metabolite_ml_wb, "DecisionTree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "DecisionTree_importance")
writeData(microbe_metabolite_ml_wb, "DecisionTree_importance", importance, rowNames = TRUE)

#25 bagging: 筛选特征微生物与代谢物-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(microbe_metabolite_ml_wb, "Bagging_accuracy")
writeData(microbe_metabolite_ml_wb, "Bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "Bagging_importance")
writeData(microbe_metabolite_ml_wb, "Bagging_importance", importance, rowNames = TRUE)

#26 lasso: 筛选特征微生物与代谢物-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(microbe_metabolite_ml_wb, "Lasso_accuracy")
writeData(microbe_metabolite_ml_wb, "Lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "Lasso_importance")
writeData(microbe_metabolite_ml_wb, "Lasso_importance", importance, rowNames = TRUE)

#27 naivebayes: 筛选特征微生物与代谢物------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(microbe_metabolite_ml_wb, "NaiveBayes_accuracy")
writeData(microbe_metabolite_ml_wb, "NaiveBayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "NaiveBayes_importance")
writeData(microbe_metabolite_ml_wb, "NaiveBayes_importance", importance, rowNames = TRUE)

#28 nnet神经网络: 筛选特征微生物与代谢物------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(microbe_metabolite_ml_wb, "NeuralNet_accuracy")
writeData(microbe_metabolite_ml_wb, "NeuralNet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(microbe_metabolite_ml_wb, "NeuralNet_importance")
writeData(microbe_metabolite_ml_wb, "NeuralNet_importance", importance, rowNames = TRUE)

# 保存机器学习分析总表
saveWorkbook(microbe_metabolite_ml_wb, file.path(microbe_metabolite_ml_path, "machine_learning.xlsx"), overwrite = TRUE)

#29 微生物与代谢物跨域网络-----
# 创建网络分析目录
microbe_metabolite_network_path <- file.path(microbe_metabolite_path, "network_analysis")
dir.create(microbe_metabolite_network_path, recursive = TRUE)

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
ggsave(file.path(microbe_metabolite_network_path, "cross_domain_network.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(microbe_metabolite_network_path, "cross_domain_network.pdf"), plot = p, width = 12, height = 10)

# 创建网络分析总表
microbe_metabolite_network_wb <- createWorkbook()
addWorksheet(microbe_metabolite_network_wb, "network_correlation_data")
writeData(microbe_metabolite_network_wb, "network_correlation_data", dat$cortab, rowNames = TRUE)

#30 hub.network.omics-------
res <- hub.network.omics(cor = cor, top = 10)
p = res [[1]]
dat = res [[2]]

# 保存Hub网络结果
ggsave(file.path(microbe_metabolite_network_path, "hub_network.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(microbe_metabolite_network_path, "hub_network.pdf"), plot = p, width = 10, height = 8)
addWorksheet(microbe_metabolite_network_wb, "hub_network_data")
writeData(microbe_metabolite_network_wb, "hub_network_data", dat, rowNames = TRUE)

#31 module.compare.net.pip -----
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

# 保存网络模块比较结果
addWorksheet(microbe_metabolite_network_wb, "module_comparison")
writeData(microbe_metabolite_network_wb, "module_comparison", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(microbe_metabolite_network_wb, file.path(microbe_metabolite_network_path, "network_analysis.xlsx"), overwrite = TRUE)
