library(tidyverse)
library(data.table)
library(phyloseq)
library(ggClusterNet)
library(EasyStat)
library(tidyfst)
library(fs)

# 宏基因组成与功能-------
rm(list=ls())

ps.metm= EasyMultiOmics::ps.micro%>% filter_OTU_ps(Top = 1000)
ps.metf = EasyMultiOmics::ps.kegg%>% filter_OTU_ps(Top = 1000)

sample_data(ps.metf)
sample_data(ps.metm)
#-主题--颜色等
package.amp()

res = theme_my(ps.metf)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")

# 创建宏基因组分析目录
metagenome_path =  "./result/metf_metm/"
dir.create(metagenome_path, recursive = TRUE)

# 创建宏基因组成分与功能分析目录
metf_path <- file.path(metagenome_path, "composition_function")
dir.create(metf_path, recursive = TRUE)

# 创建宏基因组成分与功能分析总表
metf_wb <- createWorkbook()

#0 微生物物筛选-------
res = loadingPCA.ms(ps = ps.metm)
# p = res[[1]]
# p
dat = res[[2]]
dat$id = row.names(dat)
id = dat %>% arrange(desc(PCone)) %>% head(15) %>% .$id
id

ftab = ps.metm %>% scale_micro() %>%
  subset_taxa(
    row.names(tax_table(ps.metm)) %in% c(id)
  ) %>% vegan_otu() %>%
  as.data.frame()

head(ftab)

# 保存微生物筛选结果
addWorksheet(metf_wb, "microbe_selection_PCA")
writeData(metf_wb, "microbe_selection_PCA", dat, rowNames = TRUE)
addWorksheet(metf_wb, "filtered_microbes")
writeData(metf_wb, "filtered_microbes", ftab, rowNames = TRUE)

#1 RDA_CCA:对群落功能影响的微生物共排序-----
result = RDA_CCA(ps = ps.metf,
                 env = ftab,
                 # path = RDApath,
                 chose.env = FALSE
)

p1 = result[[1]]
p1+theme_classic()+scale_fill_manual(values = col.g)
# 提取作图数据
dataplot = result[[2]]
dataplot

# 提取带有标记的图片
p2 = result[[3]]

p2+theme_classic()+scale_fill_manual(values = col.g)
aov = result[[4]]
aov

# 保存RDA/CCA结果
ggsave(file.path(metf_path, "RDA_CCA_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(metf_path, "RDA_CCA_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(metf_path, "RDA_CCA_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(metf_path, "RDA_CCA_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(metf_wb, "RDA_CCA_plotdata")
writeData(metf_wb, "RDA_CCA_plotdata", dataplot, rowNames = TRUE)
addWorksheet(metf_wb, "RDA_CCA_anova")
writeData(metf_wb, "RDA_CCA_anova", aov, rowNames = TRUE)

#2 RDA_CCA_explain_percent: 对微生物群落功能的解释比例-----
result = RDA_CCA_explain_percent(ps = ps.metf,
                                 env.dat = ftab)
sam = result[[1]]
sam
all = result[[2]]
all

# 保存解释比例结果
addWorksheet(metf_wb, "RDA_CCA_explain_sample")
writeData(metf_wb, "RDA_CCA_explain_sample", sam, rowNames = TRUE)
addWorksheet(metf_wb, "RDA_CCA_explain_all")
writeData(metf_wb, "RDA_CCA_explain_all", all, rowNames = TRUE)

#3 rdacca.hp.micro: 层次分割对群落功能影响的微生物---------
library(rdacca.hp)
library(vegan)
library(ggClusterNet)
library(tidyverse)
res = rdacca.hp.micro(OTU = ps.metf %>% filter_OTU_ps(500) %>%vegan_otu() %>% as.data.frame(),
                      env = ftab[,1:5],
                      cca = FALSE)
p3_rda = res[[1]]
p3_rda
dat1 = res[[2]]
dat1
p3_db_rda  = res[[3]]
p3_db_rda
dat2 = res[[4]]
dat2

# 保存层次分割结果
ggsave(file.path(metf_path, "rdacca_hp_plot1.png"), plot = p3_rda, width = 10, height = 8, dpi = 300)
ggsave(file.path(metf_path, "rdacca_hp_plot1.pdf"), plot = p3_rda, width = 10, height = 8)
ggsave(file.path(metf_path, "rdacca_hp_plot2.png"), plot = p3_db_rda, width = 10, height = 8, dpi = 300)
ggsave(file.path(metf_path, "rdacca_hp_plot2.pdf"), plot = p3_db_rda, width = 10, height = 8)
addWorksheet(metf_wb, "rdacca_hp_data1")
writeData(metf_wb, "rdacca_hp_data1", dat1, rowNames = TRUE)
addWorksheet(metf_wb, "rdacca_hp_data2")
writeData(metf_wb, "rdacca_hp_data2", dat2, rowNames = TRUE)

#4 rdacca.hp.micro.p: 层次分割(显著)对群落功能影响的微生物---------
res <-rdacca.hp.micro.p(
  OTU = ps.metf %>% filter_OTU_ps(200) %>%vegan_otu() %>% as.data.frame(),
  env = ftab[,1:5],
  cca = FALSE,
  dbRDA = FALSE)

p3_rda = res[[1]]
p3_rda
dat1 = res[[2]]
dat1

# 保存显著层次分割结果
ggsave(file.path(metf_path, "rdacca_hp_sig.png"), plot = p3_rda, width = 10, height = 8, dpi = 300)
ggsave(file.path(metf_path, "rdacca_hp_sig.pdf"), plot = p3_rda, width = 10, height = 8)
addWorksheet(metf_wb, "rdacca_hp_sig_data")
writeData(metf_wb, "rdacca_hp_sig_data", dat1, rowNames = TRUE)

# 相关性分析-------
# 创建相关性分析目录
correlation_path <- file.path(metagenome_path, "correlation")
dir.create(correlation_path, recursive = TRUE)

# 创建相关性分析总表
correlation_wb <- createWorkbook()

#7 MetalTast:  单一功能与微生物群落的相关性计算----
otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps.metf)))
head(otu1)
tabOTU1 = list(bac = otu1 )
rep = MetalTast (env.dat = ftab, tabOTU = tabOTU1,distance = "bray",method = "metal")
rep

# 保存MetalTast结果
addWorksheet(correlation_wb, "metalTast_results")
writeData(correlation_wb, "metalTast_results", rep, rowNames = TRUE)

#8 cor_env_ggcorplot_top:高丰度功能和微生物的关系探索----
result = cor_env_ggcorplot_top(ps =  ps.metf,
                               env1 = ftab,
                               label =  TRUE,
                               col_cluster = TRUE,
                               row_cluster = TRUE,
                               method = "spearman",
                               r.threshold= 0,
                               p.threshold= 0 )
p1 <- result[[1]]
p1
p2 <- result[[2]]
p2
topotu =  result[[3]]
topotu
cordata = result[[4]]
cordata

# 保存高丰度相关性结果
ggsave(file.path(correlation_path, "cor_top_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(correlation_path, "cor_top_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(correlation_path, "cor_top_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(correlation_path, "cor_top_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(correlation_wb, "cor_top_otu")
writeData(correlation_wb, "cor_top_otu", topotu, rowNames = TRUE)
addWorksheet(correlation_wb, "cor_top_data")
writeData(correlation_wb, "cor_top_data", cordata, rowNames = TRUE)

#9 cor_env_ggcorplot_top:丰度差异大的功能和微生物关系探索----
result = cor_env_ggcorplot_rcv(ps =  ps.metf,
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

# 保存差异大的相关性结果
ggsave(file.path(correlation_path, "cor_rcv_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(correlation_path, "cor_rcv_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(correlation_path, "cor_rcv_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(correlation_path, "cor_rcv_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(correlation_wb, "cor_rcv_otu")
writeData(correlation_wb, "cor_rcv_otu", rcvotu, rowNames = TRUE)
addWorksheet(correlation_wb, "cor_rcv_data")
writeData(correlation_wb, "cor_rcv_data", cordata, rowNames = TRUE)

#10 MatCorPlot:微生物群落和代谢Science组合图表------
#detach("package:mia", unload = TRUE)
otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps.metf)))
head(otu1)
tabOTU1 = list(bac = otu1 )
# tabOTU1 = list(b = otu)
# tabOTU1 = list(f = otu2)
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
                 onlysig = FALSE)

p

# 保存MatCorPlot结果
ggsave(file.path(correlation_path, "MatCorPlot.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(correlation_path, "MatCorPlot.pdf"), plot = p, width = 12, height = 10)

# 保存相关性分析总表
saveWorkbook(correlation_wb, file.path(correlation_path, "metagenome_correlation_analysis.xlsx"), overwrite = TRUE)

# 创建多组学整合分析目录
omics_path <- file.path(metagenome_path, "multi_omics")
dir.create(omics_path, recursive = TRUE)

# 创建多组学整合分析总表
omics_wb <- createWorkbook()

#11 heatmap.line.omics---------
res=heatmap.line.omics(ps01=ps.metf,
                       ps02=ps.metm,
                       lab.1 = "metf",
                       lab.2 = "metm")
grid.draw(res[[1]])
res[[2]]
res[[3]]
res[[4]]

# 保存热图线性组学结果
ggsave(file.path(omics_path, "heatmap_line_omics.png"), plot = res[[1]], width = 12, height = 10, dpi = 300)
ggsave(file.path(omics_path, "heatmap_line_omics.pdf"), plot = res[[1]], width = 12, height = 10)
addWorksheet(omics_wb, "heatmap_line_data1")
writeData(omics_wb, "heatmap_line_data1", res[[2]], rowNames = TRUE)
addWorksheet(omics_wb, "heatmap_line_data2")
writeData(omics_wb, "heatmap_line_data2", res[[3]], rowNames = TRUE)
addWorksheet(omics_wb, "heatmap_line_data3")
writeData(omics_wb, "heatmap_line_data3", res[[4]], rowNames = TRUE)

#12 volcano.line2.omics: 火山互连 ----
res= volcano.line2.omics(ps01 = ps.metm,
                         ps02 = ps.metf,
                         group1= "micro",
                         group2= "metf",
                         r.threshold = 0.6,
                         p.threshold = 0.1,
                         method = "spearman",
                         top = 500)

res[[1]]$KO_OE
res[[2]]$KO_OE

# 保存火山互连结果
if (!is.null(res[[1]]$KO_OE)) {
  ggsave(file.path(omics_path, "volcano_line_KO_OE_plot1.png"), plot = res[[1]]$KO_OE, width = 10, height = 8, dpi = 300)
  ggsave(file.path(omics_path, "volcano_line_KO_OE_plot1.pdf"), plot = res[[1]]$KO_OE, width = 10, height = 8)
}
if (!is.null(res[[2]]$KO_OE)) {
  addWorksheet(omics_wb, "volcano_line_KO_OE_data")
  writeData(omics_wb, "volcano_line_KO_OE_data", res[[2]]$KO_OE, rowNames = TRUE)
}

#13 quare.line.omics: 九象限火山图 ----
results <- quare.line.omics(ps.16s = ps.metm,
                            ps.ms =  ps.metf,
                            group1 = "community",
                            group2 = "function",
                            r.threshold = 0.7,
                            p.threshold = 0.05,
                            method = "spearman",
                            top = 500)
results$OE_WT
results$OE_WT$data

# 保存九象限火山图结果
if (!is.null(results$OE_WT)) {
  ggsave(file.path(omics_path, "quare_line_OE_WT.png"), plot = results$OE_WT, width = 10, height = 8, dpi = 300)
  ggsave(file.path(omics_path, "quare_line_OE_WT.pdf"), plot = results$OE_WT, width = 10, height = 8)
}
if (!is.null(results$OE_WT$data)) {
  addWorksheet(omics_wb, "quare_line_OE_WT_data")
  writeData(omics_wb, "quare_line_OE_WT_data", results$OE_WT$data, rowNames = TRUE)
}

#14 随机森林寻找对微生物群落功能影响较大的微生物--------
result <- ms_micro.rf.omics(ps  = ps.metf,env = ftab )
p <- result[[2]]
p
data = result[[1]]
head(data)

# 保存随机森林结果
ggsave(file.path(omics_path, "rf_importance.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(omics_path, "rf_importance.pdf"), plot = p, width = 10, height = 8)
addWorksheet(omics_wb, "rf_importance_data")
writeData(omics_wb, "rf_importance_data", data, rowNames = TRUE)

# 保存成分与功能分析总表
saveWorkbook(metf_wb, file.path(metf_path, "metagenome_composition_function.xlsx"), overwrite = TRUE)


# 机器学习----
# 创建机器学习分析目录
ml_path <- file.path(metagenome_path, "machine_learning")
dir.create(ml_path, recursive = TRUE)

# 创建机器学习分析总表
ml_wb <- createWorkbook()

library(ggClusterNet)
library(phyloseq)
library(tidyverse)
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

#15 合并微生物与功能-----
ps03 = merge.ps(ps1 = ps.metm,
                ps2 = ps.metf,
                N1 = 100,
                N2 = 100,
                scale = TRUE,
                onlygroup = TRUE,#不进行列合并，只用于区分不同域
                dat1.lab = "metm",
                dat2.lab = "metf")

map= sample_data(ps03)
head(map)
i=1
id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
pst = ps03 %>% subset_samples.wt("Group" ,id.g[,i])

# 保存合并数据信息
addWorksheet(ml_wb, "merged_data_info")
writeData(ml_wb, "merged_data_info", as.data.frame(map), rowNames = TRUE)

#16 loadingPCA.omics:载荷矩阵筛选特征微生物与代谢物------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(ml_path, "loadingPCA_omics.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "loadingPCA_omics.pdf"), plot = p34.1, width = 10, height = 8)
addWorksheet(ml_wb, "loadingPCA_results")
writeData(ml_wb, "loadingPCA_results", dat, rowNames = TRUE)

#17 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps = ps03 ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv+theme_classic()
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存交叉验证结果
ggsave(file.path(ml_path, "rfcv_omics.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "rfcv_omics.pdf"), plot = prfcv, width = 10, height = 8)
addWorksheet(ml_wb, "rfcv_results")
writeData(ml_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#18 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = pst %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
AUC =  res[[2]]
AUC
dat =  res[[3]]
dat

# 保存ROC结果
ggsave(file.path(ml_path, "ROC_omics_plot1.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "ROC_omics_plot1.pdf"), plot = p33.1, width = 10, height = 8)
addWorksheet(ml_wb, "ROC_AUC")
writeData(ml_wb, "ROC_AUC", AUC, rowNames = TRUE)
addWorksheet(ml_wb, "ROC_results")
writeData(ml_wb, "ROC_results", dat, rowNames = TRUE)

#19 randomforest.omics: 随机森林筛选----
res = randomforest.omics( ps = ps03 %>% filter_OTU_ps(100),
                          group  = "Group",
                          optimal = 50)

p42.1 = res[[1]]
p42.1+mytheme1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
colnames(dat)
p42.4 =res[[4]]
p42.4

# 保存随机森林结果
ggsave(file.path(ml_path, "randomforest_omics_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "randomforest_omics_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(ml_path, "randomforest_omics_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "randomforest_omics_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(ml_path, "randomforest_omics_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(ml_path, "randomforest_omics_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(ml_wb, "randomforest_results")
writeData(ml_wb, "randomforest_results", dat, rowNames = TRUE)

#20 svm: 筛选特征微生物与代谢物------
res <- svm.omics(ps = ps03 %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(ml_wb, "svm_AUC")
writeData(ml_wb, "svm_AUC", AUC, rowNames = TRUE)
addWorksheet(ml_wb, "svm_importance")
writeData(ml_wb, "svm_importance", importance, rowNames = TRUE)

#21 GLM: 筛选特征微生物与代谢物---------
res <- glm.omics(ps = ps03 %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(ml_wb, "glm_AUC")
writeData(ml_wb, "glm_AUC", AUC, rowNames = TRUE)
addWorksheet(ml_wb, "glm_importance")
writeData(ml_wb, "glm_importance", importance, rowNames = TRUE)

#22 xgboost: 筛选特征微生物与代谢物-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance = importance$importance

# 保存XGBoost结果
addWorksheet(ml_wb, "xgboost_accuracy")
writeData(ml_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(ml_wb, "xgboost_importance")
writeData(ml_wb, "xgboost_importance", importance, rowNames = TRUE)

#23 kNN: 筛选特征微生物与代谢物--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存kNN结果
addWorksheet(ml_wb, "knn_accuracy")
writeData(ml_wb, "knn_accuracy", accuracy, rowNames = TRUE)
addWorksheet(ml_wb, "knn_importance")
writeData(ml_wb, "knn_importance", importance, rowNames = TRUE)

#24 decisiontree: 筛选特征微生物与代谢物 有问题------
res =decisiontree.omics(ps=ps03, top =500,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(ml_wb, "decisiontree_accuracy")
writeData(ml_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(ml_wb, "decisiontree_importance")
writeData(ml_wb, "decisiontree_importance", importance, rowNames = TRUE)

#25 bagging: 筛选特征微生物与代谢物-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(ml_wb, "bagging_accuracy")
writeData(ml_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(ml_wb, "bagging_importance")
writeData(ml_wb, "bagging_importance", importance, rowNames = TRUE)

#26 lasso: 筛选特征微生物与代谢物-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(ml_wb, "lasso_accuracy")
writeData(ml_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(ml_wb, "lasso_importance")
writeData(ml_wb, "lasso_importance", importance, rowNames = TRUE)

#27 naivebayes: 筛选特征微生物与代谢物------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(ml_wb, "naivebayes_accuracy")
writeData(ml_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(ml_wb, "naivebayes_importance")
writeData(ml_wb, "naivebayes_importance", importance, rowNames = TRUE)

#28 nnet神经网络: 筛选特征微生物与代谢物------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(ml_wb, "nnet_accuracy")
writeData(ml_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(ml_wb, "nnet_importance")
writeData(ml_wb, "nnet_importance", importance, rowNames = TRUE)

# 保存机器学习总表
saveWorkbook(ml_wb, file.path(ml_path, "metagenome_machine_learning.xlsx"), overwrite = TRUE)

# 网络分析----
# 创建网络分析目录
network_path <- file.path(metagenome_path, "network")
dir.create(network_path, recursive = TRUE)

# 创建网络分析总表
network_wb <- createWorkbook()

#29 corBionetwork.st：微生物与代谢物跨域网络-----
library(sna)
library(igraph)
otu1 =  otu_table(ps03)
row.names(otu1) = gsub("-","_", row.names(otu1))
otu_table(ps03)  = otu1
detach("package:mia", unload = TRUE)

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
  select_layout = TRUE,layout_net = "model_maptree2",
  r.threshold=0.6,
  p.threshold=0.05,
  maxnode = 3,
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
ggsave(file.path(network_path, "corBionetwork.png"), plot = p, width = 12, height = 10, dpi = 300)
ggsave(file.path(network_path, "corBionetwork.pdf"), plot = p, width = 12, height = 10)
addWorksheet(network_wb, "corBionetwork_data")
writeData(network_wb, "corBionetwork_data", dat$cortab, rowNames = TRUE)

#30 hub.network.omics-------
res <- hub.network.omics(cor = cor, top = 10)
p = res [[1]]
p
dat = res [[2]]

# 保存Hub网络结果
ggsave(file.path(network_path, "hub_network.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(network_path, "hub_network.pdf"), plot = p, width = 10, height = 8)
addWorksheet(network_wb, "hub_network_data")
writeData(network_wb, "hub_network_data", dat, rowNames = TRUE)

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
res

# 保存模块比较结果
addWorksheet(network_wb, "module_compare_results")
writeData(network_wb, "module_compare_results", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(network_wb, file.path(network_path, "metagenome_network_analysis.xlsx"), overwrite = TRUE)


