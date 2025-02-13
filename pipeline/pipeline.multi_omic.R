# 多组学联合----
rm(list=ls())
library(vegan)
library(phyloseq)
library(dplyr)
library(ggClusterNet)
# 导入数据-
ps.trans=ps.trans
ps.ms=ps.ms
ps.micro=ps.micro
ps.meta=ps.kegg

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
otu = ps.all3 %>% vegan_otu() %>% t() %>% as.data.frame()
head(otu)

#转录与宏基因分组map文件要完全一致
map = ps03 %>% sample_data()
head(map)

#1 ordinate.omics: 排序分析#----------
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

#2 ordinateTest.omics:群落水平差异检测#-------
dat1 = ordinateTest.omics(ps = ps03, Micromet = "adonis", dist = "bray")
dat1

#3 pairordinateTest.omics:两两分组群落水平差异检测#-------
dat2 = pairordinateTest.omics(ps = ps03, Micromet = "MRPP", dist = "bray")
dat2

#4 cluster.omics:样品聚类#-----
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

map= sample_data(ps03)
head(map)
i=1
id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
ps.cs = ps03 %>% subset_samples.wt("Group" ,id.g[,i])

#5 loadingPCA.omics------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#6 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps = ps03 ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#7 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = ps.cs %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

#8 randomforest.omics: 随机森林筛选多组学标志物----
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

#9 svm: 筛选多组学标志物------
res <- svm.omics(ps = ps03 %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#10 GLM: 筛选多组学标志物---------
res <- glm.omics(ps = ps.cs %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#11 xgboost: 筛选多组学标志物-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#12 kNN: 筛选多组学标志物--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#13 decisiontree: 筛选多组学标志物------
res =decisiontree.omics(ps=ps03, top =1000,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#14 bagging: 筛选多组学标志物-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#15 lasso: 筛选多组学标志物-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#16 naivebayes: 筛选多组学标志物------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#17 nnet神经网络: 筛选多组学标志物------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance



# lavaan 结构方程模型----
res <- lavaan.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
                    ps.meta=ps.kegg,filter=0.05)
res

# piecewiseSEM结构方程模型----
library(piecewiseSEM)
res <- piecewiseSEM.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
                    ps.meta=ps.kegg,filter=0.05)
r_squared= res[[1]]
aic = res[[2]]

# plspm结构方程模型----
res <- plspm.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
                          ps.meta=ps.kegg,filter=0.05)
dat = res [[1]]
cor =  res [[2]]
effect =  res [[3]]
