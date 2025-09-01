# 转录与宏基因组联合------
rm(list=ls())
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(EasyStat)

#转录与宏基因联合------
ps02 = EasyMultiOmics::ps.trans
ps01 = EasyMultiOmics::ps.kegg
sample_data(ps01)
sample_data(ps02)
#0 转录与宏基因数据合并----
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

#4 mantal.omics ：mantal相似性检测#------
dat = mantal.omics(ps01= ps01,ps02= ps02)
dat

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

#6 heatmap.line.omics---------
res=heatmap.line.omics(ps01=ps.kegg,
                       ps02=ps.trans,
                       lab.1 = "meta",
                       lab.2 = "trans")
grid.draw(res[[1]])
res[[2]]
res[[3]]
res[[4]]

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

#9 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps =pst ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#10 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = pst %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

#11 loadingPCA.omics:载荷矩阵筛选特征功能------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

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

#13 svm: 筛选转录与宏基因组特征基因------
res <- svm.omics(ps = ps.cs %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#14 GLM: 筛选转录与宏基因组特征基因---------
res <- glm.omics(ps = ps.cs %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#15 xgboost: 筛选转录与宏基因组特征基因-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#16 kNN: 筛选转录与宏基因组特征基因--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#17 decisiontree: 筛选转录与宏基因组特征基因------
res =decisiontree.omics(ps=ps03, top =1000,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#18 bagging: 筛选转录与宏基因组特征基因-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#19 lasso: 筛选转录与宏基因组特征基因-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#20 naivebayes: 筛选转录与宏基因组特征基因------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#21 nnet神经网络: 筛选转录与宏基因组特征基因------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#22 宏基因与转录跨域网络-----
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

#23 hub.network.omics-------
res <- hub.network.omics(cor = cor, top = 10)
p = res [[1]]
dat = res [[2]]

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
