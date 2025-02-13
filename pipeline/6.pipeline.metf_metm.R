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

#1 RDA_CCA:对群落功能影响的微生物共排序-----
result = RDA_CCA(ps = ps.metf,
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

#2 RDA_CCA_explain_percent: 对微生物群落功能的解释比例-----
result = RDA_CCA_explain_percent(ps = ps.metf,
                                 env.dat = ftab)
sam = result[[1]]
sam
all = result[[2]]

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

# 相关性分析-------
#7 MetalTast:  单一功能与微生物群落的相关性计算----
otu1 = as.data.frame(t(ggClusterNet::vegan_otu(ps.metf)))
head(otu1)
tabOTU1 = list(bac = otu1 )
rep = MetalTast (env.dat = ftab, tabOTU = tabOTU1,distance = "bray",method = "metal")
rep

#8 cor_env_ggcorplot_top:高丰度功能和微生物的关系探索----
result = cor_env_ggcorplot_top(ps =  ps.metf,
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

#10 MatCorPlot:微生物群落和代谢Science组合图表------
detach("package:mia", unload = TRUE)
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
                 onlysig = TRUE)
p

#11 heatmap.line.omics---------
res=heatmap.line.omics(ps01=ps.metf,
                       ps02=ps.metm,
                       lab.1 = "metf",
                       lab.2 = "metm")
grid.draw(res[[1]])
res[[2]]
res[[3]]
res[[4]]

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

#13 quare.line.omics: 九象限火山图 ----
results <- quare.line.omics(ps.16s = ps.metm,
                            ps.ms =  ps.metf,
                            group1 = "community",
                            group2 = "function",
                            r.threshold = 0.7,
                            p.threshold = 0.05,
                            method = "spearman",
                            top = 500)
results$KO_WT
results$OE_WT$data

#14 随机森林寻找对微生物群落功能影响较大的微生物--------
result <- ms_micro.rf.omics(ps  = ps.metf,env = ftab )
p <- result[[2]]
p
data = result[[1]]
head(data)

# 机器学习----
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

#16 loadingPCA.omics:载荷矩阵筛选特征微生物与代谢物------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#17 rfcv.omics :交叉验证结果-------
result =rfcv.omics(ps = ps03 ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#18 Roc.omics:ROC 曲线绘制----
res = Roc.omics( ps = pst %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

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

#20 svm: 筛选特征微生物与代谢物------
res <- svm.omics(ps = ps03 %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#21 GLM: 筛选特征微生物与代谢物---------
res <- glm.omics(ps = ps03 %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#22 xgboost: 筛选特征微生物与代谢物-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#23 kNN: 筛选特征微生物与代谢物--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#24 decisiontree: 筛选特征微生物与代谢物 有问题------
res =decisiontree.omics(ps=ps03, top =500,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#25 bagging: 筛选特征微生物与代谢物-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#26 lasso: 筛选特征微生物与代谢物-----
res =lasso.omics (ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#27 naivebayes: 筛选特征微生物与代谢物------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#28 nnet神经网络: 筛选特征微生物与代谢物------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#29 corBionetwork.st：微生物与代谢物跨域网络-----
library(sna)
library(igraph)
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
  select_layout = TRUE,layout_net = "model_maptree2",
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

#30 hub.network.omics-------
res <- hub.network.omics(cor = cor, top = 10)
p = res [[1]]
dat = res [[2]]

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

#32 wgcna.omics-----
