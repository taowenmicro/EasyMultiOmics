#细菌-真菌联合------
rm(list=ls())
# library(EasyMultiOmics)
ps02 = EasyMultiOmics::ps.16s
ps01 = EasyMultiOmics::ps.its

#0 细菌/真菌数据合并----
ps03 =  merge16S_ITS(ps16s = ps01,
                                   psITS = ps02,
                                   NITS = 500,
                                   N16s = 500,
                                   scale = FALSE,
                                   onlygroup =FALSE)

ps03

# 注意tax文件所产生的变化，OTU的名字加上了前缀，增加了filed列，定义细菌还是真菌数据
tax= ps03%>% vegan_tax() %>% as.data.frame()
head(tax)

# otu文件的变化
otu = ps03 %>% vegan_otu() %>% t() %>% as.data.frame()
head(otu)

#细菌真菌分组map文件要完全一致
map = ps03 %>% sample_data()
head(map)

# 定义后续参数-----
map= sample_data(ps03)
head(map)
# 提取分组因子数量
gnum = phyloseq::sample_data(ps03)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.16s对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps03)$Group %>%unique();axis_order

#-主题--颜色等
package.amp()
res = theme_my(ps03)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

# alpha-diversity-----
#1 alpha.omics: 6中alpha多样性计算 ----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
tab = alpha.omics(ps = ps03,group = "Group")
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

result = EasyStat::MuiKwWlx2(data = data,num =3:6)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:6,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+mytheme1

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(6)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+mytheme1
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(6)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+mytheme1

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0+mytheme1

# beta diversity-----
#2 ordinate.omics: 排序分析#----------
result = ordinate.omics(ps = ps03, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                        pair = F)
p3_1 = result[[1]]
p3_1+mytheme1

#带标签图形出图
p3_2 = result[[3]]
p3_2+mytheme1

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
p3_3+mytheme1

#3 ordinateTest.omics:群落水平差异检测#-------
dat1 = ordinateTest.omics(ps = ps03, Micromet = "adonis", dist = "bray")
dat1

#4 pairordinateTest.omics:两两分组群落水平差异检测#-------
dat2 = pairordinateTest.omics(ps = ps03, Micromet = "MRPP", dist = "bray")
dat2

#5 mantal.omics ：群落差异检测普鲁士分析#------
dat = mantal.omics(ps01= ps01,ps02= ps02)

#6 cluster.omics:样品聚类#-----
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

# constraint analysis----
ps.tem = ps.16s %>% tax_glom_wt("Genus")
# 真菌筛选
res = loadingPCA.micro(ps = ps.its%>% tax_glom_wt("Genus"))
# p = res[[1]]
# p
dat = res[[2]]
dat$id = row.names(dat)
id = dat %>% arrange(desc(PCone)) %>% head(15) %>% .$id
id

ftab = ps.its%>% tax_glom_wt("Genus") %>% scale_micro() %>%
  phyloseq::subset_taxa(
    row.names(phyloseq::tax_table(ps.its%>% tax_glom_wt("Genus"))) %in% c(id)
  ) %>% vegan_otu() %>%
  as.data.frame()

head(ftab)
#7 RDA_CCA:细菌真菌 共排序-----
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

#8 RDA_CCA_explain_percent: 对微生物群落的解释比例-----
result = RDA_CCA_explain_percent(ps = ps.tem, env.dat = ftab)

sam = result[[1]]
sam

all = result[[2]]

p= result[[3]]
p
#9 rdacca.hp.micro:  层次分割对细菌群落影响的真菌---------
library(rdacca.hp)
library(vegan)
library(ggClusterNet)
library(tidyverse)
res = rdacca.hp.micro(OTU = ps.tem %>% filter_OTU_ps(500) %>%vegan_otu() %>% as.data.frame(), env = ftab[,1:5], cca = FALSE)
p3_rda = res[[1]]
p3_rda
dat1 = res[[2]]
dat1
p3_db_rda  = res[[3]]
p3_db_rda
dat2 = res[[4]]
dat2

#10 rdacca.hp.micro.p: 层次分割对微生物影响的代谢物#---------
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
p3_db_rda  = res[[3]]
p3_db_rda
dat2 = res[[4]]
dat2

# 网络分析---------
#11 corBionetwork.st-----
pst = ps03 %>% tax_glom_wt("Genus")
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
  minsize = 4,maxsize = 14)

dat = res[[2]]
cor = dat$cortab
p = res[[1]]
p

#12 net_properties.4-------
id = names(cor)
id

for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  dat = net_properties.4(igraph,n.hub = FALSE)
  head(dat,n = 16)
  colnames(dat) = id[i]

  if (i == 1) {
    dat2 = dat
  } else{
    dat2 = cbind(dat2,dat)
  }
}
head(dat2)

#13 netproperties.sample-----

id = pst %>% sample_data() %>% .$Group %>% unique()
id
id0 = names(cor)
for (i in 1:length(id)) {
  pst1 = pst %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(pst = pst1,cor = cor[[id0[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}


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

#14 node_properties-------

for (i in 1:length(id0)) {

  igraph= cor[[id0[i]]] %>% make_igraph()
  nodepro = node_properties(igraph) %>% as.data.frame()
  nodepro$Group = id[i]
  head(nodepro)
  colnames(nodepro) = paste0(colnames(nodepro),".",id[i])
  nodepro = nodepro %>%
    as.data.frame() %>%
    rownames_to_column("ASV.name")

  # head(dat.f)
  if (i == 1) {
    nodepro2 = nodepro
  } else{

    nodepro2 = nodepro2 %>% full_join(nodepro,by = "ASV.name")

  }
}
head(nodepro2)

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

#rfcv.omics :交叉验证结果-------
ps = ps03 %>% subset_taxa.wt("Family","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Order","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Genus","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Phylum","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("class","Unassigned",TRUE)
ps = ps %>% tax_glom_wt("Genus")

result =rfcv.omics(ps = ps ,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv

rfcvtable = result[[3]]
rfcvtable

#8 Roc.omics:ROC 曲线绘制----
id = sample_data(ps03)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)
sample_data(ps03)

pst = ps03 %>% subset_samples.wt("Group",group)

res = Roc.omics( ps = pst %>% filter_OTU_ps(100),repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

#9 loadingPCA.omics:载荷矩阵筛选特征功能------
res = loadingPCA.omics(ps = ps03,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#10 randomforest.omics: 随机森林筛选特征功能----
res = randomforest.omics( ps = ps03 %>% filter_OTU_ps(100),
                          group  = "Group",
                          optimal = 50)
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

#11 svm: 筛选特征细菌与真菌------
res <- svm.omics(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#12 GLM: 筛选特征细菌与真菌---------
res <- glm.omics(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#13 xgboost: 筛选特征细菌与真菌-----
res =xgboost.omics(ps = ps03, top = 5,k =5  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#14 kNN: 筛选特征细菌与真菌--------
res = knn.omics(ps =ps03, top = 20,k =5 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# #15 decisiontree: 筛选特征细菌与真菌------
res =decisiontree.omics(ps=ps03, top =500,  k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#16 bagging: 筛选特征细菌与真菌-----
res =bagging.omics(ps =  ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#17 lasso: 筛选特征细菌与真菌-----
res = lasso.omics (ps =  ps03, top = 20, seed = 10, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#18 naivebayes: 筛选特征细菌与真菌------
res = naivebayes.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# #19 nnet神经网络: 筛选特征细菌与真菌------
res =nnet.omics(ps=ps03, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance
