rm(list=ls())
# 转录组分析流程------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(EasyMultiOmics)
library(EasyMultiOmics.db)

ps.trans = EasyMultiOmics::ps.trans
ps.trans = ps.trans   %>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
otu_table(ps.trans )  =  round( otu_table(ps.trans ))
#--提取有多少个分组
gnum = phyloseq::sample_data(ps.trans)$Group %>% unique() %>% length()
gnum

# 设定排序顺序
axis_order = phyloseq::sample_data(ps.trans)$Group %>% unique()


#-主题--颜色等
res = theme_my(ps.trans)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

# data preprocess------
#1 按照通路（KEGG）合并基因#------
detach("package:mia", unload = TRUE)
ps.trans2 =trans_kegg ( ps.trans)
phyloseq::tax_table(ps.trans2) %>% head()
ps1

#2 基于 Mkegg 通路合并#-------
ps.trans2 =trans_mkegg ( ps.trans)

#3 基于reaction 层面合并#----
ps.trans3 =trans_rekegg ( ps.trans)
ps.trans3

# ordinate analysis------
#4 ordinate.trans: 排序分析#----------
result = ordinate.trans(ps =ps.trans, group = "Group", dist = "bray",
                       method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                       pair = F)
p3_1 = result[[1]]
p3_1+mytheme1

#带标签图mytheme1#带标签图形出图
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

#5 transTest.metf:群落功能差异检测#-------
dat1 = transtest.trans(ps =ps.trans, Micromet = "adonis", dist = "bray")
dat1

#6 pairtranstest.trans:两两分组功能差异检测#-------
dat2 = pairtranstest.trans(ps =ps.trans, Micromet = "MRPP", dist = "bray")
dat2

#7 mantal.trans：功能差异检测普鲁士分析#------
result <- mantal.trans(ps = ps.trans,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)
data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7+mytheme1

#8 cluster.trans:样品聚类#-----
res = cluster.trans(ps= ps.trans,
                   hcluter_method = "complete",
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

#  gene classification-------
#9 Ven.Upset.trans: 用于展示共有、特有的功能----
# 分组小于6时使用
res = Ven.Upset.trans(ps =   ps.trans,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

#10 ggflower.trans:花瓣图展示共有特有功能------
res <- ggflower.trans (ps =  ps.trans ,
                      # rep = 1,
                      group = "ID",
                      start = 5, # 风车效果
                      m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                      a = 0.2, # 花瓣胖瘦
                      b = 1, # 花瓣距离花心的距离
                      lab.leaf = 3, # 花瓣标签到圆心的距离
                      col.cir = "#FF0000",
                      N = 0.5 )

p13.1 = res[[1]]

p13.1

dat = res[[2]]
dat

#11 ven.network.trans: ven网络展示共有特有功能----
rank.names( ps.trans)
result =ven.network.trans(
  ps =  ps.trans,
  N = 0.5,
  fill = "KO_id")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat

#12 barMainplot.trans: 堆积柱状图展示组成----
phyloseq::tax_table(ps.trans) %>% colnames()
result = barMainplot.trans(ps = ps.trans,
                          j = "KO_id",
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)

p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","level_2","Abundance(%)")
head(databar)

#13 cluMicro.bar.trans: 聚类堆积柱状图展示组成 -----
phyloseq::tax_table(ps.trans2) %>% colnames()
result <-  cluMicro.bar.trans(dist = "bray",
                               ps =  ps.trans2,
                               j = "MDESCRPTION",
                               Top = 10, # 提取丰度前十的物种注释
                               tran = TRUE, # 转化为相对丰度值
                               hcluter_method = "complete",
                               Group = "Group",
                               cuttree = length(unique(phyloseq::sample_data(ps.trans2)$Group)))

result[[1]]
p5_2 <- result[[2]]
p5_2
p5_3 <- result[[3]]
p5_3
p5_4 <- result[[4]]
p5_4
clubardata <- result[[5]]
clubardata

# gene difference analysis
#14 aldex2.trans-------
library(ALDEx2)
?aldex2.trans
dat = aldex2.trans(ps = ps.trans %>% filter_OTU_ps(50),
                   group = "Group",
                   test = "t")
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#15 ancom.trans:-------------
dat = ancom.trans(ps =  ps.trans %>% filter_OTU_ps(50),
                  group = "Group",
                  p_adj_method = "BH",
                  ANCOM.value = "detected_0.6",
                  alpha=0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#16 corncob.trans:----------
library(corncob)
dat = corncob.trans(
  ps =  ps.trans %>% filter_OTU_ps(50),
  group =  "Group",
  alpha = 0.5)

dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#17 lefse.trans: ----------
dat = lefse.trans(ps = ps.trans %>% filter_OTU_ps(50),group =  "Group",alpha = 0.05)

#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)


#18 limma.v.TMM.trans:-------
library(limma)
dat = limma.v.TMM.trans (ps =  ps.trans %>% filter_OTU_ps(50),
                         group =  "Group",alpha = 0.05,method = "TMMwsp")

#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#19 maaslin2.trans: 输出文件----------
library(Maaslin2)
dat = maaslin2.trans(ps =  ps.trans %>% filter_OTU_ps(50),group =  "Group",alpha = 0.05,
                     rare = FALSE)

#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#20 metaSeq.trans:---------
library(metagenomeSeq)
dat = metaSeq.trans(ps =  ps.trans %>% filter_OTU_ps(50),group = "Group",alpha = 0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#21 wilcox.sampl.trans:-------
dat = wilcox.sampl.trans(ps =  ps.trans %>% filter_OTU_ps(50),
                         group =  "Group",
                         alpha = 0.05)

#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#22 wilcox.clr.trans:----------
# 未经稀释情况下使用CLR转换的丰度(应用伪计数1后)
dat = wilcox.clr.trans(ps =  ps.trans %>% filter_OTU_ps(50),
                       group =  "Group",
                       alpha = 0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#23 maaslin2.trans:----------
library(Maaslin2)
dat = maaslin2.trans(ps =  ps.trans %>% filter_OTU_ps(50),group =  "Group",alpha = 0.05,
                     rare = TRUE)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#24 ancombc.trans:-------
library("ANCOMBC")
library(DT)
dat = ancombc.trans(ps =  ps.trans %>% filter_OTU_ps(50),
                    group = "Group",
                    alpha = 0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#25 ancombc2.trans-------
#速度比较慢
dat = ancombc2.trans(ps =  ps.trans %>% filter_OTU_ps(20),
                     group =  "Group",
                     alpha = 0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#26 zicoseq.trans-------
library(GUniFrac)
dat = zicoseq.trans(ps =  ps.trans %>% filter_OTU_ps(1000) ,
                    group =  "Group",
                    alpha = 0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#27 linda.trans------
library(mia)
library(MicrobiomeStat)
dat = linda.trans(ps =  ps.trans %>% filter_OTU_ps(50),
                  group =  "Group",
                  alpha = 0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#28 dacomp.trans------
library(dacomp)
dat = dacomp.trans (ps =  ps.trans %>% filter_OTU_ps(50) ,
                    sd = 0.6,
                    q_BH =  0.1,
                    q_DSFDR =0.1)

#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

#29 EdgerSuper.trans:EdgeR计算差异功能基因----
res =EdgerSuper.trans (ps =  ps.trans %>% filter_OTU_ps(1000),
                       group  = "Group",
                       artGroup = NULL)
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][2]
p16.2
p16.3 = res[[1]][3]
p16.3
dat = res[[2]]
dat


#30 DESep2Super.trans:DESep2计算差异功能基因 ----
res = DESep2Super.trans (ps = ps.trans %>% filter_OTU_ps(1000),
                        group  = "Group",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][2]
p15.2
p15.3 = res[[1]][3]
p15.3
dat = res[[2]]
dat


#31 t.trans: 差异分析t检验----
res = t.trans(ps = ps.trans%>%filter_OTU_ps(50),group  = "Group",artGroup =NULL,pvalue=0.05,FCvalue=2,padjust=TRUE,padjust_method="BH",
              log2FC_threshold=10,ncol=3,nrow=1,outpath=NULL)
p17=res[[1]]
p17
data=res[[2]]
head(data)
inter_union=res[[3]]
head(inter_union)
p17.1 = res [[4]][1]
p17.1
p17.2 = res [[4]][2]
p17.2
p17.3 = res [[4]][3]
p17.3

#32 wlx.trans:非参数检验-----
res = wlx.trans(ps = ps.trans%>%filter_OTU_ps(50),
                group  = "Group",
                artGroup =NULL,
                pvalue=0.05,FCvalue=2,
                padjust=TRUE,
                padjust_method="BH",
                log2FC_threshold=10,
                ncol=3,nrow=1,
                outpath=NULL)
p18=res[[1]]
p18
data=res[[2]]
head(data)
inter_union=res[[3]]
head(inter_union)
p18.1 = res [[4]][1]
p18.1
p18.2 = res [[4]][2]
p18.2
p18.3 = res [[4]][3]
p18.3

#33 Mui.Group.volcano.trans: 聚类火山图------
res =  EdgerSuper2.trans(ps = ps.trans,group  = "Group",artGroup =NULL,
                         j = "gene")
res2 = Mui.Group.volcano.trans(res = res)

p29.1 = res2[[1]]
p29.1
p29.2 = res2[[2]]
p29.2
dat = res2[[3]]
dat


# biomarker identification----
id = sample_data( ps.trans)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])


pst = ps.trans %>% subset_samples.wt("Group",group) %>%
      filter_taxa(function(x) sum(x ) > 10, TRUE)
#33 rfcv.trans :交叉验证结果-------
library(ROCR)
result =rfcv.trans(ps = ps.trans %>%filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#34 randomforest.trans: 随机森林筛选特征功能----
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
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
res = randomforest.trans( ps = ps.trans %>% filter_OTU_ps(500),group  = "Group",
                          optimal = 50,
                          fill = "KO_id")

p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

#35 loadingPCA.trans:载荷矩阵筛选特征功能------
res = loadingPCA.trans(ps =  ps.trans,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#36 svm.trans------
res <- svm.trans(ps = ps.trans %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#37 glm.trans---------
res <- glm.trans(ps = ps.trans %>% filter_OTU_ps(500), k = 5)
accuracy= res[[1]]
accuracy
importance = res[[2]]
importance

#38 xgboost.trans-----
res = xgboost.trans(ps = ps.trans, seed = 200, top = 200 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#39 decisiontree ------
res = decisiontree.trans(ps=ps.trans, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#40 bagging.trans-----
res = bagging.trans(ps =  ps.trans , top = 20, seed = 110, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#41 lasso.trans-----
res = lasso.trans (ps =  pst, top = 500, seed = 1, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#42 LDA.trans--------
res= LDA.trans(ps=ps.trans,group = "Group", Top = 100)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#43 naivebayes.trans------
res = naivebayes.trans(ps=pst, top = 200, seed = 100, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#44 nnet.trans 神经网络------
res = nnet.trans  (ps=ps.trans, top = 200, seed = 10, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# enrich analysis --------
#45 KEGG_enrich.trans-------
library(clusterProfiler)
res = KEGG_enrich.trans(ps = ps.trans %>% filter_OTU_ps(1000))
res$plots[[1]]$`KO-OE`

#46 gsva.trans  有错----
library(GSVA)
res = gsva.trans(ps = ps.trans2 %>% filter_OTU_ps(1000),
                 FC = FALSE,
                 adj = FALSE,lg2FC = 0.3,
                 padj = 0.05)
#总的热图
res$plots[[1]][[1]]
#两两组合
res$plots[[2]]$`KO OE`

#47 gsea.trans  缺-----


#48 pathway_enrich.trans：通路富集-------
phyloseq::tax_table(ps.trans3) %>% colnames()
ps.trans3 = ps.trans %>% tax_glom_wt("KO_id")

res2 = pathway_enrich.trans(ps = ps.trans3,
                         dif.method = "wilcox")

res2$plots$OE.WT.plot
res2$plots$KO.OE.plot

#49 reaction.show.trans：反应展示-------
res3 = reaction.show.trans(ps= ps.trans,dif.method = "wilcox")

res3$plots$OE.WT.plot

res3$plotdata$OE.WT

#50 enrich.module.trans  缺-----

#51 enrich.reaction.trans  缺-----


Go_enrich

#network analysis------

#52 network.pip:网络分析主#--------

library(igraph)
library(sna)
tab.r = network.pip(
  ps = ps.trans,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = TRUE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label =  FALSE,
  lab = "elements",
  group = "Group",
  fill = "Super_class",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1
)

dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]
# 提取网络图可视化结果
p0 = plot[[1]]
p0

#53 net_properties.4:网络属性计算#---------
i = 1
cor = cortab
id = names(cor)
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  dat = net_properties.4(igraph,n.hub = F)
  head(dat,n = 16)
  colnames(dat) = id[i]

  if (i == 1) {
    dat2 = dat
  } else{
    dat2 = cbind(dat2,dat)
  }
}
head(dat2)

#54 netproperties.sample:单个样本的网络属性#-------
for (i in 1:length(id)) {
  ps.mst = ps.trans %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(ps.mst = ps.mst,cor = cor[[id[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}

map= sample_data(ps.trans)
map$ID = row.names(map)
map = map %>% as.tibble()
dat3 = dat.f2 %>% rownames_to_column("ID") %>% inner_join(map,by = "ID")

head(dat3)

#55 node_properties:计算节点属性#---------
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
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

#56 module.compare.net.pip:网络显著性比较#-----
dat = module.compare.net.pip(
  ps.ms = NULL,
  corg = cor,
  degree = TRUE,
  zipi = FALSE,
  r.threshold= 0.8,
  p.threshold=0.05,
  method = "spearman",
  padj = F,
  n = 3)
res = dat[[1]]
head(res)

#57 CNPS.network: CNPS 基因功能网络 -----
library(data.table)
library(EasyStat)
library(tidyfst)
library(sna)
dat = db.cnps
head(dat)
library(sna)
res = CNPS.network(ps = ps.trans,dat = dat,id.0 = "C")
p= res[[1]]
p

dat1= res[[2]]
dat1

dat2= res[[3]]
dat2


wgcna.hub.trans

hub.gene.trans
