# EasyMultiOmics


## EasyMultiOmics: An R package for comprehensive multi-omics data mining in systematic microbiome research

TaoWen Penghao Xie JunYuan et al



### Abstract

Advancements in high-throughput sequencing technology have significantly increased omics data, necessitating a shift from single-omics to multi-omics joint analysis to address systematic microbiome research requirements. However, tools that can comprehensively and efficiently perform multi-omics analysis are largely deficient. Here, we design and develop an R package EasyMultiOmics to solve this short slab. Built on the phyloseq class, EasyMultiOmics is a comprehensive tool to integrate statistical analysis, multi-omics analysis, and visualization . EasyMultiOmics contains over 300 functions and provides 14 workflows for single-omics data mining, including microbiome, metabolome, transcriptome, and metagenome analyses. For multi-omics joint analysis, it offers seven workflows that enable combined analysis of bacteria and fungi, composition and function of metagenome, microbiome and metabolome, transcriptome and metabolome, transcriptome and metagenome, as well as integrated analysis of four major omics (transcriptome, metabolome, metagenome, and microbiome). Compared to other multi-omics related tools, EasyMultiOmics is a fast, flexible, and modular, providing powerful and convenient tool for researchers.


![1718185672438.jpg](![1718629217418](https://github.com/taowenmicro/EasyMultiOmics/assets/48119869/afd855b8-240d-41c3-a2c5-7017b5599a0b)


## pipeline for microbiome data

```

# 基于整理好的函数，重新写流程
# rm(list=ls())


library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)

data("ps16s")

# 定义后续参数
map= sample_data(ps)
head(map)
# 提取分组因子数量
gnum = phyloseq::sample_data(ps)$Group %>% unique() %>% length()
gnum
# 扩增子微生物组学分析#-----
#1 alpha.micro: 6中alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.micro(ps = ps,group = "Group",Plot = TRUE )
head(tab)


data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)
# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = EasyStat::MuiKwWlx2(data = data,num = 3:9)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:9,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(9)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(9)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0



#2 alpha.pd :用于计算pd多样性#-------
library(ape)
library(picante)
tab2 = alpha.pd(ps)
head(tab2)
result = EasyStat::MuiKwWlx2(data = tab2,num = 3)
result1 = EasyStat::FacetMuiPlotresultBox(data = tab2,num = 3,
                                          result = result,
                                          sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]]
p1_1



#3 alpha.rare.line:alpha多样性稀释曲线#---------
rare <- mean(phyloseq::sample_sums(ps))/10
result = alpha.rare.line(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)
#--提供单个样本溪稀释曲线的绘制
p2_1 <- result[[1]]
p2_1
## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]]
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]]
p2_3



#4 ordinate.micro: 排序分析#----------
result = ordinate.micro(ps = ps, group = "Group", dist = "bray",
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


#5 MicroTest:群落水平差异检测#-------
dat1 = MicroTest.micro(ps = ps, Micromet = "adonis", dist = "bray")
dat1

#6 pairMicroTest:两两分组群落水平差异检测#-------
dat2 = pairMicroTest.micro(ps = ps, Micromet = "MRPP", dist = "bray")
dat2


#7 mantal.micro ：群落差异检测普鲁士分析#------
map= sample_data(ps)
head(map)

result <- mantal.micro(ps = ps,
                       method =  "spearman",
                       group = "Group",
                       ncol = gnum,
                       nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7


#8 cluster.micro:样品聚类#-----
cluster.micro
sample_data(ps)
res = cluster.micro (ps= ps,hcluter_method = "complete",
                     dist = "bray",
                     cuttree = 3,
                     row_cluster = T,
                     col_cluster =  T)

p4 = res[[1]]
p4
p4_1 = res[[2]]
p4_1
p4_2 = res[[3]]
p4_2
dat = res[[4]]# 聚类矩阵
head(dat)


#9 distance.micro:分组之间距离比对#-----
res = distance.micro(ps,group = "Group")
p5.1 = res[[1]]
p5.1
p5.2 = res[[2]]
p5.2
p5.3 = res[[3]]
p5.3
dat = res[[4]]
head(dat)

#10 Ven.Upset.micro: 用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = Ven.Upset.micro(ps =  ps,
                      group = "Group",
                      N = 0.5,
                      size = 3)
p10.1 = res[[1]]
p10.1
p10.2 = res[[2]]
p10.2
p10.3 = res[[4]]
p10.3
dat = res[[3]]
head(dat)

#11 ggVen.Upset.micro:用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = ggVen.Upset.micro(ps = ps,group = "Group")
p11.1 = grid::grid.draw(res[[1]])

dat = res[[2]]



#12 VenSeper.micro: 错误 详细展示每一组中OTU的物种组成 -------
#---每个部分
library("ggpubr")
library(agricolae)
library(reshape2)

result = VenSuper.micro(ps = ps,
                        group = group,
                        num = 6
                         )
# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]
p7_1
# 每部分的otu门类冲积图
p7_2 <- result[[3]]
p7_2
#每个部分序列的数量占比，并作差异
p8 <- result[[2]]
p8




#13 ggflower.micro:花瓣图展示共有特有微生物------
res <- ggflower.micro (ps = ps ,
                       # rep = 1,
                       group = "ID",
                       start = 1, # 风车效果
                       m1 = 2, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                       a = 0.3, # 花瓣胖瘦
                       b = 1, # 花瓣距离花心的距离
                       lab.leaf = 1, # 花瓣标签到圆心的距离
                       col.cir = "yellow",
                       N = 0.5 )

p13.1 = res[[1]]
p13.1

dat = res[[2]]
head(dat)



#14ven.network.micro:venn 网络----
result = ven.network.micro(
  ps = ps,
  N = 0.5,
  fill = "Phylum")
p14  = result[[1]]
p14

dat = result[[2]]
head(dat)


#15 Micro_tern.micro: 三元图展示组成----
ps1 = ps %>% filter_OTU_ps(500)
res = Micro_tern.micro(ps1)
p15 = res[[1]]
p15[[1]] +theme_bw()


dat =  res[[2]]
head(dat)

res$groups

#16 barMainplot.micro: 堆积柱状图展示组成----
library(ggalluvial)
pst = ps %>% subset_taxa.wt("Species","Unassigned",TRUE)
pst = pst %>% subset_taxa.wt("Genus","Unassigned",TRUE)
result = barMainplot.micro(ps = pst,
                           j = "Genus",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           Top =5)
p4_1 <- result[[1]]
# scale_fill_manual(values = colset2) +
# scale_x_discrete(limits = axis_order) +
# mytheme1
p4_1

p4_2  <- result[[3]]
# # scale_fill_brewer(palette = "Paired") +
# scale_fill_manual(values = colset2) +
# scale_x_discrete(limits = axis_order) +
# mytheme1
p4_2

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group",j,"Abundance(%)")
head(databar)




#17 cluMicro.bar.micro: 聚类堆积柱状图展示组成 -----
result <-  cluMicro.bar.micro (dist = "bray",
                               ps = ps16s,
                               j = "Genus",
                               Top = 10, # 提取丰度前十的物种注释
                               tran = TRUE, # 转化为相对丰度值
                               hcluter_method = "complete",
                               Group = "Group",
                               cuttree = length(unique(phyloseq::sample_data(ps)$Group)))
result[[1]]
p5_2 <- result[[2]]
p5_2
p5_3 <- result[[3]]
p5_4 <- result[[4]]
clubardata <- result[[5]]


#18 cir_barplot.micro:环状堆积柱状图 -----

library(ggtree) # j = "Phylum"
res = cir_barplot.micro(
  ps = ps,
  Top = 15,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete")

p17 = res[[1]]
p17
dat= res[[2]]
head(dat)

#19 cir_plot.micro:和弦图展示物种组成-----
res = cir_plot.micro(ps  = ps,Top = 12,rank = 6)




#20 micro_maptree.micro；圈图展示物种组成-----

library(ggraph)
library(data.tree)
library(igraph)

tax = ps %>% vegan_tax() %>%
  as.data.frame()
head(tax)
ps.bac <- ps %>% subset_taxa.wt("Kingdom", "Bacteria")
ps.bac

res = micro.maptree (ps = ps.bac,
                     Top = 100,
                     labtab =  NULL,
                     seed = 11)
p20 = res[[1]]
p20
dat = res[[2]]
head(dat)

#21 phy_tree.micro；树状图展示物种组成以及相关性-----
library(ggstar)
result <- phy_tree.micro(ps = ps,Top = 100)

p0 = result[[1]]
p0
p1 = result[[2]]
p2 = result[[3]]
p3 = result[[4]]
p4 = result[[5]]
p5 = result[[6]]
p6 = result[[7]]
p7 = result[[8]]
p7
dat = result[[9]]

#22 sankey.micro: 桑基图展示物种组成------

res = sankey.micro (ps = ps,
                    rank = 6,
                    Top = 50)

p22 = res[[1]]
dat = res[[2]]
saveNetwork(p22,paste("./sankey_Group.html", sep = ""))


#23 sankey.m.Group.micro: 按照分组绘制  -----
result = sankey.m.Group.micro(
  ps = ps  %>% subset_taxa.wt("Species","Unassigned",TRUE),
  rank = 6,
  Top = 50

)

p = result[[1]]
p

saveNetwork(p,paste("./sankey_Group.html", sep = ""))

dat = result[[2]]
dat


#24 Microheatmap.micro: 热图展示物种相对丰度差异-----

heatnum = 30
map = phyloseq::sample_data(ps)
map$ID = row.names(map)
phyloseq::sample_data(ps) = map

ps_tem = ps %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus")
id <- ps %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus") %>%
  ggClusterNet::filter_OTU_ps(100) %>%
  ggClusterNet::vegan_otu() %>%
  t() %>% as.data.frame() %>%rowCV %>%
  sort(decreasing = TRUE) %>%
  head(heatnum) %>%
  names()

result <- Microheatmap.micro(ps_rela = ps_tem,id = id ,col_cluster = FALSE,row_cluster = FALSE)

p24.1 <- result[[1]]
p24.1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p24.2 <- result[[2]]
p24.2
dat = result[[3]]
head(dat)

#25 EdgerSuper.micro:EdgeR计算差异微生物----
res = EdgerSuper.micro(ps = ps %>% ggClusterNet::filter_OTU_ps(500),group  = "Group",artGroup = NULL, j = "OTU")


p25.1 =  res[[1]][1]
p25.1
p25.2 =  res[[1]][2]
p25.2
p25.3 =  res[[1]][3]
p25.3
dat =  res[[2]]
head(dat)

#25.1 EdgerSuper2.micro:EdgeR计算差异微生物-----
res =  EdgerSuper2.micro (ps = ps,group  = "Group",artGroup =NULL, j = "OTU")
head(res)



#26 DESep2Super.micro:DESep2计算差异微生物-----
res = DESep2Super.micro(ps = ps %>% ggClusterNet::filter_OTU_ps(500),
                        group  = "Group",
                        artGroup =NULL,
                        j = "OTU"
                        # path = diffpath.1
)


p26.1 =  res[[1]][1]
p26.1
p26.2 =  res[[1]][2]
p26.2
p26.3 =  res[[1]][3]
p26.3
dat =  res[[2]]
dat

#27 edge_Manhattan.micro: 曼哈顿图展示差异微生物------

res = edge_Manhattan.micro(
  ps = ps%>% ggClusterNet::filter_OTU_ps(500),
  pvalue = 0.05,
  lfc = 0
)
p27.1= res[[1]]
p27.1
p27.2= res[[2]]
p27.2
p27.3= res[[3]]
p27.3

#28 stamp_diff.micro: stamp展示差异微生物----
map = phyloseq::sample_data(ps)
map$ID = row.names(map)

# map$Group = as.factor(map$Group)
sample_data(ps) = map

allgroup <- combn(unique(map$Group),2)
plot_list <- list()
for (i in 1:dim(allgroup)[2]) {
  ps_sub <- phyloseq::subset_samples(ps,Group %in% allgroup[,i]);ps_sub
  p <- stemp_diff.micro(ps = ps_sub,Top = 20,ranks = 6)
  plot_list[[i]] <- p

}

p28.1 = plot_list[[1]]
p28.1
p28.2 = plot_list[[2]]
p28.2
p28.3 = plot_list[[3]]
p28.3


#29 Mui.Group.volcano.micro: 聚类火山图------

res =  EdgerSuper2.micro (ps = ps,group  = "Group",artGroup =NULL, j = "OTU")

res2 = Mui.Group.volcano.micro(res = res)

p29.1 = res2[[1]]
p29.1
p29.2 = res2[[2]]
p29.2
dat = res2[[3]]
dat


#30 Mui.cluster.volcano.micro: 指定分组绘制聚类火山图------
library(ggrepel)
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
# 设置分组
group2 = c(aaa[1,1],aaa[2,1])
b= data.frame(group2)

res =  EdgerSuper2.micro (ps = ps,group  = "Group",artGroup =b, j = "OTU")

res2 = Mui.cluster.volcano.micro (res = res,rs.k = 6)
p30.1 = res2[[1]]+ggtitle(group2)
p30.1
p30.2 = res2[[2]]+ggtitle(group2)
p30.2
dat= res2[[3]]
dat




#31 Plot.CompareWithCK.micro: 相对对比柱状图-----
library(parallel)
rm(list=ls())

res = EdgerSuper.micro(ps = ps,group  = "Group",artGroup = NULL,
                       j = "Genus")
result = res [[2]]
res = Plot.CompareWithCK.micro(ps = ps,
                               CK = "WT",
                               j = "Genus",
                               abun = 0.001,
                               cpu = 6)

p = res[[1]]
p
data = res[[2]]
data


#32 rfcv.micro :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.Micro(ps = ps %>% filter_OTU_ps(500),
                   group  = "Group",optimal = 20,nrfcvnum = 6)

prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable




#33 Roc.micro:ROC 曲线绘制----
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)

ps = ps %>% subset_taxa.wt("Family","Unassigned",T)
ps = ps %>% subset_taxa.wt("Order","Unassigned",T)
ps = ps %>% subset_taxa.wt("Genus","Unassigned",T)
ps = ps %>% subset_taxa.wt("Phylum","Unassigned",T)
ps = ps %>% subset_taxa.wt("class","Unassigned",T)

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)
res = Roc.micro( ps = pst %>% filter_OTU_ps(1000),group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat


#34 loadingPCA.micro: 载荷矩阵筛选特征微生物------
res = loadingPCA.micro(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat



#35 LDA.micro: LDA筛选特征微生物-----

p1 <- p_base.micro(ps,Top = 100)
p1$data

tablda = LDA.micro(ps = ps,
                   Top = 100,
                   p.lvl = 0.05,
                   lda.lvl = 1,
                   seed = 11,
                   adjust.p = F)
tablda[[1]]

p35 <- lefse_bar(taxtree = tablda[[2]])
p35
dat = tablda[[2]]
dat


#36 svm.micro: 空----


#37 glm.micro: 空----


#38 xgboost.micro: 空----


#39 lr.micro: 空----


#40 decisiontree.micro: 空----


#41 naivebayes.micro: 空----


#42 randomforest.micro: 筛选特征微生物----
res = randomforest.micro( ps = ps,group  = "Group", optimal = 50)
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4


#43 neutralModel: 中性模型----
library(picante)
library(ape)
# library(vegan)
library(FSA)
library(eulerr)
# library(grid)
# library(gridExtra)
require(minpack.lm)
require(Hmisc)
require(stats4)
# library(parallel)
psphy = filter_taxa(ps, function(x) sum(x ) > 100, TRUE);psphy
map = sample_data(psphy)
n = map$Group %>% unique() %>%
  length()
n

result = neutralModel(ps = psphy,group  = "Group",ncol = n)
#--合并图表
p43 =  result[[1]]
p43
dat = result[[3]]
dat
dat2 = result[[4]]
dat2


# #44 phyloSignal: 系统发育信号 空----
# data("env")
# env
#
# phyloSignal(ps = psphy,
#             group  = "Group",
#             env = envRDA[,1:2],
#             path = phypath2)
#
# result = phySigPlot(ps = ps,group  = "Group",env = env[,1:2],path = phypath2)
# #
# #提取图片
# p2 = result[[1]] + mytheme1
# p2
# #-提取作图数据
# data = result[[2]]
# head(data)


#45 phySigPlot: 空 ----

#46 nullModel: 零模型----
result <- nullModel(ps = psphy,
                    group="Group",
                    dist.method =  "bray",
                    gamma.method = "total",
                    transfer = "none",
                    null.model = "ecosphere"
)

#--分组零模型运行结果
nullModeltab <- result[[1]]

# 比例
ratiotab <- result[[2]]
#-统计量统计差异
aovtab <- result[[3]]


#47 bNTICul: bNTI----
result = bNTICul(ps = psphy,
                 group  = "Group",
                 num = 10,
                 thread = 1
)
bNTI = result[[1]]
head(bNTI)


#47 RCbary: ----
result = RCbary(ps = psphy ,group  = "Group",num = 10,thread = 1)

RCbary = result[[1]]
head(RCbary)

#48 bNTIRCPlot: BetaNTI和RCbray联合出图--------
RCb = RCbary %>%
  dplyr::mutate(Sample_1 = Site2, Sample_2 = Site1)

result = bNTIRCPlot(ps = psphy ,RCb  =Rcb, bNTI = bNTI,group  = "Group")

#--bNTI出图片
p3 <- result[[1]]
p3

#RCbary可视化
p4 <- result[[2]]
p4

#组合图片BNTI，RCbray
p5 <- result[[3]]
p5
plotdata = result[[4]]
head(plotdata)

dat = result[[5]]
head(dat)








#49 buplotall.micro: 筛选特征微生物----
#50 KEGG_enrich.micro: 筛选特征微生物----
#51 FEAST.micro 错: 筛选特征微生物----


#43 randomforest.micro: 筛选特征微生物----
#43 randomforest.micro: 筛选特征微生物----
#43 randomforest.micro: 筛选特征微生物----

```
