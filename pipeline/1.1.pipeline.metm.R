# 宏基因组物种分析#-----
rm(list=ls())

# BiocManager::install("MicrobiotaProcess")

library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)


# 有效基因数量
sample_sums(ps.micro)



## 参数设定-----
map= sample_data(ps.micro)
head(map)

phyloseq::tax_table(ps.micro) %>% head()
# 提取分组因子数量
gnum = phyloseq::sample_data(ps.micro)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.metm对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps.micro)$Group %>%unique();axis_order
col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")

#-主题--
package.amp()

res = theme_my(ps.micro)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

# alpha diversity -----
#1 alpha.metm: 6中alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metm(ps = ps.micro,group = "Group" )
head(tab)


data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)
# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)






result = EasySigR::MuiKwWlx2(data = data,num = 3:6)

result1 = EasySigR::FacetMuiPlotresultBox(data = data,num = 3:6,
                                          result = result,
                                          sig_show ="abc",ncol = 4,width = 0.4 )

p1_1 = result1[[1]] +scale_fill_manual(values = col.g)

p1_1+
  ggplot2::scale_x_discrete(limits = axis_order) +
  # theme_cell()+
  theme_nature()+

  ggplot2::guides(fill = guide_legend(title = none))



res = EasySigR::FacetMuiPlotresultBar(data = data,num = c(3:6),
                                      result = result,sig_show ="abc",ncol = 4,
                            mult.y = 0.3

                            )
p1_2 = res[[1]]
p1_2 +
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme_nature +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)

p = p1_2 +
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme_cell +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)

# #  推荐保存参数
# ggsave(
#   filename = "Figure1.pdf",
#   plot = p,
#   width = 10.5,
#   height = 4.5,
#   units = "in"
# )


res = FacetMuiPlotReBoxBar(data = data,num = c(3:6),result = result,sig_show ="abc",ncol = 4,
                           mult.y = 0.3,
                           lab.yloc = 1.1
                           )
p1_3 = res[[1]]
p1_3+
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme_nature +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)

#基于输出数据使用ggplot出图
lab = result1[[2]] %>% distinct(group ,name,.keep_all = TRUE)
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(data = lab,aes(x=group , y=y ,label=stat),alpha = 0.6)
p1_0+
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme_nature +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)+
  ggplot2::scale_color_manual(values = col.g)

#3 alpha_rare.metm:alpha多样性稀释曲线#---------
  rare <- mean(phyloseq::sample_sums(ps.metm))/10

  library(microbiome)
  library(vegan)

  result = alpha_rare.metm(ps = ps.micro ,
                          group = "Group", method = "Shannon", start = 100, step = rare)

  #--提供单个样本溪稀释曲线的绘制
  p2_1 <- result[[1]] +scale_color_manual(values = col.g)
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



# beta diversity  -----
#4 ordinate.metm: 排序分析#----------
result = ordinate.metm(ps = ps.metm, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05)
p3_1 = result[[1]]
p3_1 +
  scale_fill_manual(values = colset1)+
  scale_color_manual(values = colset1,guide = F) +
  mytheme1

#带标签图形出图
p3_2 = result[[3]]
p3_2 +
  scale_fill_manual(values = colset1)+
  scale_color_manual(values = colset1,guide = F) +
  mytheme1

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
p3_3 +
  scale_fill_manual(values = colset1)+
  scale_color_manual(values = colset1,guide = F) +
  mytheme1

#5 MicroTest.metm:群落水平差异检测-------
dat1 = MicroTest.metm(ps = ps.metm, Micromet = "adonis", dist = "bray")
dat1

#6 pairMicroTest.metm:两两分组群落水平差异检测-------
dat2 = pairMicroTest.metm(ps = ps.metm, Micromet = "MRPP", dist = "bray")

dat2


#7 mantal.metm：群落差异检测普鲁士分析#------
result <- mantal.metm(ps = ps.metm,
                       method =  "spearman",
                       group = "Group",
                       ncol = gnum,
                       nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7 +
  # scale_fill_manual(values = colset1)+
  # scale_color_manual(values = colset1,guide = F) +
  mytheme1


#8 cluster.metm:样品聚类-----
res = cluster.metm (ps= ps.metm,
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
dat = res[[4]]# 聚类矩阵
head(dat)

# compositipn -----
#9 Ven.Upset.metm: 用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = Ven.Upset.metm(ps =  ps.metm,
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

#10 ggVen.Upset.metm:用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = ggVen.Upset.metm(ps = ps.metm,group = "Group")
grid::grid.draw(res[[1]])
dat = res[[2]]
dat
#11 VenSeper.micro: 详细展示每一组中OTU的物种组成 错-------
#---每个部分
library("ggpubr")
library(agricolae)
library(reshape2)

result = VenSuper.metm(ps = ps.metm,
                        group = "Group",
                        num = 6
)
# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]
p7_1+
  # scale_fill_manual(values = colset1)+
  # scale_color_manual(values = colset1,guide = F) +
  mytheme1
# 每部分的otu门类冲积图
p7_2 <- result[[3]]
p7_2+
  # scale_fill_manual(values = colset1)+
  # scale_color_manual(values = colset1,guide = F) +
  theme_bw()
#每个部分序列的数量占比，并作差异
p8 <- result[[2]]
p8


#12 ggflower.metm:花瓣图展示共有特有微生物------
res <- ggflower.metm(ps = ps.metm ,
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

#13 ven.network.metm:venn 网络----
result = ven.network.metm(
  ps = ps.metm,
  N = 0.5,
  fill = "Phylum")
p14  = result[[1]]
p14

dat = result[[2]]
head(dat)

#14 Micro_tern.metm: 三元图展示组成----
ps1 = ps.metm %>% filter_OTU_ps(500)
res = Micro_tern.metm(ps1)
p15 = res[[1]]
p15[[1]] +theme_bw()

dat =  res[[2]]
head(dat)

#15 barMainplot.metm: 堆积柱状图展示组成----
library(ggalluvial)
pst = ps.metm %>% subset_taxa.wt("Species","Unassigned",TRUE)
pst = pst %>% subset_taxa.wt("Genus","Unassigned",TRUE)

result = barMainplot.metm(ps = pst,
                           j = "Genus",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           Top =5)

 p4_1 <- result[[1]]
# scale_fill_manual(values = colset2) +
# scale_x_discrete(limits = axis_order) +
# mytheme1
p4_1+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  mytheme1

p4_2  <- result[[3]]
# # scale_fill_brewer(palette = "Paired") +
# scale_fill_manual(values = colset2) +
# scale_x_discrete(limits = axis_order) +
# mytheme1
p4_2+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group",j,"Abundance(%)")
head(databar)


#16 cluMicro.bar.metm: 聚类堆积柱状图展示组成-----
result <-  cluMicro.bar.metm (dist = "bray",
                               ps = ps.metm,
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

#17 cir_barplot.metf:环状堆积柱状图 -----
library(ggtree) # j = "Phylum"
res = cir_barplot.metm(
  ps = ps.metm,
  Top = 15,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete")

p17 = res[[1]]
p17
dat= res[[2]]
head(dat)

#18 cir_plot.metm:和弦图展示物种组成-----
res = cir_plot.metm(ps  = ps.metm,Top = 12,rank = 6)


#19 Microheatmap.metm: 热图展示物种相对丰度差异-----

heatnum = 30
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
# phyloseq::sample_data(ps) = map

ps_tem = ps.metm %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus")
id <- ps.metm %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus") %>%
  ggClusterNet::filter_OTU_ps(100) %>%
  ggClusterNet::vegan_otu() %>%
  t() %>% as.data.frame() %>%rowCV %>%
  sort(decreasing = TRUE) %>%
  head(heatnum) %>%
  names()
? Microheatmap.metm
result <- Microheatmap.metm(ps_rela = ps_tem,id = id ,col_cluster = FALSE,row_cluster = FALSE)

p24.1 <- result[[1]]
p24.1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p24.2 <- result[[2]]
p24.2
dat = result[[3]]
head(dat)

# difference analysis -----
#20 EdgerSuper.metm:EdgeR计算差异微生物----
res = EdgerSuper.metm(ps = ps.metm %>% ggClusterNet::filter_OTU_ps(500),group  = "Group",artGroup = NULL, j = "Species")

p25.1 =  res[[1]][1]
p25.1
p25.2 =  res[[1]][2]
p25.2
p25.3 =  res[[1]][3]
p25.3
dat =  res[[2]]
head(dat)

#21 EdgerSuper2.metm:EdgeR计算差异微生物-----
res =  EdgerSuper2.metm (ps = ps.metm,group  = "Group",artGroup =NULL, j = "Species")
head(res)

#22 DESep2Super.metm:DESep2计算差异微生物-----
res = DESep2Super.metm(ps = ps.metm %>% ggClusterNet::filter_OTU_ps(500),
                        group  = "Group",
                        artGroup =NULL,
                        j = "Species"
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

#23 edge_Manhattan.metm: 曼哈顿图展示差异微生物------
res = edge_Manhattan.metm(
  ps = ps.metm%>% ggClusterNet::filter_OTU_ps(500),
  pvalue = 0.05,
  lfc = 0
)
p27.1= res[[1]]
p27.1
p27.2= res[[2]]
p27.2
p27.3= res[[3]]
p27.3

#24 stemp_diff.metm: stamp展示差异微生物----
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
#
# # map$Group = as.factor(map$Group)
# sample_data(ps) = map
allgroup <- combn(unique(map$Group),2)
plot_list <- list()
for (i in 1:dim(allgroup)[2]) {
  ps_sub <- phyloseq::subset_samples(ps.metm,Group %in% allgroup[,i]);ps_sub
  p <- stemp_diff.metm(ps = ps_sub,Top = 20,ranks = 6)
  plot_list[[i]] <- p

}

p28.1 = plot_list[[1]]
p28.1
p28.2 = plot_list[[2]]
p28.2
p28.3 = plot_list[[3]]
p28.3

#25 Mui.Group.volcano.metm: 聚类火山图------
res =  EdgerSuper2.metm (ps = ps.metm,group  = "Group",artGroup =NULL, j = "OTU")

res2 = Mui.Group.volcano.metm(res = res)

p29.1 = res2[[1]]
p29.1
p29.2 = res2[[2]]
p29.2
dat = res2[[3]]
dat

# biomarker identification -----
#26 rfcv.metm :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
result = rfcv.metm(ps = ps.metm %>% filter_OTU_ps(100),
                   group  = "Group",optimal = 20,nrfcvnum = 6)

prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#27 Roc.metm:ROC 曲线绘制----
id = sample_data(ps.metm)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)

ps = ps.metm %>% subset_taxa.wt("Family","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Order","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Genus","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Phylum","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("class","Unassigned",TRUE)

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)
res = Roc.metm( ps = pst %>% filter_OTU_ps(1000),group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

#28 loadingPCA.metm: 载荷矩阵筛选特征微生物------
res = loadingPCA.metm(ps = ps.metm,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#29 LDA.metm: LDA筛选特征微生物-----
tablda = LDA.metm(ps = ps.metm,
                   Top = 100,
                   p.lvl = 0.05,
                   lda.lvl = 1,
                   seed = 11,
                   adjust.p = F)

p35 <- lefse_bar(taxtree = tablda[[2]])
p35
dat = tablda[[2]]
dat

#30 svm.metm:svm筛选特征微生物 ----
res <- svm.metm(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#31 glm.metm :glm筛选特征微生物----
res <- glm.metm (ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#32 xgboost.metm: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)
res = xgboost.metm(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#33 lasso.metm: lasso筛选特征微生物----
library(glmnet)
res =lasso.metm(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#34 decisiontree.micro: 错----
library(rpart)
res =decisiontree.metm(ps=pst, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#35 naivebayes.metm: bayes筛选特征微生物----
res = naivebayes.metm(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#36 randomforest.metm: 随机森林筛选特征微生物----
res = randomforest.metm( ps = pst,group  = "Group", optimal = 50)
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

#37 bagging.metm: Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metm(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#38 nnet.metm : 神经网络筛选特征微生物 ------
res =nnet.metm (ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#network analysis -----
#39 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
tab.r = network.pip(
  ps = ps.metm,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label = TRUE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
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
#40 net_properties.4:网络属性计算#---------
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

#41 netproperties.sample:单个样本的网络属性#-------
for (i in 1:length(id)) {
  pst = ps.metm %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(pst = pst,cor = cor[[id[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}

map = map %>% as.tibble()
dat3 = dat.f2 %>% rownames_to_column("ID") %>% inner_join(map,by = "ID")

head(dat3)

#42 node_properties:计算节点属性#---------
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



#43 module.compare.net.pip:网络显著性比较#-----
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
head(res)


