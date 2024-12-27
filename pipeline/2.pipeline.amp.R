# 扩增子微生物组学分析#-----
rm(list=ls())
# library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)

# 定义后续参数-----
map= sample_data(ps.16s)
head(map)
# 提取分组因子数量
gnum = phyloseq::sample_data(ps.16s)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.16s对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps.16s)$Group %>%unique();axis_order

#-主题--颜色等
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

# alpha diversity -----
#1 alpha.micro: 6中alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.micro(ps = ps.16s,group = "Group" )
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
p1_1+
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)


res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(9)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)

res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(9)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0+
  ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)

#2 alpha.pd :用于计算pd多样性#-------
library(ape)
library(picante)
tab2 = alpha.pd(ps16s)
head(tab2)
result = EasyStat::MuiKwWlx2(data = tab2,num = 3)
result1 = EasyStat::FacetMuiPlotresultBox(data = tab2,num = 3,
                                          result = result,
                                          sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]]
p1_1+
  #ggplot2::scale_x_discrete(limits = axis_order) +
  mytheme2 +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)




#3 alpha.rare.line:alpha多样性稀释曲线#---------
rare <- mean(phyloseq::sample_sums(ps.16s))/10
result = alpha.rare.line(ps = ps.16s, group = "Group", method = "Richness", start = 100, step = rare)
#-- Plot the rarefaction curve for a single sample
p2_1 <- result[[1]]
p2_1
## Provide a data table for convenient output
raretab <- result[[2]]
head(raretab)
#-- Display rarefaction curves grouped by categories
p2_2 <- result[[3]]
p2_2
#-- Plot rarefaction curves with standard deviations by groups
p2_3 <- result[[4]]
p2_3

# beta diversity-----
#4 ordinate.micro: 排序分析#----------
result = ordinate.micro(ps = ps.16s, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                        pair = F)
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

#5 MicroTest.micro:群落水平差异检测-------
dat1 = MicroTest.micro(ps = ps.16s, Micromet = "adonis", dist = "bray")
dat1

#6 pairMicroTest:两两分组群落水平差异检测-------
dat2 = pairMicroTest.micro(ps = ps.16s, Micromet = "MRPP", dist = "bray")
dat2

#7 mantal.micro ：群落差异检测普鲁士分析#------
result <- mantal.micro(ps = ps.16s,
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

#8 cluster.micro:样品聚类#-----
res = cluster.micro (ps= ps.16s,
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


#9 distance.micro:分组之间距离比对#-----
res = distance.micro(ps.16s,group = "Group")
p5.1 = res[[1]]
p5.1 +
  # scale_fill_manual(values = colset1)+
  # scale_color_manual(values = colset1,guide = F) +
  mytheme1
p5.2 = res[[2]]
p5.2 +
   scale_fill_manual(values = colset1)+
   scale_color_manual(values = colset1,guide = F) +
  mytheme1
p5.3 = res[[3]]
p5.3 +
   scale_fill_manual(values = colset1)+
   scale_color_manual(values = colset1,guide = F) +
  mytheme1
dat = res[[4]]
head(dat)


# compositipn -----
#10 Ven.Upset.micro: 用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = Ven.Upset.micro(ps =  ps.16s,
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
res = ggVen.Upset.micro(ps = ps.16s,group = "Group")
grid::grid.draw(res[[1]])
dat = res[[2]]
dat


#12 VenSeper.micro: 详细展示每一组中OTU的物种组成 -------
#---每个部分
library("ggpubr")
library(agricolae)
library(reshape2)

result = VenSuper.micro(ps = ps.16s,
                        group =  "Group",
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


#13 ggflower.micro:花瓣图展示共有特有微生物------
res <- ggflower.micro (ps = ps.16s ,
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



#14 ven.network.micro:venn 网络----
result = ven.network.micro(
  ps = ps.16s,
  N = 0.5,
  fill = "Phylum")
p14  = result[[1]]
p14

dat = result[[2]]
head(dat)


#15 Micro_tern.micro: 三元图展示组成----
ps1 = ps.16s %>% filter_OTU_ps(500)
res = Micro_tern.micro(ps1)
p15 = res[[1]]
p15[[1]] +theme_bw()

dat =  res[[2]]
head(dat)

#16 barMainplot.micro: 堆积柱状图展示组成----
library(ggalluvial)
pst = ps.16s %>% subset_taxa.wt("Species","Unassigned",TRUE)
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




#17 cluMicro.bar.micro: 聚类堆积柱状图展示组成 问题-----
result <-  cluMicro.bar.micro (dist = "bray",
                               ps = ps.16s,
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
  ps = ps.16s,
  Top = 15,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete")

p17 = res[[1]]
p17
dat= res[[2]]
head(dat)

#19 cir_plot.micro:和弦图展示物种组成-----
res = cir_plot.micro(ps  = ps.16s,Top = 12,rank = 6)




#20 maptree.micro；圈图展示物种组成-----

library(ggraph)
library(data.tree)
library(igraph)

tax = ps.16s %>% vegan_tax() %>%
  as.data.frame()
head(tax)
ps.bac <- ps.16s %>% subset_taxa.wt("Kingdom", "Bacteria")
ps.bac

res =maptree.micro (ps = ps.bac,
                     Top = 100,
                     labtab =  NULL,
                     seed = 11)
p20 = res[[1]]
p20
dat = res[[2]]
head(dat)

#21 phy_tree.micro；树状图展示物种组成以及相关性-----
library(ggstar)
result <- phy_tree.micro(ps = ps.16s,Top = 100)

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
head(dat)
#22 sankey.micro: 桑基图展示物种组成------
res = sankey.micro (ps = ps.16s,
                    rank = 6,
                    Top = 50)
p22 = res[[1]]
p22
dat = res[[2]]
dat

#23 sankey.m.Group.micro: 按照分组绘制-----
result = sankey.m.Group.micro(
  ps = ps.16s  %>% subset_taxa.wt("Species","Unassigned",TRUE),
  rank = 6,
  Top = 50)

p = result[[1]]
p

saveNetwork(p,paste("./sankey_Group.html", sep = ""))

dat = result[[2]]
dat


#24 Microheatmap.micro: 热图展示物种相对丰度差异-----

heatnum = 30
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
# phyloseq::sample_data(ps) = map

ps_tem = ps.16s %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus")
id <- ps.16s %>%
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


# difference analysis -----
#25 EdgerSuper.micro:EdgeR计算差异微生物----
res = EdgerSuper.micro(ps = ps.16s %>% ggClusterNet::filter_OTU_ps(500),group  = "Group",artGroup = NULL, j = "OTU")


p25.1 =  res[[1]][1]
p25.1
p25.2 =  res[[1]][2]
p25.2
p25.3 =  res[[1]][3]
p25.3
dat =  res[[2]]
head(dat)

#25.1 EdgerSuper2.micro:EdgeR计算差异微生物-----
res =  EdgerSuper2.micro (ps = ps.16s,group  = "Group",artGroup =NULL, j = "OTU")
head(res)


#26 DESep2Super.micro:DESep2计算差异微生物-----
res = DESep2Super.micro(ps = ps.16s %>% ggClusterNet::filter_OTU_ps(500),
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
  ps = ps.16s%>% ggClusterNet::filter_OTU_ps(500),
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
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
#
# # map$Group = as.factor(map$Group)
# sample_data(ps) = map

allgroup <- combn(unique(map$Group),2)
plot_list <- list()
for (i in 1:dim(allgroup)[2]) {
  ps_sub <- phyloseq::subset_samples(ps.16s,Group %in% allgroup[,i]);ps_sub
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

res =  EdgerSuper2.micro (ps = ps.16s,group  = "Group",artGroup =NULL, j = "OTU")

res2 = Mui.Group.volcano.micro(res = res)

p29.1 = res2[[1]]
p29.1
p29.2 = res2[[2]]
p29.2
dat = res2[[3]]
dat

#30 Mui.cluster.volcano.micro: 指定分组绘制聚类火山图------
library(ggrepel)
id = sample_data(ps.16s)$Group %>% unique()
aaa = combn(id,2)
# 设置分组
group2 = c(aaa[1,1],aaa[2,1])
b= data.frame(group2)

res =  EdgerSuper2.micro (ps = ps.16s,group  = "Group",artGroup =b, j = "OTU")

res2 = Mui.cluster.volcano.micro (res = res,rs.k = 6)
p30.1 = res2[[1]]+ggtitle(group2)
p30.1
p30.2 = res2[[2]]+ggtitle(group2)
p30.2
dat= res2[[3]]
dat





# biomarker identification -----
#32 rfcv.micro :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.Micro(ps = ps.16s %>% filter_OTU_ps(100),
                   group  = "Group",optimal = 20,nrfcvnum = 6)

prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable


#33 Roc.micro:ROC 曲线绘制----
id = sample_data(ps.16s)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)

ps = ps.16s %>% subset_taxa.wt("Family","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Order","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Genus","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Phylum","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("class","Unassigned",TRUE)

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
res = loadingPCA.micro(ps = ps.16s,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat


#35 LDA.micro: LDA筛选特征微生物-----
p1 <- p_base.micro(ps.16s,Top = 100)
tablda = LDA.micro(ps = ps.16s,
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

#36 svm.micro:svm筛选特征微生物 ----
res <- svm.micro(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#37 glm.micro:glm筛选特征微生物----
res <- glm.micro(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#38 xgboost.micro: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.micro(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#39 lasso.micro: lasso筛选特征微生物----
library(glmnet)
res =lasso.micro (ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance




#40 decisiontree.micro: 错----

top = 20
seed = 6358
set.seed(seed)

ps.cs <- pst %>% filter_OTU_ps(top)

map <- as.data.frame(phyloseq::sample_data(ps.cs))
otutab <- as.data.frame(t(ggClusterNet::vegan_otu(ps.cs %>% scale_micro()) ))
colnames(otutab) <- gsub("-", "_", colnames(otutab))
test <- as.data.frame(t(otutab))
test$group <- factor(map$Group)


#  group

data <-  test
# 分割数据为训练集和测试集


split <- caTools::sample.split(data$group , SplitRatio = 0.7)  # 将数据按照指定比例分割
train_data <- subset(data, split == TRUE)  # 训练集
test_data <- subset(data, split == FALSE)
str(train_data)
head(train_data)
# 训练决策树模型
a_rpart <- rpart(group ~ ., data =train_data, method = 'class',
                 parms = list(split = 'information'))


plot(a_rpart, margin = 0.1)
text(a_rpart, cex = 0.5)
a_rpart$cptable

#决策树划分细节概要，各个分支结构等
summary( a_rpart)

# 得到测试集的预测值
pred <- predict(a_rpart, newdata =test_data , type = 'class')




#41 naivebayes.micro: bayes筛选特征微生物----
res = naivebayes.micro(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#42 randomforest.micro: 随机森林筛选特征微生物----
res = randomforest.micro( ps = pst,group  = "Group", optimal = 50)
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4


#43 bagging.micro : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.micro(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance
#44 nnet.micro 神经网络筛选特征微生物  ------

res =nnet.micro(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#network analysis -----
#45 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
tab.r = network.pip(
  ps = ps.16s,
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
#46 net_properties.4:网络属性计算#---------
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

#47 netproperties.sample:单个样本的网络属性#-------
for (i in 1:length(id)) {
  pst = ps.16s %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
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

#48 node_properties:计算节点属性#---------
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



#49 negative.correlation.ratio:网络稳定性-计算负相关的比例#----
res4 = negative.correlation.ratio(ps = ps.16s,
                                  corg = cortab,
                                  # Top = 500,
                                  degree = TRUE,
                                  zipi = FALSE)

p5 = res4[[1]]
p5
dat6 = res4[[2]]


#50 community.stability:网络稳定性-群落稳定性-只有pair样本使用-----
treat = ps.16s %>% sample_data()
treat$pair = paste( "A",c(rep(1:10,3)),sep = "")
# head(treat)
sample_data(ps.16s) = treat

#一般性的处理，没有时间梯度的，这里设定time为F，意味着每两个群落的组合都进行比较
res5 = community.stability( ps = ps.16s,
                            corg = cor,
                            time = FALSE)
p6 = res5[[1]]
p6
dat7 = res5[[2]]
dat7

#51 natural.con.microp:网络稳定性-网络抗毁性#------
library("pulsar")
res6 = natural.con.microp (
  ps = ps.16s,
  corg = cor,
  norm = TRUE,
  end = 150,# 小于网络包含的节点数量
  start = 0
)
p7 = res6[[1]]
p7
dat8  = res6[[2]]
dat8


#52 module.compare.net.pip:网络显著性比较#-----
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


#53 module.compare.m:网络稳定性-模块相似性#--------
library(tidyfst)

res1 = module.compare.m(
  ps = NULL,
  corg = cor,
  zipi = FALSE,
  zoom = 0.6,
  padj = FALSE,
  n = 3)

#不同分组使用一个圆圈展示，圆圈内一个点代表一个模块，相连接的模块代表了相似的模块。
p1 = res1[[1]]
p1
#--提取模块的OTU，分组等的对应信息
dat1 = res1[[2]]
head(dat1)
#模块相似度结果表格
dat2 = res1[[3]]
head(dat2)

dat2$m1 = dat2$module1 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$m2 = dat2$module2 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$cross = paste(dat2$m1,dat2$m2,sep = "_Vs_")
# head(dat2)
dat2 = dat2 %>% filter(module1 != "none")

p2 = ggplot(dat2) + geom_bar(aes(x = cross,fill = cross)) +
  labs(x = "",
       y = "numbers.of.similar.modules"
  )+ theme_classic()

p2




#54 Robustness.Targeted.remova:l网络稳定性-去除关键节点-网络鲁棒性#-----
# 鲁棒性计算需要物种丰富，所以即使计算好了相关矩阵，也需要输入ps对象
library(patchwork)
res2= Robustness.Targeted.removal(ps = ps.16s,
                                  corg = cor,
                                  degree = TRUE,
                                  zipi = FALSE
)
p3 = res2[[1]]
p3
#提取数据
dat4 = res2[[2]]
dat4

#55 Robustness.Random.removal网络稳定性-随即取出任意比例节点-网络鲁棒性#---------
res3 = Robustness.Random.removal(ps = ps.16s,
                                 corg = cortab,
                                 Top = 0
)
p4 = res3[[1]]
p4
#提取数据
dat5 = res3[[2]]
dat5


# community assemble-----
#56 neutralModel: 中性模型----
library(picante)
library(ape)
library(FSA)
library(eulerr)
require(minpack.lm)
require(Hmisc)
require(stats4)
sample_data(ps16s)
psphy = filter_taxa(ps16s, function(x) sum(x ) > 100, TRUE);psphy
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



#57 nullModel: 零模型----
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


#58 bNTICul:β最近分类单元指数计算----
result = bNTICul(ps = psphy,
                 group  = "Group",
                 num = 10,
                 thread = 1
)
bNTI = result[[1]]
head(bNTI)


#59 RCbary: RCbary 计算----
result = RCbary(ps = psphy ,group  = "Group",num = 10,thread = 1)

RCbary = result[[1]]
head(RCbary)

#60 bNTIRCPlot: BetaNTI和RCbray联合出图--------
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





# other------
#61 FEAST.micro:溯源分析----
result = FEAST.micro(ps = ps.16s,
                     group = "Group",
                     sinkG = "OE",
                     sourceG = c("WT","KO")
)

#62 Plot_FEAST:溯源分析可视化分组----
p <- Plot_FEAST(data = result)
p
#63 FEAST.micro:溯源分析可视化样品----
p2 = MuiPlot_FEAST(data = result)
p2



# function prediction-----
library(GO.db)
library(DOSE)
library(GO.db)
library(GSEABase)
library(ggtree)
library(aplot)
library(clusterProfiler)
library("GSVA")
library(dplyr)
ps.kegg = ps.kegg %>% filter_OTU_ps(Top = 1000)
tax =ps.kegg %>% tax_table()
colnames(tax)[3] = "KOid"
tax_table(ps.kegg) =tax
res = EdgerSuper.metf (ps = ps.kegg,
                       group  = "Group",
                       artGroup = NULL)
dat = res[[2]]
dat
KEGG_enrich.metf
#64 KEGG_enrich.micro: taxfun2功能富集分析----
res2 = KEGG_enrich.micro(ps =  ps.kegg,
                        #  diffpath = diffpath,
                        dif = dat
)
dat1= res2$`KO-OE`
dat2= res2$`KO-WT`
dat3= res2$`OE-WT`

#65 buplotall.micro: taxfun2功能富集分析气泡图----

Desep_group <-ps.kegg %>% sample_data() %>%
  .$Group %>%
  as.factor() %>%
  levels() %>%
  as.character()
Desep_group
cbtab = combn(Desep_group,2)
cbtab
Desep_group = cbtab[,1]
group = paste(Desep_group[1],Desep_group[2],sep = "-")
id = paste(group,"level",sep = "")

result = buplot.micro(dt  =dat1,dif = dif,id = id)
p1 = result[[1]]
p1
p2 = result[[2]]
p2



