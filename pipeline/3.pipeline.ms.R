# 代谢组学数据挖掘-----
rm(list=ls())
# library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(EasyMultiOmics.db)

ps.ms = EasyMultiOmics::ps.ms

#--提取有多少个分组
gnum = phyloseq::sample_data(ps.ms)$Group %>% unique() %>% length()
gnum

# 设定排序顺序
axis_order = phyloseq::sample_data(ps.ms)$Group %>% unique()

#-主题--
package.amp()
res = theme_my(ps.ms)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

# annoation & normalization
#1 ann.HMDB: 代谢物注释------
tax= ps.ms %>% vegan_tax() %>%
  as.data.frame()
head(tax)
id = tax$Metabolite

#-HMDB数据库注释
detach("mia")
tax1 = ann.HMDB (id = id)
head(tax1)
colnames(tax1)
tax1 = tax1 %>% distinct(id,.keep_all = TRUE) %>%
  column_to_rownames("id")
head(tax1)

tax$ID = row.names(tax)
colnames(tax)
tax3 = tax %>% left_join(tax1,by = c("Metabolite" = "id.org"))
head(tax3)
row.names(tax3) = tax3$metab_id
tax_table(ps.ms) = phyloseq::tax_table(as.matrix(tax3))


#2 ann.kegg2: 代谢物注释KEGG数据库导入文件 ----
tax2 = ann.kegg(id)

head(tax2)
#--注释kegg数据库算法2
tax2 = ann.kegg2(id)
head(tax2)


# tax0 = ps.ms %>% vegan_tax() %>% as.data.frame() #%>% rownames_to_column("ID")
# head(tax0)
# tax2 = tax0 %>% left_join(tax,by = "ID") %>% column_to_rownames("ID")
# head(tax2)
# phyloseq::tax_table(ps.ms) = as.matrix(tax1)
# tax_table(ps.ms)
# saveRDS(ps.ms,".//ps.ms_GC_upper.rds")


#3 zone.fill.ms ----
?zone.fill.ms
ps.ms2 = zone.fill.ms(ps = ps.ms,method = "repeat")

#4 scale_IS.ms ----
ps.ms2 = scale_IS.ms(ps = ps.ms,IS = "metab_8497")

#5 scale_QC.ms  ----
ps.ms2 = scale_QC.ms(ps = ps.ms,
                     QC = c("OE36_2","OE36_3","OE2_1"))

#6 normalize.ms  ----
?normalize.ms
ps.ms2 = normalize.ms(ps = ps.ms,method = "rela")


# ordinate analysis
#7 ordinate.ms:  代谢物排序分析----
?ordinate.ms
result = ordinate.ms(ps = ps.ms, group = "Group", dist = "bray",
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
                          mapping = aes(xend = oNMDS1,
                                        yend = oNMDS2,color = Group),show.legend=F,
                          alpha = 0.5
                          ) + # spiders
  geom_point(mapping = aes(x = x, y = y,fill = Group ),
             data = cent, size = 5,pch = 23,
             color = "black")
p3_3 +
  scale_fill_manual(values = colset1)+
  scale_color_manual(values = colset1,guide = F) +
  mytheme1


#8 MicroTest.ms:代谢组总体差异检测#-------
?MicroTest.ms
dat1 = MicroTest.ms(ps = ps.ms, Micromet = "adonis", dist = "bray")
dat1

#9 pairMicroTest.ms:两两分组代谢总体水平差异检测#-------
?pairMicroTest.ms
dat2 = pairMicroTest.ms(ps= ps.ms, Micromet = "MRPP", dist = "bray")
dat2


#10 mantal.ms:代谢群落差异检测普鲁士分析-----
map= sample_data(ps.ms)
head(map)
?mantal.ms
result <- mantal.ms(ps = ps.ms,
                       method =  "spearman",
                       group = "Group",
                       ncol = gnum,
                       nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7


#11 plsda.ms:plsda排序分析-----
library(ggalt)
library(vegan)
library(mixOmics)
library(ggplot2)
library(ggalt)
library(ggforce)
library(caret)
?plsda.ms
res = plsda.ms(ps=ps.ms,Group = "Group")
p11 = res [[1]]
p11
dat = res[[2]]
dat

#12 oplsda.ms:oplsda排序分析 ----
library(pacman)
library(ggsci)
library(ropls)
# BiocManager::install("ropls")
?oplsda.ms
res = oplsda.ms(ps = ps.ms, ncol=3,nrow = 1)

p12 = res[1]
dat = res[2]


#  metabolite classification
#13 Ven.Ups.mset.ms: 用于展示共有、特有的代谢物----
# 分组小于6时使用
?Ven.Upset.ms
res = Ven.Upset.ms(ps =  ps.ms,
                        group = "Group",
                        N = 0.5,
                        size = 3)
p10.1 = res[[1]]
p10.1
dat = res[[2]]
dat

#14 ggflower.ms: 花瓣图展示共有特有代谢物------
?ggflower.ms
res <- ggflower.ms(ps= ps.ms ,
                       # rep = 1,
                       group = "Group",
                       start = 1, # 风车效果
                       m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                       a = 0.2, # 花瓣胖瘦
                       b = 1, # 花瓣距离花心的距离
                       lab.leaf = 1, # 花瓣标签到圆心的距离
                       col.cir = "yellow",
                       N = 0.5 )

p14 = res[[1]]
p14

dat = res[[2]]
dat



#15 Ms_tern.ms: 三元图展示组成----
tax_table(ps.ms)
ps1 = ps.ms %>% filter_OTU_ps(500)
?Ms_tern.ms
res = Ms_tern.ms(ps1, color = "Mode")
p15 = res[[1]]
p15[[1]] +theme_bw()

dat =  res[[2]]
head(dat)

#16 barMainplot.ms 代谢物分类堆叠柱状图#--------

j = "Class"
rank.names(ps.ms)
strbar = c("Super_class","Class" , "Sub_class"  )
?barMainplot.ms
 result = barMainplot.ms(ps = ps.ms,
                          j = "Class",
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top = 12)
  p4_1 <- result[[1]] +scale_fill_brewer(palette = "Paired")
  p4_1
  p4_2  <- result[[3]] +
     scale_fill_manual(values = colset3) +
   scale_x_discrete(limits = axis_order)
    # mytheme1
  p4_2

  databar <- result[[2]] %>%
    dplyr::group_by(Group,aa) %>%
    dplyr::summarise(sum(Abundance)) %>% as.data.frame()
  colnames(databar) = c("Group",j,"Abundance(%)")
  head(databar)



#17 cluMicro.bar.micro: 聚类堆积柱状图展示组成 -----
  ?cluMicro.bar.ms
  res <-  cluMicro.bar.ms (dist = "bray",
                                 ps= ps.ms,
                                 j = "Class",
                                 Top = 10, # 提取丰度前十的物种注释
                                 tran = TRUE, # 转化为相对丰度值
                                 hcluter_method = "complete",
                                 Group = "Group",
                                 cuttree = length(unique(phyloseq::sample_data(ps.ms)$Group)))


  p17.1 = res[[1]]
  p17.1
  p17.2 <- res[[2]]
  p17.2
  p17.3 <- res[[3]]
  p17.3
  p17.4 <- res[[4]]
  p17.4
  clubardata <- result[[5]]
  clubardata





# difference analysis
#18 cluster_plot.ms:  代谢物 层次聚类--------
  ?cluster_plot.ms
  res = cluster_plot.ms (ps= ps.ms,
                         hcluter_method = "complete",
                         dist = "bray",cuttree = gnum,
                         row_cluster = TRUE,
                         col_cluster =  TRUE)

  p0 = res[[1]]
  p0
  p1 = res[[2]]
  p1
  p2 = res[[3]]
  p2
  dat = res[4]
  dat

#19 heatmap.ms:  热图展示代谢物差异----
ps.ms_rela <- ps.ms %>% scale_micro(method = "rela") %>%
      tax_glom_wt(ranks = "Class")
?heatmap.ms

result <- heatmap.ms (ps_rela= ps.ms_rela,
                          label =  TRUE,
                          col_cluster = TRUE,
                          row_cluster =TRUE)
    p19 <- result[[1]]
    p19
    # p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
    p19.1<- result[[2]]
    p19.1
    dat = result[[3]]
    dat

#20 statSuper:  差异代谢物#----------
 #--非参数检验
    ?statSuper
   result1 = statSuper(ps = ps.ms,group  = "Group",artGroup = NULL,method = "wilcox")
    head(result1)

    #--t检验检验--建议四个重复以上
    result2 = statSuper(ps = ps.ms,group  = "Group",artGroup = NULL,method = "ttext")
    head(result2)

#21 MuiKwWlx2: 分类化合物分组差异-------
library(EasyStat)
library(ggClusterNet)

map = sample_data(ps.ms)
head(map)
map = map[,1:2]
sample_data(ps.ms) = map


dat <- ps.ms %>% scale_micro(method = "rela") %>%
    tax_glom_wt(ranks = "Class") %>%
    vegan_otu() %>%
    as.data.frame()
  head(dat)

  dat$id = row.names(dat)

  dat2 = dat %>%
    dplyr::left_join(as.tibble(sample_data(ps.ms)),by = c("id" = "ID")) %>%
    dplyr::rename(group = Group) %>%
    dplyr::select(id,group,everything())
  # dat2 %>%
  #   dim()

  dat2$group = as.factor(dat2$group)
  head(dat2)

  result = MuiKwWlx2(data = dat2,num = c(3:dim(dat2)[2]))

  result1 = EasyStat::FacetMuiPlotresultBox(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",
                                            ncol = 4 )
  p1_1 = result1[[1]] +
    # scale_x_discrete(limits = axis_order) +
    mytheme2 +
    guides(fill = guide_legend(title = NULL)) +
    scale_fill_manual(values = colset1)
  p1_1

  res = FacetMuiPlotresultBar(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",
                              ncol = 4)
  p1_2 = res[[1]]+
    # scale_x_discrete(limits = axis_order) +
    guides(color = FALSE) +
    mytheme2 +
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_2

  res = FacetMuiPlotReBoxBar(data = dat2,num = c(3:dim(dat2)[2]),result = result,sig_show ="abc",ncol = 4)
  p1_3 = res[[1]]+
    # scale_x_discrete(limits = axis_order) +
    mytheme2 +
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_3



#22 FacetMuiPlotresultBox :单变量统计分析-箱线图可视化--------

#--提取差异代谢物标签
dat = ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()
head(dat)
# colnames(map)
map = sample_data(ps.ms) %>% as.tibble() %>%
  select(ID,Group)
data = cbind(map[,c(1,2)],dat)
head(data)
colnames(data)[2] = "group"

# num = c(3:ncol(data))

# num = 18#--为了减少运行压力，修改为18个化合物
num

#--分割数据
n.fac = length(num)/ 25
n.fac2 = ceiling(n.fac)
A = list()
# j  =2
for (j in 1:n.fac2) {

  if (j == 1) {
    A[[j]] = num[1:25]
  } else if(j != n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:y]
  }else if (j == n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:length(num)]
  }

}


plot_list1 = list()

# i=1

for (i in 1:n.fac2) {
  result = EasyStat::MuiaovMcomper2(data = data,num = A[[i]])

  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = A[[i]],
                                            result = result,
                                            sig_show ="abc",ncol = 5 )
  p1_1 = result1[[1]] +
    ggplot2::scale_x_discrete(limits = axis_order) +
    mytheme2 +
    ggplot2::guides(fill = guide_legend(title = NULL)) +
    ggplot2::scale_fill_manual(values = colset1)
  p1_1

  plot_list1[[i]] =  p1_1


  return(plot_list1)

}
res = plot_list1


#23 FacetMuiPlotresultBar :单变量统计分析-柱状图可视化--------

#--提取差异代谢物标签
dat = ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()
head(dat)
# colnames(map)
map = sample_data(ps.ms) %>% as.tibble() %>%
  select(ID,Group)
data = cbind(map[,c(1,2)],dat)
head(data)
colnames(data)[2] = "group"

# num = c(3:ncol(data))

# num = 18#--为了减少运行压力，修改为18个化合物
num

#--分割数据
n.fac = length(num)/ 25
n.fac2 = ceiling(n.fac)
A = list()
# j  =2
for (j in 1:n.fac2) {

  if (j == 1) {
    A[[j]] = num[1:25]
  } else if(j != n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:y]
  }else if (j == n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:length(num)]
  }

}


plot_list2 = list()

# i=1

for (i in 1:n.fac2) {


  res = EasyStat::FacetMuiPlotresultBar(data = data,num = A[[i]],
                                        result = result,sig_show ="abc",ncol = 5)
  p1_2 = res[[1]]+ scale_x_discrete(limits = axis_order) + guides(color = FALSE) +
    mytheme2+
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_2

  plot_list2[[i]] =  p1_2

  return(plot_list2)

}



#24 FacetMuiPlotReBoxBar: 单变量统计分析-柱状图结合散点图可视化--------

#--提取差异代谢物标签
dat = ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()
head(dat)
# colnames(map)
map = sample_data(ps.ms) %>% as.tibble() %>%
  select(ID,Group)
data = cbind(map[,c(1,2)],dat)
head(data)
colnames(data)[2] = "group"

# num = c(3:ncol(data))

# num = 18#--为了减少运行压力，修改为18个化合物
num

#--分割数据
n.fac = length(num)/ 25
n.fac2 = ceiling(n.fac)
A = list()
# j  =2
for (j in 1:n.fac2) {

  if (j == 1) {
    A[[j]] = num[1:25]
  } else if(j != n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:y]
  }else if (j == n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:length(num)]
  }

}

plot_list3 = list()

# i=1

for (i in 1:n.fac2) {



  res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = A[[i]],
                                       result = result,sig_show ="abc",ncol = 5)
  p1_3 = res[[1]]+ scale_x_discrete(limits = axis_order) +
    mytheme2 +
    guides(fill = guide_legend(title = NULL))+
    scale_fill_manual(values = colset1)
  p1_3

  plot_list3[[i]] =  p1_3



  return(plot_list3)

}





#25 MuiHeatmapBubplot: 单变量统计分析-气泡图可视化--------

#--提取差异代谢物标签
dat = ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()
head(dat)
# colnames(map)
map = sample_data(ps.ms) %>% as.tibble() %>%
  select(ID,Group)
data = cbind(map[,c(1,2)],dat)
head(data)
colnames(data)[2] = "group"

# num = c(3:ncol(data))

# num = 18#--为了减少运行压力，修改为18个化合物
num

#--分割数据
n.fac = length(num)/ 25
n.fac2 = ceiling(n.fac)
A = list()
# j  =2
for (j in 1:n.fac2) {

  if (j == 1) {
    A[[j]] = num[1:25]
  } else if(j != n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:y]
  }else if (j == n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:length(num)]
  }

}


plot_list4 = list()
plot_list5 = list()
plot_list6 = list()
plot_list7 = list()
# i=1

for (i in 1:n.fac2) {



  res = EasyStat::MuiHeatmapBubplot(
    data = data,
    i =A[[i]],
    col_cluster = F,
    row_cluster = F,
    label = TRUE,
    result = result,
    sample = TRUE,
    scale = TRUE
  )
  p1 = res[[1]]
  p1

  plot_list4[[i]] =  p1

  p2 = res[[2]]
  p2

  plot_list5[[i]] =  p2

  res = EasyStat::MuiHeatmapBubplot(
    data = data,
    i =A[[i]],
    result = result,
    col_cluster = F,
    row_cluster = F,
    label = TRUE,
    sample = FALSE,
    scale = TRUE


  )

  p1 = res[[1]]
  p1

  plot_list6[[i]] =  p1


  p2 = res[[2]]
  p2
  plot_list7[[i]] =  p2



  return(plot_list4,plot_list5,plot_list6,plot_list7)

}



#26 MuiHeatmapBubplot 单变量统计分析-热图可视化--------

#--提取差异代谢物标签
dat = ps.ms %>%
  ggClusterNet::vegan_otu() %>%
  as.data.frame()
head(dat)
# colnames(map)
map = sample_data(ps.ms) %>% as.tibble() %>%
  select(ID,Group)
data = cbind(map[,c(1,2)],dat)
head(data)
colnames(data)[2] = "group"

# num = c(3:ncol(data))

# num = 18#--为了减少运行压力，修改为18个化合物
num

#--分割数据
n.fac = length(num)/ 25
n.fac2 = ceiling(n.fac)
A = list()
# j  =2
for (j in 1:n.fac2) {

  if (j == 1) {
    A[[j]] = num[1:25]
  } else if(j != n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:y]
  }else if (j == n.fac2){
    x = (25*(j - 1) + 1)
    y = 25*j
    A[[j]] = num[x:length(num)]
  }

}



plot_list8 = list()
# i=1

for (i in 1:n.fac2) {

  res = EasyStat::value_stackBar(
    data = data,
    i =A[[i]],
    result = result,
    add_abc = TRUE)

  p1 = res[[1]]
  p1
  plot_list8[[i]] =  p1

  return(plot_list8)

}




















# biomarker identification-----
id = sample_data(ps.ms)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])


pst = ps.ms %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)



#27 rfcv.ms :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
result =rfcv.ms(ps = ps.ms,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#28 randomforest.ms: 筛选特征代谢物 ----
?randomforest.ms
res <- randomforest.ms( ps= ps.ms, group  = "Group", optimal = 50)

p25 = res[[1]]
p25
dat =res[[2]]
dat



#29 loadingPCA.ms 载荷矩阵挑选重要代谢物#------
?loadingPCA.ms
res = loadingPCA.ms(ps = ps.ms,Top = 20)
p = res[[1]]
p
dat = res[[2]]
dat

#30 LDA.ms: LDA筛选特征代谢物 -----
?LDA.ms

tablda = LDA.ms(ps = ps.ms,
                   Top = 100,
                   p.lvl = 0.05,
                   lda.lvl = 1,
                   seed = 11,
                   adjust.p = F)

p35 <- lefse_bar(taxtree = tablda[[2]])
p35
dat = tablda[[2]]
dat


#31 svm.ms:svm筛选特征微生物 ----
?svm.ms
res <- svm.ms(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#32 xgboost.ms: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
?xgboost.ms
 res = xgboost.ms(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#33 glm.ms:glm筛选特征微生物----
?glm.ms
res <- glm.ms(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#34 lasso.ms: lasso筛选特征微生物----
library(glmnet)
?lasso.ms
res =lasso.ms (ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#35 decisiontree.ms----
library(rpart)
res =decisiontree.ms(ps=ps.ms, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#37 nnet.ms 神经网络筛选特征微生物  ------
?nnet.ms
res =nnet.ms(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance



#38 bagging.micro : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
?bagging.ms
res =bagging.ms(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance




# network analysis-------
#39 network.pip:网络分析主#--------
tab.r = network.pip(
  ps = ps.ms,
  N = 400,
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
  ps.mst = ps.ms %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(ps.mst = ps.mst,cor = cor[[id[i]]])
  # head(dat.f)
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}

map= sample_data(ps.ms)
map$ID = row.names(map)
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



# pathway enrich------
#44 pathway_enrich.ms：通路富集-------

ps.ms3 = ps.ms %>% tax_glom_wt("KEGGID")
?pathway_enrich.ms

res2 = pathway_enrich.ms(ps = ps.ms3, dif.method = "wilcox")

res2$plots$OE.WT.plot
res2$plots$KO.OE.plot

res2$plotdata$OE.WT

#45 reaction.show.ms：反应展示-------
EasyMultiOmics::reaction.show.ms
res3 = reaction.show.ms(ps= ps.ms3,dif.method = "wilcox")
res3$plots$OE.WT.plot
res3$plotdata$OE.WT

#46 buplotall.ms：气泡图展示富集分析结果-------

res= buplotall.ms(ps= ps.ms3,dif.method = "wilcox")

res$plots$OE.WT.plot

res3$plotdata$OE.WT


