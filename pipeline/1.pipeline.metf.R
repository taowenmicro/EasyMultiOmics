# 宏基因组功能
library(EasyMultiOmics)
library(tidyverse)
library(data.table)
library(phyloseq)
library(ggClusterNet)
library(EasyStat)
library(tidyfst)
library(fs)

#-主题--颜色----
res = theme_my()
mytheme1 = res[[1]]
mytheme2 = res[[2]]; 
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]


# 按照功能数据库分类-------
# kegg数据库--------

data("ps.kegg")

ps =ps.kegg %>% filter_OTU_ps(Top = 1000)

tax = ps %>% tax_table() 
colnames(tax)[3] = "KOid"
tax_table(ps) =tax 

head(tax)
# function diversity -----
#1 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps,group = "Group",Plot = TRUE )
head(tab)


data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)
# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = EasyStat::MuiKwWlx2(data = data,num = 3:5)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:5,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0


#2 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps))/10
result = alpha_rare.metf(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]]
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

#3 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps,
                       group = "Group",
                       dist = "bray",
                       method = "PCA",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                       pair = FALSE)
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




#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2


#6 mantal.metf：群落功能差异检测普鲁士分析#------

result <- mantal.metf(ps = ps,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

#7 cluster.metf:样品聚类#-----

res = cluster.metf(ps= ps,
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


#8 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)

res = Micro_tern.metf(ps=ps1,color = "Pathway"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# function classfication 

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps.t,
                          j = "Pathway",
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Pathway","Abundance(%)")
head(databar)

#10 Ven.Upset.metf: 用于展示共有、特有的功能----
# 分组小于6时使用
res = Ven.Upset.metf(ps =  ps,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

#11 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps)$Group %>% table()
num =10
result = VenSuper.metf(ps = ps,
                       group = "Group",
                       num =10 )

# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p7_1
#每个部分序列的数量占比，并作差异
dat<- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]+scale_fill_brewer(palette = "Paired")
p7_2



#12 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps ,
                      # rep = 1,
                      group = "ID",
                      start = 1, # 风车效果
                      m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                      a = 0.2, # 花瓣胖瘦
                      b = 1, # 花瓣距离花心的距离
                      lab.leaf = 1, # 花瓣标签到圆心的距离
                      col.cir = "yellow",
                      N = 0.5 )

p13.1 = res[[1]]
p13.1

dat = res[[2]]
dat

#13 ven.network.metf: ven网络展示共有特有功能----
map =  sample_data(ps)
result =ven.network.metf(
  ps = ps,
  N = 0.5,
  fill = "Pathway")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat


# function differential analysis----
#14 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps,
                        group  = "Group",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

#15 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps,
                       group  = "Group",
                       artGroup = NULL)
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][1]
p16.2
p16.3 = res[[1]][1]
p16.3
dat = res[[2]]
dat

#16 t.metf: 差异分析t检验----
res = t.metf(ps = ps,group  = "Group",artGroup =NULL)
p17.1 = res [[1]][1]
p17.1
p17.2 = res [[1]][2]
p17.2
p17.3 = res [[1]][3]
p17.3

dat =  res [[2]]
dat

#17 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps,group  = "Group",artGroup =NULL)
p18.1 = res [[1]][1]
p18.1$`KO-OE_plot`
p18.2 = res [[1]][2]
p18.2
p18.3 = res [[1]][3]
p18.3

dat =  res [[2]]
dat %>% head()

#18 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2








# biobiomarker identification-----
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)
#19 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable


#20 Roc.metf:ROC 曲线绘制----


res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat


#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#22 svm.metf:svm筛选特征微生物 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#23 glm.metf:glm筛选特征微生物----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance


#24 xgboost.metf: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance



#25 lasso.metf: lasso筛选特征微生物----

library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#26 decisiontree.micro: 错----
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


#27 naivebayes.metf: bayes筛选特征微生物----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#28 LDA.metf: LDA筛选特征微生物-----
tablda = LDA.metf(ps = pst,
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

#29 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps,group  = "Group", optimal = 50 ,fill ="Pathway"  )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4


#30 nnet.metf: 神经网络筛选特征微生物  ------
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#31 bagging.metf : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#network analysis -----
#32 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
tab.r = network.pip(
  ps = ps,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label = FALSE,
  lab = "elements",
  group = "Group",
  #fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1,
  fill = "Pathway" 
)

dat = tab.r[[2]]

cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]
p0
#33 net_properties.4:网络属性计算#---------
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

#34 netproperties.sample:单个样本的网络属性#-------
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

#35 node_properties:计算节点属性#---------
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



#36 module.compare.net.pip:网络显著性比较#-----
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









# pathway enrich------
#37 KEGG_enrich2:KEGG富集分析#---------
library(GO.db)
library(DOSE)
library(GO.db)
library(GSEABase)
library(ggtree)
library(aplot)
library(clusterProfiler)
library("GSVA")
library(dplyr)

# 调用差异分析的结果
res = EdgerSuper.metf (ps = ps,
                       group  = "Group",
                       artGroup = NULL)
dat = res[[2]]
dat
res2 = KEGG_enrich.metf(ps = ps,
            #  diffpath = diffpath,
            dif = dat
            )
dat1= res2$`KO-OE`
dat2= res2$`KO-WT`
dat3= res2$`OE-WT`


#38 buplot.metf：富集分析气泡图-----
Desep_group <- ps %>% sample_data() %>%
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
id
result = buplot.metf(dt  =dat1,dif = dif,id = id)
p1 = result[[1]]
p1
p2 = result[[2]]
p2

#39 GSVA.metf:GSVA富集分析--------
library("GSVA")
library(limma)

res= GSVA.metf(ps = ps, 
           # GSVApath = GSVApath,
           lg2FC = 0.5,
           padj = 0.05)
dat1 = res[[1]][[1]]
head(dat1)
p39.1= res[[1]][[2]]
p39.1

dat3 = res[[2]]$`KO-WT`
dat3 = res[[2]]$`KO-OE`
dat3 = res[[2]]$`OE-WT`
p39.2 = res[[3]]$`KO-WT`
p39.2 = res[[3]]$`OE-WT`

#40 kegg_function:按照kegg通路合并基因-------
ps.kegg.function = kegg_function(ps = ps)

#41 基于Mkegg通路合并基因-------
ps.mkegg.function =mkegg_function(ps = ps)

#42 基于reaction合并基因-------
library(EasyMultiOmics.db)
ps_kegg_function.reaction = reaction_function(ps = ps)




# cnps 循环数据库-----
# data prepration------
#1 cnps_gene: KEGG数据库中过滤CNPS 元素循环相关基因------
res  = cnps_gene(ps=ps.kegg)
res$C_gene    
res$N_gene 
res$P_gene  
res$S_gene 

#2 cnps_ps: KEGG数据库中过滤CNPS 元素循环相关基因构建ps -----
ps.cnps =cnps_ps(psko = ps.kegg)
ps = ps.cnps
tax_table(ps) %>% colnames()
otu_table(ps)

# function diversity -----
#3 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps,group = "Group",Plot = TRUE )
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = EasyStat::MuiKwWlx2(data = data,num = 3:5)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:5,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0


#4 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps))/10
result = alpha_rare.metf(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]]
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

#5 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps,
                       group = "Group",
                       dist = "bray",
                       method = "NMDS",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                       pair = FALSE)
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




#6 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

#7 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2


#8 mantal.metf：群落功能差异检测普鲁士分析#------

result <- mantal.metf(ps = ps,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

#9 cluster.metf:样品聚类#-----

res = cluster.metf(ps= ps,
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


#10 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)
tax_table(ps1)
res = Micro_tern.metf(ps=ps1,color = "Group"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# function classfication 

#11 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps,
                          j =  "module" ,
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Function","Abundance(%)")
head(databar)

#12 Ven.Upset.metf: 用于展示共有、特有的功能 ----
# 分组小于6时使用
library(dplyr)
res = Ven.Upset.metf(ps =  ps,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

#13 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps)$Group %>% table()
num =10
result = VenSuper.metf(ps = ps,
                       group = "Group",
                       num =10, 
                       j= "module" )

# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p7_1
#每个部分序列的数量占比，并作差异
dat<- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]+scale_fill_brewer(palette = "Paired")
p7_2



#14 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps ,
                      # rep = 1,
                      group = "ID",
                      start = 1, # 风车效果
                      m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                      a = 0.2, # 花瓣胖瘦
                      b = 1, # 花瓣距离花心的距离
                      lab.leaf = 1, # 花瓣标签到圆心的距离
                      col.cir = "yellow",
                      N = 0.5 )

p13.1 = res[[1]]
p13.1

dat = res[[2]]
dat

#15 ven.network.metf: ven网络展示共有特有功能----
map =  sample_data(ps)
result =ven.network.metf(
  ps = ps,
  N = 0.5,
  fill = "module")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat


# function differential analysis----
#16 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps,
                        group  = "Group",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

#17 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps,
                       group  = "Group",
                       artGroup = NULL)
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][1]
p16.2
p16.3 = res[[1]][1]
p16.3
dat = res[[2]]
dat

#18 t.metf: 差异分析t检验----
res = t.metf(ps = ps,group  = "Group",artGroup =NULL)
p17.1 = res [[1]][1]
p17.1
p17.2 = res [[1]][2]
p17.2
p17.3 = res [[1]][3]
p17.3

dat =  res [[2]]
dat

#19 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps,group  = "Group",artGroup =NULL)
p18.1 = res [[1]][1]
p18.1$`KO-OE_plot`
p18.2 = res [[1]][2]
p18.2
p18.3 = res [[1]][3]
p18.3

dat =  res [[2]]
dat %>% head()

#20 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2

# biobiomarker identification-----
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

#21 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#22 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat


#23 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#24 svm.metf:svm筛选特征微生物 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#25 glm.metf:glm筛选特征微生物----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance


#26 xgboost.metf: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#27 lasso.metf: lasso筛选特征微生物----
library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#28 decisiontree.micro: 错----
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


#29 naivebayes.metf: bayes筛选特征微生物----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#30 LDA.metf: LDA筛选特征微生物-----
tablda = LDA.metf(ps = pst,
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

#31 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps,group  = "Group", optimal = 50 ,fill = "module" )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

#32 nnet.metf: 神经网络筛选特征微生物  ------
library(nnet)
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#33 bagging.metf : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#34 CNPS.network:CNPS基因网络-------
library(data.table)
library(phyloseq)
library(ggClusterNet)
library(EasyStat)
library(tidyfst)
library(sna)
res = CNPS.network2(ps = ps.kegg,id.0 = "C")
p = res[[1]]
dat = res[[2]]
dat2 = res[[3]]

res = CNPS.network2(ps = ps.kegg,id.0 = "N")
p = res[[1]]
dat = res[[2]]
dat2 = res[[3]]

res = CNPS.network2(ps = ps.kegg,id.0 = "P")
p = res[[1]]
dat = res[[2]]
dat2 = res[[3]]
res = CNPS.network2(ps = ps.kegg,id.0 = "S")
p = res[[1]]
dat = res[[2]]
dat2 = res[[3]]









# card数据库------
data("ps.card")

# 数据过滤-----
#  通过比例和长度过滤--常规操作-文章常用
tax = ps.card %>% vegan_tax() %>% as.data.frame() %>% rownames_to_column("ID")
head(tax)

colnames(tax)[is.na(colnames(tax))] = "others"
tax$`Identity(%)` = as.numeric(tax$`Identity(%)`)
tax$Align_len = as.numeric(tax$Align_len)

tax = tax %>% 
  dplyr::filter(`Identity(%)` > 0.6,Align_len > 25)
dim(tax )

tax = tax %>% column_to_rownames("ID")

# tax$Drug.Class = tax$Drug_Class
# tax$Resistance.Mechanism = tax$Resistance_Mechanism
tax$ARO.Name1 = tax$ARO.Name

head(tax)
tax$Drug.Class[str_detect(tax$Drug_Class, "[;]")] = "muti-type"
tax_table(ps.card) = tax_table(as.matrix(tax))
ps = ps.card %>% filter_OTU_ps(Top = 1000)

ps

#--提取有多少个分组
Top = 20
gnum = sample_data(ps)$Group %>%unique() %>% length()
gnum
# 设定排序顺序--Group

axis_order = sample_data(ps)$Group %>%unique()

# 设定排序顺序--ID
map = sample_data(ps)
map$ID = row.names(map)
map = map %>% 
  as.tibble() %>%
  dplyr::arrange(desc(Group))
axis_order.s = map$ID


# function diversity -----

#1 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps.card,group = "Group",Plot = TRUE )
head(tab)


data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)
# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = EasyStat::MuiKwWlx2(data = data,num = 3:5)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:5,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0





#2 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps))/10
result = alpha_rare.metf(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)

#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3


#3 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps,
                       group = "Group",
                       dist = "bray",
                      method = "PCA",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                        pair = FALSE)
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




#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2



#6 mantal.metf：群落功能差异检测普鲁士分析#------

result <- mantal.meta(ps = ps,
                       method =  "spearman",
                       group = "Group",
                       ncol = gnum,
                       nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7




#7 cluster.metf:样品聚类#-----

res = cluster.metf(ps= ps,
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




#8 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)
res = Micro_tern.metf(ps1, fill = "Drug_Class" )
p15 = res[[1]]
p15
dat =  res[[2]]
dat


# function classfication 

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
ps.t = ps
otu = vegan_otu(ps.t) %>% t()
otu[otu > 0] = 1
otu_table(ps.t) = otu_table(otu,taxa_are_rows = TRUE)


result = barMainplot.metf(ps = ps.t,
                           j = "Drug_Class",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Drug_Class","Abundance(%)")
head(databar)


#10 Ven.Upset.metf: 用于展示共有、特有的功能----
# 分组小于6时使用
res = Ven.Upset.metf(ps =  ps,
                      group = "Group",
                      N = 0.5,
                      size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat


#11 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps)$Group %>% table()
num =10
result = VenSuper.metf(ps = ps,
                        group = "Group",
                       num =10 )

# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p7_1
#每个部分序列的数量占比，并作差异
dat<- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]+scale_fill_brewer(palette = "Paired")
p7_2



#13 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps ,
                       # rep = 1,
                       group = "ID",
                       start = 1, # 风车效果
                       m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                       a = 0.2, # 花瓣胖瘦
                       b = 1, # 花瓣距离花心的距离
                       lab.leaf = 1, # 花瓣标签到圆心的距离
                       col.cir = "yellow",
                       N = 0.5 )

p13.1 = res[[1]]

p13.1

dat = res[[2]]
dat






#14 ven.network.metf: ven网络展示共有特有功能----
map =  sample_data(ps)
result =ven.network.metf(
  ps = ps,
  N = 0.5,
  fill = "Drug_Class")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat





# function differential analysis----
# 数据转换为"Drug_Class" 或其他
ps.g = ps %>%
  tax_glom_meta(ranks = "Drug_Class")

#15 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps.g,
                        group  = "Drug_Class",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

#16 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps.g,
                 group  = "Group",
                 artGroup = NULL)
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][1]
p16.2
p16.3 = res[[1]][1]
p16.3
dat = res[[2]]
dat

#17 t.metf: 差异分析t检验----

res = t.metf(ps = ps.g,group  = "Group",artGroup =NULL)
p17.1 = res [[1]][1]
p17.1
p17.2 = res [[1]][2]
p17.2
p17.3 = res [[1]][3]
p17.3

dat =  res [[2]]
dat

#18 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps.g,group  = "Group",artGroup =NULL)
p18.1 = res [[1]][1]
p18.1$`KO-OE_plot`
p18.2 = res [[1]][2]
p18.2
p18.3 = res [[1]][3]
p18.3

dat =  res [[2]]
dat %>% head()

#19 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2


# biomarker identification-----
#20 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable


#21 Roc.metf:ROC 曲线绘制----
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat


#22 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#23 svm.metf:svm筛选特征微生物 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#24 glm.metf:glm筛选特征微生物----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance






#25 xgboost.metf: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance



#26 lasso.metf: lasso筛选特征微生物----

library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#27 decisiontree.micro: 错----
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


#28 naivebayes.metf: bayes筛选特征微生物----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#29 LDA.metf: LDA筛选特征微生物-----
tablda = LDA.metf(ps = pst,
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

#30 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps,group  = "Group", optimal = 50 ,fill ="Drug_Class" )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4


#31 nnet.metf: 神经网络筛选特征微生物  ------
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#32 bagging.metf : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance



#network analysis -----
#33 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
tab.r = network.pip(
  ps = ps,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label = FALSE,
  lab = "elements",
  group = "Group",
  #fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1,
  fill = "Resistance_Mechanism"
)

dat = tab.r[[2]]

cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]
p0
#34 net_properties.4:网络属性计算#---------
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

#35 netproperties.sample:单个样本的网络属性#-------
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

#36 node_properties:计算节点属性#---------
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



#37 module.compare.net.pip:网络显著性比较#-----
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


# 







# cazy数据库-----
data("ps.cazy")

ps =ps.cazy %>% filter_OTU_ps(Top = 1000)



# function diversity -----
#1 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps,group = "Group",Plot = TRUE )
head(tab)


data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)
# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = EasyStat::MuiKwWlx2(data = data,num = 3:5)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:5,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0


#2 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps))/10
result = alpha_rare.metf(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]]
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

#3 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps,
                       group = "Group",
                       dist = "bray",
                       method = "PCA",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                       pair = FALSE)
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




#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2


#6 mantal.metf：群落功能差异检测普鲁士分析#------

result <- mantal.metf(ps = ps,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

#7 cluster.metf:样品聚类#-----

res = cluster.metf(ps= ps,
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


#8 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)
tax_table(ps1)
res = Micro_tern.metf(ps=ps1,color = "Class_description"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# function classfication 

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps,
                          j =  "Class_description" ,
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Pathway","Abundance(%)")
head(databar)

#10 Ven.Upset.metf: 用于展示共有、特有的功能 ----
# 分组小于6时使用
library(dplyr)
res = Ven.Upset.metf(ps =  ps,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

#11 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps)$Group %>% table()
num =10
result = VenSuper.metf(ps = ps,
                       group = "Group",
                       num =10, 
                       j= "Class" )

# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p7_1
#每个部分序列的数量占比，并作差异
dat<- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]+scale_fill_brewer(palette = "Paired")
p7_2



#12 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps ,
                      # rep = 1,
                      group = "ID",
                      start = 1, # 风车效果
                      m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                      a = 0.2, # 花瓣胖瘦
                      b = 1, # 花瓣距离花心的距离
                      lab.leaf = 1, # 花瓣标签到圆心的距离
                      col.cir = "yellow",
                      N = 0.5 )

p13.1 = res[[1]]
p13.1

dat = res[[2]]
dat

#13 ven.network.metf: ven网络展示共有特有功能----
map =  sample_data(ps)
result =ven.network.metf(
  ps = ps,
  N = 0.5,
  fill = "Class")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat


# function differential analysis----
#14 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps,
                        group  = "Group",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

#15 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps,
                       group  = "Group",
                       artGroup = NULL)
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][1]
p16.2
p16.3 = res[[1]][1]
p16.3
dat = res[[2]]
dat

#16 t.metf: 差异分析t检验----
res = t.metf(ps = ps,group  = "Group",artGroup =NULL)
p17.1 = res [[1]][1]
p17.1
p17.2 = res [[1]][2]
p17.2
p17.3 = res [[1]][3]
p17.3

dat =  res [[2]]
dat

#17 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps,group  = "Group",artGroup =NULL)
p18.1 = res [[1]][1]
p18.1$`KO-OE_plot`
p18.2 = res [[1]][2]
p18.2
p18.3 = res [[1]][3]
p18.3

dat =  res [[2]]
dat %>% head()

#18 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2








# biobiomarker identification-----
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)
#19 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable


#20 Roc.metf:ROC 曲线绘制----


res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat


#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#22 svm.metf:svm筛选特征微生物 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#23 glm.metf:glm筛选特征微生物----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance


#24 xgboost.metf: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance



#25 lasso.metf: lasso筛选特征微生物----

library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#26 decisiontree.micro: 错----
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


#27 naivebayes.metf: bayes筛选特征微生物----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#28 LDA.metf: LDA筛选特征微生物-----
tablda = LDA.metf(ps = pst,
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

#29 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps,group  = "Group", optimal = 50 ,fill ="Class" )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4


#30 nnet.metf: 神经网络筛选特征微生物  ------
library(nnet)
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#31 bagging.metf : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#network analysis -----
#32 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
tab.r = network.pip(
  ps = ps,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label = FALSE,
  lab = "elements",
  group = "Group",
  #fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1,
  fill = "Class" 
)

dat = tab.r[[2]]

cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]
p0
#33 net_properties.4:网络属性计算#---------
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

#34 netproperties.sample:单个样本的网络属性#-------
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

#35 node_properties:计算节点属性#---------
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



#36 module.compare.net.pip:网络显著性比较#-----
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










# cog数据库-----
data("ps.cog")

ps = ps.cog %>% filter_OTU_ps(Top = 1000)

tax_table(ps)  
# Function
# function diversity -----
#1 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps,group = "Group",Plot = TRUE )
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = EasyStat::MuiKwWlx2(data = data,num = 3:5)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:5,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0


#2 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps))/10
result = alpha_rare.metf(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]]
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

#3 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps,
                       group = "Group",
                       dist = "bray",
                       method = "PCA",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                       pair = FALSE)
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




#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2


#6 mantal.metf：群落功能差异检测普鲁士分析#------

result <- mantal.metf(ps = ps,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

#7 cluster.metf:样品聚类#-----

res = cluster.metf(ps= ps,
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


#8 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)
tax_table(ps1)
res = Micro_tern.metf(ps=ps1,color = "Function"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# function classfication 

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps,
                          j =  "Function" ,
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Function","Abundance(%)")
head(databar)

#10 Ven.Upset.metf: 用于展示共有、特有的功能 ----
# 分组小于6时使用
library(dplyr)
res = Ven.Upset.metf(ps =  ps,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

#11 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps)$Group %>% table()
num =10
result = VenSuper.metf(ps = ps,
                       group = "Group",
                       num =10, 
                       j= "Function" )

# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p7_1
#每个部分序列的数量占比，并作差异
dat<- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]+scale_fill_brewer(palette = "Paired")
p7_2



#12 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps ,
                      # rep = 1,
                      group = "ID",
                      start = 1, # 风车效果
                      m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                      a = 0.2, # 花瓣胖瘦
                      b = 1, # 花瓣距离花心的距离
                      lab.leaf = 1, # 花瓣标签到圆心的距离
                      col.cir = "yellow",
                      N = 0.5 )

p13.1 = res[[1]]
p13.1

dat = res[[2]]
dat

#13 ven.network.metf: ven网络展示共有特有功能----
map =  sample_data(ps)
result =ven.network.metf(
  ps = ps,
  N = 0.5,
  fill = "Function")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat


# function differential analysis----
#14 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps,
                        group  = "Group",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

#15 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps,
                       group  = "Group",
                       artGroup = NULL)
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][1]
p16.2
p16.3 = res[[1]][1]
p16.3
dat = res[[2]]
dat

#16 t.metf: 差异分析t检验----
res = t.metf(ps = ps,group  = "Group",artGroup =NULL)
p17.1 = res [[1]][1]
p17.1
p17.2 = res [[1]][2]
p17.2
p17.3 = res [[1]][3]
p17.3

dat =  res [[2]]
dat

#17 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps,group  = "Group",artGroup =NULL)
p18.1 = res [[1]][1]
p18.1$`KO-OE_plot`
p18.2 = res [[1]][2]
p18.2
p18.3 = res [[1]][3]
p18.3

dat =  res [[2]]
dat %>% head()

#18 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2

# biobiomarker identification-----
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

#19 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#20 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat


#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#22 svm.metf:svm筛选特征微生物 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#23 glm.metf:glm筛选特征微生物----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance


#24 xgboost.metf: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#25 lasso.metf: lasso筛选特征微生物----
library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#26 decisiontree.micro: 错----
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


#27 naivebayes.metf: bayes筛选特征微生物----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#28 LDA.metf: LDA筛选特征微生物-----
tablda = LDA.metf(ps = pst,
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

#29 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps,group  = "Group", optimal = 50 ,fill ="Function"  )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

#30 nnet.metf: 神经网络筛选特征微生物  ------
library(nnet)
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#31 bagging.metf : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#network analysis -----
#32 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
tab.r = network.pip(
  ps = ps,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label = FALSE,
  lab = "elements",
  group = "Group",
  #fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1,
  fill = "Function" 
)

dat = tab.r[[2]]

cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]
p0
#33 net_properties.4:网络属性计算#---------
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

#34 netproperties.sample:单个样本的网络属性#-------
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

#35 node_properties:计算节点属性#---------
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



#36 module.compare.net.pip:网络显著性比较#-----
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


# vfdb 数据库-----
data("ps.vfdb")

ps = ps.vfdb %>% filter_OTU_ps(Top = 1000)

tax_table(ps)  
# Function
# function diversity -----
#1 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps,group = "Group",Plot = TRUE )
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = EasyStat::MuiKwWlx2(data = data,num = 3:5)
result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = 3:5,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)
res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+ 
  scale_x_discrete(limits = axis_order) + 
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0


#2 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps))/10
result = alpha_rare.metf(ps = ps, group = "Group", method = "Richness", start = 100, step = rare)
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]]
#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

#3 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps,
                       group = "Group",
                       dist = "bray",
                       method = "PCA",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                       pair = FALSE)
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




#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2


#6 mantal.metf：群落功能差异检测普鲁士分析#------

result <- mantal.metf(ps = ps,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

#7 cluster.metf:样品聚类#-----

res = cluster.metf(ps= ps,
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


#8 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)
tax_table(ps1)
res = Micro_tern.metf(ps=ps1,color = "Function"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# function classfication 

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps,
                          j =  "Function" ,
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Function","Abundance(%)")
head(databar)

#10 Ven.Upset.metf: 用于展示共有、特有的功能 ----
# 分组小于6时使用
library(dplyr)
res = Ven.Upset.metf(ps =  ps,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

#11 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps)$Group %>% table()
num =10
result = VenSuper.metf(ps = ps,
                       group = "Group",
                       num =10, 
                       j= "Function" )

# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p7_1
#每个部分序列的数量占比，并作差异
dat<- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]+scale_fill_brewer(palette = "Paired")
p7_2



#12 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps ,
                      # rep = 1,
                      group = "ID",
                      start = 1, # 风车效果
                      m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                      a = 0.2, # 花瓣胖瘦
                      b = 1, # 花瓣距离花心的距离
                      lab.leaf = 1, # 花瓣标签到圆心的距离
                      col.cir = "yellow",
                      N = 0.5 )

p13.1 = res[[1]]
p13.1

dat = res[[2]]
dat

#13 ven.network.metf: ven网络展示共有特有功能----
map =  sample_data(ps)
result =ven.network.metf(
  ps = ps,
  N = 0.5,
  fill = "Function")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat


# function differential analysis----
#14 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps,
                        group  = "Group",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

#15 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps,
                       group  = "Group",
                       artGroup = NULL)
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][1]
p16.2
p16.3 = res[[1]][1]
p16.3
dat = res[[2]]
dat

#16 t.metf: 差异分析t检验----
res = t.metf(ps = ps,group  = "Group",artGroup =NULL)
p17.1 = res [[1]][1]
p17.1
p17.2 = res [[1]][2]
p17.2
p17.3 = res [[1]][3]
p17.3

dat =  res [[2]]
dat

#17 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps,group  = "Group",artGroup =NULL)
p18.1 = res [[1]][1]
p18.1$`KO-OE_plot`
p18.2 = res [[1]][2]
p18.2
p18.3 = res [[1]][3]
p18.3

dat =  res [[2]]
dat %>% head()

#18 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2

# biobiomarker identification-----
id = sample_data(ps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

#19 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

#20 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat


#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

#22 svm.metf:svm筛选特征微生物 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

#23 glm.metf:glm筛选特征微生物----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance


#24 xgboost.metf: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#25 lasso.metf: lasso筛选特征微生物----
library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance


#26 decisiontree.micro: 错----
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


#27 naivebayes.metf: bayes筛选特征微生物----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#28 LDA.metf: LDA筛选特征微生物-----
tablda = LDA.metf(ps = pst,
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

#29 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps,group  = "Group", optimal = 50 ,fill ="Function"  )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

#30 nnet.metf: 神经网络筛选特征微生物  ------
library(nnet)
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#31 bagging.metf : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

#network analysis -----
#32 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
tab.r = network.pip(
  ps = ps,
  N = 200,
  # ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  # method = "sparcc",
  label = FALSE,
  lab = "elements",
  group = "Group",
  #fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1,
  fill = "Function" 
)

dat = tab.r[[2]]

cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]
p0
#33 net_properties.4:网络属性计算#---------
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

#34 netproperties.sample:单个样本的网络属性#-------
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

#35 node_properties:计算节点属性#---------
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



#36 module.compare.net.pip:网络显著性比较#-----
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














