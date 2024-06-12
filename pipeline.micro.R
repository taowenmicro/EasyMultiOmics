

# 基于整理好的函数，重新写流程

library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)


# 扩增子微生物组学分析#-----
# alpha.micro: 6中alpha多样性计算#----
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

#alpha.pd :用于计算pd多样性#-------
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


#  alpha.rare.line:alpha多样性稀释曲线#---------
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


# ordinate.micro: 排序分析#----------
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

# MicroTest:群落水平差异检测#-------
dat1 = MicroTest(ps = ps, Micromet = "adonis", dist = "bray")
dat1
# pairMicroTest:两两分组群落水平差异检测#-------
dat2 = pairMicroTest(ps = ps, Micromet = "MRPP", dist = "bray")
dat2

# mantal.micro ：群落差异检测普鲁士分析#------
map= sample_data(ps)
head(map)

result <- mantal.micro(ps = ps,
                       method =  "spearman",
                       group = "Group",
                       ncol = 3,
                       nrow = 1)

data <- result[[1]]

p3_7 <- result[[2]]
p3_7

# distance.micro:分组之间距离比对#-----
res = distance.micro(ps,group = "Group")
p4.1 = res[[1]]
p4.1
p4.2 = res[[2]]
p4.2
p4.3 = res[[3]]
p4.3



