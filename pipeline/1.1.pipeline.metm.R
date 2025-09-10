# 宏基因组物种分析#-----
rm(list=ls())

# BiocManager::install("MicrobiotaProcess")
library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(openxlsx)

# 有效基因数量
sample_sums(ps.micro)

## 参数设定-----
map= sample_data(ps.micro)
head(map)

phyloseq::tax_table(ps.micro) %>% head()
# 提取分组因子数量
gnum = phyloseq::sample_data(ps.micro)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.micro对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps.micro)$Group %>%unique();axis_order
col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")

#-主题--
package.amp()

res = theme_my(ps.micro)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]

# 保存路径设置---
path =  "../result/metm/"
fs::dir_create(path)

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

result = MuiKwWlx2(data = data,num = 3:6)

result1 = FacetMuiPlotresultBox(data = data,num = 3:6,
                                          result = result,
                                          sig_show ="abc",ncol = 4,width = 0.4 )


p1_1 = result1[[1]] +scale_fill_manual(values = col.g)+
  ggplot2::scale_x_discrete(limits = axis_order) +
  # theme_cell()+
  theme_nature()+
  ggplot2::guides(fill = guide_legend(title = none))

p1_1

res = FacetMuiPlotresultBar(data = data,num = c(3:6),
                                      result = result,sig_show ="abc",ncol = 4,
                            mult.y = 0.3
                            )
p1_2 = res[[1]]+
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)

p = res[[1]] +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_cell() +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)


res = FacetMuiPlotReBoxBar(data = data,num = c(3:6),result = result,sig_show ="abc",ncol = 4,
                           mult.y = 0.3,
                           lab.yloc = 1.1
                           )
p1_3 = res[[1]]+
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
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
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = none)) +
  ggplot2::scale_fill_manual(values = col.g)+
  ggplot2::scale_color_manual(values = col.g)

# 保存路径
alppath <- file.path(path, "alpha")
dir.create(alppath, showWarnings = FALSE)

# 保存Alpha多样性图片
ggsave(file.path(alppath, "alpha_diversity_boxplot.png"), plot = p1_1, width = 12, height = 6)
ggsave(file.path(alppath, "alpha_diversity_boxplot.pdf"), plot = p1_1, width = 12, height = 6)
ggsave(file.path(alppath, "alpha_diversity_barplot.png"), plot = p1_2, width = 12, height = 6)
ggsave(file.path(alppath, "alpha_diversity_barplot.pdf"), plot = p1_2, width = 12, height = 6)
ggsave(file.path(alppath, "alpha_diversity_reboxbar.png"), plot = p1_3, width = 12, height = 6)
ggsave(file.path(alppath, "alpha_diversity_reboxbar.pdf"), plot = p1_3, width = 12, height = 6)
ggsave(file.path(alppath, "alpha_diversity_violinplot.png"), plot = p1_0, width = 12, height = 6)
ggsave(file.path(alppath, "alpha_diversity_violinplot.pdf"), plot = p1_0, width = 12, height = 6)


#3 alpha_rare.metm:alpha多样性稀释曲线#---------
library(microbiome)
library(vegan)

rare <- mean(phyloseq::sample_sums(ps.micro))/10

result = alpha_rare.metm(ps = ps.micro ,
                          group = "Group", method = "Richness", start = 100, step = rare)

#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +scale_color_manual(values = col.g)+
  theme_nature()
p2_1
## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]]+
  theme_nature()
p2_2
#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]]+
  theme_nature()
p2_3

# 保存路径
alppath <- file.path(path, "alpha")
dir.create(alppath, showWarnings = FALSE)

# 保存稀释曲线图片
ggsave(file.path(alppath, "alpha_rarefaction_curve.png"), plot = p2_1, width = 10, height = 8)
ggsave(file.path(alppath, "alpha_rarefaction_curve.pdf"), plot = p2_1, width = 10, height = 8)
ggsave(file.path(alppath, "alpha_rarefaction_group_curve.png"), plot = p2_2, width = 10, height = 8)
ggsave(file.path(alppath, "alpha_rarefaction_group_curve.pdf"), plot = p2_2, width = 10, height = 8)
ggsave(file.path(alppath, "alpha_rarefaction_group_sd_curve.png"), plot = p2_3, width = 10, height = 8)
ggsave(file.path(alppath, "alpha_rarefaction_group_sd_curve.pdf"), plot = p2_3, width = 10, height = 8)


# beta diversity  -----
#4 ordinate.metm: 排序分析#----------
result = ordinate.metm(ps = ps.micro, group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05)
p3_1 = result[[1]]
p3_1 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = "none") +
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

#带标签图形出图
p3_2 = result[[3]]
p3_2 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


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
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F)+
  theme_cell()+
  theme(axis.title.y = element_text(angle = 90))

# 创建beta多样性主目录
betapath <- file.path(path, "beta")
dir.create(betapath, showWarnings = FALSE, recursive = TRUE)

# 保存排序分析图片
ggsave(file.path(betapath, "pcoa_basic.png"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(betapath, "pcoa_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(betapath, "pcoa_labeled.png"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(betapath, "pcoa_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(betapath, "pcoa_refined.png"), plot = p3_3, width = 10, height = 8)
ggsave(file.path(betapath, "pcoa_refined.pdf"), plot = p3_3, width = 10, height = 8)

# 保存排序分析数据
beta_wb <- createWorkbook()
addWorksheet(beta_wb, "plotdata")
addWorksheet(beta_wb, "cent")
addWorksheet(beta_wb, "segs")
writeData(beta_wb, "plotdata", plotdata)
writeData(beta_wb, "cent", cent)
writeData(beta_wb, "segs", segs)


#5 MicroTest.metm:群落水平差异检测-------
dat1 = MicroTest.metm(ps = ps.micro, Micromet = "adonis", dist = "bray")
dat1

# 保存数据
addWorksheet(beta_wb, "adonis_results")
writeData(beta_wb, "adonis_results", dat1)

#6 pairMicroTest.metm:两两分组群落水平差异检测-------
dat2 = pairMicroTest.metm(ps = ps.micro, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较数据
addWorksheet(beta_wb, "pairwise_MRPP")
writeData(beta_wb, "pairwise_MRPP", dat2)


#7 mantal.metm：群落差异检测普鲁士分析#------
result <- mantal.metm(ps = ps.micro,
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
  theme_cell

# 保存Mantel分析图片
ggsave(file.path(mantelpath, "mantel_plot.png"), plot = p3_7_final, width = 10, height = 8)
ggsave(file.path(mantelpath, "mantel_plot.pdf"), plot = p3_7_final, width = 10, height = 8)

# 保存Mantel分析数据
addWorksheet(beta_wb, "mantel_results")
writeData(beta_wb, "mantel_results", mantel_data)
saveWorkbook(beta_wb, file.path(betapath, "beta_results.xlsx"), overwrite = TRUE)


#8 cluster.metm:样品聚类-----
res = cluster.metm (ps= ps.micro,
                     hcluter_method = "complete",
                     dist = "bray",
                     cuttree = 3,
                     row_cluster = TRUE,
                     col_cluster =  TRUE)

p4 = res[[1]] +
  scale_fill_manual(values = col.g)

p4
p4_1 = res[[2]]
p4_1
p4_2 = res[[3]]
p4_2
dat = res[[4]] # 聚类矩阵
head(dat)

# 创建聚类分析保存目录
clusterpath <- file.path(path, "cluster")
dir.create(clusterpath, showWarnings = FALSE, recursive = TRUE)

# 保存聚类分析图片
ggsave(file.path(clusterpath, "cluster_heatmap.png"), plot = p4, width = 10, height = 8)
ggsave(file.path(clusterpath, "cluster_heatmap.pdf"), plot = p4, width = 10, height = 8)
ggsave(file.path(clusterpath, "cluster_dendrogram.png"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(clusterpath, "cluster_dendrogram.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(clusterpath, "cluster_scatter.png"), plot = p4_2, width = 10, height = 8)
ggsave(file.path(clusterpath, "cluster_scatter.pdf"), plot = p4_2, width = 10, height = 8)

# 保存聚类数据
cluster_wb <- createWorkbook()
addWorksheet(cluster_wb, "cluster_matrix")
writeData(cluster_wb, "cluster_matrix", dat, rowNames = TRUE)
saveWorkbook(cluster_wb, file.path(clusterpath, "cluster_results.xlsx"), overwrite = TRUE)

# compositipn -----
# 创建组成分析保存目录
comppath <- file.path(path, "composition")
dir.create(comppath, showWarnings = FALSE, recursive = TRUE)

#9 Ven.Upset.metm: 用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = Ven.Upset.metm(ps =  ps.micro,
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

# 保存Venn/Upset图片
ggsave(file.path(comppath, "venn_diagram.png"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(comppath, "venn_diagram.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(comppath, "upset_plot.png"), plot = p10.2, width = 10, height = 8)
ggsave(file.path(comppath, "upset_plot.pdf"), plot = p10.2, width = 10, height = 8)
ggsave(file.path(comppath, "venn_additional.png"), plot = p10.3, width = 10, height = 8)
ggsave(file.path(comppath, "venn_additional.pdf"), plot = p10.3, width = 10, height = 8)

# 保存Venn/Upset数据
comp_wb <- createWorkbook()
addWorksheet(comp_wb, "venn_upset_data")
writeData(comp_wb, "venn_upset_data", dat, rowNames = TRUE)

#10 ggVen.Upset.metm:用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = ggVen.Upset.metm(ps = ps.micro,group = "Group")

# 保存ggVenn图片 (grid对象需要特殊处理)
png(file.path(comppath, "ggvenn_plot.png"), width = 10, height = 8, units = "in", res = 300)
grid::grid.draw(res[[1]])
dev.off()

pdf(file.path(comppath, "ggvenn_plot.pdf"), width = 10, height = 8)
grid::grid.draw(res[[1]])
dev.off()

dat = res[[2]]
dat

# 保存ggVenn数据
addWorksheet(comp_wb, "ggvenn_data")
writeData(comp_wb, "ggvenn_data", dat, rowNames = TRUE)

#11 VenSeper.micro: 详细展示每一组中OTU的物种组成 错-------
#---每个部分
library("ggpubr")
library(agricolae)
library(reshape2)

result = VenSuper.metm(ps = ps.micro,
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

# 保存VenSuper图片
ggsave(file.path(comppath, "vensuper_barplot.png"), plot = p7_1, width = 10, height = 8)
ggsave(file.path(comppath, "vensuper_barplot.pdf"), plot = p7_1, width = 10, height = 8)
ggsave(file.path(comppath, "vensuper_alluvial.png"), plot = p7_2, width = 10, height = 8)
ggsave(file.path(comppath, "vensuper_alluvial.pdf"), plot = p7_2, width = 10, height = 8)


#12 ggflower.metm:花瓣图展示共有特有微生物------
res <- ggflower.metm(ps = ps.micro ,
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

# 保存花瓣图图片
ggsave(file.path(comppath, "flower_plot.png"), plot = p13.1, width = 10, height = 8)
ggsave(file.path(comppath, "flower_plot.pdf"), plot = p13.1, width = 10, height = 8)

# 保存花瓣图数据
addWorksheet(comp_wb, "flower_data")
writeData(comp_wb, "flower_data", dat, rowNames = TRUE)

#13 ven.network.metm:venn 网络----
result = ven.network.metm(
  ps = ps.micro ,#%>%  filter_OTU_ps(10000),
  N = 0.5,
  fill = "Phylum")
p14  = result[[1]]
p14

dat = result[[2]]
head(dat)

# 保存Venn网络图片
ggsave(file.path(comppath, "venn_network.png"), plot = p14, width = 10, height = 8)
ggsave(file.path(comppath, "venn_network.pdf"), plot = p14, width = 10, height = 8)

# 保存Venn网络数据（保留行名）
addWorksheet(comp_wb, "venn_network_data")
writeData(comp_wb, "venn_network_data", dat, rowNames = TRUE)

#14 Micro_tern.metm: 三元图展示组成----
ps1 = ps.micro %>% filter_OTU_ps(500)
res = Micro_tern.metm(ps1)
p15 = res[[1]]
p15[[1]] +theme_classic()

dat =  res[[2]]
head(dat)

# 保存三元图图片
ggsave(file.path(comppath, "ternary_plot.png"), plot = p15[[1]] + theme_classic(), width = 10, height = 8)
ggsave(file.path(comppath, "ternary_plot.pdf"), plot = p15[[1]] + theme_classic(), width = 10, height = 8)

# 保存三元图数据（保留行名）
addWorksheet(comp_wb, "ternary_data")
writeData(comp_wb, "ternary_data", dat, rowNames = TRUE)

#15 barMainplot.metm: 堆积柱状图展示组成----
library(ggalluvial)
# pst = ps.micro %>% subset_taxa.wt("Species","Unassigned",TRUE)
pst = ps.micro %>% subset_taxa.wt("Genus","Unassigned",TRUE)

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
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))


p4_2  <- result[[3]]
# # scale_fill_brewer(palette = "Paired") +
# scale_fill_manual(values = colset2) +
# scale_x_discrete(limits = axis_order) +
# mytheme1
p4_2+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","genus","Abundance(%)")
head(databar)

# 保存堆积柱状图图片
ggsave(file.path(comppath, "barplot_samples.png"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(comppath, "barplot_samples.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(comppath, "barplot_groups.png"), plot = p4_2, width = 10, height = 8)
ggsave(file.path(comppath, "barplot_groups.pdf"), plot = p4_2, width = 10, height = 8)

# 保存堆积柱状图数据
addWorksheet(comp_wb, "barplot_raw_data")
addWorksheet(comp_wb, "barplot_summary_data")
writeData(comp_wb, "barplot_raw_data", result[[2]], rowNames = TRUE)
writeData(comp_wb, "barplot_summary_data", databar)  # 汇总数据通常不需要行名


#16 cluMicro.bar.metm: 聚类堆积柱状图展示组成-----
result <-  cluMicro.bar.metm (dist = "bray",
                               ps = ps.micro,
                               j = "Genus",
                               Top = 7, # 提取丰度前十的物种注释
                               tran = TRUE, # 转化为相对丰度值
                               hcluter_method = "complete",
                               Group = "Group",
                               cuttree = length(unique(phyloseq::sample_data(ps.micro)$Group)))

p5_2 <- result[[2]]
p5_2

p5_3 <- result[[3]]
p5_3
p5_4 <- result[[4]]
p5_4
clubardata <- result[[5]]

# 保存聚类堆积柱状图图片
ggsave(file.path(comppath, "cluster_barplot_1.png"), plot = p5_2, width = 10, height = 8)
ggsave(file.path(comppath, "cluster_barplot_1.pdf"), plot = p5_2, width = 10, height = 8)
ggsave(file.path(comppath, "cluster_barplot_2.png"), plot = p5_3, width = 10, height = 8)
ggsave(file.path(comppath, "cluster_barplot_2.pdf"), plot = p5_3, width = 10, height = 8)
ggsave(file.path(comppath, "cluster_barplot_3.png"), plot = p5_4, width = 10, height = 8)
ggsave(file.path(comppath, "cluster_barplot_3.pdf"), plot = p5_4, width = 10, height = 8)

# 保存聚类堆积柱状图数据
addWorksheet(comp_wb, "cluster_barplot_data")
writeData(comp_wb, "cluster_barplot_data", clubardata, rowNames = TRUE)

#17 cir_barplot.metf:环状堆积柱状图 -----
library(ggtree) # j = "Phylum"
res = cir_barplot.metm(
  ps = ps.micro,
  Top = 7,
  dist = "bray",
  cuttree = 3,
  hcluter_method = "complete")

p17 = res[[1]]
p17
dat= res[[2]]
head(dat)

# 保存环状堆积柱状图图片
ggsave(file.path(comppath, "circular_barplot.png"), plot = p17, width = 10, height = 8)
ggsave(file.path(comppath, "circular_barplot.pdf"), plot = p17, width = 10, height = 8)

# 保存环状堆积柱状图数据（保留行名）
addWorksheet(comp_wb, "circular_barplot_data")
writeData(comp_wb, "circular_barplot_data", dat, rowNames = TRUE)

#18 cir_plot.metm:和弦图展示物种组成-----
res = cir_plot.metm(ps  = ps.micro,Top = 12,rank = 6)
p <- recordPlot()

# 保存为 PNG 文件
png(file.path(comppath, "cir_plot.png"), width = 10, height = 8, units = "in", res = 300)
replayPlot(p)
dev.off()

# 保存为 PDF 文件
pdf(file.path(comppath, "cir_plot.pdf"), width = 10, height = 8)
replayPlot(p)
dev.off()

# 保存和弦图数据
addWorksheet(comp_wb, "cir_plot_data")
writeData(comp_wb, "cir_plot_data", res[[1]], rowNames = TRUE)

#19 Microheatmap.metm: 热图展示物种相对丰度差异-----
heatnum = 30
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
# phyloseq::sample_data(ps) = map

ps_tem = ps.micro %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus")
id <- ps.micro %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus") %>%
  ggClusterNet::filter_OTU_ps(100) %>%
  ggClusterNet::vegan_otu() %>%
  t() %>% as.data.frame() %>%rowCV %>%
  sort(decreasing = TRUE) %>%
  head(heatnum) %>%
  names()
#? Microheatmap.metm
result <- Microheatmap.metm(ps_rela = ps_tem,id = id ,col_cluster = FALSE,
                            row_cluster = FALSE,
                            col1 = (ggsci::pal_gsea(alpha = 1))(12))

p24.1 <- result[[1]]

p24.1
# p1 +  scale_fill_gradientn(colours =colorRampPalette(RColorBrewer::brewer.pal(11,"Set3"))(60))
p24.2 <- result[[2]]
p24.2
dat = result[[3]]
head(dat)

# 保存微生物热图图片
ggsave(file.path(comppath, "microbiome_heatmap_1.png"), plot = p24.1, width = 10, height = 8)
ggsave(file.path(comppath, "microbiome_heatmap_1.pdf"), plot = p24.1, width = 10, height = 8)
ggsave(file.path(comppath, "microbiome_heatmap_2.png"), plot = p24.2, width = 10, height = 8)
ggsave(file.path(comppath, "microbiome_heatmap_2.pdf"), plot = p24.2, width = 10, height = 8)

# 保存微生物热图数据（保留行名）
addWorksheet(comp_wb, "heatmap_data")
addWorksheet(comp_wb, "selected_otus")
writeData(comp_wb, "heatmap_data", dat, rowNames = TRUE)
writeData(comp_wb, "selected_otus", data.frame(OTU_ID = id))
saveWorkbook(comp_wb, file.path(comppath, "comp_results.xlsx"), overwrite = TRUE)

# difference analysis -----
# 创建差异分析主目录
diffpath <- file.path(path, "differential")
dir.create(diffpath, showWarnings = FALSE, recursive = TRUE)
#20 EdgerSuper.metm:EdgeR计算差异微生物----
res = EdgerSuper.metm(ps = ps.micro %>% ggClusterNet::filter_OTU_ps(500),group  = "Group",artGroup = NULL, j = "Species")

p25.1 =  res[[1]][1]
p25.1
p25.2 =  res[[1]][2]
p25.2
p25.3 =  res[[1]][3]
p25.3
dat =  res[[2]]
head(dat)

# 保存EdgeR图片
ggsave(file.path(diffpath, "edger_volcano_1.png"), plot = p25.1, width = 10, height = 8)
ggsave(file.path(diffpath, "edger_volcano_1.pdf"), plot = p25.1, width = 10, height = 8)
ggsave(file.path(diffpath, "edger_volcano_2.png"), plot = p25.2, width = 10, height = 8)
ggsave(file.path(diffpath, "edger_volcano_2.pdf"), plot = p25.2, width = 10, height = 8)
ggsave(file.path(diffpath, "edger_volcano_3.png"), plot = p25.3, width = 10, height = 8)
ggsave(file.path(diffpath, "edger_volcano_3.pdf"), plot = p25.3, width = 10, height = 8)

# 保存EdgeR数据
diff_wb <- createWorkbook()
addWorksheet(diff_wb, "edger_results")
writeData(diff_wb, "edger_results", dat, rowNames = TRUE)

#21 EdgerSuper2.metm:EdgeR计算差异微生物-----
res =  EdgerSuper2.metm (ps = ps.micro,group  = "Group",artGroup =NULL, j = "Species")
head(res)

# 保存EdgeR2数据
addWorksheet(diff_wb, "edger2_results")
writeData(diff_wb, "edger2_results", res, rowNames = TRUE)

#22 DESep2Super.metm:DESep2计算差异微生物-----
res = DESep2Super.metm(ps = ps.micro %>% ggClusterNet::filter_OTU_ps(500),
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

# 保存DESeq2图片
ggsave(file.path(diffpath, "deseq2_volcano_1.png"), plot = p26.1[[1]], width = 10, height = 8)
ggsave(file.path(diffpath, "deseq2_volcano_1.pdf"), plot = p26.1[[1]], width = 10, height = 8)
ggsave(file.path(diffpath, "deseq2_volcano_2.png"), plot = p26.2[[1]], width = 10, height = 8)
ggsave(file.path(diffpath, "deseq2_volcano_2.pdf"), plot = p26.2[[1]], width = 10, height = 8)
ggsave(file.path(diffpath, "deseq2_volcano_3.png"), plot = p26.3[[1]], width = 10, height = 8)
ggsave(file.path(diffpath, "deseq2_volcano_3.pdf"), plot = p26.3[[1]], width = 10, height = 8)

# 保存DESeq2数据
addWorksheet(diff_wb, "deseq2_results")
writeData(diff_wb, "deseq2_results", dat, rowNames = TRUE)

#23 edge_Manhattan.metm: 曼哈顿图展示差异微生物------

res = edge_Manhattan.metm(
  ps = ps.micro%>% ggClusterNet::filter_OTU_ps(500),
  pvalue = 0.05,
  lfc = 0
)
p27.1= res[[1]]
p27.1
p27.2= res[[2]]
p27.2
p27.3= res[[3]]
p27.3

# 保存曼哈顿图图片
ggsave(file.path(diffpath, "manhattan_plot_1.png"), plot = p27.1, width = 10, height = 8)
ggsave(file.path(diffpath, "manhattan_plot_1.pdf"), plot = p27.1, width = 10, height = 8)
ggsave(file.path(diffpath, "manhattan_plot_2.png"), plot = p27.2, width = 10, height = 8)
ggsave(file.path(diffpath, "manhattan_plot_2.pdf"), plot = p27.2, width = 10, height = 8)
ggsave(file.path(diffpath, "manhattan_plot_3.png"), plot = p27.3, width = 10, height = 8)
ggsave(file.path(diffpath, "manhattan_plot_3.pdf"), plot = p27.3, width = 10, height = 8)

#24 stemp_diff.metm: stamp展示差异微生物----
# map = phyloseq::sample_data(ps)
# map$ID = row.names(map)
#
# # map$Group = as.factor(map$Group)
# sample_data(ps) = map
allgroup <- combn(unique(map$Group),2)
plot_list <- list()
for (i in 1:dim(allgroup)[2]) {
  ps_sub <- phyloseq::subset_samples(ps.micro,Group %in% allgroup[,i]);ps_sub
  p <- stemp_diff.metm(ps = ps_sub,Top = 20,ranks = 6)
  plot_list[[i]] <- p

}

p28.1 = plot_list[[1]]
p28.1
p28.2 = plot_list[[2]]
p28.2
p28.3 = plot_list[[3]]
p28.3

# 保存图片
ggsave(file.path(diffpath, "stamp_plot_1.png"), plot = p28.1, width = 10, height = 8)
ggsave(file.path(diffpath, "stamp_plot_1.pdf"), plot = p28.1, width = 10, height = 8)
ggsave(file.path(diffpath, "stamp_plot_2.png"), plot = p28.2, width = 10, height = 8)
ggsave(file.path(diffpath, "stamp_plot_2.pdf"), plot = p28.2, width = 10, height = 8)
ggsave(file.path(diffpath, "stamp_plot_3.png"), plot = p28.3, width = 10, height = 8)
ggsave(file.path(diffpath, "stamp_plot_3.pdf"), plot = p28.3, width = 10, height = 8)


#25 Mui.Group.volcano.metm: 聚类火山图------
res =  EdgerSuper2.metm (ps = ps.micro,group  = "Group",artGroup =NULL, j = "OTU")
res2 = Mui.Group.volcano.metm(res = res)

p29.1 = res2[[1]]
p29.1
p29.2 = res2[[2]]
p29.2
dat = res2[[3]]
dat

# 保存聚类火山图图片
ggsave(file.path(diffpath, "cluster_volcano_1.png"), plot = p29.1, width = 10, height = 8)
ggsave(file.path(diffpath, "cluster_volcano_1.pdf"), plot = p29.1, width = 10, height = 8)
ggsave(file.path(diffpath, "cluster_volcano_2.png"), plot = p29.2, width = 10, height = 8)
ggsave(file.path(diffpath, "cluster_volcano_2.pdf"), plot = p29.2, width = 10, height = 8)

# 保存聚类火山图数据
addWorksheet(diff_wb, "cluster_volcano_results")
writeData(diff_wb, "cluster_volcano_results", dat, rowNames = TRUE)
saveWorkbook(diff_wb, file.path(diffpath, "diff_results.xlsx"), overwrite = TRUE)

# biomarker identification -----
# 创建生物标记物分析主目录
biomarkerpath <- file.path(path, "biomarker")
dir.create(biomarkerpath, showWarnings = FALSE, recursive = TRUE)
#26 rfcv.metm :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
result = rfcv.metm(ps = ps.micro %>% filter_OTU_ps(100),
                   group  = "Group",optimal = 20,nrfcvnum = 6)

prfcv = result[[1]]

prfcv+theme_classic()

# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存随机森林交叉验证图片
ggsave(file.path(biomarkerpath, "rfcv_plot.png"), plot = prfcv, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "rfcv_plot.pdf"), plot = prfcv, width = 10, height = 8)

# 保存随机森林交叉验证数据
biomarker_wb <- createWorkbook()
addWorksheet(biomarker_wb, "rfcv_plot_data")
addWorksheet(biomarker_wb, "rfcv_summary_table")
writeData(biomarker_wb, "rfcv_plot_data", result[[2]], rowNames = TRUE)
writeData(biomarker_wb, "rfcv_summary_table", rfcvtable, rowNames = TRUE)

#27 Roc.metm:ROC 曲线绘制----
id = sample_data(ps.micro)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])
b= data.frame(group)

ps = ps.micro %>% subset_taxa.wt("Family","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Order","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Genus","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("Phylum","Unassigned",TRUE)
ps = ps %>% subset_taxa.wt("class","Unassigned",TRUE)

pst = ps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)
res = Roc.metm( ps = pst %>% filter_OTU_ps(1000),group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1+theme_classic()
AUC =  res[[2]]
AUC
dat =  res[[3]]
dat

# 保存ROC分析图片
ggsave(file.path(biomarkerpath, "roc_curve.png"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "roc_curve.pdf"), plot = p33.1, width = 10, height = 8)

# 保存ROC分析数据
addWorksheet(biomarker_wb, "roc_results")
addWorksheet(biomarker_wb, "auc_summary")
writeData(biomarker_wb, "roc_results", dat, rowNames = TRUE)
writeData(biomarker_wb, "auc_summary", data.frame(AUC = AUC), rowNames = TRUE)

#28 loadingPCA.metm: 载荷矩阵筛选特征微生物------
res = loadingPCA.metm(ps = ps.micro,Top = 20)
p34.1 = res[[1]]
p34.1+theme_classic()
dat = res[[2]]
dat

# 保存PCA载荷分析图片
ggsave(file.path(biomarkerpath, "pca_loading.png"), plot = p34.1, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "pca_loading.pdf"), plot = p34.1, width = 10, height = 8)

# 保存PCA载荷分析数据
addWorksheet(biomarker_wb, "pca_loading_data")
writeData(biomarker_wb, "pca_loading_data", dat, rowNames = TRUE)

#29 LDA.metm: LDA筛选特征微生物-----
tablda = LDA.metm(ps = ps.micro,
                   Top = 10,
                   p.lvl = 0.05,
                   lda.lvl = 4,
                   seed = 11,
                   adjust.p = F)

p35 <- lefse_bar(taxtree = tablda[[2]])
p35+theme_classic()
dat = tablda[[2]]
dat

# 保存LDA分析图片
ggsave(file.path(biomarkerpath, "lda_barplot.png"), plot = p35, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "lda_barplot.pdf"), plot = p35, width = 10, height = 8)

# 保存LDA分析数据
addWorksheet(biomarker_wb, "lda_results")
addWorksheet(biomarker_wb, "lda_parameters")
writeData(biomarker_wb, "lda_results", dat, rowNames = TRUE)
writeData(biomarker_wb, "lda_parameters", data.frame(
  Parameter = c("Top", "p.lvl", "lda.lvl", "seed", "adjust.p"),
  Value = c("10", "0.05", "4", "11", "FALSE")
))

#30 svm.metm:svm筛选特征微生物 ----
res <- svm.metm(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM分析数据
addWorksheet(biomarker_wb, "svm_auc")
addWorksheet(biomarker_wb, "svm_feature_importance")
writeData(biomarker_wb, "svm_auc", data.frame(AUC = AUC))
writeData(biomarker_wb, "svm_feature_importance", importance, rowNames = TRUE)

#31 glm.metm :glm筛选特征微生物----
res <- glm.metm (ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM分析数据
addWorksheet(biomarker_wb, "glm_auc")
addWorksheet(biomarker_wb, "glm_feature_importance")
writeData(biomarker_wb, "glm_auc", data.frame(AUC = AUC))
writeData(biomarker_wb, "glm_feature_importance", importance, rowNames = TRUE)

#32 xgboost.metm: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)
res = xgboost.metm(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance
importance_data <- as.data.frame(importance_data)

# 保存XGBoost分析数据
addWorksheet(biomarker_wb, "xgboost_accuracy")
addWorksheet(biomarker_wb, "xgboost_feature_importance")
writeData(biomarker_wb, "xgboost_accuracy", data.frame(Accuracy = accuracy))
writeData(biomarker_wb, "xgboost_feature_importance", importance_data, rowNames = TRUE)

#33 lasso.metm: lasso筛选特征微生物----
library(glmnet)
res =lasso.metm(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso分析数据
addWorksheet(biomarker_wb, "lasso_accuracy")
addWorksheet(biomarker_wb, "lasso_feature_importance")
writeData(biomarker_wb, "lasso_accuracy", data.frame(Accuracy = accuracy))
writeData(biomarker_wb, "lasso_feature_importance", importance, rowNames = TRUE)

#34 decisiontree.micro: 错----
library(rpart)
res =decisiontree.metm(ps=pst, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树分析数据
addWorksheet(biomarker_wb, "tree_accuracy")
addWorksheet(biomarker_wb, "tree_feature_importance")
writeData(biomarker_wb, "tree_accuracy", data.frame(Accuracy = accuracy))
writeData(biomarker_wb, "tree_feature_importance", importance, rowNames = TRUE)

#35 naivebayes.metm: bayes筛选特征微生物----
res = naivebayes.metm(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯分析数据
addWorksheet(biomarker_wb, "naivebayes_accuracy")
addWorksheet(biomarker_wb, "naivebayes_feature_importance")
writeData(biomarker_wb, "naivebayes_accuracy", data.frame(Accuracy = accuracy))
writeData(biomarker_wb, "naivebayes_feature_importance", importance, rowNames = TRUE)

#36 randomforest.metm: 随机森林筛选特征微生物----
res = randomforest.metm( ps = pst%>% filter_OTU_ps(20),group  = "Group", optimal = 20)
p42.1 = res[[1]]
p42.1+theme_classic()
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林图片
ggsave(file.path(biomarkerpath, "rf_importance.png"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "rf_importance.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "rf_error.png"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "rf_error.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "rf_additional.png"), plot = p42.4, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "rf_additional.pdf"), plot = p42.4, width = 10, height = 8)

# 保存随机森林数据
addWorksheet(biomarker_wb, "rf_feature_importance")
writeData(biomarker_wb, "rf_feature_importance", dat, rowNames = TRUE)

#37 bagging.metm: Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res =bagging.metm(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging分析数据
addWorksheet(biomarker_wb, "bagging_accuracy")
addWorksheet(biomarker_wb, "bagging_feature_importance")
writeData(biomarker_wb, "bagging_accuracy", data.frame(Accuracy = accuracy))
writeData(biomarker_wb, "bagging_feature_importance", importance, rowNames = TRUE)

#38 nnet.metm : 神经网络筛选特征微生物 ------
res =nnet.metm (ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络分析数据
addWorksheet(biomarker_wb, "nnet_accuracy")
addWorksheet(biomarker_wb, "nnet_feature_importance")
writeData(biomarker_wb, "nnet_accuracy", data.frame(Accuracy = accuracy))
writeData(biomarker_wb, "nnet_feature_importance", importance, rowNames = TRUE)
saveWorkbook(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"), overwrite = TRUE)

#network analysis -----
# 创建网络分析主目录
networkpath <- file.path(path, "network")
dir.create(networkpath, showWarnings = FALSE, recursive = TRUE)
#39 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
detach("package:mia", unload = TRUE)
tab.r = network.pip(
  ps = ps.micro,
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

# 保存主网络分析图片
ggsave(file.path(networkpath, "main_network_plot.png"), plot = p0, width = 12, height = 10)
ggsave(file.path(networkpath, "main_network_plot.pdf"), plot = p0, width = 12, height = 10)

# 保存主网络分析数据
network_wb <- createWorkbook()

# 保存相关性矩阵信息
if (!is.null(cortab)) {
  for (i in 1:length(cortab)) {
    group_name <- names(cortab)[i]
    addWorksheet(main_network_wb, paste0("correlation_", group_name))

    # 转换相关性矩阵为数据框
    if (is.matrix(cortab[[i]])) {
      cor_df <- as.data.frame(cortab[[i]])
      writeData(main_network_wb, paste0("correlation_", group_name), cor_df, rowNames = TRUE)
    }
  }
}

# 保存网络基本信息
network_info <- data.frame(
  Parameter = c("N_features", "r_threshold", "p_threshold", "maxnode", "layout", "cluster_method"),
  Value = c("200", "0.6", "0.05", "2", "model_maptree2", "cluster_fast_greedy")
)
addWorksheet(network_wb, "network_parameters")
writeData(network_wb, "network_parameters", network_info)


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

# 保存网络属性数据
addWorksheet(network_wb, "network_properties_summary")
writeData(network_wb, "network_properties_summary", dat2, rowNames = TRUE)

#41 netproperties.sample:单个样本的网络属性#-------
for (i in 1:length(id)) {
  pst = ps.micro %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
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

# 保存样本网络属性数据
addWorksheet(network_wb, "sample_network_properties")
addWorksheet(network_wb, "sample_metadata_combined")
writeData(network_wb, "sample_network_properties", dat.f2, rowNames = TRUE)
writeData(network_wb, "sample_metadata_combined", dat3)

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

# 保存节点属性数据
addWorksheet(network_wb, "node_properties_combined")
writeData(network_wb, "node_properties_combined", nodepro2, rowNames = FALSE)


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

# 保存网络比较数据
addWorksheet(network_wb, "network_comparison_results")
writeData(network_wb, "network_comparison_results", res)
saveWorkbook(network_wb, file.path(networkpath, "network_results.xlsx"), overwrite = TRUE)


