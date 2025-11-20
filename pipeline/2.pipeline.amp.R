# 扩增子微生物组学分析#-----
rm(list=ls())
library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(ggsci)

# 定义后续参数-----
# # 弹出文件选择窗口，选择 RDS 文件
# file_path <- file.choose()
# # 读入数据
# ps.16s <- readRDS(file_path)
library(tcltk)
file_path <- tcltk::tk_choose.files(
  caption = "请选择 ps_ITS.rds 文件",
  multi = FALSE,
  filters = matrix(c("RDS files", ".rds",
                     "All files", "*"), ncol = 2, byrow = TRUE)
)


ps.16s <- readRDS(file_path)
ps.16s



map= sample_data(ps.16s)
head(map)


# 提取分组因子数量
gnum = phyloseq::sample_data(ps.16s)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.16s对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps.16s)$Group %>%unique();axis_order

# col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")
col.g <- gen_colors(axis_order, palette = "Dark2")
col.g

# 创建扩增子微生物组分析目录
amplicon_path =  "../2.pipeline.amp.amplicon.16S"
fs::dir_create(amplicon_path)

#-主题--颜色等
package.amp()
res = theme_my(ps.16s)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]


# alpha diversity -----
# 创建Alpha多样性分析目录
amplicon_alpha_path <- file.path(amplicon_path, "alpha_diversity")
dir.create(amplicon_alpha_path, recursive = TRUE)

# 创建Alpha多样性分析总表
amplicon_alpha_wb <- createWorkbook()

#1 alpha.micro: 6中alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.micro(ps = ps.16s,group = "Group" )
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

result = MuiaovMcomper2(data = data,num = 3:9)
result1 = FacetMuiPlotresultBox2(data = data,num = 3:9,
                                result = result,
                                sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1+
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

res = FacetMuiPlotresultBar(data = data,num = c(3:(9)),result = result,sig_show ="abc",ncol = 4)

p1_2 = res[[1]]
p1_2+
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

res = FacetMuiPlotReBoxBar(data = data,num = c(3:(9)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

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
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

# 保存Alpha多样性结果
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_box.png"), plot = p1_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_box.pdf"), plot = p1_1, width = 12, height = 8)
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_bar.png"), plot = p1_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_bar.pdf"), plot = p1_2, width = 12, height = 8)
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_boxbar.png"), plot = p1_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_boxbar.pdf"), plot = p1_3, width = 12, height = 8)
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_violin.png"), plot = p1_0, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "alpha_diversity_violin.pdf"), plot = p1_0, width = 12, height = 8)
addWorksheet(amplicon_alpha_wb, "alpha_diversity_data")
writeData(amplicon_alpha_wb, "alpha_diversity_data", data, rowNames = TRUE)

#2 alpha.pd :用于计算pd多样性#-------
library(ape)
library(picante)

map = ps.16s %>% sample_data()
map$ID = row.names(map)
sample_data(ps.16s) = map

tab2 = alpha.pd.micro(ps.16s)
head(tab2)
result = MuiKwWlx2(data = tab2,num = 3)
result1 = FacetMuiPlotresultBox(data = tab2,num = 3,
                                result = result,
                                sig_show ="abc",ncol = 1 )
p1_1 = result1[[1]]
p1_1+
  #ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)

# 保存PD多样性结果
ggsave(file.path(amplicon_alpha_path, "pd_diversity.png"), plot = p1_1, width = 8, height = 6, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "pd_diversity.pdf"), plot = p1_1, width = 8, height = 6)
addWorksheet(amplicon_alpha_wb, "pd_diversity_data")
writeData(amplicon_alpha_wb, "pd_diversity_data", tab2, rowNames = TRUE)

#3 alpha.rare.line:alpha多样性稀释曲线#---------
rare <- mean(phyloseq::sample_sums(ps.16s))/10
result = alpha.rare.line.micro(ps = ps.16s, group = "Group", method = "Richness", start = 100, step = rare)
#-- Plot the rarefaction curve for a single sample
p2_1 <- result[[1]]
p2_1+theme_nature()+theme(axis.title.y = element_text(angle = 90))
## Provide a data table for convenient output
raretab <- result[[2]]
head(raretab)
#-- Display rarefaction curves grouped by categories
p2_2 <- result[[3]]
p2_2+theme_nature()+theme(axis.title.y = element_text(angle = 90))
#-- Plot rarefaction curves with standard deviations by groups
p2_3 <- result[[4]]
p2_3+theme_nature()+theme(axis.title.y = element_text(angle = 90))

# 保存稀释曲线结果
ggsave(file.path(amplicon_alpha_path, "rarefaction_individual.png"), plot = p2_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "rarefaction_individual.pdf"), plot = p2_1, width = 10, height = 8)
ggsave(file.path(amplicon_alpha_path, "rarefaction_group.png"), plot = p2_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "rarefaction_group.pdf"), plot = p2_2, width = 10, height = 8)
ggsave(file.path(amplicon_alpha_path, "rarefaction_group_sd.png"), plot = p2_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_alpha_path, "rarefaction_group_sd.pdf"), plot = p2_3, width = 10, height = 8)
addWorksheet(amplicon_alpha_wb, "rarefaction_data")
writeData(amplicon_alpha_wb, "rarefaction_data", raretab, rowNames = TRUE)

# 保存Alpha多样性总表
saveWorkbook(amplicon_alpha_wb, file.path(amplicon_alpha_path, "alpha_diversity.xlsx"), overwrite = TRUE)

# beta diversity-----
# 创建Beta多样性分析目录
amplicon_beta_path <- file.path(amplicon_path, "beta_diversity")
dir.create(amplicon_beta_path, recursive = TRUE)

# 创建Beta多样性分析总表
amplicon_beta_wb <- createWorkbook()

#4 ordinate.micro: 排序分析#----------
result = ordinate.micro(ps = ps.16s , group = "Group", dist = "bray",
                        method = "PCoA", Micromet = "anosim", pvalue.cutoff = 0.05,
                        pair = F)
p3_1 = result[[1]]
p3_1 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

#带标签图形出图
p3_2 = result[[3]]
p3_2+
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +theme_nature()+
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

p3_3 = p3_1 +geom_segment(data = segs,
                          mapping = aes(xend = oNMDS1, yend = oNMDS2,color = Group),show.legend=F) + # spiders
  geom_point(mapping = aes(x = x, y = y),data = cent, size = 5,pch = 24,color = "black",fill = "yellow")
p3_3 +
  scale_fill_manual(values = col.g)+
  scale_color_manual(values = col.g,guide = F) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

# 保存排序分析结果
ggsave(file.path(amplicon_beta_path, "ordination_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "ordination_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(amplicon_beta_path, "ordination_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "ordination_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(amplicon_beta_path, "ordination_refined.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "ordination_refined.pdf"), plot = p3_3, width = 10, height = 8)
addWorksheet(amplicon_beta_wb, "ordination_data")
writeData(amplicon_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

#5 MicroTest.micro:群落水平差异检测-------
dat1 = MicroTest.micro(ps = ps.16s, Micromet = "adonis", dist = "bray")
dat1

# 保存群落差异检测结果
addWorksheet(amplicon_beta_wb, "microtest_results")
writeData(amplicon_beta_wb, "microtest_results", dat1, rowNames = TRUE)

#6 pairMicroTest:两两分组群落水平差异检测-------
dat2 = pairMicroTest.micro(ps = ps.16s, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两分组差异检测结果
addWorksheet(amplicon_beta_wb, "pair_microtest_results")
writeData(amplicon_beta_wb, "pair_microtest_results", dat2, rowNames = TRUE)

#7 mantal.micro ：群落差异检测普鲁士分析#------
result <- mantal.micro(ps = ps.16s,
                       method =  "spearman",
                       group = "Group",
                       ncol = gnum,
                       nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

# 保存Mantel分析结果
ggsave(file.path(amplicon_beta_path, "mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "mantel_test.pdf"), plot = p3_7, width = 10, height = 8)
addWorksheet(amplicon_beta_wb, "mantel_results")
writeData(amplicon_beta_wb, "mantel_results", data, rowNames = TRUE)

#8 cluster.micro:样品聚类#-----
res = cluster_micro(ps= ps.16s,
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

# 保存聚类分析结果
ggsave(file.path(amplicon_beta_path, "cluster_heatmap.png"), plot = p4, width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_beta_path, "cluster_heatmap.pdf"), plot = p4, width = 12, height = 10)
ggsave(file.path(amplicon_beta_path, "cluster_dendrogram1.png"), plot = p4_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "cluster_dendrogram1.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(amplicon_beta_path, "cluster_dendrogram2.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "cluster_dendrogram2.pdf"), plot = p4_2, width = 10, height = 8)
addWorksheet(amplicon_beta_wb, "cluster_results")
writeData(amplicon_beta_wb, "cluster_results", dat, rowNames = TRUE)

#9 distance.micro:分组之间距离比对#-----
res = distance.micro2(ps = ps.16s,group = "Group")
p5.1 = res[[1]]
p5.1
p5.2 = res[[2]]
p5.2 +
  scale_fill_manual(values = colset2)+
  scale_color_manual(values = colset1,guide = F)+theme_nature()+
  theme(axis.title.y = element_text(angle = 90))
p5.3 = res[[3]]
p5.3 +
  scale_fill_manual(values = colset2)+
  theme(axis.title.x = element_text(angle = 90))
dat = res[[4]]
head(dat)

# 保存距离比对结果
ggsave(file.path(amplicon_beta_path, "distance_plot1.png"), plot = p5.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "distance_plot1.pdf"), plot = p5.1, width = 10, height = 8)
ggsave(file.path(amplicon_beta_path, "distance_plot2.png"), plot = p5.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "distance_plot2.pdf"), plot = p5.2, width = 10, height = 8)
ggsave(file.path(amplicon_beta_path, "distance_plot3.png"), plot = p5.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_beta_path, "distance_plot3.pdf"), plot = p5.3, width = 10, height = 8)
addWorksheet(amplicon_beta_wb, "distance_data")
writeData(amplicon_beta_wb, "distance_data", dat, rowNames = TRUE)

# 保存Beta多样性总表
saveWorkbook(amplicon_beta_wb, file.path(amplicon_beta_path, "beta_diversity.xlsx"), overwrite = TRUE)

# composition -----
# 创建微生物组成分析目录
amplicon_composition_path <- file.path(amplicon_path, "composition")
dir.create(amplicon_composition_path, recursive = TRUE)

# 创建微生物组成分析总表
amplicon_composition_wb <- createWorkbook()

#10 Ven.Upset.micro: 用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = Ven.Upset.metm(ps =  ps.16s,
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

# 保存韦恩图和upset图结果
ggsave(file.path(amplicon_composition_path, "venn_diagram.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "venn_diagram.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(amplicon_composition_path, "upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "upset_plot.pdf"), plot = p10.2, width = 12, height = 8)
ggsave(file.path(amplicon_composition_path, "upset_plot3.png"), plot = p10.3, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "upset_plot3.pdf"), plot = p10.3, width = 12, height = 8)
addWorksheet(amplicon_composition_wb, "venn_upset_data")
writeData(amplicon_composition_wb, "venn_upset_data", dat, rowNames = TRUE)

#11 ggVen.Upset.micro:用于展示共有、特有的OTU/ASV----
# 分组小于6时使用
res = ggVen.Upset.micro(ps = ps.16s,group = "Group")
grid::grid.draw(res[[1]])
dat = res[[2]]
dat

# 保存ggVen结果
ggsave(file.path(amplicon_composition_path, "ggVen_upset.png"), plot = res[[1]], width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "ggVen_upset.pdf"), plot = res[[1]], width = 12, height = 10)
addWorksheet(amplicon_composition_wb, "ggVen_data")
writeData(amplicon_composition_wb, "ggVen_data", dat, rowNames = TRUE)

#12 VenSeper.micro: 详细展示每一组中OTU的物种组成 -------
#---每个部分
library("ggpubr")
library(agricolae)
library(reshape2)

result = VenSuper.metm(ps = ps.16s,
                        group =  "Group",
                        num = 3
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

# 保存韦恩详细分析结果
ggsave(file.path(amplicon_composition_path, "venn_detail_bar.png"), plot = p7_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "venn_detail_bar.pdf"), plot = p7_1, width = 12, height = 8)
ggsave(file.path(amplicon_composition_path, "venn_detail_alluvial.png"), plot = p7_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "venn_detail_alluvial.pdf"), plot = p7_2, width = 12, height = 8)
ggsave(file.path(amplicon_composition_path, "venn_detail_proportion.png"), plot = p8, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "venn_detail_proportion.pdf"), plot = p8, width = 10, height = 8)
addWorksheet(amplicon_composition_wb, "venn_detail_data")
writeData(amplicon_composition_wb, "venn_detail_data", result[[4]], rowNames = TRUE)

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

# 保存花瓣图结果
ggsave(file.path(amplicon_composition_path, "flower_plot.png"), plot = p13.1, width = 10, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "flower_plot.pdf"), plot = p13.1, width = 10, height = 10)
addWorksheet(amplicon_composition_wb, "flower_data")
writeData(amplicon_composition_wb, "flower_data", dat, rowNames = TRUE)

#14 ven.network.micro:venn 网络----
result = ven.network.micro(
  ps = ps.16s,
  N = 0.5,
  fill = "Phylum")
p14  = result[[1]]
p14

dat = result[[2]]
head(dat)

# 保存韦恩网络图结果
ggsave(file.path(amplicon_composition_path, "venn_network.png"), plot = p14, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "venn_network.pdf"), plot = p14, width = 10, height = 8)
addWorksheet(amplicon_composition_wb, "venn_network_data")
writeData(amplicon_composition_wb, "venn_network_data", dat, rowNames = TRUE)

#15 Micro_tern.micro: 三元图展示组成----
p <- ps_polygon_plot(ps.16s %>% filter_OTU_ps(10000), group = "Group",taxrank = "Phylum")
print(p)


# res = Micro_tern.micro(ps.16s %>% filter_OTU_ps(100))
# p15 = res[[1]]
# p15[[1]] +theme_bw()
#
# dat =  res[[2]]
# head(dat)

# 保存三元图结果
ggsave(file.path(amplicon_composition_path, "ternary_plot.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "ternary_plot.pdf"), plot = p, width = 10, height = 8)
# addWorksheet(amplicon_composition_wb, "ternary_data")
# writeData(amplicon_composition_wb, "ternary_data", dat, rowNames = TRUE)

#16 barMainplot.micro: 堆积柱状图展示组成----
library(ggalluvial)
pst = ps.16s %>% subset_taxa.wt("Species","Unassigned",TRUE)
pst = pst %>% subset_taxa.wt("Genus","Unassigned",TRUE)
rank.names(pst)
result = barMainplot.micro(ps = pst,
                           j = "Genus",
                           # axis_ord = axis_order,
                           label = FALSE,
                           sd = FALSE,
                           Top =10)
p4_1 <- result[[1]]
p4_1+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

p4_2  <- result[[3]]
p4_2+
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +theme_nature()+
  theme(axis.title.y = element_text(angle = 90))

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Genus","Abundance(%)")
head(databar)

# 保存堆积柱状图结果
ggsave(file.path(amplicon_composition_path, "barplot_main.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "barplot_main.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(amplicon_composition_path, "barplot_secondary.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "barplot_secondary.pdf"), plot = p4_2, width = 12, height = 8)
addWorksheet(amplicon_composition_wb, "barplot_data")
writeData(amplicon_composition_wb, "barplot_data", databar, rowNames = TRUE)

#17 cluMicro.bar.micro: 聚类堆积柱状图展示组成 问题-----
result <-  cluMicro.bar.metm (dist = "bray",
                               ps = ps.16s,
                               j = "Genus",
                               Top = 10, # 提取丰度前十的物种注释
                               tran = TRUE, # 转化为相对丰度值
                               hcluter_method = "complete",
                               Group = "Group",
                               cuttree = length(unique(phyloseq::sample_data(ps.16s)$Group)))
result[[1]]

p5_2 <- result[[2]]
p5_2
p5_3 <- result[[3]]
p5_3
p5_4 <- result[[4]]
p5_4
clubardata <- result[[5]]

# 保存聚类堆积柱状图结果
ggsave(file.path(amplicon_composition_path, "cluster_bar1.png"), plot = result[[1]], width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "cluster_bar1.pdf"), plot = result[[1]], width = 12, height = 8)
ggsave(file.path(amplicon_composition_path, "cluster_bar2.png"), plot = p5_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "cluster_bar2.pdf"), plot = p5_2, width = 12, height = 8)
ggsave(file.path(amplicon_composition_path, "cluster_bar3.png"), plot = p5_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "cluster_bar3.pdf"), plot = p5_3, width = 12, height = 8)
ggsave(file.path(amplicon_composition_path, "cluster_bar4.png"), plot = p5_4, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "cluster_bar4.pdf"), plot = p5_4, width = 12, height = 8)
addWorksheet(amplicon_composition_wb, "cluster_bar_data")
writeData(amplicon_composition_wb, "cluster_bar_data", clubardata, rowNames = TRUE)

#18 cir_barplot.micro:环状堆积柱状图 -----
library(ggtree)
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

# 保存环状堆积柱状图结果
ggsave(file.path(amplicon_composition_path, "circular_barplot.png"), plot = p17, width = 10, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "circular_barplot.pdf"), plot = p17, width = 10, height = 10)
addWorksheet(amplicon_composition_wb, "circular_barplot_data")
writeData(amplicon_composition_wb, "circular_barplot_data", dat, rowNames = TRUE)








#19 cir_plot.micro:和弦图展示物种组成-----
res = cir_plot.micro(ps  = ps.16s,Top = 12,rank = 6)

# 保存和弦图结果
ggsave(file.path(amplicon_composition_path, "chord_plot.png"), plot = res, width = 10, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "chord_plot.pdf"), plot = res, width = 10, height = 10)

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

# 保存圈图结果
ggsave(file.path(amplicon_composition_path, "maptree.png"), plot = p20, width = 12, height = 12, dpi = 300)
ggsave(file.path(amplicon_composition_path, "maptree.pdf"), plot = p20, width = 12, height = 12)
addWorksheet(amplicon_composition_wb, "maptree_data")
writeData(amplicon_composition_wb, "maptree_data", dat, rowNames = TRUE)

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

# 保存系统发育树结果
ggsave(file.path(amplicon_composition_path, "phy_tree_main.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "phy_tree_main.pdf"), plot = p0, width = 12, height = 10)
ggsave(file.path(amplicon_composition_path, "phy_tree_p1.png"), plot = p1, width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "phy_tree_p1.pdf"), plot = p1, width = 12, height = 10)
ggsave(file.path(amplicon_composition_path, "phy_tree_p7.png"), plot = p7, width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "phy_tree_p7.pdf"), plot = p7, width = 12, height = 10)
addWorksheet(amplicon_composition_wb, "phy_tree_data")
writeData(amplicon_composition_wb, "phy_tree_data", dat, rowNames = TRUE)

#22 sankey.micro: 桑基图展示物种组成------
res = sankey.micro (ps = ps.16s,
                    rank = 6,
                    Top = 100)
p22 = res[[1]]
p22
dat = res[[2]]
dat

# 保存桑基图结果
ggsave(file.path(amplicon_composition_path, "sankey.png"), plot = p22, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_composition_path, "sankey.pdf"), plot = p22, width = 12, height = 8)
addWorksheet(amplicon_composition_wb, "sankey_data")
writeData(amplicon_composition_wb, "sankey_data", dat, rowNames = TRUE)

#23 sankey.m.Group.micro: 按照分组绘制-----
result = sankey.m.Group.micro(
  ps = ps.16s  %>% subset_taxa.wt("Species","Unassigned",TRUE),
  rank = 6,
  Top = 50)

p = result[[1]]
p

saveNetwork(p,paste(file.path(amplicon_composition_path, "sankey_Group.html"), sep = ""))

dat = result[[2]]
dat

# 保存分组桑基图数据
addWorksheet(amplicon_composition_wb, "sankey_group_data")
writeData(amplicon_composition_wb, "sankey_group_data", dat, rowNames = TRUE)

#24 Microheatmap.micro: 热图展示物种相对丰度差异-----
heatnum = 30
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

result <- Microheatmap.micro(ps_rela = ps_tem,id = id ,col_cluster = TRUE,
                             row_cluster = TRUE)

p24.1 <- result[[1]]
p24.1
p24.2 <- result[[2]]
p24.2
dat = result[[3]]
head(dat)

# 保存微生物热图结果
ggsave(file.path(amplicon_composition_path, "heatmap1.png"), plot = p24.1, width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "heatmap1.pdf"), plot = p24.1, width = 12, height = 10)
ggsave(file.path(amplicon_composition_path, "heatmap2.png"), plot = p24.2, width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_composition_path, "heatmap2.pdf"), plot = p24.2, width = 12, height = 10)
addWorksheet(amplicon_composition_wb, "heatmap_data")
writeData(amplicon_composition_wb, "heatmap_data", dat, rowNames = TRUE)

# 保存微生物组成总表
saveWorkbook(amplicon_composition_wb, file.path(amplicon_composition_path, "composition_analysis.xlsx"), overwrite = TRUE)

# difference analysis -----
# 创建差异分析目录
amplicon_diff_path <- file.path(amplicon_path, "differential")
dir.create(amplicon_diff_path, recursive = TRUE)

# 创建差异分析总表
amplicon_diff_wb <- createWorkbook()

#25 EdgerSuper.micro:EdgeR计算差异微生物----
res = EdgerSuper.micro(ps = ps.16s %>%
                         tax_glom_wt(6) %>%
                         ggClusterNet::filter_OTU_ps(500),group  = "Group",artGroup = NULL, j = "OTU")

p25.1 =  res[[1]][[1]]
p25.1
p25.2 =  res[[1]][[2]]
p25.2
p25.3 =  res[[1]][[3]]
p25.3
dat =  res[[2]]
head(dat)

# 保存EdgeR结果
ggsave(file.path(amplicon_diff_path, "EdgeR_plot1.png"), plot = p25.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "EdgeR_plot1.pdf"), plot = p25.1, width = 10, height = 8)
ggsave(file.path(amplicon_diff_path, "EdgeR_plot2.png"), plot = p25.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "EdgeR_plot2.pdf"), plot = p25.2, width = 10, height = 8)
ggsave(file.path(amplicon_diff_path, "EdgeR_plot3.png"), plot = p25.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "EdgeR_plot3.pdf"), plot = p25.3, width = 10, height = 8)
addWorksheet(amplicon_diff_wb, "EdgeR_results")
writeData(amplicon_diff_wb, "EdgeR_results", dat, rowNames = TRUE)

#25.1 EdgerSuper2.micro:EdgeR计算差异微生物-----
res =  EdgerSuper2.micro (ps = ps.16s,group  = "Group",artGroup =NULL, j = "OTU")
head(res)

res %>% filter(level != "nosig") %>% .$group %>% table()





# 保存EdgeR2结果
addWorksheet(amplicon_diff_wb, "EdgeR2_results")
writeData(amplicon_diff_wb, "EdgeR2_results", res, rowNames = TRUE)

#26 DESep2Super.micro:DESep2计算差异微生物-----
res = DESep2Super.micro(ps = ps.16s %>% ggClusterNet::filter_OTU_ps(500),
                        group  = "Group",
                        artGroup =NULL,
                        j = "OTU")
p26.1 =  res[[1]][[1]]
p26.1
p26.2 =  res[[1]][[2]]
p26.2
p26.3 =  res[[1]][[3]]
p26.3
dat =  res[[2]]
dat

# 保存DESep2结果
ggsave(file.path(amplicon_diff_path, "DESep2_plot1.png"), plot = p26.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "DESep2_plot1.pdf"), plot = p26.1, width = 10, height = 8)
ggsave(file.path(amplicon_diff_path, "DESep2_plot2.png"), plot = p26.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "DESep2_plot2.pdf"), plot = p26.2, width = 10, height = 8)
ggsave(file.path(amplicon_diff_path, "DESep2_plot3.png"), plot = p26.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "DESep2_plot3.pdf"), plot = p26.3, width = 10, height = 8)
addWorksheet(amplicon_diff_wb, "DESep2_results")
writeData(amplicon_diff_wb, "DESep2_results", dat, rowNames = TRUE)

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

# 保存曼哈顿图结果
ggsave(file.path(amplicon_diff_path, "Manhattan_plot1.png"), plot = p27.1, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "Manhattan_plot1.pdf"), plot = p27.1, width = 12, height = 8)
ggsave(file.path(amplicon_diff_path, "Manhattan_plot2.png"), plot = p27.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "Manhattan_plot2.pdf"), plot = p27.2, width = 12, height = 8)
ggsave(file.path(amplicon_diff_path, "Manhattan_plot3.png"), plot = p27.3, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "Manhattan_plot3.pdf"), plot = p27.3, width = 12, height = 8)

#28 stamp_diff.micro: stamp展示差异微生物----
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

# 保存STAMP结果
ggsave(file.path(amplicon_diff_path, "stamp_plot1.png"), plot = p28.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "stamp_plot1.pdf"), plot = p28.1, width = 10, height = 8)
ggsave(file.path(amplicon_diff_path, "stamp_plot2.png"), plot = p28.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "stamp_plot2.pdf"), plot = p28.2, width = 10, height = 8)
ggsave(file.path(amplicon_diff_path, "stamp_plot3.png"), plot = p28.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "stamp_plot3.pdf"), plot = p28.3, width = 10, height = 8)

#29 Mui.Group.volcano.micro: 聚类火山图------
res =  EdgerSuper2.micro (ps = ps.16s,group  = "Group",artGroup =NULL, j = "OTU")
res2 = Mui.Group.volcano.micro(res = res)

p29.1 = res2[[1]]
p29.1
p29.2 = res2[[2]]
p29.2
dat = res2[[3]]
dat

# 保存聚类火山图结果
ggsave(file.path(amplicon_diff_path, "volcano_cluster1.png"), plot = p29.1, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "volcano_cluster1.pdf"), plot = p29.1, width = 12, height = 8)
ggsave(file.path(amplicon_diff_path, "volcano_cluster2.png"), plot = p29.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "volcano_cluster2.pdf"), plot = p29.2, width = 12, height = 8)
addWorksheet(amplicon_diff_wb, "volcano_cluster_data")
writeData(amplicon_diff_wb, "volcano_cluster_data", dat, rowNames = TRUE)

#30 Mui.cluster.volcano.micro: 指定分组绘制聚类火山图------
library(ggrepel)
id = sample_data(ps.16s)$Group %>% unique()
aaa = combn(id,2)
# 设置分组
group2 = c(aaa[1,1],aaa[2,1])
b= data.frame(group2)

res =  EdgerSuper2.micro (ps = ps.16s,group  = "Group",artGroup =b, j = "OTU")

res2 = Mui.cluster.volcano.micro (res = res,rs.k = 6)
p30.1 = res2[[1]]+ggtitle(paste(group2, collapse = "-"))
p30.1
p30.2 = res2[[2]]+ggtitle(paste(group2, collapse = "-"))
p30.2
dat= res2[[3]]
dat

# 保存指定分组聚类火山图结果
ggsave(file.path(amplicon_diff_path, "volcano_specific1.png"), plot = p30.1, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "volcano_specific1.pdf"), plot = p30.1, width = 12, height = 8)
ggsave(file.path(amplicon_diff_path, "volcano_specific2.png"), plot = p30.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_diff_path, "volcano_specific2.pdf"), plot = p30.2, width = 12, height = 8)
addWorksheet(amplicon_diff_wb, "volcano_specific_data")
writeData(amplicon_diff_wb, "volcano_specific_data", dat, rowNames = TRUE)

# 保存差异分析总表
saveWorkbook(amplicon_diff_wb, file.path(amplicon_diff_path, "differential_analysis.xlsx"), overwrite = TRUE)

# biomarker identification -----
# 创建生物标志物识别目录
amplicon_biomarker_path <- file.path(amplicon_path, "biomarker")
dir.create(amplicon_biomarker_path, recursive = TRUE)
# 创建生物标志物识别总表
amplicon_biomarker_wb <- createWorkbook()

library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

#32 rfcv.micro :交叉验证结果-------
result =rfcv.Micro(ps = ps.16s %>% filter_OTU_ps(100),
                   group  = "Group",optimal = 20,nrfcvnum = 6)

prfcv = result[[1]]+theme_classic()
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存交叉验证结果
ggsave(file.path(amplicon_biomarker_path, "rfcv.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "rfcv.pdf"), plot = prfcv, width = 10, height = 8)
addWorksheet(amplicon_biomarker_wb, "rfcv_results")
writeData(amplicon_biomarker_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

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
p33.1+theme_classic()
p33.2 =  res[[2]]

dat =  res[[3]]
dat

# 保存ROC结果
ggsave(file.path(amplicon_biomarker_path, "ROC_plot1.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "ROC_plot1.pdf"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(amplicon_biomarker_path, "ROC_plot2.png"), plot = p33.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "ROC_plot2.pdf"), plot = p33.2, width = 10, height = 8)
addWorksheet(amplicon_biomarker_wb, "ROC_results")
writeData(amplicon_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#34 loadingPCA.micro: 载荷矩阵筛选特征微生物------
res = loadingPCA.micro(ps = ps.16s,Top = 20)
p34.1 = res[[1]]
p34.1+theme_classic()
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(amplicon_biomarker_path, "loadingPCA.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "loadingPCA.pdf"), plot = p34.1, width = 10, height = 8)
addWorksheet(amplicon_biomarker_wb, "loadingPCA_results")
writeData(amplicon_biomarker_wb, "loadingPCA_results", dat, rowNames = TRUE)

#35 LDA.micro: LDA筛选特征微生物-----
p1 <- p_base.micro(ps.16s,Top = 20)
tablda = LDA.micro(ps = ps.16s,
                   Top = 20,
                   p.lvl = 0.05,
                   lda.lvl = 2,
                   seed = 11,
                   adjust.p = F)
tablda[[1]]

p35 <- lefse_bar(taxtree = tablda[[2]])
p35+theme_classic()
dat = tablda[[2]]
dat

# 保存LDA结果
ggsave(file.path(amplicon_biomarker_path, "LDA.png"), plot = p35, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "LDA.pdf"), plot = p35, width = 10, height = 8)
addWorksheet(amplicon_biomarker_wb, "LDA_results")
writeData(amplicon_biomarker_wb, "LDA_results", dat, rowNames = TRUE)

#36 svm.micro:svm筛选特征微生物 ----
res <- svm.micro(ps = ps.16s %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(amplicon_biomarker_wb, "svm_AUC")
writeData(amplicon_biomarker_wb, "svm_AUC", AUC, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "svm_importance")
writeData(amplicon_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#37 glm.micro:glm筛选特征微生物----
res <- glm.micro(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(amplicon_biomarker_wb, "glm_AUC")
writeData(amplicon_biomarker_wb, "glm_AUC", AUC, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "glm_importance")
writeData(amplicon_biomarker_wb, "glm_importance", importance, rowNames = TRUE)

#38 xgboost.micro: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)

res = xgboost.micro(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(amplicon_biomarker_wb, "xgboost_accuracy")
writeData(amplicon_biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "xgboost_importance")
writeData(amplicon_biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#39 lasso.micro: lasso筛选特征微生物----
library(glmnet)
res =lasso.micro (ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(amplicon_biomarker_wb, "lasso_accuracy")
writeData(amplicon_biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "lasso_importance")
writeData(amplicon_biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#40 decisiontree.micro----
library(rpart)
res =decisiontree.micro(ps=pst, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(amplicon_biomarker_wb, "decisiontree_accuracy")
writeData(amplicon_biomarker_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "decisiontree_importance")
writeData(amplicon_biomarker_wb, "decisiontree_importance", importance, rowNames = TRUE)

#41 naivebayes.micro: bayes筛选特征微生物----
res = naivebayes.micro(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(amplicon_biomarker_wb, "naivebayes_accuracy")
writeData(amplicon_biomarker_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "naivebayes_importance")
writeData(amplicon_biomarker_wb, "naivebayes_importance", importance, rowNames = TRUE)

#42 randomforest.micro: 随机森林筛选特征微生物----
res = randomforest.micro( ps = pst,group  = "Group", optimal = 50)
p42.1 = res[[1]]

p42.1 +theme_classic()
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林结果
ggsave(file.path(amplicon_biomarker_path, "randomforest_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "randomforest_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(amplicon_biomarker_path, "randomforest_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "randomforest_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(amplicon_biomarker_path, "randomforest_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_biomarker_path, "randomforest_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(amplicon_biomarker_wb, "randomforest_results")
writeData(amplicon_biomarker_wb, "randomforest_results", dat, rowNames = TRUE)

#43 bagging.micro : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
ps = pst
res = bagging.micro(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(amplicon_biomarker_wb, "bagging_accuracy")
writeData(amplicon_biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "bagging_importance")
writeData(amplicon_biomarker_wb, "bagging_importance", importance, rowNames = TRUE)

#44 nnet.micro 神经网络筛选特征微生物  ------
res =nnet.micro(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(amplicon_biomarker_wb, "nnet_accuracy")
writeData(amplicon_biomarker_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(amplicon_biomarker_wb, "nnet_importance")
writeData(amplicon_biomarker_wb, "nnet_importance", importance, rowNames = TRUE)

# 保存生物标志物识别总表
saveWorkbook(amplicon_biomarker_wb, file.path(amplicon_biomarker_path, "biomarker_analysis.xlsx"), overwrite = TRUE)

#network analysis -----
# 创建网络分析目录
amplicon_network_path <- file.path(amplicon_path, "network")
dir.create(amplicon_network_path, recursive = TRUE)

# 创建网络分析总表
amplicon_network_wb <- createWorkbook()

#45 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
detach("package:mia", unload = TRUE)
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
  label =FALSE,
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

# 保存网络分析主图
ggsave(file.path(amplicon_network_path, "network_main.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(amplicon_network_path, "network_main.pdf"), plot = p0, width = 12, height = 10)

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

# 保存网络属性数据
addWorksheet(amplicon_network_wb, "network_properties")
writeData(amplicon_network_wb, "network_properties", dat2, rowNames = TRUE)

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

# 保存样本网络属性数据
addWorksheet(amplicon_network_wb, "sample_network_properties")
writeData(amplicon_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

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

# 保存节点属性数据
addWorksheet(amplicon_network_wb, "node_properties")
writeData(amplicon_network_wb, "node_properties", nodepro2, rowNames = TRUE)

#49 negative.correlation.ratio:网络稳定性-计算负相关的比例----
res4 = negative.correlation.ratio(ps = ps.16s,
                                  corg = cortab,
                                  # Top = 500,
                                  degree = TRUE,
                                  zipi = FALSE)

p5 = res4[[1]]
p5+theme_classic()
dat6 = res4[[2]]
dat6

# 保存负相关比例结果
ggsave(file.path(amplicon_network_path, "negative_correlation_ratio.png"), plot = p5, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_network_path, "negative_correlation_ratio.pdf"), plot = p5, width = 10, height = 8)
addWorksheet(amplicon_network_wb, "negative_correlation_ratio")
writeData(amplicon_network_wb, "negative_correlation_ratio", dat6, rowNames = TRUE)

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

# 保存群落稳定性结果
ggsave(file.path(amplicon_network_path, "community_stability.png"), plot = p6, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_network_path, "community_stability.pdf"), plot = p6, width = 10, height = 8)
addWorksheet(amplicon_network_wb, "community_stability")
writeData(amplicon_network_wb, "community_stability", dat7, rowNames = TRUE)

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

# 保存网络抗毁性结果
ggsave(file.path(amplicon_network_path, "network_robustness.png"), plot = p7, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_network_path, "network_robustness.pdf"), plot = p7, width = 10, height = 8)
addWorksheet(amplicon_network_wb, "network_robustness")
writeData(amplicon_network_wb, "network_robustness", dat8, rowNames = TRUE)

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

# 保存网络显著性比较结果
addWorksheet(amplicon_network_wb, "module_compare_results")
writeData(amplicon_network_wb, "module_compare_results", res, rowNames = TRUE)

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

# 保存模块相似性结果
ggsave(file.path(amplicon_network_path, "module_similarity1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_network_path, "module_similarity1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(amplicon_network_path, "module_similarity2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_network_path, "module_similarity2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(amplicon_network_wb, "module_similarity_otu")
writeData(amplicon_network_wb, "module_similarity_otu", dat1, rowNames = TRUE)
addWorksheet(amplicon_network_wb, "module_similarity_results")
writeData(amplicon_network_wb, "module_similarity_results", dat2, rowNames = TRUE)

#54 Robustness.Targeted.removal:网络稳定性-去除关键节点-网络鲁棒性#-----
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

# 保存目标移除鲁棒性结果
ggsave(file.path(amplicon_network_path, "robustness_targeted.png"), plot = p3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_network_path, "robustness_targeted.pdf"), plot = p3, width = 10, height = 8)
addWorksheet(amplicon_network_wb, "robustness_targeted")
writeData(amplicon_network_wb, "robustness_targeted", dat4, rowNames = TRUE)

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

# 保存随机移除鲁棒性结果
ggsave(file.path(amplicon_network_path, "robustness_random.png"), plot = p4, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_network_path, "robustness_random.pdf"), plot = p4, width = 10, height = 8)
addWorksheet(amplicon_network_wb, "robustness_random")
writeData(amplicon_network_wb, "robustness_random", dat5, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(amplicon_network_wb, file.path(amplicon_network_path, "network_analysis.xlsx"), overwrite = TRUE)

# community assemble-----
# 创建群落组装分析目录
amplicon_assembly_path <- file.path(amplicon_path, "community_assembly")
dir.create(amplicon_assembly_path, recursive = TRUE)

# 创建群落组装分析总表
amplicon_assembly_wb <- createWorkbook()

library(picante)
library(ape)
library(FSA)
library(eulerr)
library(minpack.lm)
library(Hmisc)
data("ps16s")
ps16s
# 使用绝对丰度
otu_table(ps16s) <- round(otu_table(ps16s) * 1000000)
psphy = filter_taxa(ps16s, function(x) sum(x ) > 1000, TRUE);psphy
#psphy =filter_taxa(ps16s, function(x) mean(x) > 0.001, TRUE)
otu_table(ps16s)
map = sample_data(psphy)
n = map$Group %>% unique() %>%
  length()
n

#56 neutralModel: 中性模型----
result = neutralModel(ps = psphy,group  = "Group",ncol = n)
#--合并图表
p43 =  result[[1]]
p43
dat = result[[3]]
dat
dat2 = result[[4]]
dat2

# 保存中性模型结果
ggsave(file.path(amplicon_assembly_path, "neutral_model.png"), plot = p43, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_assembly_path, "neutral_model.pdf"), plot = p43, width = 12, height = 8)
addWorksheet(amplicon_assembly_wb, "neutral_model_data1")
writeData(amplicon_assembly_wb, "neutral_model_data1", dat, rowNames = TRUE)
addWorksheet(amplicon_assembly_wb, "neutral_model_data2")
writeData(amplicon_assembly_wb, "neutral_model_data2", dat2, rowNames = TRUE)

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

# 保存零模型结果
addWorksheet(amplicon_assembly_wb, "null_model_results")
writeData(amplicon_assembly_wb, "null_model_results", nullModeltab, rowNames = TRUE)
addWorksheet(amplicon_assembly_wb, "null_model_ratio")
writeData(amplicon_assembly_wb, "null_model_ratio", ratiotab, rowNames = TRUE)
addWorksheet(amplicon_assembly_wb, "null_model_aov")
writeData(amplicon_assembly_wb, "null_model_aov", aovtab, rowNames = TRUE)

#58 bNTICul:β最近分类单元指数计算----
library(parallel)
result = bNTICul(ps = psphy,
                 group  = "Group",
                 num = 10,
                 thread = 1
)
bNTI = result[[1]]
head(bNTI)

# 保存βNTI结果
addWorksheet(amplicon_assembly_wb, "bNTI_results")
writeData(amplicon_assembly_wb, "bNTI_results", bNTI, rowNames = TRUE)

#59 RCbary: RCbary 计算----
result = RCbary(ps = psphy ,group  = "Group",num = 10,thread = 1)

RCbary = result[[1]]
head(RCbary)

# 保存RCbray结果
addWorksheet(amplicon_assembly_wb, "RCbray_results")
writeData(amplicon_assembly_wb, "RCbray_results", RCbary, rowNames = TRUE)

#60 bNTIRCPlot: BetaNTI和RCbray联合出图--------
RCb = RCbary %>%
  dplyr::mutate(Sample_1 = Site2, Sample_2 = Site1)

result = bNTIRCPlot(ps = psphy ,RCb  =RCb, bNTI = bNTI,group  = "Group")

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

# 保存βNTI和RCbray联合分析结果
ggsave(file.path(amplicon_assembly_path, "bNTI_plot.png"), plot = p3, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_assembly_path, "bNTI_plot.pdf"), plot = p3, width = 10, height = 8)
ggsave(file.path(amplicon_assembly_path, "RCbray_plot.png"), plot = p4, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_assembly_path, "RCbray_plot.pdf"), plot = p4, width = 10, height = 8)
ggsave(file.path(amplicon_assembly_path, "bNTI_RCbray_combined.png"), plot = p5, width = 12, height = 8, dpi = 300)
ggsave(file.path(amplicon_assembly_path, "bNTI_RCbray_combined.pdf"), plot = p5, width = 12, height = 8)
addWorksheet(amplicon_assembly_wb, "bNTI_RCbray_plotdata")
writeData(amplicon_assembly_wb, "bNTI_RCbray_plotdata", plotdata, rowNames = TRUE)
addWorksheet(amplicon_assembly_wb, "bNTI_RCbray_summary")
writeData(amplicon_assembly_wb, "bNTI_RCbray_summary", dat, rowNames = TRUE)

# 保存群落组装分析总表
saveWorkbook(amplicon_assembly_wb, file.path(amplicon_assembly_path, "community_assembly.xlsx"), overwrite = TRUE)

# other------
# 创建其他分析目录
amplicon_other_path <- file.path(amplicon_path, "other_analysis")
dir.create(amplicon_other_path, recursive = TRUE)

# 创建其他分析总表
amplicon_other_wb <- createWorkbook()

#61 FEAST.micro:溯源分析----
result = FEAST.micro(ps = ps.16s,
                     group = "Group",
                     sinkG = "OE",
                     sourceG = c("WT","KO")
)

# 保存FEAST溯源分析数据
addWorksheet(amplicon_other_wb, "FEAST_results")
writeData(amplicon_other_wb, "FEAST_results", result, rowNames = TRUE)

#62 Plot_FEAST:溯源分析可视化分组----
p <- Plot_FEAST(data = result)
p

# 保存FEAST分组可视化结果
ggsave(file.path(amplicon_other_path, "FEAST_group.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_other_path, "FEAST_group.pdf"), plot = p, width = 10, height = 8)

#63 MuiPlot_FEAST:溯源分析可视化样品----
p2 = MuiPlot_FEAST(data = result)
p2

# 保存FEAST样品可视化结果
ggsave(file.path(amplicon_other_path, "FEAST_sample.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_other_path, "FEAST_sample.pdf"), plot = p2, width = 10, height = 8)

# 保存其他分析总表
saveWorkbook(amplicon_other_wb, file.path(amplicon_other_path, "other_analysis.xlsx"), overwrite = TRUE)

# function prediction-----
# 创建功能预测分析目录
amplicon_function_path <- file.path(amplicon_path, "function_prediction")
dir.create(amplicon_function_path, recursive = TRUE)

# 创建功能预测分析总表
amplicon_function_wb <- createWorkbook()

library(DOSE)
library(GO.db)
library(GSEABase)
library(ggtree)
library(aplot)
library(clusterProfiler)
library("GSVA")
library(dplyr)
ps.kegg = ps.kegg %>% filter_OTU_ps(Top = 1000)

tax =ps.kegg %>% phyloseq::tax_table()

colnames(tax)[3] = "KOid"
tax_table(ps.kegg) =tax
res = EdgerSuper.metf (ps = ps.kegg,
                       group  = "Group",
                       artGroup = NULL)
dat = res[[2]]
dat

# 保存差异分析结果
addWorksheet(amplicon_function_wb, "differential_analysis")
writeData(amplicon_function_wb, "differential_analysis", dat, rowNames = TRUE)

#64 KEGG_enrich.micro: taxfun2功能富集分析----
res2 = KEGG_enrich.micro(ps =  ps.kegg,
                         #  diffpath = diffpath,
                         dif = dat
)
dat1= res2$`KO-OE`
dat2= res2$`KO-WT`
dat3= res2$`OE-WT`

# 保存KEGG富集分析结果
addWorksheet(amplicon_function_wb, "KEGG_enrich_KO_OE")
writeData(amplicon_function_wb, "KEGG_enrich_KO_OE", dat1, rowNames = TRUE)
addWorksheet(amplicon_function_wb, "KEGG_enrich_KO_WT")
writeData(amplicon_function_wb, "KEGG_enrich_KO_WT", dat2, rowNames = TRUE)
addWorksheet(amplicon_function_wb, "KEGG_enrich_OE_WT")
writeData(amplicon_function_wb, "KEGG_enrich_OE_WT", dat3, rowNames = TRUE)

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

result = buplot.micro(dt=dat1,id = id)
p1 = result[[1]]
p1
p2 = result[[2]]
p2

# 保存功能富集气泡图结果
ggsave(file.path(amplicon_function_path, "function_bubble_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_function_path, "function_bubble_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(amplicon_function_path, "function_bubble_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(amplicon_function_path, "function_bubble_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(amplicon_function_wb, "bubble_plot_data")
writeData(amplicon_function_wb, "bubble_plot_data", result[[3]], rowNames = TRUE)

# 保存功能预测分析总表
saveWorkbook(amplicon_function_wb, file.path(amplicon_function_path, "function_prediction.xlsx"), overwrite = TRUE)


