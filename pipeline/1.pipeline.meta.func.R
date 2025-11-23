# 宏基因组功能
rm(list=ls())
library(tidyverse)
library(data.table)
library(ggClusterNet)
library(EasyStat)
library(tidyfst)
library(fs)
library(EasyMultiOmics)
library(openxlsx)
#-主题--颜色----
package.amp()
ps.kegg =  EasyMultiOmics::ps.kegg
# sample_data(ps.kegg)
# phyloseq::tax_table(ps.kegg) %>% head()
res = theme_my(ps.kegg)
mytheme1 = res[[1]]
mytheme2 = res[[2]];
colset1 = res[[3]];colset2 = res[[4]];colset3 = res[[5]];colset4 = res[[6]]
col.g =  c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")

# 提取分组因子数量
gnum = phyloseq::sample_data(ps.kegg)$Group %>% unique() %>% length()
gnum
#--设定排序顺序1：按照ps.16s对象中map文件顺序进行
axis_order =  phyloseq::sample_data(ps.kegg)$Group %>%unique();axis_order

# 保存路径设置---
path =  "./result/metagenome/"
dir.create(path)



# kegg数据库--------
keggpath = file.path(path, "kegg")
ps.metf =ps.kegg %>% filter_OTU_ps(Top = 1000)
tax = ps.metf %>% phyloseq::tax_table()
colnames(tax)[3] = "KOid"
tax_table(ps.metf) =tax
head(tax)

# function diversity -----
#1 alpha.metf: 6种alpha多样性计算#----
# 创建alpha多样性分析目录
alphapath <- file.path(keggpath, "alpha")
dir.create(alphapath, recursive = TRUE)
# alppath = paste(otupath,"/diversity_of_gene/",sep = "")
# dir.create(alppath)

all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps.metf,group = "Group",Plot = TRUE )
head(tab)
data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)
# data$Inv_Simpson[is.na(data$Inv_Simpson)]
# data$Inv_Simpson %>% tail(1000)

result = MuiKwWlx2(data = data,num = 3:5)
result1 = FacetMuiPlotresultBox(data = data,num = 3:5,
                                          result = result,
                                          sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]+
  # scale_x_discrete(limits = axis_order) +
  # mytheme1 +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic() +
#  ggplot2::guides(fill = guide_legend(title = none))+
  theme(axis.title.y = element_text(angle = 90))


p1_1

res = FacetMuiPlotresultBar(data = data,num = c(3:(5)),
                                      result = result,
                                      sig_show ="abc",ncol = 4)
p1_2 = res[[1]]+
  # scale_x_discrete(limits = axis_order) +
  # mytheme1 +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  #  ggplot2::guides(fill = guide_legend(title = none))+
  theme(axis.title.y = element_text(angle = 90))
p1_2


res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:5),result = result,
                                     sig_show ="abc",ncol = 3)
p1_3 = res[[1]]+
  # scale_x_discrete(limits = axis_order) +
  # mytheme1 +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  #  ggplot2::guides(fill = guide_legend(title = none))+
  theme(axis.title.y = element_text(angle = 90))
p1_3

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=0.75, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.9)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))+
  # scale_x_discrete(limits = axis_order) +
  # mytheme1 +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  #  ggplot2::guides(fill = guide_legend(title = none))+
  theme(axis.title.y = element_text(angle = 90))

p1_0

# 保存alpha多样性分析结果
ggsave(file.path(alphapath, "alpha_diversity_box.png"), plot = p1_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(alphapath, "alpha_diversity_box.pdf"), plot = p1_1, width = 12, height = 8)
ggsave(file.path(alphapath, "alpha_diversity_bar.png"), plot = p1_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(alphapath, "alpha_diversity_bar.pdf"), plot = p1_2, width = 12, height = 8)
ggsave(file.path(alphapath, "alpha_diversity_boxbar.png"), plot = p1_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(alphapath, "alpha_diversity_boxbar.pdf"), plot = p1_3, width = 12, height = 8)
ggsave(file.path(alphapath, "alpha_diversity_violin.png"), plot = p1_0, width = 14, height = 8, dpi = 300)
ggsave(file.path(alphapath, "alpha_diversity_violin.pdf"), plot = p1_0, width = 14, height = 8)

# 保存alpha多样性表格
alpha_wb <- createWorkbook()
addWorksheet(alpha_wb, "alpha_diversity_data")
addWorksheet(alpha_wb, "alpha_diversity_stats")
writeData(alpha_wb, "alpha_diversity_data", data)
writeData(alpha_wb, "alpha_diversity_stats", result$sta)
saveWorkbook(alpha_wb, file.path(alphapath, "alpha_diversity_results.xlsx"), overwrite = TRUE)


#3 ordinate.metf: 排序分析#----------
# 创建beta多样性分析目录
betapath <- file.path(keggpath, "beta")
dir.create(betapath, recursive = TRUE)

result = ordinate.metf(ps = ps.metf,
                       group = "Group",
                       dist = "bray",
                       method = "PCA",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                       pair = FALSE)
p3_1 = result[[1]]
p3_1+
  scale_fill_manual(values = col.g)+
  theme_classic()+
  #  ggplot2::guides(fill = guide_legend(title = none))+
  theme(axis.title.y = element_text(angle = 90))

#带标签图形出图
p3_2 = result[[3]]
p3_2+
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_nature()+
  #  ggplot2::guides(fill = guide_legend(title = none))+
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
p3_3+
  scale_fill_manual(values = col.g)+
  theme_classic()+
  #  ggplot2::guides(fill = guide_legend(title = none))+
  theme(axis.title.y = element_text(angle = 90))

# 保存PCA分析结果
ggsave(file.path(betapath, "pca_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(betapath, "pca_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(betapath, "pca_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(betapath, "pca_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(betapath, "pca_enhanced.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(betapath, "pca_enhanced.pdf"), plot = p3_3, width = 10, height = 8)

# 保存beta多样性表格
beta_wb <- createWorkbook()
addWorksheet(beta_wb, "plotdata")
addWorksheet(beta_wb, "centroids")
addWorksheet(beta_wb, "segments")
writeData(beta_wb, "plotdata", plotdata)
writeData(beta_wb, "centroids", cent)
writeData(beta_wb, "segments", segs)


#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps.metf, Micromet = "adonis", dist = "bray")
dat1

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps.metf, Micromet = "MRPP", dist = "bray")
dat2

# 保存群落差异检测结果
addWorksheet(beta_wb, "adonis_test")
addWorksheet(beta_wb, "pairwise_test")
writeData(beta_wb, "adonis_test", dat1)
writeData(beta_wb, "pairwise_test", dat2)

#6 mantal.metf：群落功能差异检测普鲁士分析#------
result <- mantal.metf(ps = ps.metf,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

#保存图片
ggsave(file.path(betapath, "mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(betapath, "mantel_test.pdf"), plot = p3_7, width = 10, height = 8)

# mantal分析
addWorksheet(beta_wb, "mantel_test")
writeData(beta_wb, "mantel_test", data)
saveWorkbook(beta_wb, file.path(betapath, "beta_results.xlsx"), overwrite = TRUE)

#7 cluster.metf:样品聚类#-----
# 创建组成分析保存目录
clusterpath <- file.path(keggpath, "cluster")
dir.create(clusterpath, recursive = TRUE)

res = cluster.metf(ps= ps.metf,
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

#保存图片
ggsave(file.path(clusterpath, "cluster_heatmap.png"), plot = p4, width = 12, height = 10, dpi = 300)
ggsave(file.path(clusterpath, "cluster_heatmap.pdf"), plot = p4, width = 12, height = 10)
ggsave(file.path(clusterpath, "cluster_dendrogram.png"), plot = p4_1, width = 10, height = 6, dpi = 300)
ggsave(file.path(clusterpath, "cluster_dendrogram.pdf"), plot = p4_1, width = 10, height = 6)
ggsave(file.path(clusterpath, "cluster_groups.png"), plot = p4_2, width = 8, height = 6, dpi = 300)
ggsave(file.path(clusterpath, "cluster_groups.pdf"), plot = p4_2, width = 8, height = 6)

# 聚类表格保存
cluster_wb <- createWorkbook()
addWorksheet(cluster_wb, "cluster_data")
writeData(cluster_wb, "cluster_data", dat, rowNames = TRUE)
saveWorkbook(cluster_wb, file.path(clusterpath, "cluster_results.xlsx"), overwrite = TRUE)

#8 Micro_tern.metf: 三元图展示功能----
# 创建组成分析保存目录
comppath <- file.path(keggpath, "composition")
dir.create(comppath, recursive = TRUE)

ps.metf2 =  ps.metf %>% tax_glom("Level3")


res = Micro_tern.metf(ps= ps.metf2 %>% filter_OTU_ps(500),color = "Pathway"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图图片
ggsave(file.path(comppath, "ternary_plot.png"), plot = p15[[1]] + theme_classic(), width = 10, height = 8)
ggsave(file.path(comppath, "ternary_plot.pdf"), plot = p15[[1]] + theme_classic(), width = 10, height = 8)

# 保存三元图数据（保留行名）
comp_wb <- createWorkbook()
addWorksheet(comp_wb, "ternary_data")
writeData(comp_wb, "ternary_data", dat, rowNames = TRUE)


# function classfication
#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
rank_names(ps.metf)

result = barMainplot.metf(ps =ps.metf2%>%filter_OTU_ps(Top = 500),
                          j = "Pathway",
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1 <- p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2 <- p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Pathway","Abundance(%)")
head(databar)

ggsave(file.path(comppath, "barplot_main.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(comppath, "barplot_main.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(comppath, "barplot_summary.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(comppath, "barplot_summary.pdf"), plot = p4_2, width = 10, height = 8)

# 保存柱状图数据
addWorksheet(comp_wb, "barplot_data")
writeData(comp_wb, "barplot_data", databar, rowNames = TRUE)

#10 Ven.Upset.metf: 用于展示共有、特有的功能----
# 分组小于6时使用
res = Ven.Upset.metf(ps = ps.metf,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

ggsave(file.path(comppath, "venn_diagram.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(comppath, "venn_diagram.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(comppath, "upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(comppath, "upset_plot.pdf"), plot = p10.2, width = 12, height = 8)

# 保存Venn图数据
addWorksheet(comp_wb, "venn_data")
writeData(comp_wb, "venn_data", dat, rowNames = TRUE)


#12 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps.metf,
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

# 保存花瓣图图片
ggsave(file.path(comppath, "flower_plot.png"), plot = p13.1, width = 10, height = 8)
ggsave(file.path(comppath, "flower_plot.pdf"), plot = p13.1, width = 10, height = 8)

# 保存花瓣图数据
addWorksheet(comp_wb, "flower_data")
writeData(comp_wb, "flower_data", dat, rowNames = TRUE)


#13 ven.network.metf: ven网络展示共有特有功能----
map =  sample_data(ps.metf)
ps_ven=ps.metf
tax=data.frame(phyloseq:: tax_table(ps.metf))
colnames(tax)[[3]]="KOnumber"

tax_table(ps_ven) <- as.matrix(tax)

result =ven.network.metf(
  ps = ps_ven,
  N = 0.5,
  fill = "Pathway")

p14  = result[[1]] +
  theme(legend.position = "none")

p14
dat = result[[2]]
dat

# 保存Venn网络图片
ggsave(file.path(comppath, "venn_network.png"), plot = p14, width = 10, height = 8)
ggsave(file.path(comppath, "venn_network.pdf"), plot = p14, width = 10, height = 8)

# 保存Venn网络数据（保留行名）
addWorksheet(comp_wb, "venn_network_data")
writeData(comp_wb, "venn_network_data", dat, rowNames = TRUE)
saveWorkbook(comp_wb, file.path(comppath, "composition_results.xlsx"), overwrite = TRUE)

# function differential analysis----
# 创建差异分析主目录
diffpath <- file.path(keggpath, "differential")
dir.create(diffpath, recursive = TRUE)

# 创建差异分析总表
diff_wb <- createWorkbook()

#14 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps.metf,
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

# 保存DESeq2图片
ggsave(file.path(diffpath, "deseq2_volcano_1.png"), plot = p15.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "deseq2_volcano_1.pdf"), plot = p15.1, width = 10, height = 8)
ggsave(file.path(diffpath, "deseq2_volcano_2.png"), plot = p15.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "deseq2_volcano_2.pdf"), plot = p15.2, width = 10, height = 8)
ggsave(file.path(diffpath, "deseq2_volcano_3.png"), plot = p15.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "deseq2_volcano_3.pdf"), plot = p15.3, width = 10, height = 8)

# 保存DESeq2数据到总表
addWorksheet(diff_wb, "deseq2_results")
writeData(diff_wb, "deseq2_results", dat, rowNames = TRUE)

#15 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps.metf,
                       group  = "Group",
                       artGroup = NULL,
                       j="meta")
p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][2]
p16.2
p16.3 = res[[1]][3]
p16.3
dat = res[[2]]
dat

# 保存EdgeR图片
ggsave(file.path(diffpath, "edger_volcano_1.png"), plot = p16.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "edger_volcano_1.pdf"), plot = p16.1, width = 10, height = 8)
ggsave(file.path(diffpath, "edger_volcano_2.png"), plot = p16.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "edger_volcano_2.pdf"), plot = p16.2, width = 10, height = 8)
ggsave(file.path(diffpath, "edger_volcano_3.png"), plot = p16.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "edger_volcano_3.pdf"), plot = p16.3, width = 10, height = 8)

# 保存EdgeR数据到总表
addWorksheet(diff_wb, "edger_results")
writeData(diff_wb, "edger_results", dat, rowNames = TRUE)

#16 t.metf: 差异分析t检验----
res = t.metf(ps = ps.metf,group  = "Group",artGroup =NULL)
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

# 保存t检验图片
ggsave(file.path(diffpath, "ttest_main.png"), plot = p17, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "ttest_main.pdf"), plot = p17, width = 10, height = 8)
ggsave(file.path(diffpath, "ttest_volcano_1.png"), plot = p17.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "ttest_volcano_1.pdf"), plot = p17.1, width = 10, height = 8)
ggsave(file.path(diffpath, "ttest_volcano_2.png"), plot = p17.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "ttest_volcano_2.pdf"), plot = p17.2, width = 10, height = 8)
ggsave(file.path(diffpath, "ttest_volcano_3.png"), plot = p17.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "ttest_volcano_3.pdf"), plot = p17.3, width = 10, height = 8)

# 保存t检验数据到总表
addWorksheet(diff_wb, "ttest_results")
addWorksheet(diff_wb, "ttest_intersections")
writeData(diff_wb, "ttest_results", data, rowNames = TRUE)
writeData(diff_wb, "ttest_intersections", inter_union, rowNames = TRUE)

#17 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps.metf,group  = "Group",artGroup =NULL)
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

# 保存wlx非参数检验图片
ggsave(file.path(diffpath, "wilcox_main.png"), plot = p18, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "wilcox_main.pdf"), plot = p18, width = 10, height = 8)
ggsave(file.path(diffpath, "wilcox_volcano_1.png"), plot = p18.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "wilcox_volcano_1.pdf"), plot = p18.1, width = 10, height = 8)
ggsave(file.path(diffpath, "wilcox_volcano_2.png"), plot = p18.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "wilcox_volcano_2.pdf"), plot = p18.2, width = 10, height = 8)
ggsave(file.path(diffpath, "wilcox_volcano_3.png"), plot = p18.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(diffpath, "wilcox_volcano_3.pdf"), plot = p18.3, width = 10, height = 8)

# 保存wlx非参数检验数据到总表
addWorksheet(diff_wb, "wilcox_results")
addWorksheet(diff_wb, "wilcox_intersections")
writeData(diff_wb, "wilcox_results", data, rowNames = TRUE)
writeData(diff_wb, "wilcox_intersections", inter_union, rowNames = TRUE)

#18 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps.metf)$Group),2)
ps_sub <- subset_samples(ps.metf,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2

# 保存stamp差异分析图片
ggsave(file.path(diffpath, "stamp_analysis.png"), plot = p19, width = 12, height = 8, dpi = 300)
ggsave(file.path(diffpath, "stamp_analysis.pdf"), plot = p19, width = 12, height = 8)

# 保存stamp差异分析数据到总表
addWorksheet(diff_wb, "stamp_plot_data")
addWorksheet(diff_wb, "stamp_results")
writeData(diff_wb, "stamp_plot_data", dat1, rowNames = TRUE)
writeData(diff_wb, "stamp_results", dat2, rowNames = TRUE)

# 保存差异分析总表
saveWorkbook(diff_wb, file.path(diffpath, "differential_analysis_results.xlsx"), overwrite = TRUE)

# biobiomarker identification-----
# 创建生物标记物识别目录
biomarkerpath <- file.path(keggpath, "biomarker")
dir.create(biomarkerpath, recursive = TRUE)

# 创建生物标记物总表
biomarker_wb <- createWorkbook()

id = sample_data(ps.metf)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps.metf %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

#19 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)
result =rfcv.metf(ps = ps.metf %>% filter_OTU_ps(200),group  = "Group",optimal = 30,nrfcvnum = 6)
prfcv = result[[1]]
prfcv+theme_classic()
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存随机森林交叉验证图片
ggsave(file.path(biomarkerpath, "random_forest_cv.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "random_forest_cv.pdf"), plot = prfcv, width = 10, height = 8)

# 保存随机森林交叉验证数据到总表
addWorksheet(biomarker_wb, "rfcv_plot_data")
addWorksheet(biomarker_wb, "rfcv_table")
writeData(biomarker_wb, "rfcv_plot_data", result[[2]], rowNames = TRUE)
writeData(biomarker_wb, "rfcv_table", rfcvtable, rowNames = TRUE)

#20 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1+theme_classic()
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

# 保存ROC曲线图片
ggsave(file.path(biomarkerpath, "roc_curve.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "roc_curve.pdf"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "roc_boxplot.png"), plot = p33.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "roc_boxplot.pdf"), plot = p33.2, width = 10, height = 8)

# 保存ROC数据到总表
addWorksheet(biomarker_wb, "roc_data")
writeData(biomarker_wb, "roc_data", dat, rowNames = TRUE)

#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.metf,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存载荷矩阵PCA图片
ggsave(file.path(biomarkerpath, "loading_pca.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "loading_pca.pdf"), plot = p34.1, width = 10, height = 8)

# 保存载荷矩阵PCA数据到总表
addWorksheet(biomarker_wb, "loading_pca_data")
writeData(biomarker_wb, "loading_pca_data", dat, rowNames = TRUE)

#22 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = ps.metf %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果数据到总表
addWorksheet(biomarker_wb, "svm_auc")
addWorksheet(biomarker_wb, "svm_importance")
writeData(biomarker_wb, "svm_auc", AUC, rowNames = TRUE)
writeData(biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#23 glm.metf:glm筛选特征功能----
pst = subset_samples(ps.metf,Group %in% c("KO" ,"OE"))
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)

AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果数据到总表
addWorksheet(biomarker_wb, "glm_auc")
addWorksheet(biomarker_wb, "glm_importance")
writeData(biomarker_wb, "glm_auc", AUC, rowNames = TRUE)
writeData(biomarker_wb, "glm_importance", importance, rowNames = TRUE)

#24 xgboost.metf: xgboost筛选特征功能----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)
res = xgboost.metf(ps =ps.metf, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果数据到总表
addWorksheet(biomarker_wb, "xgboost_accuracy")
addWorksheet(biomarker_wb, "xgboost_importance")
writeData(biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
writeData(biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#25 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果数据到总表
addWorksheet(biomarker_wb, "lasso_accuracy")
addWorksheet(biomarker_wb, "lasso_importance")
writeData(biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
writeData(biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#26 decisiontree.micro: ----
library(dplyr)
library(caret)
library(rpart)
res = decisiontree.metf(ps=ps.metf, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果数据到总表
addWorksheet(biomarker_wb, "decision_tree_accuracy")
addWorksheet(biomarker_wb, "decision_tree_importance")
writeData(biomarker_wb, "decision_tree_accuracy", accuracy, rowNames = TRUE)
writeData(biomarker_wb, "decision_tree_importance", importance, rowNames = TRUE)

#27 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps.metf, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果数据到总表
addWorksheet(biomarker_wb, "naive_bayes_accuracy")
addWorksheet(biomarker_wb, "naive_bayes_importance")
writeData(biomarker_wb, "naive_bayes_accuracy", accuracy, rowNames = TRUE)
writeData(biomarker_wb, "naive_bayes_importance", importance, rowNames = TRUE)

#28 LDA.metf: LDA筛选特征功能-----
tablda = LDA.metf(ps = ps.metf,
                  Top = 100,
                  p.lvl = 0.05,
                  lda.lvl = 1,
                  seed = 11,
                  adjust.p = F)
tablda[[1]]

p35 <- lefse_bar(taxtree = tablda[[2]])+ theme_classic()+ ylab("")
p35
dat = tablda[[2]]
dat

# 保存LDA图片
ggsave(file.path(biomarkerpath, "lda_lefse.png"), plot = p35, width = 12, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "lda_lefse.pdf"), plot = p35, width = 12, height = 8)

# 保存LDA数据到总表
addWorksheet(biomarker_wb, "lda_summary")
addWorksheet(biomarker_wb, "lda_detailed_results")
writeData(biomarker_wb, "lda_summary", tablda[[1]], rowNames = TRUE)
writeData(biomarker_wb, "lda_detailed_results", dat, rowNames = TRUE)

#29 randomforest.metf: 随机森林筛选特征功能----
rank_names(ps.metf)
res = randomforest.metf( ps = ps.metf,group  = "Group", optimal = 50 ,fill ="Level1")
p42.1 = res[[1]]
p42.1+theme_classic()+ ylab("")
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林图片
ggsave(file.path(biomarkerpath, "random_forest_importance.png"), plot = p42.1, width = 12, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "random_forest_importance.pdf"), plot = p42.1, width = 12, height = 8)
ggsave(file.path(biomarkerpath, "random_forest_error.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "random_forest_error.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(biomarkerpath, "random_forest_prediction.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(biomarkerpath, "random_forest_prediction.pdf"), plot = p42.4, width = 10, height = 8)

# 保存随机森林数据到总表
addWorksheet(biomarker_wb, "random_forest_data")
writeData(biomarker_wb, "random_forest_data", dat, rowNames = TRUE)

#30 nnet.metf: 神经网络筛选特征功能  ------
res =nnet.metf(ps=ps.metf, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果数据到总表
addWorksheet(biomarker_wb, "neural_network_accuracy")
addWorksheet(biomarker_wb, "neural_network_importance")
writeData(biomarker_wb, "neural_network_accuracy", accuracy, rowNames = TRUE)
writeData(biomarker_wb, "neural_network_importance", importance, rowNames = TRUE)

#31 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps =  ps.metf, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果数据到总表
addWorksheet(biomarker_wb, "bagging_accuracy")
addWorksheet(biomarker_wb, "bagging_importance")
writeData(biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
writeData(biomarker_wb, "bagging_importance", importance, rowNames = TRUE)

# 保存生物标记物总表
saveWorkbook(biomarker_wb, file.path(biomarkerpath, "biomarker_results.xlsx"), overwrite = TRUE)

#network analysis -----
# 创建网络分析目录
networkpath <- file.path(keggpath, "network")
dir.create(networkpath, recursive = TRUE)

# 创建网络分析总表
network_wb <- createWorkbook()

#32 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
detach("package:mia", unload = TRUE)
tab.r = network.pip.metf(
  ps = ps.metf,
  N = 50,
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
  ram.net =FALSE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R=10,
  ncpus = 1,
  fill = "Level1"
)
dat = tab.r[[2]]

cortab = dat$net.cor.matrix$cortab

#-提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]
p0

# 保存网络图
ggsave(file.path(networkpath, "network_plot.png"), plot = p0, width = 12, height = 12, dpi = 300)
ggsave(file.path(networkpath, "network_plot.pdf"), plot = p0, width = 12, height = 12)

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

# 保存网络属性数据到总表
addWorksheet(network_wb, "network_properties")
writeData(network_wb, "network_properties", dat2, rowNames = TRUE)

#35 node_properties:计算节点属性#---------
for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  nodepro = node_properties(igraph) %>% as.data.frame()
  nodepro$Group = id[i]
  head(nodepro)
  colnames(nodepro) = paste0(colnames(nodepro),".",id[i])
  nodepro = nodepro %>%
    as.data.frame() %>%
    rownames_to_column("KO.name")


  # head(dat.f)
  if (i == 1) {
    nodepro2 = nodepro
  } else{
    nodepro2 = nodepro2 %>% full_join(nodepro,by = "KO.name")
  }
}
head(nodepro2)

# 保存节点属性数据到总表
addWorksheet(network_wb, "node_properties")
writeData(network_wb, "node_properties", nodepro2, rowNames = TRUE)

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

# 保存网络比较数据到总表
addWorksheet(network_wb, "network_comparison")
writeData(network_wb, "network_comparison", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(network_wb, file.path(networkpath, "network_analysis_results.xlsx"), overwrite = TRUE)

# pathway enrich------
# 创建通路富集分析目录
enrichpath <- file.path(keggpath, "enrichment")
dir.create(enrichpath, recursive = TRUE)

# 创建富集分析总表
enrich_wb <- createWorkbook()

#37 KEGG_enrich2:KEGG富集分析#---------
library(GO.db)
library(DOSE)
library(GO.db)
library(GSEABase)
library(ggtree)
library(aplot)
library(clusterProfiler)
library("GSVA")

# 调用差异分析的结果
res = EdgerSuper.metf(ps = ps.metf,
                      group  = "Group",
                      artGroup = NULL)
dat = res[[2]]
dat

res2 = KEGG_enrich.metf(ps = ps.metf,
                        #  diffpath = diffpath,
                        dif = dat
)
dat1= res2$`KO-OE`
dat2= res2$`KO-WT`
dat3= res2$`OE-WT`

# 保存KEGG富集分析数据到总表
addWorksheet(enrich_wb, "KO_vs_OE")
addWorksheet(enrich_wb, "KO_vs_WT")
addWorksheet(enrich_wb, "OE_vs_WT")
writeData(enrich_wb, "KO_vs_OE", dat1, rowNames = TRUE)
writeData(enrich_wb, "KO_vs_WT", dat2, rowNames = TRUE)
writeData(enrich_wb, "OE_vs_WT", dat3, rowNames = TRUE)

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

result = buplot.micro(dt  =dat1,id = id)
p1 = result[[1]]
p1
p2 = result[[2]]
p2

# 保存富集分析气泡图
ggsave(file.path(enrichpath, "enrichment_bubble_plot1.png"), plot = p1, width = 12, height = 10, dpi = 300)
ggsave(file.path(enrichpath, "enrichment_bubble_plot1.pdf"), plot = p1, width = 12, height = 10)
ggsave(file.path(enrichpath, "enrichment_bubble_plot2.png"), plot = p2, width = 12, height = 10, dpi = 300)
ggsave(file.path(enrichpath, "enrichment_bubble_plot2.pdf"), plot = p2, width = 12, height = 10)

# 保存气泡图数据到总表
addWorksheet(enrich_wb, "bubble_plot_data1")
addWorksheet(enrich_wb, "bubble_plot_data2")
writeData(enrich_wb, "bubble_plot_data1", result[[3]], rowNames = TRUE)
writeData(enrich_wb, "bubble_plot_data2", result[[4]], rowNames = TRUE)

#39 GSVA.metf:GSVA富集分析--------
library(limma)
res= GSVA.metf(ps = ps.metf,
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

# 保存GSVA图片
ggsave(file.path(enrichpath, "gsva_heatmap.png"), plot = p39.1, width = 12, height = 10, dpi = 300)
ggsave(file.path(enrichpath, "gsva_heatmap.pdf"), plot = p39.1, width = 12, height = 10)
ggsave(file.path(enrichpath, "gsva_comparison1.png"), plot = p39.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(enrichpath, "gsva_comparison1.pdf"), plot = p39.2, width = 10, height = 8)

#40 kegg_function:按照kegg通路合并基因-------
detach(package:mia, unload = TRUE)
detach(package:ANCOMBC, unload = TRUE)
ps.kegg.function = kegg_function(ps = ps.metf)

# 提取OTU表和税分类表
otu_kegg_func <- as.data.frame(phyloseq::otu_table(ps.kegg.function))
tax_kegg_func <- as.data.frame(phyloseq::tax_table(ps.kegg.function))

# 保存KEGG通路合并基因数据到总表
addWorksheet(enrich_wb, "kegg_function_data")
addWorksheet(enrich_wb, "kegg_function_taxonomy")
writeData(enrich_wb, "kegg_function_data", otu_kegg_func, rowNames = TRUE)
writeData(enrich_wb, "kegg_function_taxonomy", tax_kegg_func, rowNames = TRUE)

#41 基于Mkegg通路合并基因-------
ps.mkegg.function =mkegg_function(ps = ps.metf)

# 提取OTU表和税分类表
otu_mkegg_func <- as.data.frame(phyloseq::otu_table(ps.mkegg.function))
tax_mkegg_func <- as.data.frame(phyloseq::tax_table(ps.mkegg.function))

# 保存Mkegg通路合并基因数据到总表
addWorksheet(enrich_wb, "mkegg_function_data")
addWorksheet(enrich_wb, "mkegg_function_taxonomy")
writeData(enrich_wb, "mkegg_function_data", otu_mkegg_func, rowNames = TRUE)
writeData(enrich_wb, "mkegg_function_taxonomy", tax_mkegg_func, rowNames = TRUE)

#42 基于reaction合并基因-------
ps_kegg_function.reaction = reaction_function(ps = ps.metf)

# 提取OTU表和税分类表
otu_reaction_func <- as.data.frame(phyloseq::otu_table(ps_kegg_function.reaction))
tax_reaction_func <- as.data.frame(phyloseq::tax_table(ps_kegg_function.reaction))

# 保存reaction合并基因数据到总表
addWorksheet(enrich_wb, "reaction_function_data")
addWorksheet(enrich_wb, "reaction_function_taxonomy")
writeData(enrich_wb, "reaction_function_data", otu_reaction_func, rowNames = TRUE)
writeData(enrich_wb, "reaction_function_taxonomy", tax_reaction_func, rowNames = TRUE)

# 保存富集分析总表
saveWorkbook(enrich_wb, file.path(enrichpath, "enrichment_results.xlsx"), overwrite = TRUE)

# cnps 循环数据库-----
# 创建CNPS循环分析目录
cnpspath <- file.path(path, "cnps_cycling")
dir.create(cnpspath, recursive = TRUE)

# 创建CNPS总表
cnps_wb <- createWorkbook()

# data prepration------
#1 cnps_gene: KEGG数据库中过滤CNPS 元素循环相关基因------
res  = cnps_gene(ps=ps.metf)
C_gene <- res$C_gene
N_gene <- res$N_gene
P_gene <- res$P_gene
S_gene <- res$S_gene

# 保存CNPS基因数据到总表
addWorksheet(cnps_wb, "C_genes")
addWorksheet(cnps_wb, "N_genes")
addWorksheet(cnps_wb, "P_genes")
addWorksheet(cnps_wb, "S_genes")
writeData(cnps_wb, "C_genes", C_gene, rowNames = TRUE)
writeData(cnps_wb, "N_genes", N_gene, rowNames = TRUE)
writeData(cnps_wb, "P_genes", P_gene, rowNames = TRUE)
writeData(cnps_wb, "S_genes", S_gene, rowNames = TRUE)

#2 cnps_ps: KEGG数据库中过滤CNPS 元素循环相关基因构建ps -----
ps.cnps = cnps_ps(psko = ps.metf)
ps.cnps
phyloseq::tax_table(ps.cnps) %>% colnames()
otu_table(ps.cnps)

# 保存CNPS phyloseq对象数据到总表
otu_cnps <- as.data.frame(phyloseq::otu_table(ps.cnps))
tax_cnps <- as.data.frame(phyloseq::tax_table(ps.cnps))
sample_cnps <- as.data.frame(phyloseq::sample_data(ps.cnps))
addWorksheet(cnps_wb, "cnps_otu_table")
addWorksheet(cnps_wb, "cnps_taxonomy")
addWorksheet(cnps_wb, "cnps_sample_data")
writeData(cnps_wb, "cnps_otu_table", otu_cnps, rowNames = TRUE)
writeData(cnps_wb, "cnps_taxonomy", tax_cnps, rowNames = TRUE)
writeData(cnps_wb, "cnps_sample_data", sample_cnps, rowNames = TRUE)

# 保存CNPS总表
saveWorkbook(cnps_wb, file.path(cnpspath, "cnps_data.xlsx"), overwrite = TRUE)

# function diversity -----
# 创建CNPS alpha多样性目录
cnps_alphapath <- file.path(cnpspath, "alpha")
dir.create(cnps_alphapath, recursive = TRUE)

# 创建CNPS alpha多样性总表
cnps_alpha_wb <- createWorkbook()

#3 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE" )
#--alpha多样性指标运算
tab = alpha.metf(ps = ps.cnps,group = "Group",Plot = TRUE )
head(tab)

data = cbind(data.frame(ID = 1:length(tab$Group),group = tab$Group),tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

result =MuiKwWlx2(data = data,num = 3:5)
result1 = FacetMuiPlotresultBox(data = data,num = 3:5,
                                result = result,
                                sig_show ="abc",ncol = 4 )
p1_1 = result1[[1]]
p1_1 = p1_1+
  scale_x_discrete(limits = axis_order) +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

# 保存CNPS alpha多样性图片
ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_box.png"), plot = p1_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_box.pdf"), plot = p1_1, width = 12, height = 8)

res = FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2 = p1_2+
  scale_x_discrete(limits = axis_order) +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_bar.png"), plot = p1_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_bar.pdf"), plot = p1_2, width = 12, height = 8)

res = FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3 = p1_3+
  scale_x_discrete(limits = axis_order) +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_boxbar.png"), plot = p1_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_boxbar.pdf"), plot = p1_3, width = 12, height = 8)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  geom_text(aes(x=group , y=y ,label=stat))
p1_0 = p1_0+guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_violin.png"), plot = p1_0, width = 14, height = 8, dpi = 300)
ggsave(file.path(cnps_alphapath, "cnps_alpha_diversity_violin.pdf"), plot = p1_0, width = 14, height = 8)

# 保存CNPS alpha多样性数据到总表
addWorksheet(cnps_alpha_wb, "cnps_alpha_data")
addWorksheet(cnps_alpha_wb, "cnps_alpha_stats")
writeData(cnps_alpha_wb, "cnps_alpha_data", data, rowNames = TRUE)
writeData(cnps_alpha_wb, "cnps_alpha_stats", result$sta, rowNames = TRUE)

#4 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps.cnps))/10
result = alpha_rare.metf(ps = ps.cnps, group = "Group", method = "Richness", start = 100, step = rare)

p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)

ggsave(file.path(cnps_alphapath, "cnps_rarefaction_individual.png"), plot = p2_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_alphapath, "cnps_rarefaction_individual.pdf"), plot = p2_1, width = 10, height = 8)

raretab <- result[[2]]
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)

ggsave(file.path(cnps_alphapath, "cnps_rarefaction_group.png"), plot = p2_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_alphapath, "cnps_rarefaction_group.pdf"), plot = p2_2, width = 10, height = 8)

# 保存稀释曲线数据到总表
addWorksheet(cnps_alpha_wb, "rarefaction_data")
writeData(cnps_alpha_wb, "rarefaction_data", raretab, rowNames = TRUE)

# 保存CNPS alpha多样性总表
saveWorkbook(cnps_alpha_wb, file.path(cnps_alphapath, "cnps_alpha_diversity.xlsx"), overwrite = TRUE)

# 创建CNPS beta多样性目录
cnps_betapath <- file.path(cnpspath, "beta")
dir.create(cnps_betapath, recursive = TRUE)

# 创建CNPS beta多样性总表
cnps_beta_wb <- createWorkbook()

#5 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps.cnps,
                       group = "Group",
                       dist = "bray",
                       method = "NMDS",
                       Micromet = "anosim",
                       pvalue.cutoff = 0.05,
                       pair = FALSE)
p3_1 = result[[1]]
p3_1 = p3_1 +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

ggsave(file.path(cnps_betapath, "cnps_nmds_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_betapath, "cnps_nmds_basic.pdf"), plot = p3_1, width = 10, height = 8)

#带标签图形出图
p3_2 = result[[3]]
p3_2 = p3_2+
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

ggsave(file.path(cnps_betapath, "cnps_nmds_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_betapath, "cnps_nmds_labeled.pdf"), plot = p3_2, width = 10, height = 8)

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
p3_3 = p3_3+
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

ggsave(file.path(cnps_betapath, "cnps_nmds_enhanced.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_betapath, "cnps_nmds_enhanced.pdf"), plot = p3_3, width = 10, height = 8)

# 保存CNPS beta多样性数据到总表
addWorksheet(cnps_beta_wb, "nmds_plotdata")
addWorksheet(cnps_beta_wb, "nmds_centroids")
addWorksheet(cnps_beta_wb, "nmds_segments")
writeData(cnps_beta_wb, "nmds_plotdata", plotdata, rowNames = TRUE)
writeData(cnps_beta_wb, "nmds_centroids", cent, rowNames = TRUE)
writeData(cnps_beta_wb, "nmds_segments", segs, rowNames = TRUE)

#6 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps.cnps, Micromet = "adonis", dist = "bray")
dat1

#7 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps =  ps.cnps, Micromet = "MRPP", dist = "bray")
dat2

# 保存CNPS群落差异检测结果到总表
addWorksheet(cnps_beta_wb, "adonis_test")
addWorksheet(cnps_beta_wb, "pairwise_mrpp")
writeData(cnps_beta_wb, "adonis_test", dat1, rowNames = TRUE)
writeData(cnps_beta_wb, "pairwise_mrpp", dat2, rowNames = TRUE)

#8 mantal.metf：群落功能差异检测普鲁士分析#------
result <- mantal.metf(ps =  ps.cnps,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

ggsave(file.path(cnps_betapath, "cnps_mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_betapath, "cnps_mantel_test.pdf"), plot = p3_7, width = 10, height = 8)

# 保存Mantel检验结果到总表
addWorksheet(cnps_beta_wb, "mantel_test")
writeData(cnps_beta_wb, "mantel_test", data, rowNames = TRUE)

# 保存CNPS beta多样性总表
saveWorkbook(cnps_beta_wb, file.path(cnps_betapath, "cnps_beta_diversity.xlsx"), overwrite = TRUE)

# 创建CNPS聚类分析目录
cnps_clusterpath <- file.path(cnpspath, "cluster")
dir.create(cnps_clusterpath, recursive = TRUE)

# 创建CNPS聚类总表
cnps_cluster_wb <- createWorkbook()

#9 cluster.metf:样品聚类#-----
res = cluster.metf(ps= ps.cnps,
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

# 保存CNPS聚类图片
ggsave(file.path(cnps_clusterpath, "cnps_cluster_heatmap.png"), plot = p4, width = 12, height = 10, dpi = 300)
ggsave(file.path(cnps_clusterpath, "cnps_cluster_heatmap.pdf"), plot = p4, width = 12, height = 10)
ggsave(file.path(cnps_clusterpath, "cnps_cluster_dendrogram.png"), plot = p4_1, width = 10, height = 6, dpi = 300)
ggsave(file.path(cnps_clusterpath, "cnps_cluster_dendrogram.pdf"), plot = p4_1, width = 10, height = 6)
ggsave(file.path(cnps_clusterpath, "cnps_cluster_groups.png"), plot = p4_2, width = 8, height = 6, dpi = 300)
ggsave(file.path(cnps_clusterpath, "cnps_cluster_groups.pdf"), plot = p4_2, width = 8, height = 6)

# 保存CNPS聚类数据到总表
addWorksheet(cnps_cluster_wb, "cluster_data")
writeData(cnps_cluster_wb, "cluster_data", dat, rowNames = TRUE)

# 保存CNPS聚类总表
saveWorkbook(cnps_cluster_wb, file.path(cnps_clusterpath, "cnps_cluster_results.xlsx"), overwrite = TRUE)

# 创建CNPS组成分析目录
cnps_comppath <- file.path(cnpspath, "composition")
dir.create(cnps_comppath, recursive = TRUE)

# 创建CNPS组成分析总表
cnps_comp_wb <- createWorkbook()

#10 Micro_tern.metf: 三元图展示功能----
res = Micro_tern.metf(ps= ps.cnps %>% filter_OTU_ps(500),
                      color = "Group"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat
phyloseq::tax_table(ps.cnps)

# 保存CNPS三元图
ggsave(file.path(cnps_comppath, "cnps_ternary_plot.png"), plot = p15, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_comppath, "cnps_ternary_plot.pdf"), plot = p15, width = 10, height = 8)

# 保存三元图数据到总表
addWorksheet(cnps_comp_wb, "ternary_data")
writeData(cnps_comp_wb, "ternary_data", dat, rowNames = TRUE)

#11 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps.cnps,
                          j =  "module" ,
                          # axis_ord = axis_order,
                          label = FALSE,
                          sd = FALSE,
                          Top =10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1 <- p4_1+mytheme1

ggsave(file.path(cnps_comppath, "cnps_barplot_main.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_comppath, "cnps_barplot_main.pdf"), plot = p4_1, width = 12, height = 8)

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2 <- p4_2+mytheme1

ggsave(file.path(cnps_comppath, "cnps_barplot_summary.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_comppath, "cnps_barplot_summary.pdf"), plot = p4_2, width = 10, height = 8)

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Function","Abundance(%)")
head(databar)

# 保存柱状图数据到总表
addWorksheet(cnps_comp_wb, "barplot_data")
writeData(cnps_comp_wb, "barplot_data", databar, rowNames = TRUE)

#12 Ven.Upset.metf: 用于展示共有、特有的功能 ----
# 分组小于6时使用
library(dplyr)
res = Ven.Upset.metf(ps =  ps.cnps,
                     group = "Group",
                     N = 1,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

# 保存CNPS Venn图
ggsave(file.path(cnps_comppath, "cnps_venn_diagram.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_comppath, "cnps_venn_diagram.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(cnps_comppath, "cnps_upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_comppath, "cnps_upset_plot.pdf"), plot = p10.2, width = 12, height = 8)

# 保存Venn图数据到总表
addWorksheet(cnps_comp_wb, "venn_data")
writeData(cnps_comp_wb, "venn_data", dat, rowNames = TRUE)

#14 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps.cnps ,
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

# 保存CNPS花瓣图
ggsave(file.path(cnps_comppath, "cnps_flower_plot.png"), plot = p13.1, width = 10, height = 10, dpi = 300)
ggsave(file.path(cnps_comppath, "cnps_flower_plot.pdf"), plot = p13.1, width = 10, height = 10)

# 保存花瓣图数据到总表
addWorksheet(cnps_comp_wb, "flower_data")
writeData(cnps_comp_wb, "flower_data", dat, rowNames = TRUE)

#15 ven.network.metf: ven网络展示共有特有功能----
result =ven.network.metf(
  ps = ps.cnps,
  N = 0.5,
  fill = "Group")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat

# 保存CNPS Venn网络图
ggsave(file.path(cnps_comppath, "cnps_venn_network.png"), plot = p14, width = 12, height = 10, dpi = 300)
ggsave(file.path(cnps_comppath, "cnps_venn_network.pdf"), plot = p14, width = 12, height = 10)

# 保存Venn网络数据到总表
addWorksheet(cnps_comp_wb, "venn_network_data")
writeData(cnps_comp_wb, "venn_network_data", dat, rowNames = TRUE)

# 保存CNPS组成分析总表
saveWorkbook(cnps_comp_wb, file.path(cnps_comppath, "cnps_composition_results.xlsx"), overwrite = TRUE)

# function differential analysis----
# 创建CNPS差异分析目录
cnps_diffpath <- file.path(cnpspath, "differential")
dir.create(cnps_diffpath, recursive = TRUE)

# 创建CNPS差异分析总表
cnps_diff_wb <- createWorkbook()

#16 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps.cnps,
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

# 保存CNPS DESeq2图片
ggsave(file.path(cnps_diffpath, "cnps_deseq2_volcano1.png"), plot = p15.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_deseq2_volcano1.pdf"), plot = p15.1, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_deseq2_volcano2.png"), plot = p15.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_deseq2_volcano2.pdf"), plot = p15.2, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_deseq2_volcano3.png"), plot = p15.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_deseq2_volcano3.pdf"), plot = p15.3, width = 10, height = 8)

# 保存DESeq2数据到总表
addWorksheet(cnps_diff_wb, "cnps_deseq2_results")
writeData(cnps_diff_wb, "cnps_deseq2_results", dat, rowNames = TRUE)

#17 EdgerSuper.metf:EdgeR计算差异功能基因----
res = EdgerSuper.metf (ps = ps.cnps,
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

# 保存CNPS EdgeR图片
ggsave(file.path(cnps_diffpath, "cnps_edger_volcano1.png"), plot = p16.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_edger_volcano1.pdf"), plot = p16.1, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_edger_volcano2.png"), plot = p16.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_edger_volcano2.pdf"), plot = p16.2, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_edger_volcano3.png"), plot = p16.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_edger_volcano3.pdf"), plot = p16.3, width = 10, height = 8)

# 保存EdgeR数据到总表
addWorksheet(cnps_diff_wb, "cnps_edger_results")
writeData(cnps_diff_wb, "cnps_edger_results", dat, rowNames = TRUE)

#18 t.metf: 差异分析t检验----
res = t.metf(ps = ps.cnps,group  = "Group",artGroup =NULL)
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

# 保存CNPS t检验图片
ggsave(file.path(cnps_diffpath, "cnps_ttest_main.png"), plot = p17, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_ttest_main.pdf"), plot = p17, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_ttest_volcano1.png"), plot = p17.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_ttest_volcano1.pdf"), plot = p17.1, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_ttest_volcano2.png"), plot = p17.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_ttest_volcano2.pdf"), plot = p17.2, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_ttest_volcano3.png"), plot = p17.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_ttest_volcano3.pdf"), plot = p17.3, width = 10, height = 8)

# 保存t检验数据到总表
addWorksheet(cnps_diff_wb, "cnps_ttest_results")
addWorksheet(cnps_diff_wb, "cnps_ttest_intersections")
writeData(cnps_diff_wb, "cnps_ttest_results", data, rowNames = TRUE)
writeData(cnps_diff_wb, "cnps_ttest_intersections", inter_union, rowNames = TRUE)

#19 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps.cnps,group  = "Group",artGroup =NULL)
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

# 保存CNPS Wilcoxon检验图片
ggsave(file.path(cnps_diffpath, "cnps_wilcox_main.png"), plot = p18, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_wilcox_main.pdf"), plot = p18, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_wilcox_volcano1.png"), plot = p18.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_wilcox_volcano1.pdf"), plot = p18.1, width = 10, height = 8)
ggsave(file.path(cnps_diffpath, "cnps_wilcox_volcano2.png"), plot = p18.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_wilcox_volcano2.pdf"), plot = p18.2, width = 10, height = 8, dpi = 300)

#20 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps.cnps)$Group),2)
ps_sub <- subset_samples(ps.cnps,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2

# 保存CNPS STAMP图片
ggsave(file.path(cnps_diffpath, "cnps_stamp_analysis.png"), plot = p19, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_diffpath, "cnps_stamp_analysis.pdf"), plot = p19, width = 12, height = 8)

# 保存STAMP数据到总表
addWorksheet(cnps_diff_wb, "cnps_stamp_plot_data")
addWorksheet(cnps_diff_wb, "cnps_stamp_results")
writeData(cnps_diff_wb, "cnps_stamp_plot_data", dat1, rowNames = TRUE)
writeData(cnps_diff_wb, "cnps_stamp_results", dat2, rowNames = TRUE)

# 保存CNPS差异分析总表
saveWorkbook(cnps_diff_wb, file.path(cnps_diffpath, "cnps_differential_results.xlsx"), overwrite = TRUE)

# biobiomarker identification-----
# 创建CNPS生物标记物目录
cnps_biomarkerpath <- file.path(cnpspath, "biomarker")
dir.create(cnps_biomarkerpath, recursive = TRUE)

# 创建CNPS生物标记物总表
cnps_biomarker_wb <- createWorkbook()

id = sample_data(ps.cnps)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps.cnps %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

#21 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps.cnps %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存随机森林交叉验证图片
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_cv.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_cv.pdf"), plot = prfcv, width = 10, height = 8)

# 保存随机森林交叉验证数据到总表
addWorksheet(cnps_biomarker_wb, "rfcv_plot_data")
addWorksheet(cnps_biomarker_wb, "rfcv_table")
writeData(cnps_biomarker_wb, "rfcv_plot_data", result[[2]], rowNames = TRUE)
writeData(cnps_biomarker_wb, "rfcv_table", rfcvtable, rowNames = TRUE)

#22 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

# 保存ROC曲线图片
ggsave(file.path(cnps_biomarkerpath, "cnps_roc_curve.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_roc_curve.pdf"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(cnps_biomarkerpath, "cnps_roc_boxplot.png"), plot = p33.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_roc_boxplot.pdf"), plot = p33.2, width = 10, height = 8)

# 保存ROC数据到总表
addWorksheet(cnps_biomarker_wb, "roc_data")
writeData(cnps_biomarker_wb, "roc_data", dat, rowNames = TRUE)

#23 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.cnps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存载荷矩阵PCA图片
ggsave(file.path(cnps_biomarkerpath, "cnps_loading_pca.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_loading_pca.pdf"), plot = p34.1, width = 10, height = 8)

# 保存载荷矩阵PCA数据到总表
addWorksheet(cnps_biomarker_wb, "loading_pca_data")
writeData(cnps_biomarker_wb, "loading_pca_data", dat, rowNames = TRUE)

#24 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果数据到总表
addWorksheet(cnps_biomarker_wb, "svm_auc")
addWorksheet(cnps_biomarker_wb, "svm_importance")
writeData(cnps_biomarker_wb, "svm_auc", AUC, rowNames = TRUE)
writeData(cnps_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#25 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = ps.cnps%>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果数据到总表
addWorksheet(cnps_biomarker_wb, "glm_auc")
addWorksheet(cnps_biomarker_wb, "glm_importance")
writeData(cnps_biomarker_wb, "glm_auc", AUC, rowNames = TRUE)
writeData(cnps_biomarker_wb, "glm_importance", importance, rowNames = TRUE)

#26 xgboost.metf: xgboost筛选特征功能----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)
res = xgboost.metf(ps =ps.cnps, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果数据到总表
addWorksheet(cnps_biomarker_wb, "xgboost_accuracy")
addWorksheet(cnps_biomarker_wb, "xgboost_importance")
writeData(cnps_biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
writeData(cnps_biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#27 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  ps.cnps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果数据到总表
addWorksheet(cnps_biomarker_wb, "lasso_accuracy")
addWorksheet(cnps_biomarker_wb, "lasso_importance")
writeData(cnps_biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
writeData(cnps_biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#28 decisiontree.metf: ----
library(rpart)
res =decisiontree.metf(ps=ps.cnps, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果数据到总表
addWorksheet(cnps_biomarker_wb, "decision_tree_accuracy")
addWorksheet(cnps_biomarker_wb, "decision_tree_importance")
writeData(cnps_biomarker_wb, "decision_tree_accuracy", accuracy, rowNames = TRUE)
writeData(cnps_biomarker_wb, "decision_tree_importance", importance, rowNames = TRUE)

#29 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps.cnps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果数据到总表
addWorksheet(cnps_biomarker_wb, "naive_bayes_accuracy")
addWorksheet(cnps_biomarker_wb, "naive_bayes_importance")
writeData(cnps_biomarker_wb, "naive_bayes_accuracy", accuracy, rowNames = TRUE)
writeData(cnps_biomarker_wb, "naive_bayes_importance", importance, rowNames = TRUE)

#30 LDA.metf: LDA筛选特征功能-----
tablda = LDA.metf(ps = ps.cnps,
                  Top = 100,
                  p.lvl = 0.05,
                  lda.lvl = 1,
                  seed = 11,
                  adjust.p = F)

p35 <- lefse_bar(taxtree = tablda[[2]])
p35
dat = tablda[[2]]
dat

# 保存LDA图片
ggsave(file.path(cnps_biomarkerpath, "cnps_lda_lefse.png"), plot = p35, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_lda_lefse.pdf"), plot = p35, width = 12, height = 8)

# 保存LDA数据到总表
addWorksheet(cnps_biomarker_wb, "lda_summary")
addWorksheet(cnps_biomarker_wb, "lda_detailed_results")
writeData(cnps_biomarker_wb, "lda_summary", tablda[[1]], rowNames = TRUE)
writeData(cnps_biomarker_wb, "lda_detailed_results", dat, rowNames = TRUE)

#31 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps.cnps,group  = "Group", optimal = 50 ,fill = "module" )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林图片
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_importance.png"), plot = p42.1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_importance.pdf"), plot = p42.1, width = 12, height = 8)
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_error.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_error.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_prediction.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(cnps_biomarkerpath, "cnps_random_forest_prediction.pdf"), plot = p42.4, width = 10, height = 8)

# 保存随机森林数据到总表
addWorksheet(cnps_biomarker_wb, "random_forest_data")
writeData(cnps_biomarker_wb, "random_forest_data", dat, rowNames = TRUE)

#32 nnet.metf: 神经网络筛选特征功能  ------
library(nnet)
res =nnet.metf(ps=ps.cnps, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果数据到总表
addWorksheet(cnps_biomarker_wb, "neural_network_accuracy")
addWorksheet(cnps_biomarker_wb, "neural_network_importance")
writeData(cnps_biomarker_wb, "neural_network_accuracy", accuracy, rowNames = TRUE)
writeData(cnps_biomarker_wb, "neural_network_importance", importance, rowNames = TRUE)

#33 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps =  ps.cnps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果数据到总表
addWorksheet(cnps_biomarker_wb, "bagging_accuracy")
addWorksheet(cnps_biomarker_wb, "bagging_importance")
writeData(cnps_biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
writeData(cnps_biomarker_wb, "bagging_importance", importance, rowNames = TRUE)

# 保存CNPS生物标记物总表
saveWorkbook(cnps_biomarker_wb, file.path(cnps_biomarkerpath, "cnps_biomarker_results.xlsx"), overwrite = TRUE)

# 创建CNPS网络分析目录
cnps_networkpath <- file.path(cnpspath, "network")
dir.create(cnps_networkpath, recursive = TRUE)

# 创建CNPS网络分析总表
cnps_network_wb <- createWorkbook()

#34 CNPS.network:CNPS基因网络-------
detach("package:mia", unload = TRUE)
library(data.table)
library(tidyfst)
library(sna)

# C循环网络
res_C = CNPS.network2(ps = ps.kegg,id.0 = "C")
p_C = res_C[[1]]
dat_C = res_C[[2]]
dat2_C = res_C[[3]]

# 保存C循环网络图片
ggsave(file.path(cnps_networkpath, "cnps_network_C.png"), plot = p_C, width = 12, height = 10, dpi = 300)
ggsave(file.path(cnps_networkpath, "cnps_network_C.pdf"), plot = p_C, width = 12, height = 10)

# N循环网络
res_N = CNPS.network2(ps = ps.kegg,id.0 = "N")
p_N = res_N[[1]]
dat_N = res_N[[2]]
dat2_N = res_N[[3]]

# 保存N循环网络图片
ggsave(file.path(cnps_networkpath, "cnps_network_N.png"), plot = p_N, width = 12, height = 10, dpi = 300)
ggsave(file.path(cnps_networkpath, "cnps_network_N.pdf"), plot = p_N, width = 12, height = 10)

# P循环网络
res_P = CNPS.network2(ps = ps.kegg,id.0 = "P")
p_P = res_P[[1]]
dat_P = res_P[[2]]
dat2_P = res_P[[3]]

# 保存P循环网络图片
ggsave(file.path(cnps_networkpath, "cnps_network_P.png"), plot = p_P, width = 12, height = 10, dpi = 300)
ggsave(file.path(cnps_networkpath, "cnps_network_P.pdf"), plot = p_P, width = 12, height = 10)

# S循环网络
res_S = CNPS.network2(ps = ps.kegg,id.0 = "S")
p_S = res_S[[1]]
dat_S = res_S[[2]]
dat2_S = res_S[[3]]

# 保存S循环网络图片
ggsave(file.path(cnps_networkpath, "cnps_network_S.png"), plot = p_S, width = 12, height = 10, dpi = 300)
ggsave(file.path(cnps_networkpath, "cnps_network_S.pdf"), plot = p_S, width = 12, height = 10)

# 保存CNPS网络数据到总表
addWorksheet(cnps_network_wb, "C_cycle_data1")
addWorksheet(cnps_network_wb, "C_cycle_data2")
addWorksheet(cnps_network_wb, "N_cycle_data1")
addWorksheet(cnps_network_wb, "N_cycle_data2")
addWorksheet(cnps_network_wb, "P_cycle_data1")
addWorksheet(cnps_network_wb, "P_cycle_data2")
addWorksheet(cnps_network_wb, "S_cycle_data1")
addWorksheet(cnps_network_wb, "S_cycle_data2")
writeData(cnps_network_wb, "C_cycle_data1", dat_C, rowNames = TRUE)
writeData(cnps_network_wb, "C_cycle_data2", dat2_C, rowNames = TRUE)
writeData(cnps_network_wb, "N_cycle_data1", dat_N, rowNames = TRUE)
writeData(cnps_network_wb, "N_cycle_data2", dat2_N, rowNames = TRUE)
writeData(cnps_network_wb, "P_cycle_data1", dat_P, rowNames = TRUE)
writeData(cnps_network_wb, "P_cycle_data2", dat2_P, rowNames = TRUE)
writeData(cnps_network_wb, "S_cycle_data1", dat_S, rowNames = TRUE)
writeData(cnps_network_wb, "S_cycle_data2", dat2_S, rowNames = TRUE)

# 保存CNPS网络分析总表
saveWorkbook(cnps_network_wb, file.path(cnps_networkpath, "cnps_network_results.xlsx"), overwrite = TRUE)

# card数据库------
# 创建CARD数据库目录
cardpath <- file.path(path, "card")
dir.create(cardpath, recursive = TRUE)

data("ps.card")
ps.card =  EasyMultiOmics::ps.card
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
phyloseq::tax_table(ps.card) = phyloseq:: tax_table(as.matrix(tax))
ps.card = ps.card %>% filter_OTU_ps(Top = 1000)

#--提取有多少个分组
Top = 20
gnum = sample_data(ps.card )$Group %>%unique() %>% length()
gnum
# 设定排序顺序--Group

axis_order = sample_data(ps.card )$Group %>%unique()

# 设定排序顺序--ID
map = sample_data(ps.card )
map$ID = row.names(map)
map = map %>%
  as.tibble() %>%
  dplyr::arrange(desc(Group))
axis_order.s = map$ID

# function diversity -----
# 创建CARD alpha多样性目录
card_alphapath <- file.path(cardpath, "alpha")
dir.create(card_alphapath, recursive = TRUE)

# 创建CARD alpha多样性总表
card_alpha_wb <- createWorkbook()

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

# 保存CARD alpha多样性图片
ggsave(file.path(card_alphapath, "card_alpha_diversity_box.png"), plot = p1_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_alphapath, "card_alpha_diversity_box.pdf"), plot = p1_1, width = 12, height = 8)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

ggsave(file.path(card_alphapath, "card_alpha_diversity_bar.png"), plot = p1_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_alphapath, "card_alpha_diversity_bar.pdf"), plot = p1_2, width = 12, height = 8)

res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

ggsave(file.path(card_alphapath, "card_alpha_diversity_boxbar.png"), plot = p1_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_alphapath, "card_alpha_diversity_boxbar.pdf"), plot = p1_3, width = 12, height = 8)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0

ggsave(file.path(card_alphapath, "card_alpha_diversity_violin.png"), plot = p1_0, width = 14, height = 8, dpi = 300)
ggsave(file.path(card_alphapath, "card_alpha_diversity_violin.pdf"), plot = p1_0, width = 14, height = 8)

# 保存CARD alpha多样性数据到总表
addWorksheet(card_alpha_wb, "card_alpha_data")
addWorksheet(card_alpha_wb, "card_alpha_stats")
writeData(card_alpha_wb, "card_alpha_data", data, rowNames = TRUE)
writeData(card_alpha_wb, "card_alpha_stats", result$sta, rowNames = TRUE)

#2 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps.card ))/10
result = alpha_rare.metf(ps = ps.card , group = "Group", method = "Richness", start = 100, step = rare)

#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

ggsave(file.path(card_alphapath, "card_rarefaction_individual.png"), plot = p2_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_alphapath, "card_rarefaction_individual.pdf"), plot = p2_1, width = 10, height = 8)

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2

ggsave(file.path(card_alphapath, "card_rarefaction_group.png"), plot = p2_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_alphapath, "card_rarefaction_group.pdf"), plot = p2_2, width = 10, height = 8)

#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

# 保存稀释曲线数据到alpha总表
addWorksheet(card_alpha_wb, "rarefaction_data")
writeData(card_alpha_wb, "rarefaction_data", raretab, rowNames = TRUE)

# 保存alpha多样性总表
saveWorkbook(card_alpha_wb, file.path(card_alphapath, "card_alpha_diversity.xlsx"), overwrite = TRUE)

# 创建CARD beta多样性目录
card_betapath <- file.path(cardpath, "beta")
dir.create(card_betapath, recursive = TRUE)

# 创建CARD beta多样性总表
card_beta_wb <- createWorkbook()

#3 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps.card ,
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

# 保存排序分析图片
ggsave(file.path(card_betapath, "card_ordination_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_betapath, "card_ordination_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(card_betapath, "card_ordination_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_betapath, "card_ordination_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(card_betapath, "card_ordination_spider.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_betapath, "card_ordination_spider.pdf"), plot = p3_3, width = 10, height = 8)

# 保存排序分析数据
addWorksheet(card_beta_wb, "ordination_data")
writeData(card_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps.card , Micromet = "adonis", dist = "bray")
dat1

# 保存群落差异检测结果
addWorksheet(card_beta_wb, "adonis_test")
writeData(card_beta_wb, "adonis_test", dat1, rowNames = TRUE)

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps.card , Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较结果
addWorksheet(card_beta_wb, "pairwise_test")
writeData(card_beta_wb, "pairwise_test", dat2, rowNames = TRUE)

#6 mantal.metf：群落功能差异检测普鲁士分析#------
result <- mantal.metf(ps = ps.card ,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

# 保存Mantel检验结果
ggsave(file.path(card_betapath, "card_mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_betapath, "card_mantel_test.pdf"), plot = p3_7, width = 10, height = 8)

addWorksheet(card_beta_wb, "mantel_test")
writeData(card_beta_wb, "mantel_test", data, rowNames = TRUE)

#7 cluster.metf:样品聚类#-----
res = cluster.metf(ps= ps.card ,
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

# 保存聚类分析图片
ggsave(file.path(card_betapath, "card_cluster_plot1.png"), plot = p4, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_betapath, "card_cluster_plot1.pdf"), plot = p4, width = 12, height = 8)
ggsave(file.path(card_betapath, "card_cluster_plot2.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_betapath, "card_cluster_plot2.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(card_betapath, "card_cluster_plot3.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_betapath, "card_cluster_plot3.pdf"), plot = p4_2, width = 12, height = 8)

# 保存聚类分析数据
addWorksheet(card_beta_wb, "cluster_data")
writeData(card_beta_wb, "cluster_data", dat, rowNames = TRUE)

# 保存beta多样性总表
saveWorkbook(card_beta_wb, file.path(card_betapath, "card_beta_diversity.xlsx"), overwrite = TRUE)

# 创建CARD组成分析目录
card_compositionpath <- file.path(cardpath, "composition")
dir.create(card_compositionpath, recursive = TRUE)

# 创建CARD组成分析总表
card_composition_wb <- createWorkbook()

#8 Micro_tern.metf: 三元图展示功能----
res = Micro_tern.metf(ps.card  %>% filter_OTU_ps(500), color = "Drug_Class" )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图结果
ggsave(file.path(card_compositionpath, "card_ternary_plot.png"), plot = p15, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_ternary_plot.pdf"), plot = p15, width = 10, height = 8)

addWorksheet(card_composition_wb, "ternary_data")
writeData(card_composition_wb, "ternary_data", dat, rowNames = TRUE)

# function classfication

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps.card,
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

# 保存堆积柱状图结果
ggsave(file.path(card_compositionpath, "card_barplot1.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_barplot1.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(card_compositionpath, "card_barplot2.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_barplot2.pdf"), plot = p4_2, width = 12, height = 8)

addWorksheet(card_composition_wb, "barplot_data")
writeData(card_composition_wb, "barplot_data", databar, rowNames = TRUE)

#10 Ven.Upset.metf: 用于展示共有、特有的功能----
# 分组小于6时使用
res = Ven.Upset.metf(ps =  ps.card,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

# 保存韦恩图结果
ggsave(file.path(card_compositionpath, "card_venn_plot.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_venn_plot.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(card_compositionpath, "card_upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_upset_plot.pdf"), plot = p10.2, width = 12, height = 8)

addWorksheet(card_composition_wb, "venn_upset_data")
writeData(card_composition_wb, "venn_upset_data", dat, rowNames = TRUE)

#11 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps.card)$Group %>% table()
num =10
result = VenSuper.metf(ps = ps.card,
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

# 保存韦恩详细分析结果
ggsave(file.path(card_compositionpath, "card_venn_detail1.png"), plot = p7_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_venn_detail1.pdf"), plot = p7_1, width = 12, height = 8)
ggsave(file.path(card_compositionpath, "card_venn_detail2.png"), plot = p7_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_venn_detail2.pdf"), plot = p7_2, width = 12, height = 8)

addWorksheet(card_composition_wb, "venn_detail_data")
writeData(card_composition_wb, "venn_detail_data", dat, rowNames = TRUE)

#13 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps.card ,
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

# 保存花瓣图结果
ggsave(file.path(card_compositionpath, "card_flower_plot.png"), plot = p13.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_flower_plot.pdf"), plot = p13.1, width = 10, height = 8)

addWorksheet(card_composition_wb, "flower_data")
writeData(card_composition_wb, "flower_data", dat, rowNames = TRUE)

#14 ven.network.metf: ven网络展示共有特有功能----
result =ven.network.metf(
  ps = ps.card,
  N = 0.5,
  fill = "Drug_Class")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat

# 保存韦恩网络结果
ggsave(file.path(card_compositionpath, "card_venn_network.png"), plot = p14, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_compositionpath, "card_venn_network.pdf"), plot = p14, width = 10, height = 8)

addWorksheet(card_composition_wb, "venn_network_data")
writeData(card_composition_wb, "venn_network_data", dat, rowNames = TRUE)

# 保存组成分析总表
saveWorkbook(card_composition_wb, file.path(card_compositionpath, "card_composition_analysis.xlsx"), overwrite = TRUE)

# function differential analysis----
# 创建CARD差异分析目录
card_diffpath <- file.path(cardpath, "differential")
dir.create(card_diffpath, recursive = TRUE)

# 创建CARD差异分析总表
card_diff_wb <- createWorkbook()

# 数据转换为"Drug_Class" 或其他
ps.g = ps.card %>%
  tax_glom_meta(ranks = "Drug_Class")

#15 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps.g,
                        # group  = "Drug_Class",
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

# 保存DESep2结果
ggsave(file.path(card_diffpath, "card_deseq2_plot1.png"), plot = p15.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_deseq2_plot1.pdf"), plot = p15.1, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_deseq2_plot2.png"), plot = p15.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_deseq2_plot2.pdf"), plot = p15.2, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_deseq2_plot3.png"), plot = p15.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_deseq2_plot3.pdf"), plot = p15.3, width = 10, height = 8)

addWorksheet(card_diff_wb, "DESep2_results")
writeData(card_diff_wb, "DESep2_results", dat, rowNames = TRUE)

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

# 保存EdgeR结果
ggsave(file.path(card_diffpath, "card_edger_plot1.png"), plot = p16.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_edger_plot1.pdf"), plot = p16.1, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_edger_plot2.png"), plot = p16.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_edger_plot2.pdf"), plot = p16.2, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_edger_plot3.png"), plot = p16.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_edger_plot3.pdf"), plot = p16.3, width = 10, height = 8)

addWorksheet(card_diff_wb, "EdgeR_results")
writeData(card_diff_wb, "EdgeR_results", dat, rowNames = TRUE)

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

# 保存t检验结果
ggsave(file.path(card_diffpath, "card_ttest_plot1.png"), plot = p17.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_ttest_plot1.pdf"), plot = p17.1, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_ttest_plot2.png"), plot = p17.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_ttest_plot2.pdf"), plot = p17.2, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_ttest_plot3.png"), plot = p17.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_ttest_plot3.pdf"), plot = p17.3, width = 10, height = 8)

addWorksheet(card_diff_wb, "ttest_results")
writeData(card_diff_wb, "ttest_results", dat, rowNames = TRUE)

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

# 保存非参数检验结果
ggsave(file.path(card_diffpath, "card_wilcox_plot1.png"), plot = p18.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_wilcox_plot1.pdf"), plot = p18.1, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_wilcox_plot2.png"), plot = p18.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_wilcox_plot2.pdf"), plot = p18.2, width = 10, height = 8)
ggsave(file.path(card_diffpath, "card_wilcox_plot3.png"), plot = p18.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_wilcox_plot3.pdf"), plot = p18.3, width = 10, height = 8)

addWorksheet(card_diff_wb, "wilcox_results")
writeData(card_diff_wb, "wilcox_results", dat, rowNames = TRUE)

#19 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps.card,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2

# 保存STAMP结果
ggsave(file.path(card_diffpath, "card_stamp_plot.png"), plot = p19, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_diffpath, "card_stamp_plot.pdf"), plot = p19, width = 12, height = 8)

addWorksheet(card_diff_wb, "STAMP_results1")
addWorksheet(card_diff_wb, "STAMP_results2")
writeData(card_diff_wb, "STAMP_results1", dat1, rowNames = TRUE)
writeData(card_diff_wb, "STAMP_results2", dat2, rowNames = TRUE)

# 保存差异分析总表
saveWorkbook(card_diff_wb, file.path(card_diffpath, "card_differential_analysis.xlsx"), overwrite = TRUE)

# biomarker identification-----
# 创建CARD生物标志物目录
card_biomarkerpath <- file.path(cardpath, "biomarker")
dir.create(card_biomarkerpath, recursive = TRUE)

# 创建CARD生物标志物总表
card_biomarker_wb <- createWorkbook()

#20 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps.card %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存RFCV结果
ggsave(file.path(card_biomarkerpath, "card_rfcv_plot.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_rfcv_plot.pdf"), plot = prfcv, width = 10, height = 8)

addWorksheet(card_biomarker_wb, "RFCV_results")
writeData(card_biomarker_wb, "RFCV_results", rfcvtable, rowNames = TRUE)

#21 Roc.metf:ROC 曲线绘制----
id = sample_data(ps.card)$Group %>% unique()
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

# 保存ROC结果
ggsave(file.path(card_biomarkerpath, "card_roc_plot1.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_roc_plot1.pdf"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(card_biomarkerpath, "card_roc_plot2.png"), plot = p33.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_roc_plot2.pdf"), plot = p33.2, width = 10, height = 8)

addWorksheet(card_biomarker_wb, "ROC_results")
writeData(card_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#22 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.card,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(card_biomarkerpath, "card_pca_loading.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_pca_loading.pdf"), plot = p34.1, width = 10, height = 8)

addWorksheet(card_biomarker_wb, "PCA_loading")
writeData(card_biomarker_wb, "PCA_loading", dat, rowNames = TRUE)

#23 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = ps.card %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(card_biomarker_wb, "SVM_AUC")
addWorksheet(card_biomarker_wb, "SVM_importance")
writeData(card_biomarker_wb, "SVM_AUC", AUC, rowNames = TRUE)
writeData(card_biomarker_wb, "SVM_importance", importance, rowNames = TRUE)

#24 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = ps.card %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(card_biomarker_wb, "GLM_AUC")
addWorksheet(card_biomarker_wb, "GLM_importance")
writeData(card_biomarker_wb, "GLM_AUC", AUC, rowNames = TRUE)
writeData(card_biomarker_wb, "GLM_importance", importance, rowNames = TRUE)

#25 xgboost.metf: xgboost筛选特征功能----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =ps.card, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(card_biomarker_wb, "XGBoost_accuracy")
addWorksheet(card_biomarker_wb, "XGBoost_importance")
writeData(card_biomarker_wb, "XGBoost_accuracy", accuracy, rowNames = TRUE)
writeData(card_biomarker_wb, "XGBoost_importance", importance, rowNames = TRUE)

#26 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  ps.card, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存LASSO结果
addWorksheet(card_biomarker_wb, "LASSO_accuracy")
addWorksheet(card_biomarker_wb, "LASSO_importance")
writeData(card_biomarker_wb, "LASSO_accuracy", accuracy, rowNames = TRUE)
writeData(card_biomarker_wb, "LASSO_importance", importance, rowNames = TRUE)

#27 decisiontree.metf:----
library(rpart)
res =decisiontree.metf(ps=ps.card, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(card_biomarker_wb, "DecisionTree_accuracy")
addWorksheet(card_biomarker_wb, "DecisionTree_importance")
writeData(card_biomarker_wb, "DecisionTree_accuracy", accuracy, rowNames = TRUE)
writeData(card_biomarker_wb, "DecisionTree_importance", importance, rowNames = TRUE)

#28 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps.card, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(card_biomarker_wb, "NaiveBayes_accuracy")
addWorksheet(card_biomarker_wb, "NaiveBayes_importance")
writeData(card_biomarker_wb, "NaiveBayes_accuracy", accuracy, rowNames = TRUE)
writeData(card_biomarker_wb, "NaiveBayes_importance", importance, rowNames = TRUE)

#29 LDA.metf: LDA筛选特征功能-----
tablda = LDA.metf(ps = ps.card,
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

# 保存LDA结果
ggsave(file.path(card_biomarkerpath, "card_lda_plot.png"), plot = p35, width = 12, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_lda_plot.pdf"), plot = p35, width = 12, height = 8)

addWorksheet(card_biomarker_wb, "LDA_results")
writeData(card_biomarker_wb, "LDA_results", dat, rowNames = TRUE)

#30 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps.card,group  = "Group", optimal = 50 ,fill ="Drug_Class" )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林结果
ggsave(file.path(card_biomarkerpath, "card_rf_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_rf_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(card_biomarkerpath, "card_rf_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_rf_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(card_biomarkerpath, "card_rf_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(card_biomarkerpath, "card_rf_plot4.pdf"), plot = p42.4, width = 10, height = 8)

addWorksheet(card_biomarker_wb, "RandomForest_results")
writeData(card_biomarker_wb, "RandomForest_results", dat, rowNames = TRUE)

#31 nnet.metf: 神经网络筛选特征功能  ------
res =nnet.metf(ps=ps.card, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(card_biomarker_wb, "NeuralNet_accuracy")
addWorksheet(card_biomarker_wb, "NeuralNet_importance")
writeData(card_biomarker_wb, "NeuralNet_accuracy", accuracy, rowNames = TRUE)
writeData(card_biomarker_wb, "NeuralNet_importance", importance, rowNames = TRUE)

#32 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps =  ps.card, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(card_biomarker_wb, "Bagging_accuracy")
addWorksheet(card_biomarker_wb, "Bagging_importance")
writeData(card_biomarker_wb, "Bagging_accuracy", accuracy, rowNames = TRUE)
writeData(card_biomarker_wb, "Bagging_importance", importance, rowNames = TRUE)

# 保存生物标志物总表
saveWorkbook(card_biomarker_wb, file.path(card_biomarkerpath, "card_biomarker_analysis.xlsx"), overwrite = TRUE)

#network analysis -----
# 创建CARD网络分析目录
card_networkpath <- file.path(cardpath, "network")
dir.create(card_networkpath, recursive = TRUE)

# 创建CARD网络分析总表
card_network_wb <- createWorkbook()

#33 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
phyloseq::tax_table(ps.card) %>% colnames()
tab.r = network.pip(
  ps = ps.card,
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
  lab = "AMR_Gene_Family",
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

# 保存网络分析主图
ggsave(file.path(card_networkpath, "card_network_main.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(card_networkpath, "card_network_main.pdf"), plot = p0, width = 12, height = 10)

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

# 保存网络属性数据
addWorksheet(card_network_wb, "network_properties")
writeData(card_network_wb, "network_properties", dat2, rowNames = TRUE)

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

# 保存样本网络属性数据
addWorksheet(card_network_wb, "sample_network_properties")
writeData(card_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

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

# 保存节点属性数据
addWorksheet(card_network_wb, "node_properties")
writeData(card_network_wb, "node_properties", nodepro2, rowNames = TRUE)

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

# 保存网络比较结果
addWorksheet(card_network_wb, "network_comparison")
writeData(card_network_wb, "network_comparison", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(card_network_wb, file.path(card_networkpath, "card_network_analysis.xlsx"), overwrite = TRUE)

# cazy数据库------
# 创建CAZY数据库目录
cazypath <- file.path(path, "cazy")
dir.create(cazypath, recursive = TRUE)

data("ps.cazy")
ps.cazy =  EasyMultiOmics::ps.cazy
# 数据过滤-----
#  通过比例和长度过滤--常规操作-文章常用
tax = ps.cazy %>% vegan_tax() %>% as.data.frame() %>% rownames_to_column("ID")
head(tax)

colnames(tax)[is.na(colnames(tax))] = "others"
tax$Identity = as.numeric(tax$Identity)
tax$Align_len = as.numeric(tax$Align_len)

tax = tax %>%
  dplyr::filter(Identity > 0.6, Align_len > 25)
dim(tax)
tax = tax %>% column_to_rownames("ID")

head(tax)
phyloseq::tax_table(ps.cazy) = phyloseq::tax_table(as.matrix(tax))
ps.cazy = ps.cazy %>% filter_OTU_ps(Top = 1000)

#--提取有多少个分组
Top = 20
gnum = sample_data(ps.cazy)$Group %>% unique() %>% length()
gnum
# 设定排序顺序--Group
axis_order = sample_data(ps.cazy)$Group %>% unique()

# 设定排序顺序--ID
map = sample_data(ps.cazy)
map$ID = row.names(map)
map = map %>%
  as.tibble() %>%
  dplyr::arrange(desc(Group))
axis_order.s = map$ID

# function diversity -----
# 创建CAZY alpha多样性目录
cazy_alphapath <- file.path(cazypath, "alpha")
dir.create(cazy_alphapath, recursive = TRUE)

# 创建CAZY alpha多样性总表
cazy_alpha_wb <- createWorkbook()

#1 alpha.metf: 6种alpha多样性计算#----
all.alpha = c("Shannon","Inv_Simpson","Pielou_evenness","Simpson_evenness" ,"Richness" ,"Chao1","ACE")
#--alpha多样性指标运算
tab = alpha.metf(ps = ps.cazy, group = "Group", Plot = TRUE)
head(tab)
data = cbind(data.frame(ID = 1:length(tab$Group), group = tab$Group), tab[all.alpha])
head(data)
data$ID = as.character(data$ID)

result = EasyStat::MuiKwWlx2(data = data, num = 3:5)
result1 = EasyStat::FacetMuiPlotresultBox(data = data, num = 3:5,
                                          result = result,
                                          sig_show = "abc", ncol = 4)
p1_1 = result1[[1]]
p1_1 +
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

# 保存CAZY alpha多样性图片
ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_box.png"), plot = p1_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_box.pdf"), plot = p1_1, width = 12, height = 8)

res = EasyStat::FacetMuiPlotresultBar(data = data, num = c(3:(5)), result = result, sig_show = "abc", ncol = 4)
p1_2 = res[[1]]
p1_2 +
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_bar.png"), plot = p1_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_bar.pdf"), plot = p1_2, width = 12, height = 8)

res = EasyStat::FacetMuiPlotReBoxBar(data = data, num = c(3:(5)), result = result, sig_show = "abc", ncol = 3)
p1_3 = res[[1]]
p1_3 +
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_boxbar.png"), plot = p1_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_boxbar.pdf"), plot = p1_3, width = 12, height = 8)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  geom_text(aes(x=group , y=y ,label=stat))
p1_0

ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_violin.png"), plot = p1_0, width = 14, height = 8, dpi = 300)
ggsave(file.path(cazy_alphapath, "cazy_alpha_diversity_violin.pdf"), plot = p1_0, width = 14, height = 8)

# 保存CAZY alpha多样性数据到总表
addWorksheet(cazy_alpha_wb, "cazy_alpha_data")
addWorksheet(cazy_alpha_wb, "cazy_alpha_stats")
writeData(cazy_alpha_wb, "cazy_alpha_data", data, rowNames = TRUE)
writeData(cazy_alpha_wb, "cazy_alpha_stats", result$sta, rowNames = TRUE)

#2 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps.cazy))/10
result = alpha_rare.metf(ps = ps.cazy, group = "Group", method = "Richness", start = 100, step = rare)

#--提供单个样本稀释曲线的绘制
p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_1

ggsave(file.path(cazy_alphapath, "cazy_rarefaction_individual.png"), plot = p2_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_alphapath, "cazy_rarefaction_individual.pdf"), plot = p2_1, width = 10, height = 8)

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2

ggsave(file.path(cazy_alphapath, "cazy_rarefaction_group.png"), plot = p2_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_alphapath, "cazy_rarefaction_group.pdf"), plot = p2_2, width = 10, height = 8)

#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

# 保存稀释曲线数据到alpha总表
addWorksheet(cazy_alpha_wb, "rarefaction_data")
writeData(cazy_alpha_wb, "rarefaction_data", raretab, rowNames = TRUE)

# 保存alpha多样性总表
saveWorkbook(cazy_alpha_wb, file.path(cazy_alphapath, "cazy_alpha_diversity.xlsx"), overwrite = TRUE)

# 创建CAZY beta多样性目录
cazy_betapath <- file.path(cazypath, "beta")
dir.create(cazy_betapath, recursive = TRUE)

# 创建CAZY beta多样性总表
cazy_beta_wb <- createWorkbook()

#3 ordinate.metf: 排序分析#----------
result = ordinate.metf(ps = ps.cazy,
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
plotdata = result[[2]]
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

# 保存排序分析图片
ggsave(file.path(cazy_betapath, "cazy_ordination_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_betapath, "cazy_ordination_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(cazy_betapath, "cazy_ordination_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_betapath, "cazy_ordination_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(cazy_betapath, "cazy_ordination_spider.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_betapath, "cazy_ordination_spider.pdf"), plot = p3_3, width = 10, height = 8)

# 保存排序分析数据
addWorksheet(cazy_beta_wb, "ordination_data")
writeData(cazy_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps.cazy, Micromet = "adonis", dist = "bray")
dat1

# 保存群落差异检测结果
addWorksheet(cazy_beta_wb, "adonis_test")
writeData(cazy_beta_wb, "adonis_test", dat1, rowNames = TRUE)

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps.cazy, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较结果
addWorksheet(cazy_beta_wb, "pairwise_test")
writeData(cazy_beta_wb, "pairwise_test", dat2, rowNames = TRUE)

#6 mantal.metf：群落功能差异检测普鲁士分析#------
result <- mantal.metf(ps = ps.cazy,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

# 保存Mantel检验结果
ggsave(file.path(cazy_betapath, "cazy_mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_betapath, "cazy_mantel_test.pdf"), plot = p3_7, width = 10, height = 8)

addWorksheet(cazy_beta_wb, "mantel_test")
writeData(cazy_beta_wb, "mantel_test", data, rowNames = TRUE)

#7 cluster.metf:样品聚类#-----
res = cluster.metf(ps= ps.cazy,
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

# 保存聚类分析图片
ggsave(file.path(cazy_betapath, "cazy_cluster_plot1.png"), plot = p4, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_betapath, "cazy_cluster_plot1.pdf"), plot = p4, width = 12, height = 8)
ggsave(file.path(cazy_betapath, "cazy_cluster_plot2.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_betapath, "cazy_cluster_plot2.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(cazy_betapath, "cazy_cluster_plot3.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_betapath, "cazy_cluster_plot3.pdf"), plot = p4_2, width = 12, height = 8)

# 保存聚类分析数据
addWorksheet(cazy_beta_wb, "cluster_data")
writeData(cazy_beta_wb, "cluster_data", dat, rowNames = TRUE)

# 保存beta多样性总表
saveWorkbook(cazy_beta_wb, file.path(cazy_betapath, "cazy_beta_diversity.xlsx"), overwrite = TRUE)

# 创建CAZY组成分析目录
cazy_compositionpath <- file.path(cazypath, "composition")
dir.create(cazy_compositionpath, recursive = TRUE)

# 创建CAZY组成分析总表
cazy_composition_wb <- createWorkbook()

#8 Micro_tern.metf: 三元图展示功能----
res = Micro_tern.metf(ps.cazy %>% filter_OTU_ps(500), color = "Class")
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图结果
ggsave(file.path(cazy_compositionpath, "cazy_ternary_plot.png"), plot = p15, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_ternary_plot.pdf"), plot = p15, width = 10, height = 8)

addWorksheet(cazy_composition_wb, "ternary_data")
writeData(cazy_composition_wb, "ternary_data", dat, rowNames = TRUE)

# function classfication

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
result = barMainplot.metf(ps = ps.cazy,
                          j = "Class",
                          label = FALSE,
                          sd = FALSE,
                          Top = 10)
p4_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p4_1+mytheme1

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2+mytheme1

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Class","Abundance(%)")
head(databar)

# 保存堆积柱状图结果
ggsave(file.path(cazy_compositionpath, "cazy_barplot1.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_barplot1.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(cazy_compositionpath, "cazy_barplot2.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_barplot2.pdf"), plot = p4_2, width = 12, height = 8)

addWorksheet(cazy_composition_wb, "barplot_data")
writeData(cazy_composition_wb, "barplot_data", databar, rowNames = TRUE)

#10 Ven.Upset.metf: 用于展示共有、特有的功能----
# 分组小于6时使用
res = Ven.Upset.metf(ps =  ps.cazy,
                     group = "Group",
                     N = 0.5,
                     size = 3)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat

# 保存韦恩图结果
ggsave(file.path(cazy_compositionpath, "cazy_venn_plot.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_venn_plot.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(cazy_compositionpath, "cazy_upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_upset_plot.pdf"), plot = p10.2, width = 12, height = 8)

addWorksheet(cazy_composition_wb, "venn_upset_data")
writeData(cazy_composition_wb, "venn_upset_data", dat, rowNames = TRUE)

#11 VenSeper.metf:  详细展示每一组中物种功能 小问题  分类填充-------
#---每个部分
sample_data(ps.cazy)$Group %>% table()
num = 10
result = VenSuper.metf(ps = ps.cazy,
                       group = "Group",
                       num = 10)

# 提取韦恩图中全部部分的otu极其丰度做门类柱状图
p7_1 <- result[[1]]+scale_fill_brewer(palette = "Paired")
p7_1
#每个部分序列的数量占比，并作差异
dat<- result[[2]]
# 每部分的otu门类冲积图
p7_2 <- result[[3]]+scale_fill_brewer(palette = "Paired")
p7_2

# 保存韦恩详细分析结果
ggsave(file.path(cazy_compositionpath, "cazy_venn_detail1.png"), plot = p7_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_venn_detail1.pdf"), plot = p7_1, width = 12, height = 8)
ggsave(file.path(cazy_compositionpath, "cazy_venn_detail2.png"), plot = p7_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_venn_detail2.pdf"), plot = p7_2, width = 12, height = 8)

addWorksheet(cazy_composition_wb, "venn_detail_data")
writeData(cazy_composition_wb, "venn_detail_data", dat, rowNames = TRUE)

#12 ggflower.metf:花瓣图展示共有特有功能------
res <- ggflower.metf (ps = ps.cazy,
                      group = "ID",
                      start = 1, # 风车效果
                      m1 = 1, # 花瓣形状，方形到圆形到棱形，数值逐渐减少。
                      a = 0.2, # 花瓣胖瘦
                      b = 1, # 花瓣距离花心的距离
                      lab.leaf = 1, # 花瓣标签到圆心的距离
                      col.cir = "yellow",
                      N = 0.5)

p13.1 = res[[1]]

p13.1

dat = res[[2]]
dat

# 保存花瓣图结果
ggsave(file.path(cazy_compositionpath, "cazy_flower_plot.png"), plot = p13.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_flower_plot.pdf"), plot = p13.1, width = 10, height = 8)

addWorksheet(cazy_composition_wb, "flower_data")
writeData(cazy_composition_wb, "flower_data", dat, rowNames = TRUE)

#13 ven.network.metf: ven网络展示共有特有功能----
result = ven.network.metf(
  ps = ps.cazy,
  N = 0.5,
  fill = "Class")

p14  = result[[1]] +
  theme(legend.position = "none")
p14
dat = result[[2]]
dat

# 保存韦恩网络结果
ggsave(file.path(cazy_compositionpath, "cazy_venn_network.png"), plot = p14, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_compositionpath, "cazy_venn_network.pdf"), plot = p14, width = 10, height = 8)

addWorksheet(cazy_composition_wb, "venn_network_data")
writeData(cazy_composition_wb, "venn_network_data", dat, rowNames = TRUE)

# 保存组成分析总表
saveWorkbook(cazy_composition_wb, file.path(cazy_compositionpath, "cazy_composition_analysis.xlsx"), overwrite = TRUE)

# function differential analysis----
# 创建CAZY差异分析目录
cazy_diffpath <- file.path(cazypath, "differential")
dir.create(cazy_diffpath, recursive = TRUE)

# 创建CAZY差异分析总表
cazy_diff_wb <- createWorkbook()

# 数据转换为"Class" 或其他
ps.g = ps.cazy %>%
  tax_glom_meta(ranks = "Class")

#14 DESep2Super.metf:DESep2计算差异功能基因 ----
res = DESep2Super.metf (ps = ps.g,
                        artGroup = NULL)
p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][1]
p15.2
p15.3 = res[[1]][1]
p15.3
dat = res[[2]]
dat

# 保存DESep2结果
ggsave(file.path(cazy_diffpath, "cazy_deseq2_plot1.png"), plot = p15.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_deseq2_plot1.pdf"), plot = p15.1, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_deseq2_plot2.png"), plot = p15.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_deseq2_plot2.pdf"), plot = p15.2, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_deseq2_plot3.png"), plot = p15.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_deseq2_plot3.pdf"), plot = p15.3, width = 10, height = 8)

addWorksheet(cazy_diff_wb, "DESep2_results")
writeData(cazy_diff_wb, "DESep2_results", dat, rowNames = TRUE)

#15 EdgerSuper.metf:EdgeR计算差异功能基因----
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

# 保存EdgeR结果
ggsave(file.path(cazy_diffpath, "cazy_edger_plot1.png"), plot = p16.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_edger_plot1.pdf"), plot = p16.1, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_edger_plot2.png"), plot = p16.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_edger_plot2.pdf"), plot = p16.2, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_edger_plot3.png"), plot = p16.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_edger_plot3.pdf"), plot = p16.3, width = 10, height = 8)

addWorksheet(cazy_diff_wb, "EdgeR_results")
writeData(cazy_diff_wb, "EdgeR_results", dat, rowNames = TRUE)

#16 t.metf: 差异分析t检验----
res = t.metf(ps = ps.g, group = "Group", artGroup = NULL)
p17.1 = res [[1]][1]
p17.1
p17.2 = res [[1]][2]
p17.2
p17.3 = res [[1]][3]
p17.3
dat =  res [[2]]
dat

# 保存t检验结果
ggsave(file.path(cazy_diffpath, "cazy_ttest_plot1.png"), plot = p17.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_ttest_plot1.pdf"), plot = p17.1, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_ttest_plot2.png"), plot = p17.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_ttest_plot2.pdf"), plot = p17.2, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_ttest_plot3.png"), plot = p17.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_ttest_plot3.pdf"), plot = p17.3, width = 10, height = 8)

addWorksheet(cazy_diff_wb, "ttest_results")
writeData(cazy_diff_wb, "ttest_results", dat, rowNames = TRUE)

#17 wlx.metf:非参数检验——-----
res= wlx.metf(ps = ps.g, group = "Group", artGroup = NULL)
p18.1 = res [[1]][1]
p18.1$`KO-OE_plot`
p18.2 = res [[1]][2]
p18.2
p18.3 = res [[1]][3]
p18.3
dat =  res [[2]]
dat %>% head()

# 保存非参数检验结果
ggsave(file.path(cazy_diffpath, "cazy_wilcox_plot1.png"), plot = p18.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_wilcox_plot1.pdf"), plot = p18.1, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_wilcox_plot2.png"), plot = p18.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_wilcox_plot2.pdf"), plot = p18.2, width = 10, height = 8)
ggsave(file.path(cazy_diffpath, "cazy_wilcox_plot3.png"), plot = p18.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_wilcox_plot3.pdf"), plot = p18.3, width = 10, height = 8)

addWorksheet(cazy_diff_wb, "wilcox_results")
writeData(cazy_diff_wb, "wilcox_results", dat, rowNames = TRUE)

#18 stamp.metf:stamp差异分析#------
allgroup <- combn(unique(sample_data(ps)$Group),2)
ps_sub <- subset_samples(ps.cazy,Group %in% allgroup[,1]);ps_sub
res <- stamp.metf(ps = ps_sub,Top = 20)
p19 =res[[1]]
p19
dat1= res[[1]]
dat1
dat2= res[[2]]
dat2

# 保存STAMP结果
ggsave(file.path(cazy_diffpath, "cazy_stamp_plot.png"), plot = p19, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_diffpath, "cazy_stamp_plot.pdf"), plot = p19, width = 12, height = 8)

addWorksheet(cazy_diff_wb, "STAMP_results1")
addWorksheet(cazy_diff_wb, "STAMP_results2")
writeData(cazy_diff_wb, "STAMP_results1", dat1, rowNames = TRUE)
writeData(cazy_diff_wb, "STAMP_results2", dat2, rowNames = TRUE)

# 保存差异分析总表
saveWorkbook(cazy_diff_wb, file.path(cazy_diffpath, "cazy_differential_analysis.xlsx"), overwrite = TRUE)

# biomarker identification-----
# 创建CAZY生物标志物目录
cazy_biomarkerpath <- file.path(cazypath, "biomarker")
dir.create(cazy_biomarkerpath, recursive = TRUE)

# 创建CAZY生物标志物总表
cazy_biomarker_wb <- createWorkbook()

#19 rfcv.metf :交叉验证结果-------
library(randomForest)
library(caret)
library(ROCR) ##用于计算ROC
library(e1071)

result =rfcv.metf(ps = ps.cazy %>% filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存RFCV结果
ggsave(file.path(cazy_biomarkerpath, "cazy_rfcv_plot.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_rfcv_plot.pdf"), plot = prfcv, width = 10, height = 8)

addWorksheet(cazy_biomarker_wb, "RFCV_results")
writeData(cazy_biomarker_wb, "RFCV_results", rfcvtable, rowNames = TRUE)

#20 Roc.metf:ROC 曲线绘制----
id = sample_data(ps.cazy)$Group %>% unique()
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

# 保存ROC结果
ggsave(file.path(cazy_biomarkerpath, "cazy_roc_plot1.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_roc_plot1.pdf"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(cazy_biomarkerpath, "cazy_roc_plot2.png"), plot = p33.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_roc_plot2.pdf"), plot = p33.2, width = 10, height = 8)

addWorksheet(cazy_biomarker_wb, "ROC_results")
writeData(cazy_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.cazy,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(cazy_biomarkerpath, "cazy_pca_loading.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_pca_loading.pdf"), plot = p34.1, width = 10, height = 8)

addWorksheet(cazy_biomarker_wb, "PCA_loading")
writeData(cazy_biomarker_wb, "PCA_loading", dat, rowNames = TRUE)

#22 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = ps.cazy %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(cazy_biomarker_wb, "SVM_AUC")
addWorksheet(cazy_biomarker_wb, "SVM_importance")
writeData(cazy_biomarker_wb, "SVM_AUC", AUC, rowNames = TRUE)
writeData(cazy_biomarker_wb, "SVM_importance", importance, rowNames = TRUE)

#23 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = ps.cazy %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(cazy_biomarker_wb, "GLM_AUC")
addWorksheet(cazy_biomarker_wb, "GLM_importance")
writeData(cazy_biomarker_wb, "GLM_AUC", AUC, rowNames = TRUE)
writeData(cazy_biomarker_wb, "GLM_importance", importance, rowNames = TRUE)

#24 xgboost.metf: xgboost筛选特征功能----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps = ps.cazy, top = 20)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(cazy_biomarker_wb, "XGBoost_accuracy")
addWorksheet(cazy_biomarker_wb, "XGBoost_importance")
writeData(cazy_biomarker_wb, "XGBoost_accuracy", accuracy, rowNames = TRUE)
writeData(cazy_biomarker_wb, "XGBoost_importance", importance, rowNames = TRUE)

#25 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps = ps.cazy, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存LASSO结果
addWorksheet(cazy_biomarker_wb, "LASSO_accuracy")
addWorksheet(cazy_biomarker_wb, "LASSO_importance")
writeData(cazy_biomarker_wb, "LASSO_accuracy", accuracy, rowNames = TRUE)
writeData(cazy_biomarker_wb, "LASSO_importance", importance, rowNames = TRUE)

#26 decisiontree.metf:----
library(rpart)
res =decisiontree.metf(ps = ps.cazy, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(cazy_biomarker_wb, "DecisionTree_accuracy")
addWorksheet(cazy_biomarker_wb, "DecisionTree_importance")
writeData(cazy_biomarker_wb, "DecisionTree_accuracy", accuracy, rowNames = TRUE)
writeData(cazy_biomarker_wb, "DecisionTree_importance", importance, rowNames = TRUE)

#27 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps = ps.cazy, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(cazy_biomarker_wb, "NaiveBayes_accuracy")
addWorksheet(cazy_biomarker_wb, "NaiveBayes_importance")
writeData(cazy_biomarker_wb, "NaiveBayes_accuracy", accuracy, rowNames = TRUE)
writeData(cazy_biomarker_wb, "NaiveBayes_importance", importance, rowNames = TRUE)

#28 LDA.metf: LDA筛选特征功能-----
tablda = LDA.metf(ps = ps.cazy,
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

# 保存LDA结果
ggsave(file.path(cazy_biomarkerpath, "cazy_lda_plot.png"), plot = p35, width = 12, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_lda_plot.pdf"), plot = p35, width = 12, height = 8)

addWorksheet(cazy_biomarker_wb, "LDA_results")
writeData(cazy_biomarker_wb, "LDA_results", dat, rowNames = TRUE)

#29 randomforest.metf: 随机森林筛选特征功能----
res = randomforest.metf( ps = ps.cazy,group  = "Group", optimal = 50 ,fill ="Class" )
p42.1 = res[[1]]
p42.1
p42.2 =res[[2]]
p42.2
dat =res[[3]]
dat
p42.4 =res[[4]]
p42.4

# 保存随机森林结果
ggsave(file.path(cazy_biomarkerpath, "cazy_rf_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_rf_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(cazy_biomarkerpath, "cazy_rf_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_rf_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(cazy_biomarkerpath, "cazy_rf_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(cazy_biomarkerpath, "cazy_rf_plot4.pdf"), plot = p42.4, width = 10, height = 8)

addWorksheet(cazy_biomarker_wb, "RandomForest_results")
writeData(cazy_biomarker_wb, "RandomForest_results", dat, rowNames = TRUE)

#30 nnet.metf: 神经网络筛选特征功能  ------
res =nnet.metf(ps = ps.cazy, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(cazy_biomarker_wb, "NeuralNet_accuracy")
addWorksheet(cazy_biomarker_wb, "NeuralNet_importance")
writeData(cazy_biomarker_wb, "NeuralNet_accuracy", accuracy, rowNames = TRUE)
writeData(cazy_biomarker_wb, "NeuralNet_importance", importance, rowNames = TRUE)

#31 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps = ps.cazy, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(cazy_biomarker_wb, "Bagging_accuracy")
addWorksheet(cazy_biomarker_wb, "Bagging_importance")
writeData(cazy_biomarker_wb, "Bagging_accuracy", accuracy, rowNames = TRUE)
writeData(cazy_biomarker_wb, "Bagging_importance", importance, rowNames = TRUE)

# 保存生物标志物总表
saveWorkbook(cazy_biomarker_wb, file.path(cazy_biomarkerpath, "cazy_biomarker_analysis.xlsx"), overwrite = TRUE)

#network analysis -----
# 创建CAZY网络分析目录
cazy_networkpath <- file.path(cazypath, "network")
dir.create(cazy_networkpath, recursive = TRUE)

# 创建CAZY网络分析总表
cazy_network_wb <- createWorkbook()

#32 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
phyloseq::tax_table(ps.cazy) %>% colnames()
tab.r = network.pip(
  ps = ps.cazy,
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
  lab = "Family",
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

# 保存网络分析主图
ggsave(file.path(cazy_networkpath, "cazy_network_main.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(cazy_networkpath, "cazy_network_main.pdf"), plot = p0, width = 12, height = 10)

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

# 保存网络属性数据
addWorksheet(cazy_network_wb, "network_properties")
writeData(cazy_network_wb, "network_properties", dat2, rowNames = TRUE)

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

# 保存样本网络属性数据
addWorksheet(cazy_network_wb, "sample_network_properties")
writeData(cazy_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

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

# 保存节点属性数据
addWorksheet(cazy_network_wb, "node_properties")
writeData(cazy_network_wb, "node_properties", nodepro2, rowNames = TRUE)

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

# 保存网络比较结果
addWorksheet(cazy_network_wb, "network_comparison")
writeData(cazy_network_wb, "network_comparison", res, rowNames = TRUE)
# 保存网络分析总表
saveWorkbook(cazy_network_wb, file.path(cazy_networkpath, "cazy_network_analysis.xlsx"), overwrite = TRUE)

# cog数据库-----
ps =EasyMultiOmics:: ps.cog %>% filter_OTU_ps(Top = 1000)

# Function
# function diversity -----
# 创建COG多样性分析目录
cog_diversitypath <- file.path(metapath, "cog", "diversity")
dir.create(cog_diversitypath, recursive = TRUE)

# 创建COG多样性分析总表
cog_diversity_wb <- createWorkbook()

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

# 保存alpha多样性结果
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_box.png"), plot = p1_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_box.pdf"), plot = p1_1, width = 12, height = 8)
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_bar.png"), plot = p1_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_bar.pdf"), plot = p1_2, width = 12, height = 8)
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_boxbar.png"), plot = p1_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_boxbar.pdf"), plot = p1_3, width = 12, height = 8)
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_violin.png"), plot = p1_0, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_diversitypath, "cog_alpha_diversity_violin.pdf"), plot = p1_0, width = 12, height = 8)
addWorksheet(cog_diversity_wb, "alpha_diversity_data")
writeData(cog_diversity_wb, "alpha_diversity_data", data, rowNames = TRUE)

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

# 保存稀释曲线结果
ggsave(file.path(cog_diversitypath, "cog_rarefaction_individual.png"), plot = p2_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diversitypath, "cog_rarefaction_individual.pdf"), plot = p2_1, width = 10, height = 8)
ggsave(file.path(cog_diversitypath, "cog_rarefaction_group.png"), plot = p2_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diversitypath, "cog_rarefaction_group.pdf"), plot = p2_2, width = 10, height = 8)
ggsave(file.path(cog_diversitypath, "cog_rarefaction_group_sd.png"), plot = p2_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diversitypath, "cog_rarefaction_group_sd.pdf"), plot = p2_3, width = 10, height = 8)
addWorksheet(cog_diversity_wb, "rarefaction_data")
writeData(cog_diversity_wb, "rarefaction_data", raretab, rowNames = TRUE)

# function beta diversity -----
# 创建COG beta多样性分析目录
cog_betapath <- file.path(metapath, "cog", "beta")
dir.create(cog_betapath, recursive = TRUE)

# 创建COG beta多样性分析总表
cog_beta_wb <- createWorkbook()

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

# 保存排序分析结果
ggsave(file.path(cog_betapath, "cog_ordination_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_betapath, "cog_ordination_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(cog_betapath, "cog_ordination_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_betapath, "cog_ordination_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(cog_betapath, "cog_ordination_refined.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_betapath, "cog_ordination_refined.pdf"), plot = p3_3, width = 10, height = 8)
addWorksheet(cog_beta_wb, "ordination_data")
writeData(cog_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

# 保存群落功能差异检测结果
addWorksheet(cog_beta_wb, "metatest_results")
writeData(cog_beta_wb, "metatest_results", dat1, rowNames = TRUE)

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两分组群落功能差异检测结果
addWorksheet(cog_beta_wb, "pair_metatest_results")
writeData(cog_beta_wb, "pair_metatest_results", dat2, rowNames = TRUE)

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

# 保存Mantel分析结果
ggsave(file.path(cog_betapath, "cog_mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_betapath, "cog_mantel_test.pdf"), plot = p3_7, width = 10, height = 8)
addWorksheet(cog_beta_wb, "mantel_results")
writeData(cog_beta_wb, "mantel_results", data, rowNames = TRUE)

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

# 保存聚类分析结果
ggsave(file.path(cog_betapath, "cog_cluster_heatmap.png"), plot = p4, width = 12, height = 10, dpi = 300)
ggsave(file.path(cog_betapath, "cog_cluster_heatmap.pdf"), plot = p4, width = 12, height = 10)
ggsave(file.path(cog_betapath, "cog_cluster_dendrogram1.png"), plot = p4_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_betapath, "cog_cluster_dendrogram1.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(cog_betapath, "cog_cluster_dendrogram2.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_betapath, "cog_cluster_dendrogram2.pdf"), plot = p4_2, width = 10, height = 8)
addWorksheet(cog_beta_wb, "cluster_results")
writeData(cog_beta_wb, "cluster_results", dat, rowNames = TRUE)

#8 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)
res = Micro_tern.metf(ps=ps1,color = "Function"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图结果
ggsave(file.path(cog_betapath, "cog_ternary_plot.png"), plot = p15, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_betapath, "cog_ternary_plot.pdf"), plot = p15, width = 10, height = 8)
addWorksheet(cog_beta_wb, "ternary_data")
writeData(cog_beta_wb, "ternary_data", dat, rowNames = TRUE)

# 保存beta多样性总表
saveWorkbook(cog_beta_wb, file.path(cog_betapath, "cog_beta_diversity.xlsx"), overwrite = TRUE)

# function classification -----
# 创建COG功能分类目录
cog_classificationpath <- file.path(metapath, "cog", "classification")
dir.create(cog_classificationpath, recursive = TRUE)

# 创建COG功能分类总表
cog_classification_wb <- createWorkbook()

#9 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
rank_names(ps)
result = barMainplot.metf(ps = ps,
                          j =  "Category" ,
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

# 保存堆积柱状图结果
ggsave(file.path(cog_classificationpath, "cog_barplot_main.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_barplot_main.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(cog_classificationpath, "cog_barplot_secondary.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_barplot_secondary.pdf"), plot = p4_2, width = 12, height = 8)
addWorksheet(cog_classification_wb, "barplot_data")
writeData(cog_classification_wb, "barplot_data", databar, rowNames = TRUE)

#10 Ven.Upset.metf: 用于展示共有、特有的功能 ----
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

# 保存韦恩图和upset图结果
ggsave(file.path(cog_classificationpath, "cog_venn_diagram.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_venn_diagram.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(cog_classificationpath, "cog_upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_upset_plot.pdf"), plot = p10.2, width = 12, height = 8)
addWorksheet(cog_classification_wb, "venn_upset_data")
writeData(cog_classification_wb, "venn_upset_data", dat, rowNames = TRUE)

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

# 保存韦恩详细分析结果
ggsave(file.path(cog_classificationpath, "cog_venn_detail_bar.png"), plot = p7_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_venn_detail_bar.pdf"), plot = p7_1, width = 12, height = 8)
ggsave(file.path(cog_classificationpath, "cog_venn_detail_alluvial.png"), plot = p7_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_venn_detail_alluvial.pdf"), plot = p7_2, width = 12, height = 8)
addWorksheet(cog_classification_wb, "venn_detail_data")
writeData(cog_classification_wb, "venn_detail_data", dat, rowNames = TRUE)

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

# 保存花瓣图结果
ggsave(file.path(cog_classificationpath, "cog_flower_plot.png"), plot = p13.1, width = 10, height = 10, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_flower_plot.pdf"), plot = p13.1, width = 10, height = 10)
addWorksheet(cog_classification_wb, "flower_data")
writeData(cog_classification_wb, "flower_data", dat, rowNames = TRUE)

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

# 保存韦恩网络图结果
ggsave(file.path(cog_classificationpath, "cog_venn_network.png"), plot = p14, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_classificationpath, "cog_venn_network.pdf"), plot = p14, width = 10, height = 8)
addWorksheet(cog_classification_wb, "venn_network_data")
writeData(cog_classification_wb, "venn_network_data", dat, rowNames = TRUE)

# 保存总表
saveWorkbook(cog_classification_wb, file.path(cog_classificationpath, "cog_classification_analysis.xlsx"), overwrite = TRUE)
saveWorkbook(cog_diversity_wb, file.path(cog_diversitypath, "cog_diversity_analysis.xlsx"), overwrite = TRUE)


# function differential analysis----
# 创建COG差异分析目录
cog_diffpath <- file.path(metapath, "cog", "differential")
dir.create(cog_diffpath, recursive = TRUE)

# 创建COG差异分析总表
cog_diff_wb <- createWorkbook()
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

# 保存DESep2结果
ggsave(file.path(cog_diffpath, "cog_DESep2_plot1.png"), plot = p15.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_DESep2_plot1.pdf"), plot = p15.1, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_DESep2_plot2.png"), plot = p15.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_DESep2_plot2.pdf"), plot = p15.2, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_DESep2_plot3.png"), plot = p15.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_DESep2_plot3.pdf"), plot = p15.3, width = 10, height = 8)
addWorksheet(cog_diff_wb, "DESep2_results")
writeData(cog_diff_wb, "DESep2_results", dat, rowNames = TRUE)

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

# 保存EdgeR结果
ggsave(file.path(cog_diffpath, "cog_EdgeR_plot1.png"), plot = p16.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_EdgeR_plot1.pdf"), plot = p16.1, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_EdgeR_plot2.png"), plot = p16.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_EdgeR_plot2.pdf"), plot = p16.2, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_EdgeR_plot3.png"), plot = p16.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_EdgeR_plot3.pdf"), plot = p16.3, width = 10, height = 8)
addWorksheet(cog_diff_wb, "EdgeR_results")
writeData(cog_diff_wb, "EdgeR_results", dat, rowNames = TRUE)

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

# 保存t检验结果
ggsave(file.path(cog_diffpath, "cog_ttest_plot1.png"), plot = p17.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_ttest_plot1.pdf"), plot = p17.1, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_ttest_plot2.png"), plot = p17.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_ttest_plot2.pdf"), plot = p17.2, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_ttest_plot3.png"), plot = p17.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_ttest_plot3.pdf"), plot = p17.3, width = 10, height = 8)
addWorksheet(cog_diff_wb, "ttest_results")
writeData(cog_diff_wb, "ttest_results", dat, rowNames = TRUE)

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

# 保存非参数检验结果
ggsave(file.path(cog_diffpath, "cog_wilcox_plot1.png"), plot = p18.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_wilcox_plot1.pdf"), plot = p18.1, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_wilcox_plot2.png"), plot = p18.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_wilcox_plot2.pdf"), plot = p18.2, width = 10, height = 8)
ggsave(file.path(cog_diffpath, "cog_wilcox_plot3.png"), plot = p18.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_wilcox_plot3.pdf"), plot = p18.3, width = 10, height = 8)
addWorksheet(cog_diff_wb, "wilcox_results")
writeData(cog_diff_wb, "wilcox_results", dat, rowNames = TRUE)

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

# 保存stamp结果
ggsave(file.path(cog_diffpath, "cog_stamp.png"), plot = p19, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_diffpath, "cog_stamp.pdf"), plot = p19, width = 10, height = 8)
addWorksheet(cog_diff_wb, "stamp_plot_data")
writeData(cog_diff_wb, "stamp_plot_data", dat1, rowNames = TRUE)
addWorksheet(cog_diff_wb, "stamp_results")
writeData(cog_diff_wb, "stamp_results", dat2, rowNames = TRUE)

# 保存差异分析总表
saveWorkbook(cog_diff_wb, file.path(cog_diffpath, "cog_differential_analysis.xlsx"), overwrite = TRUE)

# biomarker identification-----
# 创建COG生物标志物目录
cog_biomarkerpath <- file.path(metapath, "cog", "biomarker")
dir.create(cog_biomarkerpath, recursive = TRUE)

# 创建COG生物标志物总表
cog_biomarker_wb <- createWorkbook()

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

# 保存交叉验证结果
ggsave(file.path(cog_biomarkerpath, "cog_rfcv.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_rfcv.pdf"), plot = prfcv, width = 10, height = 8)
addWorksheet(cog_biomarker_wb, "rfcv_results")
writeData(cog_biomarker_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#20 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

# 保存ROC结果
ggsave(file.path(cog_biomarkerpath, "cog_ROC_plot1.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_ROC_plot1.pdf"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(cog_biomarkerpath, "cog_ROC_plot2.png"), plot = p33.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_ROC_plot2.pdf"), plot = p33.2, width = 10, height = 8)
addWorksheet(cog_biomarker_wb, "ROC_results")
writeData(cog_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(cog_biomarkerpath, "cog_loadingPCA.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_loadingPCA.pdf"), plot = p34.1, width = 10, height = 8)
addWorksheet(cog_biomarker_wb, "loadingPCA_results")
writeData(cog_biomarker_wb, "loadingPCA_results", dat, rowNames = TRUE)

#22 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(cog_biomarker_wb, "svm_AUC")
writeData(cog_biomarker_wb, "svm_AUC", AUC, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "svm_importance")
writeData(cog_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#23 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(cog_biomarker_wb, "glm_AUC")
writeData(cog_biomarker_wb, "glm_AUC", AUC, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "glm_importance")
writeData(cog_biomarker_wb, "glm_importance", importance, rowNames = TRUE)

#24 xgboost.metf: xgboost筛选特征功能----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =pst, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(cog_biomarker_wb, "xgboost_accuracy")
writeData(cog_biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "xgboost_importance")
writeData(cog_biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#25 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(cog_biomarker_wb, "lasso_accuracy")
writeData(cog_biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "lasso_importance")
writeData(cog_biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#26 decisiontree.metf: ----
library(rpart)
res =decisiontree.metf(ps=ps.card, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(cog_biomarker_wb, "decisiontree_accuracy")
writeData(cog_biomarker_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "decisiontree_importance")
writeData(cog_biomarker_wb, "decisiontree_importance", importance, rowNames = TRUE)

#27 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(cog_biomarker_wb, "naivebayes_accuracy")
writeData(cog_biomarker_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "naivebayes_importance")
writeData(cog_biomarker_wb, "naivebayes_importance", importance, rowNames = TRUE)

#28 LDA.metf: LDA筛选特征功能-----
tablda = LDA.metf(ps = ps,
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

# 保存LDA结果
ggsave(file.path(cog_biomarkerpath, "cog_LDA.png"), plot = p35, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_LDA.pdf"), plot = p35, width = 10, height = 8)
addWorksheet(cog_biomarker_wb, "LDA_results")
writeData(cog_biomarker_wb, "LDA_results", dat, rowNames = TRUE)

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

# 保存随机森林结果
ggsave(file.path(cog_biomarkerpath, "cog_randomforest_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_randomforest_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(cog_biomarkerpath, "cog_randomforest_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_randomforest_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(cog_biomarkerpath, "cog_randomforest_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(cog_biomarkerpath, "cog_randomforest_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(cog_biomarker_wb, "randomforest_results")
writeData(cog_biomarker_wb, "randomforest_results", dat, rowNames = TRUE)

#30 nnet.metf: 神经网络筛选特征功能  ------
library(nnet)
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(cog_biomarker_wb, "nnet_accuracy")
writeData(cog_biomarker_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "nnet_importance")
writeData(cog_biomarker_wb, "nnet_importance", importance, rowNames = TRUE)

#31 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(cog_biomarker_wb, "bagging_accuracy")
writeData(cog_biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(cog_biomarker_wb, "bagging_importance")
writeData(cog_biomarker_wb, "bagging_importance", importance, rowNames = TRUE)

# 保存生物标志物总表
saveWorkbook(cog_biomarker_wb, file.path(cog_biomarkerpath, "cog_biomarker_analysis.xlsx"), overwrite = TRUE)

#network analysis -----
# 创建COG网络分析目录
cog_networkpath <- file.path(metapath, "cog", "network")
dir.create(cog_networkpath, recursive = TRUE)

# 创建COG网络分析总表
cog_network_wb <- createWorkbook()

#32 network.pip:网络分析主函数--------
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

# 保存网络分析主图
ggsave(file.path(cog_networkpath, "cog_network_main.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(cog_networkpath, "cog_network_main.pdf"), plot = p0, width = 12, height = 10)

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

# 保存网络属性数据
addWorksheet(cog_network_wb, "network_properties")
writeData(cog_network_wb, "network_properties", dat2, rowNames = TRUE)

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

# 保存样本网络属性数据
addWorksheet(cog_network_wb, "sample_network_properties")
writeData(cog_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

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

# 保存节点属性数据
addWorksheet(cog_network_wb, "node_properties")
writeData(cog_network_wb, "node_properties", nodepro2, rowNames = TRUE)

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

# 保存网络比较结果
addWorksheet(cog_network_wb, "network_comparison")
writeData(cog_network_wb, "network_comparison", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(cog_network_wb, file.path(cog_networkpath, "cog_network_analysis.xlsx"), overwrite = TRUE)

# vfdb 数据库-----
ps =EasyMultiOmics:: ps.vfdb %>% filter_OTU_ps(Top = 1000)
phyloseq::tax_table(ps)

# Function
# function diversity -----
# 创建VFDB多样性分析目录
vfdb_diversitypath <- file.path(metapath, "vfdb", "diversity")
dir.create(vfdb_diversitypath, recursive = TRUE)

# 创建VFDB多样性分析总表
vfdb_diversity_wb <- createWorkbook()

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

# 保存alpha多样性结果
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_box.png"), plot = p1_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_box.pdf"), plot = p1_1, width = 12, height = 8)
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_bar.png"), plot = p1_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_bar.pdf"), plot = p1_2, width = 12, height = 8)
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_boxbar.png"), plot = p1_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_boxbar.pdf"), plot = p1_3, width = 12, height = 8)
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_violin.png"), plot = p1_0, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_diversitypath, "vfdb_alpha_diversity_violin.pdf"), plot = p1_0, width = 12, height = 8)
addWorksheet(vfdb_diversity_wb, "alpha_diversity_data")
writeData(vfdb_diversity_wb, "alpha_diversity_data", data, rowNames = TRUE)

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

# 保存稀释曲线结果
ggsave(file.path(vfdb_diversitypath, "vfdb_rarefaction_individual.png"), plot = p2_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diversitypath, "vfdb_rarefaction_individual.pdf"), plot = p2_1, width = 10, height = 8)
ggsave(file.path(vfdb_diversitypath, "vfdb_rarefaction_group.png"), plot = p2_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diversitypath, "vfdb_rarefaction_group.pdf"), plot = p2_2, width = 10, height = 8)
ggsave(file.path(vfdb_diversitypath, "vfdb_rarefaction_group_sd.png"), plot = p2_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diversitypath, "vfdb_rarefaction_group_sd.pdf"), plot = p2_3, width = 10, height = 8)
addWorksheet(vfdb_diversity_wb, "rarefaction_data")
writeData(vfdb_diversity_wb, "rarefaction_data", raretab, rowNames = TRUE)

# function beta diversity -----
# 创建VFDB beta多样性分析目录
vfdb_betapath <- file.path(metapath, "vfdb", "beta")
dir.create(vfdb_betapath, recursive = TRUE)

# 创建VFDB beta多样性分析总表
vfdb_beta_wb <- createWorkbook()

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

# 保存排序分析结果
ggsave(file.path(vfdb_betapath, "vfdb_ordination_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_ordination_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(vfdb_betapath, "vfdb_ordination_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_ordination_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(vfdb_betapath, "vfdb_ordination_refined.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_ordination_refined.pdf"), plot = p3_3, width = 10, height = 8)
addWorksheet(vfdb_beta_wb, "ordination_data")
writeData(vfdb_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

#4 MetaTest.metf:群落功能差异检测#-------
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

# 保存群落功能差异检测结果
addWorksheet(vfdb_beta_wb, "metatest_results")
writeData(vfdb_beta_wb, "metatest_results", dat1, rowNames = TRUE)

#5 pairMetaTest.metf:两两分组群落功能差异检测#-------
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两分组群落功能差异检测结果
addWorksheet(vfdb_beta_wb, "pair_metatest_results")
writeData(vfdb_beta_wb, "pair_metatest_results", dat2, rowNames = TRUE)

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

# 保存Mantel分析结果
ggsave(file.path(vfdb_betapath, "vfdb_mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_mantel_test.pdf"), plot = p3_7, width = 10, height = 8)
addWorksheet(vfdb_beta_wb, "mantel_results")
writeData(vfdb_beta_wb, "mantel_results", data, rowNames = TRUE)

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

# 保存聚类分析结果
ggsave(file.path(vfdb_betapath, "vfdb_cluster_heatmap.png"), plot = p4, width = 12, height = 10, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_cluster_heatmap.pdf"), plot = p4, width = 12, height = 10)
ggsave(file.path(vfdb_betapath, "vfdb_cluster_dendrogram1.png"), plot = p4_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_cluster_dendrogram1.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(vfdb_betapath, "vfdb_cluster_dendrogram2.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_cluster_dendrogram2.pdf"), plot = p4_2, width = 10, height = 8)
addWorksheet(vfdb_beta_wb, "cluster_results")
writeData(vfdb_beta_wb, "cluster_results", dat, rowNames = TRUE)

#8 Micro_tern.metf: 三元图展示功能----
res = Micro_tern.metf(ps=ps %>% filter_OTU_ps(500),color = "Function"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图结果
ggsave(file.path(vfdb_betapath, "vfdb_ternary_plot.png"), plot = p15, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_betapath, "vfdb_ternary_plot.pdf"), plot = p15, width = 10, height = 8)
addWorksheet(vfdb_beta_wb, "ternary_data")
writeData(vfdb_beta_wb, "ternary_data", dat, rowNames = TRUE)

# 保存beta多样性总表
saveWorkbook(vfdb_beta_wb, file.path(vfdb_betapath, "vfdb_beta_diversity.xlsx"), overwrite = TRUE)

# function classification -----
# 创建VFDB功能分类目录
vfdb_classificationpath <- file.path(metapath, "vfdb", "classification")
dir.create(vfdb_classificationpath, recursive = TRUE)

# 创建VFDB功能分类总表
vfdb_classification_wb <- createWorkbook()

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

# 保存堆积柱状图结果
ggsave(file.path(vfdb_classificationpath, "vfdb_barplot_main.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_barplot_main.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(vfdb_classificationpath, "vfdb_barplot_secondary.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_barplot_secondary.pdf"), plot = p4_2, width = 12, height = 8)
addWorksheet(vfdb_classification_wb, "barplot_data")
writeData(vfdb_classification_wb, "barplot_data", databar, rowNames = TRUE)

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

# 保存韦恩图和upset图结果
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_diagram.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_diagram.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(vfdb_classificationpath, "vfdb_upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_upset_plot.pdf"), plot = p10.2, width = 12, height = 8)
addWorksheet(vfdb_classification_wb, "venn_upset_data")
writeData(vfdb_classification_wb, "venn_upset_data", dat, rowNames = TRUE)

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

# 保存韦恩详细分析结果
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_detail_bar.png"), plot = p7_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_detail_bar.pdf"), plot = p7_1, width = 12, height = 8)
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_detail_alluvial.png"), plot = p7_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_detail_alluvial.pdf"), plot = p7_2, width = 12, height = 8)
addWorksheet(vfdb_classification_wb, "venn_detail_data")
writeData(vfdb_classification_wb, "venn_detail_data", dat, rowNames = TRUE)

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

# 保存花瓣图结果
ggsave(file.path(vfdb_classificationpath, "vfdb_flower_plot.png"), plot = p13.1, width = 10, height = 10, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_flower_plot.pdf"), plot = p13.1, width = 10, height = 10)
addWorksheet(vfdb_classification_wb, "flower_data")
writeData(vfdb_classification_wb, "flower_data", dat, rowNames = TRUE)

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

# 保存韦恩网络图结果
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_network.png"), plot = p14, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_classificationpath, "vfdb_venn_network.pdf"), plot = p14, width = 10, height = 8)
addWorksheet(vfdb_classification_wb, "venn_network_data")
writeData(vfdb_classification_wb, "venn_network_data", dat, rowNames = TRUE)

# 保存功能分类总表
saveWorkbook(vfdb_classification_wb, file.path(vfdb_classificationpath, "vfdb_classification_analysis.xlsx"), overwrite = TRUE)

# 保存多样性分析总表
saveWorkbook(vfdb_diversity_wb, file.path(vfdb_diversitypath, "vfdb_diversity_analysis.xlsx"), overwrite = TRUE)

# function differential analysis----
# 创建VFDB差异分析目录
vfdb_diffpath <- file.path(metapath, "vfdb", "differential")
dir.create(vfdb_diffpath, recursive = TRUE)

# 创建VFDB差异分析总表
vfdb_diff_wb <- createWorkbook()

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

# 保存DESep2结果
ggsave(file.path(vfdb_diffpath, "vfdb_DESep2_plot1.png"), plot = p15.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_DESep2_plot1.pdf"), plot = p15.1, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_DESep2_plot2.png"), plot = p15.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_DESep2_plot2.pdf"), plot = p15.2, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_DESep2_plot3.png"), plot = p15.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_DESep2_plot3.pdf"), plot = p15.3, width = 10, height = 8)
addWorksheet(vfdb_diff_wb, "DESep2_results")
writeData(vfdb_diff_wb, "DESep2_results", dat, rowNames = TRUE)

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

# 保存EdgeR结果
ggsave(file.path(vfdb_diffpath, "vfdb_EdgeR_plot1.png"), plot = p16.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_EdgeR_plot1.pdf"), plot = p16.1, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_EdgeR_plot2.png"), plot = p16.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_EdgeR_plot2.pdf"), plot = p16.2, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_EdgeR_plot3.png"), plot = p16.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_EdgeR_plot3.pdf"), plot = p16.3, width = 10, height = 8)
addWorksheet(vfdb_diff_wb, "EdgeR_results")
writeData(vfdb_diff_wb, "EdgeR_results", dat, rowNames = TRUE)

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

# 保存t检验结果
ggsave(file.path(vfdb_diffpath, "vfdb_ttest_plot1.png"), plot = p17.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_ttest_plot1.pdf"), plot = p17.1, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_ttest_plot2.png"), plot = p17.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_ttest_plot2.pdf"), plot = p17.2, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_ttest_plot3.png"), plot = p17.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_ttest_plot3.pdf"), plot = p17.3, width = 10, height = 8)
addWorksheet(vfdb_diff_wb, "ttest_results")
writeData(vfdb_diff_wb, "ttest_results", dat, rowNames = TRUE)

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

# 保存非参数检验结果
ggsave(file.path(vfdb_diffpath, "vfdb_wilcox_plot1.png"), plot = p18.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_wilcox_plot1.pdf"), plot = p18.1, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_wilcox_plot2.png"), plot = p18.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_wilcox_plot2.pdf"), plot = p18.2, width = 10, height = 8)
ggsave(file.path(vfdb_diffpath, "vfdb_wilcox_plot3.png"), plot = p18.3, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_wilcox_plot3.pdf"), plot = p18.3, width = 10, height = 8)
addWorksheet(vfdb_diff_wb, "wilcox_results")
writeData(vfdb_diff_wb, "wilcox_results", dat, rowNames = TRUE)

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

# 保存stamp结果
ggsave(file.path(vfdb_diffpath, "vfdb_stamp.png"), plot = p19, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_diffpath, "vfdb_stamp.pdf"), plot = p19, width = 10, height = 8)
addWorksheet(vfdb_diff_wb, "stamp_plot_data")
writeData(vfdb_diff_wb, "stamp_plot_data", dat1, rowNames = TRUE)
addWorksheet(vfdb_diff_wb, "stamp_results")
writeData(vfdb_diff_wb, "stamp_results", dat2, rowNames = TRUE)

# 保存差异分析总表
saveWorkbook(vfdb_diff_wb, file.path(vfdb_diffpath, "vfdb_differential_analysis.xlsx"), overwrite = TRUE)

# biomarker identification-----
# 创建VFDB生物标志物目录
vfdb_biomarkerpath <- file.path(metapath, "vfdb", "biomarker")
dir.create(vfdb_biomarkerpath, recursive = TRUE)

# 创建VFDB生物标志物总表
vfdb_biomarker_wb <- createWorkbook()

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

# 保存交叉验证结果
ggsave(file.path(vfdb_biomarkerpath, "vfdb_rfcv.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_rfcv.pdf"), plot = prfcv, width = 10, height = 8)
addWorksheet(vfdb_biomarker_wb, "rfcv_results")
writeData(vfdb_biomarker_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#20 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

# 保存ROC结果
ggsave(file.path(vfdb_biomarkerpath, "vfdb_ROC_plot1.png"), plot = p33.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_ROC_plot1.pdf"), plot = p33.1, width = 10, height = 8)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_ROC_plot2.png"), plot = p33.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_ROC_plot2.pdf"), plot = p33.2, width = 10, height = 8)
addWorksheet(vfdb_biomarker_wb, "ROC_results")
writeData(vfdb_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(vfdb_biomarkerpath, "vfdb_loadingPCA.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_loadingPCA.pdf"), plot = p34.1, width = 10, height = 8)
addWorksheet(vfdb_biomarker_wb, "loadingPCA_results")
writeData(vfdb_biomarker_wb, "loadingPCA_results", dat, rowNames = TRUE)

#22 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = ps%>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(vfdb_biomarker_wb, "svm_AUC")
writeData(vfdb_biomarker_wb, "svm_AUC", AUC, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "svm_importance")
writeData(vfdb_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#23 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(vfdb_biomarker_wb, "glm_AUC")
writeData(vfdb_biomarker_wb, "glm_AUC", AUC, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "glm_importance")
writeData(vfdb_biomarker_wb, "glm_importance", importance, rowNames = TRUE)

#24 xgboost.metf: xgboost筛选特征功能----
library(xgboost)
library(Ckmeans.1d.dp)
res = xgboost.metf(ps =ps, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存XGBoost结果
addWorksheet(vfdb_biomarker_wb, "xgboost_accuracy")
writeData(vfdb_biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "xgboost_importance")
writeData(vfdb_biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#25 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  ps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(vfdb_biomarker_wb, "lasso_accuracy")
writeData(vfdb_biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "lasso_importance")
writeData(vfdb_biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#26 decisiontree.micro: ----
library(rpart)
res =decisiontree.metf(ps=ps.card, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(vfdb_biomarker_wb, "decisiontree_accuracy")
writeData(vfdb_biomarker_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "decisiontree_importance")
writeData(vfdb_biomarker_wb, "decisiontree_importance", importance, rowNames = TRUE)

#27 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(vfdb_biomarker_wb, "naivebayes_accuracy")
writeData(vfdb_biomarker_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "naivebayes_importance")
writeData(vfdb_biomarker_wb, "naivebayes_importance", importance, rowNames = TRUE)

#28 LDA.metf: LDA筛选特征功能-----
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

# 保存LDA结果
ggsave(file.path(vfdb_biomarkerpath, "vfdb_LDA.png"), plot = p35, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_LDA.pdf"), plot = p35, width = 10, height = 8)
addWorksheet(vfdb_biomarker_wb, "LDA_results")
writeData(vfdb_biomarker_wb, "LDA_results", dat, rowNames = TRUE)

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

# 保存随机森林结果
ggsave(file.path(vfdb_biomarkerpath, "vfdb_randomforest_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_randomforest_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_randomforest_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_randomforest_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_randomforest_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(vfdb_biomarkerpath, "vfdb_randomforest_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(vfdb_biomarker_wb, "randomforest_results")
writeData(vfdb_biomarker_wb, "randomforest_results", dat, rowNames = TRUE)

#30 nnet.metf: 神经网络筛选特征功能  ------
library(nnet)
res =nnet.metf(ps=ps, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(vfdb_biomarker_wb, "nnet_accuracy")
writeData(vfdb_biomarker_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "nnet_importance")
writeData(vfdb_biomarker_wb, "nnet_importance", importance, rowNames = TRUE)

#31 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps = ps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(vfdb_biomarker_wb, "bagging_accuracy")
writeData(vfdb_biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(vfdb_biomarker_wb, "bagging_importance")
writeData(vfdb_biomarker_wb, "bagging_importance", importance, rowNames = TRUE)

# 保存生物标志物总表
saveWorkbook(vfdb_biomarker_wb, file.path(vfdb_biomarkerpath, "vfdb_biomarker_analysis.xlsx"), overwrite = TRUE)

#network analysis -----
# 创建VFDB网络分析目录
vfdb_networkpath <- file.path(metapath, "vfdb", "network")
dir.create(vfdb_networkpath, recursive = TRUE)

# 创建VFDB网络分析总表
vfdb_network_wb <- createWorkbook()

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

# 保存网络分析主图
ggsave(file.path(vfdb_networkpath, "vfdb_network_main.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(vfdb_networkpath, "vfdb_network_main.pdf"), plot = p0, width = 12, height = 10)

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

# 保存网络属性数据
addWorksheet(vfdb_network_wb, "network_properties")
writeData(vfdb_network_wb, "network_properties", dat2, rowNames = TRUE)

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

# 保存样本网络属性数据
addWorksheet(vfdb_network_wb, "sample_network_properties")
writeData(vfdb_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

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

# 保存节点属性数据
addWorksheet(vfdb_network_wb, "node_properties")
writeData(vfdb_network_wb, "node_properties", nodepro2, rowNames = TRUE)

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

# 保存网络比较结果
addWorksheet(vfdb_network_wb, "network_comparison")
writeData(vfdb_network_wb, "network_comparison", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(vfdb_network_wb, file.path(vfdb_networkpath, "vfdb_network_analysis.xlsx"), overwrite = TRUE)













