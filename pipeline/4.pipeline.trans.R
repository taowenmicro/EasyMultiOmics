rm(list=ls())
# 转录组分析流程------
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(EasyMultiOmics)
library(EasyMultiOmics.db)
library(RColorBrewer)
library(openxlsx)

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
# 创建转录组数据预处理目录
trans_preprocesspath =  "./result/transcriptome/"
dir.create(trans_preprocesspath, recursive = TRUE)

# 创建转录组数据预处理总表
trans_preprocess_wb <- createWorkbook()

#1 按照通路（KEGG）合并基因#------
detach("package:mia", unload = TRUE)
ps.trans2 = trans_kegg (ps.trans)
phyloseq::tax_table(ps.trans2) %>% head()


# 保存KEGG合并数据
addWorksheet(trans_preprocess_wb, "kegg_merged")
writeData(trans_preprocess_wb, "kegg_merged", phyloseq::tax_table(ps.trans2), rowNames = TRUE)

#2 基于 Mkegg 通路合并#-------
ps.trans2 =trans_mkegg ( ps.trans)

# 保存Mkegg合并数据
addWorksheet(trans_preprocess_wb, "mkegg_merged")
writeData(trans_preprocess_wb, "mkegg_merged", phyloseq::tax_table(ps.trans2), rowNames = TRUE)

#3 基于reaction 层面合并#----
ps.trans3 =trans_rekegg ( ps.trans)
ps.trans3

# 保存reaction合并数据
addWorksheet(trans_preprocess_wb, "reaction_merged")
writeData(trans_preprocess_wb, "reaction_merged", phyloseq::tax_table(ps.trans3), rowNames = TRUE)

# 保存数据预处理总表
saveWorkbook(trans_preprocess_wb, file.path(trans_preprocesspath, "transcriptome_preprocess.xlsx"), overwrite = TRUE)

# ordinate analysis------
# 创建转录组排序分析目录
trans_ordinatepath <- file.path(trans_preprocesspath, "transcriptome", "ordination")
dir.create(trans_ordinatepath, recursive = TRUE)

# 创建转录组排序分析总表
trans_ordinate_wb <- createWorkbook()

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

# 保存排序分析结果
ggsave(file.path(trans_ordinatepath, "ordination_basic.png"), plot = p3_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_ordinatepath, "ordination_basic.pdf"), plot = p3_1, width = 10, height = 8)
ggsave(file.path(trans_ordinatepath, "ordination_labeled.png"), plot = p3_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_ordinatepath, "ordination_labeled.pdf"), plot = p3_2, width = 10, height = 8)
ggsave(file.path(trans_ordinatepath, "ordination_refined.png"), plot = p3_3, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_ordinatepath, "ordination_refined.pdf"), plot = p3_3, width = 10, height = 8)
addWorksheet(trans_ordinate_wb, "ordination_data")
writeData(trans_ordinate_wb, "ordination_data", plotdata, rowNames = TRUE)

#5 transTest.metf:群落功能差异检测#-------
dat1 = transtest.trans(ps =ps.trans, Micromet = "adonis", dist = "bray")
dat1

# 保存群落功能差异检测结果
addWorksheet(trans_ordinate_wb, "transtest_results")
writeData(trans_ordinate_wb, "transtest_results", dat1, rowNames = TRUE)

#6 pairtranstest.trans:两两分组功能差异检测#-------
dat2 = pairtranstest.trans(ps =ps.trans, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两分组功能差异检测结果
addWorksheet(trans_ordinate_wb, "pair_transtest_results")
writeData(trans_ordinate_wb, "pair_transtest_results", dat2, rowNames = TRUE)

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

# 保存Mantel分析结果
ggsave(file.path(trans_ordinatepath, "mantel_test.png"), plot = p3_7, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_ordinatepath, "mantel_test.pdf"), plot = p3_7, width = 10, height = 8)
addWorksheet(trans_ordinate_wb, "mantel_results")
writeData(trans_ordinate_wb, "mantel_results", data, rowNames = TRUE)

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

# 保存聚类分析结果
ggsave(file.path(trans_ordinatepath, "cluster_heatmap.png"), plot = p4, width = 12, height = 10, dpi = 300)
ggsave(file.path(trans_ordinatepath, "cluster_heatmap.pdf"), plot = p4, width = 12, height = 10)
ggsave(file.path(trans_ordinatepath, "cluster_dendrogram1.png"), plot = p4_1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_ordinatepath, "cluster_dendrogram1.pdf"), plot = p4_1, width = 10, height = 8)
ggsave(file.path(trans_ordinatepath, "cluster_dendrogram2.png"), plot = p4_2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_ordinatepath, "cluster_dendrogram2.pdf"), plot = p4_2, width = 10, height = 8)
addWorksheet(trans_ordinate_wb, "cluster_results")
writeData(trans_ordinate_wb, "cluster_results", dat, rowNames = TRUE)

# 保存排序分析总表
saveWorkbook(trans_ordinate_wb, file.path(trans_ordinatepath, "transcriptome_ordination.xlsx"), overwrite = TRUE)

#  gene classification-------
# 创建转录组基因分类目录
trans_classificationpath <- file.path(trans_preprocesspath, "transcriptome", "classification")
dir.create(trans_classificationpath, recursive = TRUE)

# 创建转录组基因分类总表
trans_classification_wb <- createWorkbook()

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

# 保存韦恩图和upset图结果
ggsave(file.path(trans_classificationpath, "venn_diagram.png"), plot = p10.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "venn_diagram.pdf"), plot = p10.1, width = 10, height = 8)
ggsave(file.path(trans_classificationpath, "upset_plot.png"), plot = p10.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "upset_plot.pdf"), plot = p10.2, width = 12, height = 8)
addWorksheet(trans_classification_wb, "venn_upset_data")
writeData(trans_classification_wb, "venn_upset_data", dat, rowNames = TRUE)

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

# 保存花瓣图结果
ggsave(file.path(trans_classificationpath, "flower_plot.png"), plot = p13.1, width = 10, height = 10, dpi = 300)
ggsave(file.path(trans_classificationpath, "flower_plot.pdf"), plot = p13.1, width = 10, height = 10)
addWorksheet(trans_classification_wb, "flower_data")
writeData(trans_classification_wb, "flower_data", dat, rowNames = TRUE)

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

# 保存韦恩网络图结果
ggsave(file.path(trans_classificationpath, "venn_network.png"), plot = p14, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "venn_network.pdf"), plot = p14, width = 10, height = 8)
addWorksheet(trans_classification_wb, "venn_network_data")
writeData(trans_classification_wb, "venn_network_data", dat, rowNames = TRUE)

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

# 保存堆积柱状图结果
ggsave(file.path(trans_classificationpath, "barplot_main.png"), plot = p4_1, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "barplot_main.pdf"), plot = p4_1, width = 12, height = 8)
ggsave(file.path(trans_classificationpath, "barplot_secondary.png"), plot = p4_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "barplot_secondary.pdf"), plot = p4_2, width = 12, height = 8)
addWorksheet(trans_classification_wb, "barplot_data")
writeData(trans_classification_wb, "barplot_data", databar, rowNames = TRUE)

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

# 保存聚类堆积柱状图结果
ggsave(file.path(trans_classificationpath, "cluster_bar1.png"), plot = result[[1]], width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "cluster_bar1.pdf"), plot = result[[1]], width = 12, height = 8)
ggsave(file.path(trans_classificationpath, "cluster_bar2.png"), plot = p5_2, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "cluster_bar2.pdf"), plot = p5_2, width = 12, height = 8)
ggsave(file.path(trans_classificationpath, "cluster_bar3.png"), plot = p5_3, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "cluster_bar3.pdf"), plot = p5_3, width = 12, height = 8)
ggsave(file.path(trans_classificationpath, "cluster_bar4.png"), plot = p5_4, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_classificationpath, "cluster_bar4.pdf"), plot = p5_4, width = 12, height = 8)
addWorksheet(trans_classification_wb, "cluster_bar_data")
writeData(trans_classification_wb, "cluster_bar_data", clubardata, rowNames = TRUE)

# 保存基因分类总表
saveWorkbook(trans_classification_wb, file.path(trans_classificationpath, "transcriptome_classification.xlsx"), overwrite = TRUE)

# gene difference analysis
# 创建转录组差异分析目录
trans_diffpath <- file.path(trans_preprocesspath, "transcriptome", "differential")
dir.create(trans_diffpath, recursive = TRUE)

# 创建转录组差异分析总表
trans_diff_wb <- createWorkbook()

#14 aldex2.trans-------
library(ALDEx2)
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

# 保存ALDEx2结果
addWorksheet(trans_diff_wb, "aldex2_WT_OE")
writeData(trans_diff_wb, "aldex2_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "aldex2_WT_KO")
writeData(trans_diff_wb, "aldex2_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "aldex2_OE_KO")
writeData(trans_diff_wb, "aldex2_OE_KO", dat3, rowNames = TRUE)

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

# 保存ANCOM结果
addWorksheet(trans_diff_wb, "ancom_WT_OE")
writeData(trans_diff_wb, "ancom_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "ancom_WT_KO")
writeData(trans_diff_wb, "ancom_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "ancom_OE_KO")
writeData(trans_diff_wb, "ancom_OE_KO", dat3, rowNames = TRUE)

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

# 保存corncob结果
addWorksheet(trans_diff_wb, "corncob_WT_OE")
writeData(trans_diff_wb, "corncob_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "corncob_WT_KO")
writeData(trans_diff_wb, "corncob_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "corncob_OE_KO")
writeData(trans_diff_wb, "corncob_OE_KO", dat3, rowNames = TRUE)

#17 lefse.trans: ----------
dat = lefse.trans(ps = ps.trans %>% filter_OTU_ps(50),group =  "Group",alpha = 0.05)

#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

# 保存LEfSe结果
addWorksheet(trans_diff_wb, "lefse_WT_OE")
writeData(trans_diff_wb, "lefse_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "lefse_WT_KO")
writeData(trans_diff_wb, "lefse_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "lefse_OE_KO")
writeData(trans_diff_wb, "lefse_OE_KO", dat3, rowNames = TRUE)

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

# 保存limma结果
addWorksheet(trans_diff_wb, "limma_WT_OE")
writeData(trans_diff_wb, "limma_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "limma_WT_KO")
writeData(trans_diff_wb, "limma_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "limma_OE_KO")
writeData(trans_diff_wb, "limma_OE_KO", dat3, rowNames = TRUE)

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

# 保存MaAsLin2结果
addWorksheet(trans_diff_wb, "maaslin2_WT_OE")
writeData(trans_diff_wb, "maaslin2_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "maaslin2_WT_KO")
writeData(trans_diff_wb, "maaslin2_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "maaslin2_OE_KO")
writeData(trans_diff_wb, "maaslin2_OE_KO", dat3, rowNames = TRUE)

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

# 保存metagenomeSeq结果
addWorksheet(trans_diff_wb, "metaseq_WT_OE")
writeData(trans_diff_wb, "metaseq_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "metaseq_WT_KO")
writeData(trans_diff_wb, "metaseq_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "metaseq_OE_KO")
writeData(trans_diff_wb, "metaseq_OE_KO", dat3, rowNames = TRUE)

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

# 保存wilcox样本结果
addWorksheet(trans_diff_wb, "wilcox_sampl_WT_OE")
writeData(trans_diff_wb, "wilcox_sampl_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "wilcox_sampl_WT_KO")
writeData(trans_diff_wb, "wilcox_sampl_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "wilcox_sampl_OE_KO")
writeData(trans_diff_wb, "wilcox_sampl_OE_KO", dat3, rowNames = TRUE)

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

# 保存wilcox CLR结果
addWorksheet(trans_diff_wb, "wilcox_clr_WT_OE")
writeData(trans_diff_wb, "wilcox_clr_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "wilcox_clr_WT_KO")
writeData(trans_diff_wb, "wilcox_clr_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "wilcox_clr_OE_KO")
writeData(trans_diff_wb, "wilcox_clr_OE_KO", dat3, rowNames = TRUE)

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

# 保存MaAsLin2 (rare=TRUE)结果
addWorksheet(trans_diff_wb, "maaslin2_rare_WT_OE")
writeData(trans_diff_wb, "maaslin2_rare_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "maaslin2_rare_WT_KO")
writeData(trans_diff_wb, "maaslin2_rare_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "maaslin2_rare_OE_KO")
writeData(trans_diff_wb, "maaslin2_rare_OE_KO", dat3, rowNames = TRUE)

#24 ancombc.trans:-------
library("ANCOMBC")
library(DT)
dat = ancombc.trans(ps =  ps.trans %>% filter_OTU_ps(100),
                    group = "Group",
                    alpha = 0.05)
#--提取差异分析详细列表
dat1 = dat$WT_OE
head(dat1)

dat2 = dat$WT_KO
head(dat2)

dat3 = dat$OE_KO
head(dat3)

# 保存ANCOM-BC结果
addWorksheet(trans_diff_wb, "ancombc_WT_OE")
writeData(trans_diff_wb, "ancombc_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "ancombc_WT_KO")
writeData(trans_diff_wb, "ancombc_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "ancombc_OE_KO")
writeData(trans_diff_wb, "ancombc_OE_KO", dat3, rowNames = TRUE)

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

# 保存ANCOM-BC2结果
addWorksheet(trans_diff_wb, "ancombc2_WT_OE")
writeData(trans_diff_wb, "ancombc2_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "ancombc2_WT_KO")
writeData(trans_diff_wb, "ancombc2_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "ancombc2_OE_KO")
writeData(trans_diff_wb, "ancombc2_OE_KO", dat3, rowNames = TRUE)

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

# 保存ZicoSeq结果
addWorksheet(trans_diff_wb, "zicoseq_WT_OE")
writeData(trans_diff_wb, "zicoseq_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "zicoseq_WT_KO")
writeData(trans_diff_wb, "zicoseq_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "zicoseq_OE_KO")
writeData(trans_diff_wb, "zicoseq_OE_KO", dat3, rowNames = TRUE)

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

# 保存LinDA结果
addWorksheet(trans_diff_wb, "linda_WT_OE")
writeData(trans_diff_wb, "linda_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "linda_WT_KO")
writeData(trans_diff_wb, "linda_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "linda_OE_KO")
writeData(trans_diff_wb, "linda_OE_KO", dat3, rowNames = TRUE)

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

# 保存dacomp结果
addWorksheet(trans_diff_wb, "dacomp_WT_OE")
writeData(trans_diff_wb, "dacomp_WT_OE", dat1, rowNames = TRUE)
addWorksheet(trans_diff_wb, "dacomp_WT_KO")
writeData(trans_diff_wb, "dacomp_WT_KO", dat2, rowNames = TRUE)
addWorksheet(trans_diff_wb, "dacomp_OE_KO")
writeData(trans_diff_wb, "dacomp_OE_KO", dat3, rowNames = TRUE)

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

# 保存EdgeR结果
ggsave(file.path(trans_diffpath, "EdgeR_plot1.png"), plot = p16.1[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "EdgeR_plot1.pdf"), plot = p16.1[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "EdgeR_plot2.png"), plot = p16.2[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "EdgeR_plot2.pdf"), plot = p16.2[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "EdgeR_plot3.png"), plot = p16.3[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "EdgeR_plot3.pdf"), plot = p16.3[[1]], width = 10, height = 8)
addWorksheet(trans_diff_wb, "EdgeR_results")
writeData(trans_diff_wb, "EdgeR_results", dat, rowNames = TRUE)

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

# 保存DESep2结果
ggsave(file.path(trans_diffpath, "DESep2_plot1.png"), plot = p15.1[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "DESep2_plot1.pdf"), plot = p15.1[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "DESep2_plot2.png"), plot = p15.2[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "DESep2_plot2.pdf"), plot = p15.2[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "DESep2_plot3.png"), plot = p15.3[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "DESep2_plot3.pdf"), plot = p15.3[[1]], width = 10, height = 8)
addWorksheet(trans_diff_wb, "DESep2_results")
writeData(trans_diff_wb, "DESep2_results", dat, rowNames = TRUE)

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

# 保存t检验结果
ggsave(file.path(trans_diffpath, "ttest_overview.png"), plot = p17, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "ttest_overview.pdf"), plot = p17, width = 12, height = 8)
ggsave(file.path(trans_diffpath, "ttest_plot1.png"), plot = p17.1[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "ttest_plot1.pdf"), plot = p17.1[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "ttest_plot2.png"), plot = p17.2[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "ttest_plot2.pdf"), plot = p17.2[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "ttest_plot3.png"), plot = p17.3[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "ttest_plot3.pdf"), plot = p17.3[[1]], width = 10, height = 8)
addWorksheet(trans_diff_wb, "ttest_data")
writeData(trans_diff_wb, "ttest_data", data, rowNames = TRUE)
addWorksheet(trans_diff_wb, "ttest_inter_union")
writeData(trans_diff_wb, "ttest_inter_union", inter_union, rowNames = TRUE)

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

# 保存非参数检验结果
ggsave(file.path(trans_diffpath, "wilcox_overview.png"), plot = p18, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "wilcox_overview.pdf"), plot = p18, width = 12, height = 8)
ggsave(file.path(trans_diffpath, "wilcox_plot1.png"), plot = p18.1[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "wilcox_plot1.pdf"), plot = p18.1[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "wilcox_plot2.png"), plot = p18.2[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "wilcox_plot2.pdf"), plot = p18.2[[1]], width = 10, height = 8)
ggsave(file.path(trans_diffpath, "wilcox_plot3.png"), plot = p18.3[[1]], width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "wilcox_plot3.pdf"), plot = p18.3[[1]], width = 10, height = 8)
addWorksheet(trans_diff_wb, "wilcox_data")
writeData(trans_diff_wb, "wilcox_data", data, rowNames = TRUE)
addWorksheet(trans_diff_wb, "wilcox_inter_union")
writeData(trans_diff_wb, "wilcox_inter_union", inter_union, rowNames = TRUE)

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

# 保存聚类火山图结果
ggsave(file.path(trans_diffpath, "volcano_cluster1.png"), plot = p29.1, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "volcano_cluster1.pdf"), plot = p29.1, width = 12, height = 8)
ggsave(file.path(trans_diffpath, "volcano_cluster2.png"), plot = p29.2, width = 12, height = 8, dpi = 300)
ggsave(file.path(trans_diffpath, "volcano_cluster2.pdf"), plot = p29.2, width = 12, height = 8)
addWorksheet(trans_diff_wb, "volcano_cluster_data")
writeData(trans_diff_wb, "volcano_cluster_data", dat, rowNames = TRUE)

# 保存差异分析总表
saveWorkbook(trans_diff_wb, file.path(trans_diffpath, "transcriptome_differential_analysis.xlsx"), overwrite = TRUE)

# biomarker identification----
# 创建转录组生物标志物目录
trans_biomarkerpath <- file.path(trans_preprocesspath, "transcriptome", "biomarker")
dir.create(trans_biomarkerpath, recursive = TRUE)

# 创建转录组生物标志物总表
trans_biomarker_wb <- createWorkbook()

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


id = sample_data( ps.trans)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps.trans %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

#33 rfcv.trans :交叉验证结果-------
result =rfcv.trans(ps = ps.trans %>%filter_OTU_ps(200),group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]
prfcv
# result[[2]]# plotdata
rfcvtable = result[[3]]
rfcvtable

# 保存交叉验证结果
ggsave(file.path(trans_biomarkerpath, "rfcv.png"), plot = prfcv, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_biomarkerpath, "rfcv.pdf"), plot = prfcv, width = 10, height = 8)
addWorksheet(trans_biomarker_wb, "rfcv_results")
writeData(trans_biomarker_wb, "rfcv_results", rfcvtable, rowNames = TRUE)

#34 randomforest.trans: 随机森林筛选特征功能----
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

# 保存随机森林结果
ggsave(file.path(trans_biomarkerpath, "randomforest_plot1.png"), plot = p42.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_biomarkerpath, "randomforest_plot1.pdf"), plot = p42.1, width = 10, height = 8)
ggsave(file.path(trans_biomarkerpath, "randomforest_plot2.png"), plot = p42.2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_biomarkerpath, "randomforest_plot2.pdf"), plot = p42.2, width = 10, height = 8)
ggsave(file.path(trans_biomarkerpath, "randomforest_plot4.png"), plot = p42.4, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_biomarkerpath, "randomforest_plot4.pdf"), plot = p42.4, width = 10, height = 8)
addWorksheet(trans_biomarker_wb, "randomforest_results")
writeData(trans_biomarker_wb, "randomforest_results", dat, rowNames = TRUE)

#35 loadingPCA.trans:载荷矩阵筛选特征功能------
res = loadingPCA.trans(ps =  ps.trans,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
ggsave(file.path(trans_biomarkerpath, "loadingPCA.png"), plot = p34.1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_biomarkerpath, "loadingPCA.pdf"), plot = p34.1, width = 10, height = 8)
addWorksheet(trans_biomarker_wb, "loadingPCA_results")
writeData(trans_biomarker_wb, "loadingPCA_results", dat, rowNames = TRUE)

#36 svm.trans------
res <- svm.trans(ps = ps.trans %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
addWorksheet(trans_biomarker_wb, "svm_AUC")
writeData(trans_biomarker_wb, "svm_AUC", AUC, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "svm_importance")
writeData(trans_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#37 glm.trans---------
res <- glm.trans(ps = ps.trans %>% filter_OTU_ps(500), k = 5)
accuracy= res[[1]]
accuracy
importance = res[[2]]
importance

# 保存GLM结果
addWorksheet(trans_biomarker_wb, "glm_accuracy")
writeData(trans_biomarker_wb, "glm_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "glm_importance")
writeData(trans_biomarker_wb, "glm_importance", importance, rowNames = TRUE)

#38 xgboost.trans-----
res = xgboost.trans(ps = ps.trans, seed = 200, top = 200 )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance
importance <- as.data.frame(importance$importance)

# 保存XGBoost结果
addWorksheet(trans_biomarker_wb, "xgboost_accuracy")
writeData(trans_biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "xgboost_importance")
writeData(trans_biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#39 decisiontree ------
res = decisiontree.trans(ps=ps.trans, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
addWorksheet(trans_biomarker_wb, "decisiontree_accuracy")
writeData(trans_biomarker_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "decisiontree_importance")
writeData(trans_biomarker_wb, "decisiontree_importance", importance, rowNames = TRUE)

#40 bagging.trans-----
res = bagging.trans(ps =  ps.trans , top = 20, seed = 110, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
addWorksheet(trans_biomarker_wb, "bagging_accuracy")
writeData(trans_biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "bagging_importance")
writeData(trans_biomarker_wb, "bagging_importance", importance, rowNames = TRUE)

#41 lasso.trans-----
res = lasso.trans (ps =  pst, top = 500, seed = 1, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
addWorksheet(trans_biomarker_wb, "lasso_accuracy")
writeData(trans_biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "lasso_importance")
writeData(trans_biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#42 LDA.trans--------
res= LDA.trans(ps=ps.trans,group = "Group", Top = 100)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存LDA结果
addWorksheet(trans_biomarker_wb, "LDA_accuracy")
writeData(trans_biomarker_wb, "LDA_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "LDA_importance")
writeData(trans_biomarker_wb, "LDA_importance", importance, rowNames = TRUE)

#43 naivebayes.trans------
res = naivebayes.trans(ps=pst, top = 200, seed = 100, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
addWorksheet(trans_biomarker_wb, "naivebayes_accuracy")
writeData(trans_biomarker_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "naivebayes_importance")
writeData(trans_biomarker_wb, "naivebayes_importance", importance, rowNames = TRUE)

#44 nnet.trans 神经网络------
res = nnet.trans  (ps=ps.trans, top = 200, seed = 10, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
addWorksheet(trans_biomarker_wb, "nnet_accuracy")
writeData(trans_biomarker_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
addWorksheet(trans_biomarker_wb, "nnet_importance")
writeData(trans_biomarker_wb, "nnet_importance", importance, rowNames = TRUE)

# 保存生物标志物总表
saveWorkbook(trans_biomarker_wb, file.path(trans_biomarkerpath, "transcriptome_biomarker_analysis.xlsx"), overwrite = TRUE)

# enrich analysis --------
# 创建转录组富集分析目录
trans_enrichpath <- file.path(trans_preprocesspath, "transcriptome", "enrichment")
dir.create(trans_enrichpath, recursive = TRUE)

# 创建转录组富集分析总表
trans_enrich_wb <- createWorkbook()

#45 KEGG_enrich.trans-------
library(clusterProfiler)
res = KEGG_enrich.trans(ps = ps.trans %>% filter_OTU_ps(1000))
res$plots[[1]]$`KO-OE`

# 保存KEGG富集结果
if (!is.null(res$plots[[1]]$`KO-OE`)) {
  ggsave(file.path(trans_enrichpath, "KEGG_enrich_KO_OE.png"), plot = res$plots[[1]]$`KO-OE`, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_enrichpath, "KEGG_enrich_KO_OE.pdf"), plot = res$plots[[1]]$`KO-OE`, width = 10, height = 8)
}
addWorksheet(trans_enrich_wb, "KEGG_enrich_results")
writeData(trans_enrich_wb, "KEGG_enrich_results", res$data, rowNames = TRUE)

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

# 保存GSVA结果
if (!is.null(res$plots[[1]][[1]])) {
  ggsave(file.path(trans_enrichpath, "gsva_heatmap.png"), plot = res$plots[[1]][[1]], width = 12, height = 10, dpi = 300)
  ggsave(file.path(trans_enrichpath, "gsva_heatmap.pdf"), plot = res$plots[[1]][[1]], width = 12, height = 10)
}
if (!is.null(res$plots[[2]]$`KO OE`)) {
  ggsave(file.path(trans_enrichpath, "gsva_KO_OE.png"), plot = res$plots[[2]]$`KO OE`, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_enrichpath, "gsva_KO_OE.pdf"), plot = res$plots[[2]]$`KO OE`, width = 10, height = 8)
}
addWorksheet(trans_enrich_wb, "gsva_results")
writeData(trans_enrich_wb, "gsva_results", res$data, rowNames = TRUE)

#47 gsea.trans  缺-----

#48 pathway_enrich.trans：通路富集-------
detach(package:mia, unload = TRUE)
detach(package:ANCOMBC, unload = TRUE)
phyloseq::tax_table(ps.trans3) %>% colnames()
ps.trans3 = ps.trans %>% tax_glom_wt("KO_id")

res2 = pathway_enrich.trans(ps = ps.trans3,
                            dif.method = "wilcox")
res2$plots$OE.WT.plot
res2$plots$KO.OE.plot

# 保存通路富集结果
if (!is.null(res2$plots$OE.WT.plot)) {
  ggsave(file.path(trans_enrichpath, "pathway_enrich_OE_WT.png"), plot = res2$plots$OE.WT.plot, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_enrichpath, "pathway_enrich_OE_WT.pdf"), plot = res2$plots$OE.WT.plot, width = 10, height = 8)
}
if (!is.null(res2$plots$KO.OE.plot)) {
  ggsave(file.path(trans_enrichpath, "pathway_enrich_KO_OE.png"), plot = res2$plots$KO.OE.plot, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_enrichpath, "pathway_enrich_KO_OE.pdf"), plot = res2$plots$KO.OE.plot, width = 10, height = 8)
}
addWorksheet(trans_enrich_wb, "pathway_enrich_results")
writeData(trans_enrich_wb, "pathway_enrich_results", res2$data, rowNames = TRUE)

#49 reaction.show.trans：反应展示-------
res3 = reaction.show.trans(ps= ps.trans,dif.method = "wilcox")
res3$plots$OE.WT.plot
res3$plotdata$OE.WT

# 保存反应展示结果
if (!is.null(res3$plots$OE.WT.plot)) {
  ggsave(file.path(trans_enrichpath, "reaction_show_OE_WT.png"), plot = res3$plots$OE.WT.plot, width = 10, height = 8, dpi = 300)
  ggsave(file.path(trans_enrichpath, "reaction_show_OE_WT.pdf"), plot = res3$plots$OE.WT.plot, width = 10, height = 8)
}
addWorksheet(trans_enrich_wb, "reaction_show_results")
writeData(trans_enrich_wb, "reaction_show_results", res3$plotdata$OE.WT, rowNames = TRUE)

#50 enrich.module.trans  缺-----

# 保存富集分析总表
saveWorkbook(trans_enrich_wb, file.path(trans_enrichpath, "transcriptome_enrichment_analysis.xlsx"), overwrite = TRUE)

#network analysis------
# 创建转录组网络分析目录
trans_networkpath <- file.path(trans_preprocesspath, "transcriptome", "network")
dir.create(trans_networkpath, recursive = TRUE)

# 创建转录组网络分析总表
trans_network_wb <- createWorkbook()

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

# 保存网络分析主图
ggsave(file.path(trans_networkpath, "network_main.png"), plot = p0, width = 12, height = 10, dpi = 300)
ggsave(file.path(trans_networkpath, "network_main.pdf"), plot = p0, width = 12, height = 10)

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

# 保存网络属性数据
addWorksheet(trans_network_wb, "network_properties")
writeData(trans_network_wb, "network_properties", dat2, rowNames = TRUE)

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

# 保存样本网络属性数据
addWorksheet(trans_network_wb, "sample_network_properties")
writeData(trans_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

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

# 保存节点属性数据
addWorksheet(trans_network_wb, "node_properties")
writeData(trans_network_wb, "node_properties", nodepro2, rowNames = TRUE)

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

# 保存网络比较结果
addWorksheet(trans_network_wb, "network_comparison")
writeData(trans_network_wb, "network_comparison", res, rowNames = TRUE)

#57 CNPS.network: CNPS 基因功能网络 -----
library(data.table)
library(EasyStat)
library(tidyfst)
library(sna)
dat = db.cnps
head(dat)
res = CNPS.network(ps = ps.trans,dat = dat,id.0 = "C")
p= res[[1]]
p

dat1= res[[2]]
dat1

dat2= res[[3]]
dat2

# 保存CNPS网络结果
ggsave(file.path(trans_networkpath, "CNPS_network.png"), plot = p, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_networkpath, "CNPS_network.pdf"), plot = p, width = 10, height = 8)
addWorksheet(trans_network_wb, "CNPS_network_data1")
writeData(trans_network_wb, "CNPS_network_data1", dat1, rowNames = TRUE)
addWorksheet(trans_network_wb, "CNPS_network_data2")
writeData(trans_network_wb, "CNPS_network_data2", dat2, rowNames = TRUE)

# wgcna.hub.trans:共表达基因模块挖掘-----
library(WGCNA)
library(dplyr)
#不输出
res <- wgcna.hub.trans(ps=ps.trans, data_type="TPM",group="Group",WGCNApath=NULL,
                       alltraits=NULL,sample_filter=TRUE,clustsample_thresold=50000,
                       minModuleSize = 30,MEDissThres = 0.25,topnSelect_module=5,
                       nSelect_genes = 400)
head(res[[1]])
head(res[[2]])

# 保存WGCNA结果
addWorksheet(trans_network_wb, "wgcna_results1")
writeData(trans_network_wb, "wgcna_results1", res[[1]], rowNames = TRUE)
addWorksheet(trans_network_wb, "wgcna_results2")
writeData(trans_network_wb, "wgcna_results2", res[[2]], rowNames = TRUE)

#输出至本地
res <- wgcna.hub.trans(ps=ps.trans, data_type="TPM",group="Group",WGCNApath="./wgcna_result/",
                       alltraits=NULL,sample_filter=TRUE,clustsample_thresold=50000,
                       minModuleSize = 30,MEDissThres = 0.25,topnSelect_module=5,
                       nSelect_genes = 400)
head(res[[1]])
head(res[[2]])

# hub.gene.trans:基因模块挖掘-----
library(TCseq)
library(Mfuzz)
library(monocle)
library(Biobase)
library(clusterProfiler)  # 用于富集分析
library(DOSE)
library(AnnotationDbi)
library(ClusterGVis)
res = hub.gene.trans(ps.trans, id = "KO_id", cluster.method = "mfuzz", cluster.num = 6)

p1= res[[1]]
p1

p2= res[[2]]
p2
dat= res[[3]]

# 保存基因模块挖掘结果
ggsave(file.path(trans_networkpath, "hub_gene_plot1.png"), plot = p1, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_networkpath, "hub_gene_plot1.pdf"), plot = p1, width = 10, height = 8)
ggsave(file.path(trans_networkpath, "hub_gene_plot2.png"), plot = p2, width = 10, height = 8, dpi = 300)
ggsave(file.path(trans_networkpath, "hub_gene_plot2.pdf"), plot = p2, width = 10, height = 8)
addWorksheet(trans_network_wb, "hub_gene_results")
writeData(trans_network_wb, "hub_gene_results", dat, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(trans_network_wb, file.path(trans_networkpath, "transcriptome_network_analysis.xlsx"), overwrite = TRUE)


