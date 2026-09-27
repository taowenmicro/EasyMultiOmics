# ===================== 宏基因组物种分析（基于 EasyMultiOmics） =====================
rm(list = ls())

## ===================== 0. 基础设置 & 加载包 =====================
library(tidyverse)
library(data.table)
library(ggClusterNet)
library(tidyfst)
library(fs)
library(EasyMultiOmics)
library(openxlsx)


## ===================== 1. 读入 phyloseq 对象 & 基础信息
# 读入 phyloseq 对象
# 方式一：交互选择
# file_path <- tcltk::tk_choose.files(
#   caption = "请选择 ps_ITS.rds 文件",
#   multi   = FALSE,
#   filters = matrix(c("RDS files", ".rds",
#                      "All files", "*"), ncol = 2, byrow = TRUE)
# )
# 方式二：直接指定路径
# file_path <- "./data/ps_16s.rds"
# ps.kegg   <- readRDS(file_path)
# ps.kegg
# 读入 phyloseq 对象
ps.kegg <- EasyMultiOmics::ps.kegg
# 检查文件#
res <- infer_sequencing_type(ps = ps.kegg)  # ps 是你的 phyloseq 对象
res$label
sample_sums(ps.kegg)
map <- sample_data(ps.kegg)
head(map)

# 创建文件夹#
mg_path <- create_omics_result_dir_auto(
  ps = ps.kegg,
  base_dir     = "./result",
  include_time = FALSE
)

mg_path
# 可选——加入func与物种区分
mg_path <- gsub("metagenome", "metagenome_func", mg_path)
mg_path

# 简单检查
map <- sample_data(ps.kegg)
head(map)

sample_sums(ps.kegg)
phyloseq::tax_table(ps.kegg) %>% head()

# 分组数量与顺序
gnum       <- phyloseq::sample_data(ps.kegg)$Group %>% unique() %>% length()
gnum
axis_order <- phyloseq::sample_data(ps.kegg)$Group %>% unique()
axis_order

# 分组配色
col.g <- c("KO" = "#D55E00", "WT" = "#0072B2", "OE" = "#009E73")
col.g
col.g <- get_group_cols(axis_order, palette = "npg")
scales::show_col(col.g)

# 主题 & 颜色
package.amp()
res      <- theme_my(ps.kegg)
mytheme1 <- res[[1]]
mytheme2 <- res[[2]]
colset1  <- res[[3]]
colset2  <- res[[4]]
colset3  <- res[[5]]
colset4  <- res[[6]]

# KEGG 数据库路径
keggpath <- file.path(mg_path, "kegg")
dir.create(keggpath, recursive = TRUE, showWarnings = FALSE)

# 数据预处理
ps.metf <- ps.kegg %>% filter_OTU_ps(Top = 1000)
tax     <- ps.metf %>% phyloseq::tax_table()
colnames(tax)[3] <- "KOid"
tax_table(ps.metf) <- tax
head(tax)



## ===================== 2. Alpha 多样性分析 =====================

# 创建alpha多样性分析目录
mg_alpha_path <- file.path(mg_path, "01_alpha_diversity")
dir.create(mg_alpha_path, recursive = TRUE, showWarnings = FALSE)

alpha_xlsx_path <- file.path(mg_alpha_path, "alpha_diversity_results.xlsx")

if (file.exists(alpha_xlsx_path)) {
  mg_alpha_wb <- openxlsx::loadWorkbook(alpha_xlsx_path)
} else {
  mg_alpha_wb <- openxlsx::createWorkbook()
}



### 2.1 alpha.metf:6种alpha多样性计算 ----
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

# 盒线图
p1_1 = result1[[1]]+
  # scale_x_discrete(limits = axis_order) +
  # mytheme1 +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic() +
  #  ggplot2::guides(fill = guide_legend(title = none))+
  theme(axis.title.y = element_text(angle = 90))


p1_1

# 柱状图
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

# 盒线 + 柱状
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


# 小提琴图
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
save_plot2(p1_1, mg_alpha_path, "alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_1, mg_alpha_path, "alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_2, mg_alpha_path, "alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_2, mg_alpha_path, "alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_3, mg_alpha_path, "alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_3, mg_alpha_path, "alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_0, mg_alpha_path, "alpha_diversity_violin", width = 14, height = 8)
save_plot2(p1_0, mg_alpha_path, "alpha_diversity_violin", width = 14, height = 8)

# 保存alpha多样性表格
write_sheet2(mg_alpha_wb, "alpha_diversity_data", data)
write_sheet2(mg_alpha_wb, "alpha_diversity_stat",result$sta)
openxlsx::saveWorkbook(mg_alpha_wb, alpha_xlsx_path, overwrite = TRUE)


## ===================== 3. Beta 多样性分析 =====================

# 创建beta多样性分析目录
mg_beta_path <- file.path(mg_path, "02_beta_diversity")
dir.create(mg_beta_path, recursive = TRUE, showWarnings = FALSE)

beta_xlsx_path <- file.path(mg_beta_path, "beta_diversity_results.xlsx")

if (file.exists(beta_xlsx_path)) {
  mg_beta_wb <- openxlsx::loadWorkbook(beta_xlsx_path)
} else {
  mg_beta_wb <- openxlsx::createWorkbook()
}

### ---------- 3.1 ordinate.metf：PCoA 排序 ---------
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
save_plot2(p3_1, mg_beta_path, "pca_basic", width = 10, height = 8)
save_plot2(p3_1, mg_beta_path, "pca_basic", width = 10, height = 8)
save_plot2(p3_2, mg_beta_path, "pca_labeled", width = 10, height = 8)
save_plot2(p3_2, mg_beta_path, "pca_labeled", width = 10, height = 8)
save_plot2(p3_3, mg_beta_path, "pca_enhanced", width = 10, height = 8)
save_plot2(p3_3, mg_beta_path, "pca_enhanced", width = 10, height = 8)

# 保存排序数据
write_sheet2(mg_beta_wb, "ordination_data", plotdata)
write_sheet2(mg_beta_wb, "ordination_cent", cent)
write_sheet2(mg_beta_wb, "ordination_segs", segs)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)

### 3.2 MetaTest.metf:群落功能差异检测adonis ----
dat1 = MetaTest.metf(ps = ps.metf, Micromet = "adonis", dist = "bray")
dat1
write_sheet2(mg_beta_wb, "adonis_results", dat1)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)

### 3.3 pairMetaTest.metf:两两分组群落功能差异检测MRPP ----
dat2 = pairMetaTest.metf(ps = ps.metf, Micromet = "MRPP", dist = "bray")
dat2
write_sheet2(mg_beta_wb, "pairwise_MRPP", dat2)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)


### 3.4 mantal.metf：群落功能差异检测普鲁士分析#------
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
save_plot2(p3_7, mg_beta_path, "mantel_test", width = 10, height = 8)
save_plot2(p3_7, mg_beta_path, "mantel_test", width = 10, height = 8)

# mantal分析
# addWorksheet(beta_wb, "mantel_test")  # 已合并到 write_sheet2
write_sheet2(mg_beta_wb, "mantel_results", data)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)


## ===================== 4. 组成（Composition）分析 =====================

mg_comp_path <- file.path(mg_path, "03_composition")
dir.create(mg_comp_path, recursive = TRUE, showWarnings = FALSE)

comp_xlsx_path <- file.path(mg_comp_path, "composition_results.xlsx")

if (file.exists(comp_xlsx_path)) {
  mg_comp_wb <- openxlsx::loadWorkbook(comp_xlsx_path)
} else {
  mg_comp_wb <- openxlsx::createWorkbook()
}


### 4.1 cluster_metf:样品聚类 ----


res <- cluster_metf(
  ps = ps.metf,
  hcluter_method = "complete",
  dist = "bray",
  cuttree = 3,
  row_cluster = TRUE,
  col_cluster = TRUE,
  col.g = col.g
)
p4 = res[[1]]
p4
p4_1 = res[[2]]
p4_1
p4_2 = res[[3]]
p4_2
dat = res[[4]]
dat

#保存图片
save_plot2(p4, mg_comp_path, "cluster_heatmap", width = 12, height = 10)
save_plot2(p4, mg_comp_path, "cluster_heatmap", width = 12, height = 10)
save_plot2(p4_1, mg_comp_path, "cluster_dendrogram", width = 10, height = 6)
save_plot2(p4_1, mg_comp_path, "cluster_dendrogram", width = 10, height = 6)
save_plot2(p4_2, mg_comp_path, "cluster_groups", width = 8, height = 6)
save_plot2(p4_2, mg_comp_path, "cluster_groups", width = 8, height = 6)

write_sheet2(mg_comp_wb, "cluster_data", dat)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### 4.2 Micro_tern.metf: 三元图展示功能----
res = Micro_tern.metf(ps= ps.metf %>% filter_OTU_ps(500),color = "Pathway"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图图片

save_plot2(p15 + theme_classic(), mg_comp_path, "ternary_plot", width = 12, height = 10)
write_sheet2(mg_comp_wb, "ternary_data", dat)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### 4.3 barMainplot.metf: 堆积柱状图展示功能组成----
# 功能分组
rank_names(ps.metf)

result = barMainplot.metf(ps =ps.metf%>%filter_OTU_ps(Top = 500),
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

save_plot2(p4_1, mg_comp_path, "barplot_samples", width = 10, height = 8)
save_plot2(p4_2, mg_comp_path, "barplot_groups",  width = 10, height = 8)

write_sheet2(mg_comp_wb, "barplot_raw_data",  res_barMain[[2]])
write_sheet2(mg_comp_wb, "barplot_summary",   databar)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### 4.4 Ven.Upset.metf: 用于展示共有、特有的功能----
# 分组小于6时使用
res = Ven.Upset.metf(ps = ps.metf,
                     group = "Group",
                     N = 0.5,
                     size = 3,
                     fill_color = col.g)
p10.1 = res[[1]]
p10.1

p10.2 = res[[2]]
p10.2

dat = res[[3]]
dat


save_plot2(p10.1, mg_comp_path, "venn_diagram",  width = 10, height = 8)
save_plot2(p10.2, mg_comp_path, "upset_plot",    width = 12, height = 8)

write_sheet2(mg_comp_wb, "venn_upset_data", dat)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### 4.5 ggflower.metf:花瓣图展示共有特有功能------
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


save_plot2(p13.1, mg_comp_path, "flower_plot", width = 10, height = 8)
write_sheet2(mg_comp_wb, "flower_data", dat)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### 4.6  ven.network.metm: ven网络展示共有特有功能----
map =  sample_data(ps.metf)
ps_ven=ps.metf
tax=data.frame(phyloseq:: tax_table(ps.metf))
colnames(tax)[[3]]="KOnumber"

tax_table(ps_ven) <- as.matrix(tax)


result =ven.network.metm(
  ps = ps_ven,
  N = 0.5,
  fill = "Pathway")

p14  = result[[1]] +
  theme(legend.position = "none")

p14
dat = result[[2]]
dat


save_plot2(p14, mg_comp_path, "venn_network", width = 12, height = 10)
write_sheet2(mg_comp_wb, "venn_network_data", dat)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

## ===================== 5. 差异分析（Differential Analysis） =====================

mg_diff_path <- file.path(mg_path, "04_differential")
dir.create(mg_diff_path, recursive = TRUE, showWarnings = FALSE)

diff_xlsx_path <- file.path(mg_diff_path, "differential_results.xlsx")

if (file.exists(diff_xlsx_path)) {
  mg_diff_wb <- openxlsx::loadWorkbook(diff_xlsx_path)
} else {
  mg_diff_wb <- openxlsx::createWorkbook()
}

### 5.1 DESep2Super.metf:DESep2计算差异功能基因 ----

res= DESep2Super.metf(
  ps       = ps.metf,
  group    = "Group",
  artGroup = NULL,
  j        = "meta",
  col.g    = col.g,
  gradient = TRUE,
  top_n    = 5,
  select_by = "pvalue"
)

p15.1 = res[[1]][1]
p15.1
p15.2 = res[[1]][2]
p15.2
p15.3 = res[[1]][3]
p15.3
dat = res[[2]]
dat

# 保存DESeq2图片
save_plot2(p15.1, mg_diff_path, "deseq2_volcano_1", width = 10, height = 8)
save_plot2(p15.2, mg_diff_path, "deseq2_volcano_2", width = 10, height = 8)
save_plot2(p15.3, mg_diff_path, "deseq2_volcano_3", width = 10, height = 8)

write_sheet2(mg_diff_wb, "deseq2_results", dat)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### 5.2 EdgerSuper.metf:EdgeR计算差异功能基因----

res = EdgerSuper.metf(
  ps       = ps.metf,
  group    = "Group",
  artGroup = NULL,
  j        = "meta",
  col.g    = col.g,
  gradient = TRUE,
  top_n    = 5,
  select_by = "pvalue"  # 按 -log10(p) 选择
  # select_by = "logFC"  # 按 |logFC| 选择
)

p16.1 = res[[1]][1]
p16.1
p16.2 = res[[1]][2]
p16.2
p16.3 = res[[1]][3]
p16.3
dat = res[[2]]
dat

save_plot2(p16.1, mg_diff_path, "edger_volcano_1", width = 10, height = 8)
save_plot2(p16.2, mg_diff_path, "edger_volcano_2", width = 10, height = 8)
save_plot2(p16.3, mg_diff_path, "edger_volcano_3", width = 10, height = 8)

write_sheet2(mg_diff_wb, "edger_results", dat)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### 5.3 t_metf: 差异分析t检验----
res = t_metf(ps = ps.metf,group  = "Group",artGroup =NULL)
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
save_plot2(p17, diffpath, "ttest_main", width = 10, height = 8)
save_plot2(p17.1, diffpath, "ttest_volcano_1", width = 10, height = 8)
save_plot2(p17.2, diffpath, "ttest_volcano_2", width = 10, height = 8)
save_plot2(p17.2, diffpath, "ttest_volcano_2", width = 10, height = 8)
save_plot2(p17.3, diffpath, "ttest_volcano_3", width = 10, height = 8)

# 保存t检验数据到总表
write_sheet2(mg_diff_wb, "ttest_results",data)
write_sheet2(mg_diff_wb, "ttest_intersections", inter_union)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### 5.4 wlx.metf:非参数检验——-----
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
save_plot2(p18, diffpath, "wilcox_main", width = 10, height = 8)
save_plot2(p18.1, diffpath, "wilcox_volcano_1", width = 10, height = 8)
save_plot2(p18.2, diffpath, "wilcox_volcano_2", width = 10, height = 8)
save_plot2(p18.3, diffpath, "wilcox_volcano_3", width = 10, height = 8)

# 保存wlx非参数检验数据到总表
write_sheet2(mg_diff_wb,"wilcox_results", data)
write_sheet2(mg_diff_wb, "wilcox_intersections", inter_union)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### 5.5 stamp.metf:stamp差异分析 ----

map <- sample_data(ps.metf)
allgroup <- combn(unique(map$Group), 2)
plot_list_28 <- vector("list", ncol(allgroup))

for (i in seq_len(ncol(allgroup))) {
  ps_sub <- phyloseq::subset_samples(ps.metf, Group %in% allgroup[, i])

  sub_groups <- allgroup[, i]
  sub_col.g <- col.g[sub_groups]

  p_tmp <- stemp_diff.metm(
    ps = ps_sub,
    Top = 20,
    ranks = 6,
    col.g = sub_col.g
  )

  plot_list_28[[i]] <- p_tmp

  save_plot2(
    p_tmp,
    mg_diff_path,
    paste0("stamp_plot_", i),
    width  = 10,
    height = 8
  )
}

## ===================== 6. 生物标志物分析（Biomarker Identification） =====================
mg_biomarker_path <- file.path(mg_path, "05_biomarker")
dir.create(mg_biomarker_path, recursive = TRUE, showWarnings = FALSE)

biomarker_xlsx_path <- file.path(mg_biomarker_path, "biomarker_results.xlsx")

if (file.exists(biomarker_xlsx_path)) {
  mg_biomarker_wb <- openxlsx::loadWorkbook(biomarker_xlsx_path)
} else {
  mg_biomarker_wb <- openxlsx::createWorkbook()
}

# 选择一个二元分组子集 pst 作为 biomarker 输入
id = sample_data(ps.metf)$Group %>% unique()
aaa = combn(id,2)
i= 1
group = c(aaa[1,i],aaa[2,i])

pst = ps.metf %>% subset_samples.wt("Group",group) %>%
  filter_taxa(function(x) sum(x ) > 10, TRUE)

### 6.1  rfcv.metf :交叉验证结果-------
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
save_plot2(prfcv + theme_classic(), mg_biomarker_path, "roc_curve", width = 8, height = 6)

write_sheet2(mg_biomarker_wb, "rfcv_plot_data", result[[2]])
write_sheet2(mg_biomarker_wb, "rfcv_table", rfcvtable)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.2 Roc.metf:ROC 曲线绘制----
res = Roc.metf( ps = pst,group  = "Group",repnum = 5)
p33.1 =  res[[1]]
p33.1+theme_classic()
p33.2 =  res[[2]]
p33.2
dat =  res[[3]]
dat

# 保存ROC曲线图片
save_plot2(p33.1 + theme_classic(), mg_biomarker_path, "roc_curve", width = 8, height = 6)
save_plot2(p33.2, mg_biomarker_path, "roc_boxplot", width = 10, height = 8)

write_sheet2(mg_biomarker_wb, "roc_results", dat)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.3 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.metf,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

save_plot2(p34.1 + theme_classic(), mg_biomarker_path, "pca_loading", width = 10, height = 8)
write_sheet2(mg_biomarker_wb, "pca_loading_data", dat)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.4 svm_metf:svm筛选特征功能 ----
res <- svm_metf(ps = ps.metf %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果数据到总表
write_sheet2(mg_biomarker_wb, "svm_auc", data.frame(AUC = AUC))
write_sheet2(mg_biomarker_wb, "svm_importance", importance)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.5 glm.metf:glm筛选特征功能----
pst = subset_samples(ps.metf,Group %in% c("KO" ,"OE"))
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)

AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果数据到总表
write_sheet2(mg_biomarker_wb, "glm_auc", data.frame(AUC = AUC))
write_sheet2(mg_biomarker_wb, "glm_importance", importance)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.6 xgboost.metf: xgboost筛选特征功能----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)
res = xgboost.metf(ps =ps.metf, top = 20  )
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance
imp_xgb <- as.data.frame(res[[2]]$importance)
# 保存XGBoost结果数据到总表
write_sheet2(mg_biomarker_wb, "xgboost_accuracy", data.frame(Accuracy = accuracy))
write_sheet2(mg_biomarker_wb, "xgboost_importance", imp_xgb)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.7 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果数据到总表
write_sheet2(mg_biomarker_wb, "lasso_accuracy", data.frame(AUC = accuracy))
write_sheet2(mg_biomarker_wb, "lasso_importance", importance)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.8 decisiontree.micro: ----
library(dplyr)
library(caret)
library(rpart)
res = decisiontree.metf(ps=ps.metf, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果数据到总表
write_sheet2(mg_biomarker_wb, "tree_accuracy", data.frame(Accuracy = accuracy))
write_sheet2(mg_biomarker_wb, "tree_importance", importance)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.9 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps.metf, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果数据到总表

write_sheet2(mg_biomarker_wb, "naivebayes_accuracy", data.frame(Accuracy = accuracy))
write_sheet2(mg_biomarker_wb, "naivebayes_importance", importance)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.10 LDA.metf: LDA筛选特征功能-----
tablda = LDA.metf(
  ps = ps.metf,
  seed = 11,
  Top = 100,
  p.lvl = 0.05,
  lda.lvl = 1,
  adjust.p = F,
  col.g = col.g
)


tablda[[1]]

p35 <- lefse_bar(taxtree = tablda[[2]])+ theme_classic()+ ylab("")
p35
dat = tablda[[2]]
dat

# 保存LDA图片
save_plot2(p35, mg_biomarker_path, "lda_barplot", width = 10, height = 8)

write_sheet2(mg_biomarker_wb, "lda_results", dat)
write_sheet2(
  mg_biomarker_wb,
  "lda_parameters",
  data.frame(
    Parameter = c("Top", "p.lvl", "lda.lvl", "seed", "adjust.p"),
    Value     = c("10", "0.05", "4", "11", "FALSE")
  )
)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.11 randomforest.metf: 随机森林筛选特征功能----
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
save_plot2(p42.1 + theme_classic(), mg_biomarker_path, "rf_importance", width = 16, height = 11)
save_plot2(p42.2, mg_biomarker_path, "rf_error",      width = 10, height = 8)
save_plot2(p42.4, mg_biomarker_path, "rf_additional", width = 6, height = 4)

# 保存随机森林数据到总表

write_sheet2(mg_biomarker_wb, "rf_importance", dat)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.12 nnet.metf: 神经网络筛选特征功能  ------
res =nnet.metf(ps=ps.metf, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果数据到总表
write_sheet2(mg_biomarker_wb, "nnet_accuracy", data.frame(Accuracy = accuracy))
write_sheet2(mg_biomarker_wb, "nnet_importance", importance)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### 6.13 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging_metf(ps =  ps.metf, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果数据到总表
write_sheet2(mg_biomarker_wb, "bagging_accuracy", data.frame(Accuracy =  accuracy))
write_sheet2(mg_biomarker_wb, "bagging_importance",importance)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)


## ===================== 7. 网络分析（Network Analysis） =====================

## 如果上面你用的是 path，这里做一个映射
# mg_path <- path

mg_network_path <- file.path(mg_path, "06_network")
dir.create(mg_network_path, showWarnings = FALSE, recursive = TRUE)

network_xlsx_path <- file.path(mg_network_path, "network_results.xlsx")
if (file.exists(network_xlsx_path)) {
  mg_network_wb <- openxlsx::loadWorkbook(network_xlsx_path)
} else {
  mg_network_wb <- openxlsx::createWorkbook()
}

library(ggClusterNet)
library(igraph)

### 7.1 network.pip:网络分析主函数--------
library(ggClusterNet)
library(igraph)
detach("package:mia", unload = TRUE)


# 处理 tax_table，只保留第一个 Level1 值
ps.metf1 <- ps.metf
tax <- data.frame(phyloseq::tax_table(ps.metf1))
tax$Level1 <- sapply(strsplit(tax$Level1, ";"), function(x) x[1])
tax$Level2 <- sapply(strsplit(tax$Level2, ";"), function(x) x[1])
tax$Level3 <- sapply(strsplit(tax$Level3, ";"), function(x) x[1])
tax_table(ps.metf1) <- as.matrix(tax)

# 检查
head(phyloseq::tax_table(ps.metf1)[, "Level1"])

# 再运行函数
tab.r = network.pip.metf(
  ps = ps.metf1,
  N = 50,
  ra = 0.05,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  label = FALSE,
  lab = "elements",
  group = "Group",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R = 10,
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
save_plot2(p0, networkpath, "network_plot", width = 12, height = 12)
save_plot2(p0, networkpath, "network_plot", width = 12, height = 12)
## 主网络图（第一张）
save_plot2(p0, mg_network_path, "network_main", width = 14, height = 8)

## 相关性矩阵写表
if (!is.null(cortab)) {
  for (i in seq_along(cortab)) {
    group_name <- names(cortab)[i]
    cor_mat    <- cortab[[i]]
    if (is.matrix(cor_mat) || is.data.frame(cor_mat)) {
      cor_df <- as.data.frame(cor_mat)
      write_sheet2(
        wb       = mg_network_wb,
        sheet    = paste0("cor_", group_name),
        df     = cor_df,
        row_names = TRUE
      )
    }
  }
}

## 保存网络参数
network_info <- data.frame(
  Parameter = c("N_features", "r_threshold", "p_threshold",
                "maxnode", "layout", "cluster_method"),
  Value     = c("200", "0.6", "0.05", "2",
                "model_maptree2", "cluster_fast_greedy")
)
write_sheet2(mg_network_wb, "network_parameters", network_info)

openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)

### 7.2 net_properties.4:网络属性计算 ----
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
write_sheet2(
  wb       = mg_network_wb,
  sheet    = "network_properties_summary",
  df     = dat2,
  row_names = TRUE
)
openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


### 7.3 node_properties:计算节点属性 ----
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
write_sheet2(
  wb       = mg_network_wb,
  sheet    = "node_properties_combined",
  df    = nodepro2,
  row_names = FALSE
)
openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)

### 7.4 module.compare.net.pip:网络显著性比较 ----
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
write_sheet2(
  wb       = mg_network_wb,
  sheet    = "network_comparison_results",
  df    = res,
  row_names = FALSE
)
openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


## ===================== 8. 通路富集（pathway enrich） =====================
# 创建通路富集分析目录
mg_enrich_path <- file.path(mg_path, "07_enrich")
dir.create(mg_enrich_path, showWarnings = FALSE, recursive = TRUE)

enrich_xlsx_path <- file.path(mg_enrich_path, "enrich_results.xlsx")
if (file.exists(enrich_xlsx_path)) {
  mg_enrich_wb <- openxlsx::loadWorkbook(enrich_xlsx_path)
} else {
  mg_enrich_wb <- openxlsx::createWorkbook()
}
### 8.1 KEGG_enrich2:KEGG富集分析 ----
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
write_sheet2(mg_enrich_wb, "KO_vs_OE", dat1)
write_sheet2(mg_enrich_wb, "KO_vs_WT", dat2)
write_sheet2(mg_enrich_wb, "OE_vs_WT", dat3)
openxlsx::saveWorkbook(mg_enrich_wb, enrich_xlsx_path, overwrite = TRUE)

### 8.2  buplot.metf：富集分析气泡图-----
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


result = buplot.micro(dt = dat1, id = id, top_n = 10)
# 检查数据结构
str(dat1)

# 检查 GeneRatio 列
class(dat1$GeneRatio)
head(dat1$GeneRatio)


p1 = result[[1]]
p1
p2 = result[[2]]
p2

# 保存富集分析气泡图
save_plot2(p1, mg_enrich_path, "enrichment_bubble_plot1", width = 12, height = 10)
save_plot2(p2, mg_enrich_path, "enrichment_bubble_plot2", width = 12, height = 10)


### 8.3  GSVA.metf:GSVA富集分析--------
library(limma)
library(GSVA)


res2 <- GSVA_metf(
  ps = ps.metf,
  dif = dat,
  lg2FC = 0.5,
  padj = 0.05,
  col.g = col.g,
  heatmap_color = colorRampPalette(c("#00CED1", "#FFFFFF", "#FF4500"))(100)  # 使用自定义颜色
)
# 3. 获取结果
# GSVA 热图
dat1 <- res2[[1]][[1]]  # GSVA 分数矩阵
head(dat1)
p39.1 <- res2[[1]][[2]]  # 热图
p39.1

# 差异分析结果表
dat3 <- res2[[2]]$`KO-WT`
# dat3 <- res2[[2]]$`KO-OE`
# dat3 <- res2[[2]]$`OE-WT`

# 火山图
p39.2 <- res2[[3]]$`KO-WT`
# p39.2 <- res2[[3]]$`OE-WT`
# p39.2 <- res2[[3]]$`KO-OE`

# 保存GSVA图片
save_plot2(p39.1,  mg_enrich_path, "gsva_heatmap", width = 12, height = 10)
save_plot2(p39.2,  mg_enrich_path, "gsva_comparison1", width = 10, height = 8)

### 8.4 kegg_function:按照kegg通路合并基因-------
detach(package:mia, unload = TRUE)
detach(package:ANCOMBC, unload = TRUE)
ps.kegg.function = kegg_function(ps = ps.metf)

# 提取OTU表和税分类表
otu_kegg_func <- as.data.frame(phyloseq::otu_table(ps.kegg.function))
tax_kegg_func <- as.data.frame(phyloseq::tax_table(ps.kegg.function))

# 保存KEGG通路合并基因数据到总表
write_sheet2(mg_enrich_wb, "kegg_function_data", otu_kegg_func)
write_sheet2(mg_enrich_wb, "kegg_function_taxonomy", tax_kegg_func)
openxlsx::saveWorkbook(mg_enrich_wb, enrich_xlsx_path, overwrite = TRUE)

### 8.5 基于Mkegg通路合并基因-------
ps.mkegg.function =mkegg_function(ps = ps.metf)

# 提取OTU表和税分类表
otu_mkegg_func <- as.data.frame(phyloseq::otu_table(ps.mkegg.function))
tax_mkegg_func <- as.data.frame(phyloseq::tax_table(ps.mkegg.function))

# 保存Mkegg通路合并基因数据到总表
write_sheet2(mg_enrich_wb, "mkegg_function_data", otu_mkegg_func)
write_sheet2(mg_enrich_wb, "mkegg_function_taxonomy", tax_mkegg_func)
openxlsx::saveWorkbook(mg_enrich_wb, enrich_xlsx_path, overwrite = TRUE)

### 8.6 基于reaction合并基因-------
ps_kegg_function.reaction = reaction_function(ps = ps.metf)

# 提取OTU表和税分类表
otu_reaction_func <- as.data.frame(phyloseq::otu_table(ps_kegg_function.reaction))
tax_reaction_func <- as.data.frame(phyloseq::tax_table(ps_kegg_function.reaction))

# 保存reaction合并基因数据到总表
write_sheet2(mg_enrich_wb, "reaction_function_data", otu_reaction_func)
write_sheet2(mg_enrich_wb, "reaction_function_taxonomy", tax_reaction_func)
openxlsx::saveWorkbook(mg_enrich_wb, enrich_xlsx_path, overwrite = TRUE)



# =====================  CNPS循环数据库 =====================
# 创建CNPS循环分析目录
mg_cnps_path <- file.path(mg_path, "08_cnps")
dir.create(mg_cnps_path, showWarnings = FALSE, recursive = TRUE)

cnps_xlsx_path <- file.path(mg_cnps_path, "cnps_results.xlsx")
if (file.exists(cnps_xlsx_path)) {
  mg_cnps_wb <- openxlsx::loadWorkbook(cnps_xlsx_path)
} else {
  mg_cnps_wb <- openxlsx::createWorkbook()
}
## ===================== 1. 数据准备==============

### 1.1  cnps_gene: KEGG数据库中过滤CNPS 元素循环相关基因------
res  = cnps_gene(ps=ps.metf)
C_gene <- res$C_gene
N_gene <- res$N_gene
P_gene <- res$P_gene
S_gene <- res$S_gene

# 保存CNPS基因数据到总表
write_sheet2(mg_cnps_wb, "C_genes", C_gene)
write_sheet2(mg_cnps_wb, "N_genes", N_gene)
write_sheet2(mg_cnps_wb, "P_genes", P_gene)
write_sheet2(mg_cnps_wb, "S_genes", S_gene)
openxlsx::saveWorkbook(mg_cnps_wb, cnps_xlsx_path, overwrite = TRUE)

### 1.2 cnps_ps: KEGG数据库中过滤CNPS 元素循环相关基因构建ps -----
ps.cnps = cnps_ps(psko = ps.metf)
ps.cnps
phyloseq::tax_table(ps.cnps) %>% colnames()
otu_table(ps.cnps)

# 保存CNPS phyloseq对象数据到总表
otu_cnps <- as.data.frame(phyloseq::otu_table(ps.cnps))
tax_cnps <- as.data.frame(phyloseq::tax_table(ps.cnps))
sample_cnps <- as.data.frame(phyloseq::sample_data(ps.cnps))

write_sheet2(mg_cnps_wb, "cnps_otu_table", otu_cnps)
write_sheet2(mg_cnps_wb, "cnps_taxonomy", tax_cnps)
write_sheet2(mg_cnps_wb, "cnps_sample_data", sample_cnps)
openxlsx::saveWorkbook(mg_cnps_wb, cnps_xlsx_path, overwrite = TRUE)


## ===================== 2. Alpha 多样性分析 =====================

# 创建CNPS alpha多样性目录
cnps_alphapath <- file.path(cnpspath, "alpha")
dir.create(cnps_alphapath, recursive = TRUE)

# 创建CNPS alpha多样性总表
cnps_alpha_wb <- createWorkbook()

# 3 alpha.metf:6种alpha多样性计算 ----
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
save_plot2(p1_1, cnps_alphapath, "cnps_alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_1, cnps_alphapath, "cnps_alpha_diversity_box", width = 12, height = 8)

res = FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2 = p1_2+
  scale_x_discrete(limits = axis_order) +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

save_plot2(p1_2, cnps_alphapath, "cnps_alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_2, cnps_alphapath, "cnps_alpha_diversity_bar", width = 12, height = 8)

res = FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3 = p1_3+
  scale_x_discrete(limits = axis_order) +
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+ylab("Gene diversity") +
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

save_plot2(p1_3, cnps_alphapath, "cnps_alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_3, cnps_alphapath, "cnps_alpha_diversity_boxbar", width = 12, height = 8)

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

save_plot2(p1_0, cnps_alphapath, "cnps_alpha_diversity_violin", width = 14, height = 8)
save_plot2(p1_0, cnps_alphapath, "cnps_alpha_diversity_violin", width = 14, height = 8)

# 保存CNPS alpha多样性数据到总表
# addWorksheet(cnps_alpha_wb, "cnps_alpha_data")  # 已合并到 write_sheet2
# addWorksheet(cnps_alpha_wb, "cnps_alpha_stats")  # 已合并到 write_sheet2
writeData(cnps_alpha_wb, "cnps_alpha_data", data, rowNames = TRUE)
writeData(cnps_alpha_wb, "cnps_alpha_stats", result$sta, rowNames = TRUE)

#4 alpha_rare.metf: 稀释曲线绘制----
rare <- mean(sample_sums(ps.cnps))/10
result = alpha_rare.metf(ps = ps.cnps, group = "Group", method = "Richness", start = 100, step = rare)

p2_1 <- result[[1]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)

save_plot2(p2_1, cnps_alphapath, "cnps_rarefaction_individual", width = 10, height = 8)
save_plot2(p2_1, cnps_alphapath, "cnps_rarefaction_individual", width = 10, height = 8)

raretab <- result[[2]]
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)

save_plot2(p2_2, cnps_alphapath, "cnps_rarefaction_group", width = 10, height = 8)
save_plot2(p2_2, cnps_alphapath, "cnps_rarefaction_group", width = 10, height = 8)

# 保存稀释曲线数据到总表
# addWorksheet(cnps_alpha_wb, "rarefaction_data")  # 已合并到 write_sheet2
writeData(cnps_alpha_wb, "rarefaction_data", raretab, rowNames = TRUE)

# 保存CNPS alpha多样性总表
saveWorkbook(cnps_alpha_wb, file.path(cnps_alphapath, "cnps_alpha_diversity.xlsx"), overwrite = TRUE)

# 创建CNPS beta多样性目录
cnps_betapath <- file.path(cnpspath, "beta")
dir.create(cnps_betapath, recursive = TRUE)

# 创建CNPS beta多样性总表
cnps_beta_wb <- createWorkbook()

# 5 ordinate.metf:排序分析 ----
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

save_plot2(p3_1, cnps_betapath, "cnps_nmds_basic", width = 10, height = 8)
save_plot2(p3_1, cnps_betapath, "cnps_nmds_basic", width = 10, height = 8)

#带标签图形出图
p3_2 = result[[3]]
p3_2 = p3_2+
  guides(fill ="none") +
  scale_fill_manual(values = col.g)+
  theme_classic()+
  theme(axis.title.y = element_text(angle = 90))

save_plot2(p3_2, cnps_betapath, "cnps_nmds_labeled", width = 10, height = 8)
save_plot2(p3_2, cnps_betapath, "cnps_nmds_labeled", width = 10, height = 8)

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

save_plot2(p3_3, cnps_betapath, "cnps_nmds_enhanced", width = 10, height = 8)
save_plot2(p3_3, cnps_betapath, "cnps_nmds_enhanced", width = 10, height = 8)

# 保存CNPS beta多样性数据到总表
# addWorksheet(cnps_beta_wb, "nmds_plotdata")  # 已合并到 write_sheet2
# addWorksheet(cnps_beta_wb, "nmds_centroids")  # 已合并到 write_sheet2
# addWorksheet(cnps_beta_wb, "nmds_segments")  # 已合并到 write_sheet2
writeData(cnps_beta_wb, "nmds_plotdata", plotdata, rowNames = TRUE)
writeData(cnps_beta_wb, "nmds_centroids", cent, rowNames = TRUE)
writeData(cnps_beta_wb, "nmds_segments", segs, rowNames = TRUE)

# 6 MetaTest.metf:群落功能差异检测 ----
dat1 = MetaTest.metf(ps = ps.cnps, Micromet = "adonis", dist = "bray")
dat1

# 7 pairMetaTest.metf:两两分组群落功能差异检测 ----
dat2 = pairMetaTest.metf(ps =  ps.cnps, Micromet = "MRPP", dist = "bray")
dat2

# 保存CNPS群落差异检测结果到总表
# addWorksheet(cnps_beta_wb, "adonis_test")  # 已合并到 write_sheet2
# addWorksheet(cnps_beta_wb, "pairwise_mrpp")  # 已合并到 write_sheet2
writeData(cnps_beta_wb, "adonis_test", dat1, rowNames = TRUE)
writeData(cnps_beta_wb, "pairwise_mrpp", dat2, rowNames = TRUE)

# 8 mantal.metf：群落功能差异检测普鲁士分析#------
result <- mantal.metf(ps =  ps.cnps,
                      method =  "spearman",
                      group = "Group",
                      ncol = gnum,
                      nrow = 1)

data <- result[[1]]
data
p3_7 <- result[[2]]
p3_7

save_plot2(p3_7, cnps_betapath, "cnps_mantel_test", width = 10, height = 8)
save_plot2(p3_7, cnps_betapath, "cnps_mantel_test", width = 10, height = 8)

# 保存Mantel检验结果到总表
# addWorksheet(cnps_beta_wb, "mantel_test")  # 已合并到 write_sheet2
writeData(cnps_beta_wb, "mantel_test", data, rowNames = TRUE)

# 保存CNPS beta多样性总表
saveWorkbook(cnps_beta_wb, file.path(cnps_betapath, "cnps_beta_diversity.xlsx"), overwrite = TRUE)

# 创建CNPS聚类分析目录
cnps_clusterpath <- file.path(cnpspath, "cluster")
dir.create(cnps_clusterpath, recursive = TRUE)

# 创建CNPS聚类总表
cnps_cluster_wb <- createWorkbook()

#9 cluster.metf:样品聚类 ----
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
save_plot2(p4, cnps_clusterpath, "cnps_cluster_heatmap", width = 12, height = 10)
save_plot2(p4, cnps_clusterpath, "cnps_cluster_heatmap", width = 12, height = 10)
save_plot2(p4_1, cnps_clusterpath, "cnps_cluster_dendrogram", width = 10, height = 6)
save_plot2(p4_1, cnps_clusterpath, "cnps_cluster_dendrogram", width = 10, height = 6)
save_plot2(p4_2, cnps_clusterpath, "cnps_cluster_groups", width = 8, height = 6)
save_plot2(p4_2, cnps_clusterpath, "cnps_cluster_groups", width = 8, height = 6)

# 保存CNPS聚类数据到总表
# addWorksheet(cnps_cluster_wb, "cluster_data")  # 已合并到 write_sheet2
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
save_plot2(p15, cnps_comppath, "cnps_ternary_plot", width = 10, height = 8)
save_plot2(p15, cnps_comppath, "cnps_ternary_plot", width = 10, height = 8)

# 保存三元图数据到总表
# addWorksheet(cnps_comp_wb, "ternary_data")  # 已合并到 write_sheet2
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

save_plot2(p4_1, cnps_comppath, "cnps_barplot_main", width = 12, height = 8)
save_plot2(p4_1, cnps_comppath, "cnps_barplot_main", width = 12, height = 8)

p4_2  <- result[[3]]+scale_fill_brewer(palette = "Paired")
p4_2 <- p4_2+mytheme1

save_plot2(p4_2, cnps_comppath, "cnps_barplot_summary", width = 10, height = 8)
save_plot2(p4_2, cnps_comppath, "cnps_barplot_summary", width = 10, height = 8)

databar <- result[[2]] %>% group_by(Group,aa) %>%
  dplyr::summarise(sum(Abundance)) %>% as.data.frame()
colnames(databar) = c("Group","Function","Abundance(%)")
head(databar)

# 保存柱状图数据到总表
# addWorksheet(cnps_comp_wb, "barplot_data")  # 已合并到 write_sheet2
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
save_plot2(p10.1, cnps_comppath, "cnps_venn_diagram", width = 10, height = 8)
save_plot2(p10.1, cnps_comppath, "cnps_venn_diagram", width = 10, height = 8)
save_plot2(p10.2, cnps_comppath, "cnps_upset_plot", width = 12, height = 8)
save_plot2(p10.2, cnps_comppath, "cnps_upset_plot", width = 12, height = 8)

# 保存Venn图数据到总表
# addWorksheet(cnps_comp_wb, "venn_data")  # 已合并到 write_sheet2
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
save_plot2(p13.1, cnps_comppath, "cnps_flower_plot", width = 10, height = 10)
save_plot2(p13.1, cnps_comppath, "cnps_flower_plot", width = 10, height = 10)

# 保存花瓣图数据到总表
# addWorksheet(cnps_comp_wb, "flower_data")  # 已合并到 write_sheet2
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
save_plot2(p14, cnps_comppath, "cnps_venn_network", width = 12, height = 10)
save_plot2(p14, cnps_comppath, "cnps_venn_network", width = 12, height = 10)

# 保存Venn网络数据到总表
# addWorksheet(cnps_comp_wb, "venn_network_data")  # 已合并到 write_sheet2
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
save_plot2(p15.1, cnps_diffpath, "cnps_deseq2_volcano1", width = 10, height = 8)
save_plot2(p15.1, cnps_diffpath, "cnps_deseq2_volcano1", width = 10, height = 8)
save_plot2(p15.2, cnps_diffpath, "cnps_deseq2_volcano2", width = 10, height = 8)
save_plot2(p15.2, cnps_diffpath, "cnps_deseq2_volcano2", width = 10, height = 8)
save_plot2(p15.3, cnps_diffpath, "cnps_deseq2_volcano3", width = 10, height = 8)
save_plot2(p15.3, cnps_diffpath, "cnps_deseq2_volcano3", width = 10, height = 8)

# 保存DESeq2数据到总表
# addWorksheet(cnps_diff_wb, "cnps_deseq2_results")  # 已合并到 write_sheet2
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
save_plot2(p16.1, cnps_diffpath, "cnps_edger_volcano1", width = 10, height = 8)
save_plot2(p16.1, cnps_diffpath, "cnps_edger_volcano1", width = 10, height = 8)
save_plot2(p16.2, cnps_diffpath, "cnps_edger_volcano2", width = 10, height = 8)
save_plot2(p16.2, cnps_diffpath, "cnps_edger_volcano2", width = 10, height = 8)
save_plot2(p16.3, cnps_diffpath, "cnps_edger_volcano3", width = 10, height = 8)
save_plot2(p16.3, cnps_diffpath, "cnps_edger_volcano3", width = 10, height = 8)

# 保存EdgeR数据到总表
# addWorksheet(cnps_diff_wb, "cnps_edger_results")  # 已合并到 write_sheet2
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
save_plot2(p17, cnps_diffpath, "cnps_ttest_main", width = 10, height = 8)
save_plot2(p17, cnps_diffpath, "cnps_ttest_main", width = 10, height = 8)
save_plot2(p17.1, cnps_diffpath, "cnps_ttest_volcano1", width = 10, height = 8)
save_plot2(p17.1, cnps_diffpath, "cnps_ttest_volcano1", width = 10, height = 8)
save_plot2(p17.2, cnps_diffpath, "cnps_ttest_volcano2", width = 10, height = 8)
save_plot2(p17.2, cnps_diffpath, "cnps_ttest_volcano2", width = 10, height = 8)
save_plot2(p17.3, cnps_diffpath, "cnps_ttest_volcano3", width = 10, height = 8)
save_plot2(p17.3, cnps_diffpath, "cnps_ttest_volcano3", width = 10, height = 8)

# 保存t检验数据到总表
# addWorksheet(cnps_diff_wb, "cnps_ttest_results")  # 已合并到 write_sheet2
# addWorksheet(cnps_diff_wb, "cnps_ttest_intersections")  # 已合并到 write_sheet2
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
save_plot2(p18, cnps_diffpath, "cnps_wilcox_main", width = 10, height = 8)
save_plot2(p18, cnps_diffpath, "cnps_wilcox_main", width = 10, height = 8)
save_plot2(p18.1, cnps_diffpath, "cnps_wilcox_volcano1", width = 10, height = 8)
save_plot2(p18.1, cnps_diffpath, "cnps_wilcox_volcano1", width = 10, height = 8)
save_plot2(p18.2, cnps_diffpath, "cnps_wilcox_volcano2", width = 10, height = 8)
save_plot2(p18.2, cnps_diffpath, "cnps_wilcox_volcano2", width = 10, height = 8)

# 20 stamp.metf:stamp差异分析 ----
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
save_plot2(p19, cnps_diffpath, "cnps_stamp_analysis", width = 12, height = 8)
save_plot2(p19, cnps_diffpath, "cnps_stamp_analysis", width = 12, height = 8)

# 保存STAMP数据到总表
# addWorksheet(cnps_diff_wb, "cnps_stamp_plot_data")  # 已合并到 write_sheet2
# addWorksheet(cnps_diff_wb, "cnps_stamp_results")  # 已合并到 write_sheet2
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
save_plot2(prfcv, cnps_biomarkerpath, "cnps_random_forest_cv", width = 10, height = 8)
save_plot2(prfcv, cnps_biomarkerpath, "cnps_random_forest_cv", width = 10, height = 8)

# 保存随机森林交叉验证数据到总表
# addWorksheet(cnps_biomarker_wb, "rfcv_plot_data")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "rfcv_table")  # 已合并到 write_sheet2
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
save_plot2(p33.1, cnps_biomarkerpath, "cnps_roc_curve", width = 10, height = 8)
save_plot2(p33.1, cnps_biomarkerpath, "cnps_roc_curve", width = 10, height = 8)
save_plot2(p33.2, cnps_biomarkerpath, "cnps_roc_boxplot", width = 10, height = 8)
save_plot2(p33.2, cnps_biomarkerpath, "cnps_roc_boxplot", width = 10, height = 8)

# 保存ROC数据到总表
# addWorksheet(cnps_biomarker_wb, "roc_data")  # 已合并到 write_sheet2
writeData(cnps_biomarker_wb, "roc_data", dat, rowNames = TRUE)

#23 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.cnps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存载荷矩阵PCA图片
save_plot2(p34.1, cnps_biomarkerpath, "cnps_loading_pca", width = 10, height = 8)
save_plot2(p34.1, cnps_biomarkerpath, "cnps_loading_pca", width = 10, height = 8)

# 保存载荷矩阵PCA数据到总表
# addWorksheet(cnps_biomarker_wb, "loading_pca_data")  # 已合并到 write_sheet2
writeData(cnps_biomarker_wb, "loading_pca_data", dat, rowNames = TRUE)

#24 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果数据到总表
# addWorksheet(cnps_biomarker_wb, "svm_auc")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "svm_importance")  # 已合并到 write_sheet2
writeData(cnps_biomarker_wb, "svm_auc", AUC, rowNames = TRUE)
writeData(cnps_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#25 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = ps.cnps%>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果数据到总表
# addWorksheet(cnps_biomarker_wb, "glm_auc")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "glm_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cnps_biomarker_wb, "xgboost_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "xgboost_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cnps_biomarker_wb, "lasso_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "lasso_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cnps_biomarker_wb, "decision_tree_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "decision_tree_importance")  # 已合并到 write_sheet2
writeData(cnps_biomarker_wb, "decision_tree_accuracy", accuracy, rowNames = TRUE)
writeData(cnps_biomarker_wb, "decision_tree_importance", importance, rowNames = TRUE)

#29 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps.cnps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果数据到总表
# addWorksheet(cnps_biomarker_wb, "naive_bayes_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "naive_bayes_importance")  # 已合并到 write_sheet2
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
save_plot2(p35, cnps_biomarkerpath, "cnps_lda_lefse", width = 12, height = 8)
save_plot2(p35, cnps_biomarkerpath, "cnps_lda_lefse", width = 12, height = 8)

# 保存LDA数据到总表
# addWorksheet(cnps_biomarker_wb, "lda_summary")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "lda_detailed_results")  # 已合并到 write_sheet2
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
save_plot2(p42.1, cnps_biomarkerpath, "cnps_random_forest_importance", width = 12, height = 8)
save_plot2(p42.1, cnps_biomarkerpath, "cnps_random_forest_importance", width = 12, height = 8)
save_plot2(p42.2, cnps_biomarkerpath, "cnps_random_forest_error", width = 10, height = 8)
save_plot2(p42.2, cnps_biomarkerpath, "cnps_random_forest_error", width = 10, height = 8)
save_plot2(p42.4, cnps_biomarkerpath, "cnps_random_forest_prediction", width = 10, height = 8)
save_plot2(p42.4, cnps_biomarkerpath, "cnps_random_forest_prediction", width = 10, height = 8)

# 保存随机森林数据到总表
# addWorksheet(cnps_biomarker_wb, "random_forest_data")  # 已合并到 write_sheet2
writeData(cnps_biomarker_wb, "random_forest_data", dat, rowNames = TRUE)

#32 nnet.metf: 神经网络筛选特征功能  ------
library(nnet)
res =nnet.metf(ps=ps.cnps, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果数据到总表
# addWorksheet(cnps_biomarker_wb, "neural_network_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "neural_network_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cnps_biomarker_wb, "bagging_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cnps_biomarker_wb, "bagging_importance")  # 已合并到 write_sheet2
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
save_plot2(p_C, cnps_networkpath, "cnps_network_C", width = 12, height = 10)
save_plot2(p_C, cnps_networkpath, "cnps_network_C", width = 12, height = 10)

# N循环网络
res_N = CNPS.network2(ps = ps.kegg,id.0 = "N")
p_N = res_N[[1]]
dat_N = res_N[[2]]
dat2_N = res_N[[3]]

# 保存N循环网络图片
save_plot2(p_N, cnps_networkpath, "cnps_network_N", width = 12, height = 10)
save_plot2(p_N, cnps_networkpath, "cnps_network_N", width = 12, height = 10)

# P循环网络
res_P = CNPS.network2(ps = ps.kegg,id.0 = "P")
p_P = res_P[[1]]
dat_P = res_P[[2]]
dat2_P = res_P[[3]]

# 保存P循环网络图片
save_plot2(p_P, cnps_networkpath, "cnps_network_P", width = 12, height = 10)
save_plot2(p_P, cnps_networkpath, "cnps_network_P", width = 12, height = 10)

# S循环网络
res_S = CNPS.network2(ps = ps.kegg,id.0 = "S")
p_S = res_S[[1]]
dat_S = res_S[[2]]
dat2_S = res_S[[3]]

# 保存S循环网络图片
save_plot2(p_S, cnps_networkpath, "cnps_network_S", width = 12, height = 10)
save_plot2(p_S, cnps_networkpath, "cnps_network_S", width = 12, height = 10)

# 保存CNPS网络数据到总表
# addWorksheet(cnps_network_wb, "C_cycle_data1")  # 已合并到 write_sheet2
# addWorksheet(cnps_network_wb, "C_cycle_data2")  # 已合并到 write_sheet2
# addWorksheet(cnps_network_wb, "N_cycle_data1")  # 已合并到 write_sheet2
# addWorksheet(cnps_network_wb, "N_cycle_data2")  # 已合并到 write_sheet2
# addWorksheet(cnps_network_wb, "P_cycle_data1")  # 已合并到 write_sheet2
# addWorksheet(cnps_network_wb, "P_cycle_data2")  # 已合并到 write_sheet2
# addWorksheet(cnps_network_wb, "S_cycle_data1")  # 已合并到 write_sheet2
# addWorksheet(cnps_network_wb, "S_cycle_data2")  # 已合并到 write_sheet2
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

# 1 alpha.metf:6种alpha多样性计算 ----
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
save_plot2(p1_1, card_alphapath, "card_alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_1, card_alphapath, "card_alpha_diversity_box", width = 12, height = 8)

res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 4)
p1_2 = res[[1]]
p1_2+
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

save_plot2(p1_2, card_alphapath, "card_alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_2, card_alphapath, "card_alpha_diversity_bar", width = 12, height = 8)

res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3:(5)),result = result,sig_show ="abc",ncol = 3)
p1_3 = res[[1]]
p1_3+
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

save_plot2(p1_3, card_alphapath, "card_alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_3, card_alphapath, "card_alpha_diversity_boxbar", width = 12, height = 8)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  # theme_classic()+
  geom_text(aes(x=group , y=y ,label=stat))
p1_0

save_plot2(p1_0, card_alphapath, "card_alpha_diversity_violin", width = 14, height = 8)
save_plot2(p1_0, card_alphapath, "card_alpha_diversity_violin", width = 14, height = 8)

# 保存CARD alpha多样性数据到总表
# addWorksheet(card_alpha_wb, "card_alpha_data")  # 已合并到 write_sheet2
# addWorksheet(card_alpha_wb, "card_alpha_stats")  # 已合并到 write_sheet2
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

save_plot2(p2_1, card_alphapath, "card_rarefaction_individual", width = 10, height = 8)
save_plot2(p2_1, card_alphapath, "card_rarefaction_individual", width = 10, height = 8)

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2

save_plot2(p2_2, card_alphapath, "card_rarefaction_group", width = 10, height = 8)
save_plot2(p2_2, card_alphapath, "card_rarefaction_group", width = 10, height = 8)

#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

# 保存稀释曲线数据到alpha总表
# addWorksheet(card_alpha_wb, "rarefaction_data")  # 已合并到 write_sheet2
writeData(card_alpha_wb, "rarefaction_data", raretab, rowNames = TRUE)

# 保存alpha多样性总表
saveWorkbook(card_alpha_wb, file.path(card_alphapath, "card_alpha_diversity.xlsx"), overwrite = TRUE)

# 创建CARD beta多样性目录
card_betapath <- file.path(cardpath, "beta")
dir.create(card_betapath, recursive = TRUE)

# 创建CARD beta多样性总表
card_beta_wb <- createWorkbook()

# 3 ordinate.metf:排序分析 ----
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
save_plot2(p3_1, card_betapath, "card_ordination_basic", width = 10, height = 8)
save_plot2(p3_1, card_betapath, "card_ordination_basic", width = 10, height = 8)
save_plot2(p3_2, card_betapath, "card_ordination_labeled", width = 10, height = 8)
save_plot2(p3_2, card_betapath, "card_ordination_labeled", width = 10, height = 8)
save_plot2(p3_3, card_betapath, "card_ordination_spider", width = 10, height = 8)
save_plot2(p3_3, card_betapath, "card_ordination_spider", width = 10, height = 8)

# 保存排序分析数据
# addWorksheet(card_beta_wb, "ordination_data")  # 已合并到 write_sheet2
writeData(card_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

# 4 MetaTest.metf:群落功能差异检测 ----
dat1 = MetaTest.metf(ps = ps.card , Micromet = "adonis", dist = "bray")
dat1

# 保存群落差异检测结果
# addWorksheet(card_beta_wb, "adonis_test")  # 已合并到 write_sheet2
writeData(card_beta_wb, "adonis_test", dat1, rowNames = TRUE)

# 5 pairMetaTest.metf:两两分组群落功能差异检测 ----
dat2 = pairMetaTest.metf(ps = ps.card , Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较结果
# addWorksheet(card_beta_wb, "pairwise_test")  # 已合并到 write_sheet2
writeData(card_beta_wb, "pairwise_test", dat2, rowNames = TRUE)

# 6 mantal.metf：群落功能差异检测普鲁士分析#------
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
save_plot2(p3_7, card_betapath, "card_mantel_test", width = 10, height = 8)
save_plot2(p3_7, card_betapath, "card_mantel_test", width = 10, height = 8)

# addWorksheet(card_beta_wb, "mantel_test")  # 已合并到 write_sheet2
writeData(card_beta_wb, "mantel_test", data, rowNames = TRUE)

#7 cluster.metf:样品聚类 ----
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
save_plot2(p4, card_betapath, "card_cluster_plot1", width = 12, height = 8)
save_plot2(p4, card_betapath, "card_cluster_plot1", width = 12, height = 8)
save_plot2(p4_1, card_betapath, "card_cluster_plot2", width = 12, height = 8)
save_plot2(p4_1, card_betapath, "card_cluster_plot2", width = 12, height = 8)
save_plot2(p4_2, card_betapath, "card_cluster_plot3", width = 12, height = 8)
save_plot2(p4_2, card_betapath, "card_cluster_plot3", width = 12, height = 8)

# 保存聚类分析数据
# addWorksheet(card_beta_wb, "cluster_data")  # 已合并到 write_sheet2
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
save_plot2(p15, card_compositionpath, "card_ternary_plot", width = 10, height = 8)
save_plot2(p15, card_compositionpath, "card_ternary_plot", width = 10, height = 8)

# addWorksheet(card_composition_wb, "ternary_data")  # 已合并到 write_sheet2
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
save_plot2(p4_1, card_compositionpath, "card_barplot1", width = 12, height = 8)
save_plot2(p4_1, card_compositionpath, "card_barplot1", width = 12, height = 8)
save_plot2(p4_2, card_compositionpath, "card_barplot2", width = 12, height = 8)
save_plot2(p4_2, card_compositionpath, "card_barplot2", width = 12, height = 8)

# addWorksheet(card_composition_wb, "barplot_data")  # 已合并到 write_sheet2
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
save_plot2(p10.1, card_compositionpath, "card_venn_plot", width = 10, height = 8)
save_plot2(p10.1, card_compositionpath, "card_venn_plot", width = 10, height = 8)
save_plot2(p10.2, card_compositionpath, "card_upset_plot", width = 12, height = 8)
save_plot2(p10.2, card_compositionpath, "card_upset_plot", width = 12, height = 8)

# addWorksheet(card_composition_wb, "venn_upset_data")  # 已合并到 write_sheet2
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
save_plot2(p7_1, card_compositionpath, "card_venn_detail1", width = 12, height = 8)
save_plot2(p7_1, card_compositionpath, "card_venn_detail1", width = 12, height = 8)
save_plot2(p7_2, card_compositionpath, "card_venn_detail2", width = 12, height = 8)
save_plot2(p7_2, card_compositionpath, "card_venn_detail2", width = 12, height = 8)

# addWorksheet(card_composition_wb, "venn_detail_data")  # 已合并到 write_sheet2
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
save_plot2(p13.1, card_compositionpath, "card_flower_plot", width = 10, height = 8)
save_plot2(p13.1, card_compositionpath, "card_flower_plot", width = 10, height = 8)

# addWorksheet(card_composition_wb, "flower_data")  # 已合并到 write_sheet2
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
save_plot2(p14, card_compositionpath, "card_venn_network", width = 10, height = 8)
save_plot2(p14, card_compositionpath, "card_venn_network", width = 10, height = 8)

# addWorksheet(card_composition_wb, "venn_network_data")  # 已合并到 write_sheet2
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
save_plot2(p15.1, card_diffpath, "card_deseq2_plot1", width = 10, height = 8)
save_plot2(p15.1, card_diffpath, "card_deseq2_plot1", width = 10, height = 8)
save_plot2(p15.2, card_diffpath, "card_deseq2_plot2", width = 10, height = 8)
save_plot2(p15.2, card_diffpath, "card_deseq2_plot2", width = 10, height = 8)
save_plot2(p15.3, card_diffpath, "card_deseq2_plot3", width = 10, height = 8)
save_plot2(p15.3, card_diffpath, "card_deseq2_plot3", width = 10, height = 8)

# addWorksheet(card_diff_wb, "DESep2_results")  # 已合并到 write_sheet2
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
save_plot2(p16.1, card_diffpath, "card_edger_plot1", width = 10, height = 8)
save_plot2(p16.1, card_diffpath, "card_edger_plot1", width = 10, height = 8)
save_plot2(p16.2, card_diffpath, "card_edger_plot2", width = 10, height = 8)
save_plot2(p16.2, card_diffpath, "card_edger_plot2", width = 10, height = 8)
save_plot2(p16.3, card_diffpath, "card_edger_plot3", width = 10, height = 8)
save_plot2(p16.3, card_diffpath, "card_edger_plot3", width = 10, height = 8)

# addWorksheet(card_diff_wb, "EdgeR_results")  # 已合并到 write_sheet2
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
save_plot2(p17.1, card_diffpath, "card_ttest_plot1", width = 10, height = 8)
save_plot2(p17.1, card_diffpath, "card_ttest_plot1", width = 10, height = 8)
save_plot2(p17.2, card_diffpath, "card_ttest_plot2", width = 10, height = 8)
save_plot2(p17.2, card_diffpath, "card_ttest_plot2", width = 10, height = 8)
save_plot2(p17.3, card_diffpath, "card_ttest_plot3", width = 10, height = 8)
save_plot2(p17.3, card_diffpath, "card_ttest_plot3", width = 10, height = 8)

# addWorksheet(card_diff_wb, "ttest_results")  # 已合并到 write_sheet2
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
save_plot2(p18.1, card_diffpath, "card_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.1, card_diffpath, "card_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.2, card_diffpath, "card_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.2, card_diffpath, "card_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.3, card_diffpath, "card_wilcox_plot3", width = 10, height = 8)
save_plot2(p18.3, card_diffpath, "card_wilcox_plot3", width = 10, height = 8)

# addWorksheet(card_diff_wb, "wilcox_results")  # 已合并到 write_sheet2
writeData(card_diff_wb, "wilcox_results", dat, rowNames = TRUE)

# 19 stamp.metf:stamp差异分析 ----
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
save_plot2(p19, card_diffpath, "card_stamp_plot", width = 12, height = 8)
save_plot2(p19, card_diffpath, "card_stamp_plot", width = 12, height = 8)

# addWorksheet(card_diff_wb, "STAMP_results1")  # 已合并到 write_sheet2
# addWorksheet(card_diff_wb, "STAMP_results2")  # 已合并到 write_sheet2
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
save_plot2(prfcv, card_biomarkerpath, "card_rfcv_plot", width = 10, height = 8)
save_plot2(prfcv, card_biomarkerpath, "card_rfcv_plot", width = 10, height = 8)

# addWorksheet(card_biomarker_wb, "RFCV_results")  # 已合并到 write_sheet2
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
save_plot2(p33.1, card_biomarkerpath, "card_roc_plot1", width = 10, height = 8)
save_plot2(p33.1, card_biomarkerpath, "card_roc_plot1", width = 10, height = 8)
save_plot2(p33.2, card_biomarkerpath, "card_roc_plot2", width = 10, height = 8)
save_plot2(p33.2, card_biomarkerpath, "card_roc_plot2", width = 10, height = 8)

# addWorksheet(card_biomarker_wb, "ROC_results")  # 已合并到 write_sheet2
writeData(card_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#22 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.card,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
save_plot2(p34.1, card_biomarkerpath, "card_pca_loading", width = 10, height = 8)
save_plot2(p34.1, card_biomarkerpath, "card_pca_loading", width = 10, height = 8)

# addWorksheet(card_biomarker_wb, "PCA_loading")  # 已合并到 write_sheet2
writeData(card_biomarker_wb, "PCA_loading", dat, rowNames = TRUE)

#23 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = ps.card %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
# addWorksheet(card_biomarker_wb, "SVM_AUC")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "SVM_importance")  # 已合并到 write_sheet2
writeData(card_biomarker_wb, "SVM_AUC", AUC, rowNames = TRUE)
writeData(card_biomarker_wb, "SVM_importance", importance, rowNames = TRUE)

#24 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = ps.card %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
# addWorksheet(card_biomarker_wb, "GLM_AUC")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "GLM_importance")  # 已合并到 write_sheet2
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
# addWorksheet(card_biomarker_wb, "XGBoost_accuracy")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "XGBoost_importance")  # 已合并到 write_sheet2
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
# addWorksheet(card_biomarker_wb, "LASSO_accuracy")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "LASSO_importance")  # 已合并到 write_sheet2
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
# addWorksheet(card_biomarker_wb, "DecisionTree_accuracy")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "DecisionTree_importance")  # 已合并到 write_sheet2
writeData(card_biomarker_wb, "DecisionTree_accuracy", accuracy, rowNames = TRUE)
writeData(card_biomarker_wb, "DecisionTree_importance", importance, rowNames = TRUE)

#28 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps.card, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
# addWorksheet(card_biomarker_wb, "NaiveBayes_accuracy")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "NaiveBayes_importance")  # 已合并到 write_sheet2
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
save_plot2(p35, card_biomarkerpath, "card_lda_plot", width = 12, height = 8)
save_plot2(p35, card_biomarkerpath, "card_lda_plot", width = 12, height = 8)

# addWorksheet(card_biomarker_wb, "LDA_results")  # 已合并到 write_sheet2
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
save_plot2(p42.1, card_biomarkerpath, "card_rf_plot1", width = 10, height = 8)
save_plot2(p42.1, card_biomarkerpath, "card_rf_plot1", width = 10, height = 8)
save_plot2(p42.2, card_biomarkerpath, "card_rf_plot2", width = 10, height = 8)
save_plot2(p42.2, card_biomarkerpath, "card_rf_plot2", width = 10, height = 8)
save_plot2(p42.4, card_biomarkerpath, "card_rf_plot4", width = 10, height = 8)
save_plot2(p42.4, card_biomarkerpath, "card_rf_plot4", width = 10, height = 8)

# addWorksheet(card_biomarker_wb, "RandomForest_results")  # 已合并到 write_sheet2
writeData(card_biomarker_wb, "RandomForest_results", dat, rowNames = TRUE)

#31 nnet.metf: 神经网络筛选特征功能  ------
res =nnet.metf(ps=ps.card, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
# addWorksheet(card_biomarker_wb, "NeuralNet_accuracy")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "NeuralNet_importance")  # 已合并到 write_sheet2
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
# addWorksheet(card_biomarker_wb, "Bagging_accuracy")  # 已合并到 write_sheet2
# addWorksheet(card_biomarker_wb, "Bagging_importance")  # 已合并到 write_sheet2
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
save_plot2(p0, card_networkpath, "card_network_main", width = 12, height = 10)
save_plot2(p0, card_networkpath, "card_network_main", width = 12, height = 10)

# 34 net_properties.4:网络属性计算 ----
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
# addWorksheet(card_network_wb, "network_properties")  # 已合并到 write_sheet2
writeData(card_network_wb, "network_properties", dat2, rowNames = TRUE)

# 35 netproperties.sample:单个样本的网络属性 ----
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
# addWorksheet(card_network_wb, "sample_network_properties")  # 已合并到 write_sheet2
writeData(card_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

# 36 node_properties:计算节点属性 ----
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
# addWorksheet(card_network_wb, "node_properties")  # 已合并到 write_sheet2
writeData(card_network_wb, "node_properties", nodepro2, rowNames = TRUE)

# 37 module.compare.net.pip:网络显著性比较 ----
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
# addWorksheet(card_network_wb, "network_comparison")  # 已合并到 write_sheet2
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

# 1 alpha.metf:6种alpha多样性计算 ----
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
save_plot2(p1_1, cazy_alphapath, "cazy_alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_1, cazy_alphapath, "cazy_alpha_diversity_box", width = 12, height = 8)

res = EasyStat::FacetMuiPlotresultBar(data = data, num = c(3:(5)), result = result, sig_show = "abc", ncol = 4)
p1_2 = res[[1]]
p1_2 +
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

save_plot2(p1_2, cazy_alphapath, "cazy_alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_2, cazy_alphapath, "cazy_alpha_diversity_bar", width = 12, height = 8)

res = EasyStat::FacetMuiPlotReBoxBar(data = data, num = c(3:(5)), result = result, sig_show = "abc", ncol = 3)
p1_3 = res[[1]]
p1_3 +
  scale_x_discrete(limits = axis_order) +
  mytheme1 +
  guides(fill = guide_legend(title = NULL)) +
  scale_fill_manual(values = colset1)

save_plot2(p1_3, cazy_alphapath, "cazy_alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_3, cazy_alphapath, "cazy_alpha_diversity_boxbar", width = 12, height = 8)

#基于输出数据使用ggplot出图
p1_0 = result1[[2]] %>% ggplot(aes(x=group , y=dd )) +
  geom_violin(alpha=1, aes(fill=group)) +
  geom_jitter( aes(color = group),position=position_jitter(0.17), size=3, alpha=0.5)+
  labs(x="", y="")+
  facet_wrap(.~name,scales="free_y",ncol  = 4) +
  geom_text(aes(x=group , y=y ,label=stat))
p1_0

save_plot2(p1_0, cazy_alphapath, "cazy_alpha_diversity_violin", width = 14, height = 8)
save_plot2(p1_0, cazy_alphapath, "cazy_alpha_diversity_violin", width = 14, height = 8)

# 保存CAZY alpha多样性数据到总表
# addWorksheet(cazy_alpha_wb, "cazy_alpha_data")  # 已合并到 write_sheet2
# addWorksheet(cazy_alpha_wb, "cazy_alpha_stats")  # 已合并到 write_sheet2
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

save_plot2(p2_1, cazy_alphapath, "cazy_rarefaction_individual", width = 10, height = 8)
save_plot2(p2_1, cazy_alphapath, "cazy_rarefaction_individual", width = 10, height = 8)

## 提供数据表格，方便输出
raretab <- result[[2]]
head(raretab)
#--按照分组展示稀释曲线
p2_2 <- result[[3]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_2

save_plot2(p2_2, cazy_alphapath, "cazy_rarefaction_group", width = 10, height = 8)
save_plot2(p2_2, cazy_alphapath, "cazy_rarefaction_group", width = 10, height = 8)

#--按照分组绘制标准差稀释曲线
p2_3 <- result[[4]] +
  mytheme1 +
  guides(fill = guide_legend(title = NULL))+
  scale_fill_manual(values = colset1)
p2_3

# 保存稀释曲线数据到alpha总表
# addWorksheet(cazy_alpha_wb, "rarefaction_data")  # 已合并到 write_sheet2
writeData(cazy_alpha_wb, "rarefaction_data", raretab, rowNames = TRUE)

# 保存alpha多样性总表
saveWorkbook(cazy_alpha_wb, file.path(cazy_alphapath, "cazy_alpha_diversity.xlsx"), overwrite = TRUE)

# 创建CAZY beta多样性目录
cazy_betapath <- file.path(cazypath, "beta")
dir.create(cazy_betapath, recursive = TRUE)

# 创建CAZY beta多样性总表
cazy_beta_wb <- createWorkbook()

# 3 ordinate.metf:排序分析 ----
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
save_plot2(p3_1, cazy_betapath, "cazy_ordination_basic", width = 10, height = 8)
save_plot2(p3_1, cazy_betapath, "cazy_ordination_basic", width = 10, height = 8)
save_plot2(p3_2, cazy_betapath, "cazy_ordination_labeled", width = 10, height = 8)
save_plot2(p3_2, cazy_betapath, "cazy_ordination_labeled", width = 10, height = 8)
save_plot2(p3_3, cazy_betapath, "cazy_ordination_spider", width = 10, height = 8)
save_plot2(p3_3, cazy_betapath, "cazy_ordination_spider", width = 10, height = 8)

# 保存排序分析数据
# addWorksheet(cazy_beta_wb, "ordination_data")  # 已合并到 write_sheet2
writeData(cazy_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

# 4 MetaTest.metf:群落功能差异检测 ----
dat1 = MetaTest.metf(ps = ps.cazy, Micromet = "adonis", dist = "bray")
dat1

# 保存群落差异检测结果
# addWorksheet(cazy_beta_wb, "adonis_test")  # 已合并到 write_sheet2
writeData(cazy_beta_wb, "adonis_test", dat1, rowNames = TRUE)

# 5 pairMetaTest.metf:两两分组群落功能差异检测 ----
dat2 = pairMetaTest.metf(ps = ps.cazy, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两比较结果
# addWorksheet(cazy_beta_wb, "pairwise_test")  # 已合并到 write_sheet2
writeData(cazy_beta_wb, "pairwise_test", dat2, rowNames = TRUE)

# 6 mantal.metf：群落功能差异检测普鲁士分析#------
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
save_plot2(p3_7, cazy_betapath, "cazy_mantel_test", width = 10, height = 8)
save_plot2(p3_7, cazy_betapath, "cazy_mantel_test", width = 10, height = 8)

# addWorksheet(cazy_beta_wb, "mantel_test")  # 已合并到 write_sheet2
writeData(cazy_beta_wb, "mantel_test", data, rowNames = TRUE)

#7 cluster.metf:样品聚类 ----
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
save_plot2(p4, cazy_betapath, "cazy_cluster_plot1", width = 12, height = 8)
save_plot2(p4, cazy_betapath, "cazy_cluster_plot1", width = 12, height = 8)
save_plot2(p4_1, cazy_betapath, "cazy_cluster_plot2", width = 12, height = 8)
save_plot2(p4_1, cazy_betapath, "cazy_cluster_plot2", width = 12, height = 8)
save_plot2(p4_2, cazy_betapath, "cazy_cluster_plot3", width = 12, height = 8)
save_plot2(p4_2, cazy_betapath, "cazy_cluster_plot3", width = 12, height = 8)

# 保存聚类分析数据
# addWorksheet(cazy_beta_wb, "cluster_data")  # 已合并到 write_sheet2
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
save_plot2(p15, cazy_compositionpath, "cazy_ternary_plot", width = 10, height = 8)
save_plot2(p15, cazy_compositionpath, "cazy_ternary_plot", width = 10, height = 8)

# addWorksheet(cazy_composition_wb, "ternary_data")  # 已合并到 write_sheet2
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
save_plot2(p4_1, cazy_compositionpath, "cazy_barplot1", width = 12, height = 8)
save_plot2(p4_1, cazy_compositionpath, "cazy_barplot1", width = 12, height = 8)
save_plot2(p4_2, cazy_compositionpath, "cazy_barplot2", width = 12, height = 8)
save_plot2(p4_2, cazy_compositionpath, "cazy_barplot2", width = 12, height = 8)

# addWorksheet(cazy_composition_wb, "barplot_data")  # 已合并到 write_sheet2
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
save_plot2(p10.1, cazy_compositionpath, "cazy_venn_plot", width = 10, height = 8)
save_plot2(p10.1, cazy_compositionpath, "cazy_venn_plot", width = 10, height = 8)
save_plot2(p10.2, cazy_compositionpath, "cazy_upset_plot", width = 12, height = 8)
save_plot2(p10.2, cazy_compositionpath, "cazy_upset_plot", width = 12, height = 8)

# addWorksheet(cazy_composition_wb, "venn_upset_data")  # 已合并到 write_sheet2
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
save_plot2(p7_1, cazy_compositionpath, "cazy_venn_detail1", width = 12, height = 8)
save_plot2(p7_1, cazy_compositionpath, "cazy_venn_detail1", width = 12, height = 8)
save_plot2(p7_2, cazy_compositionpath, "cazy_venn_detail2", width = 12, height = 8)
save_plot2(p7_2, cazy_compositionpath, "cazy_venn_detail2", width = 12, height = 8)

# addWorksheet(cazy_composition_wb, "venn_detail_data")  # 已合并到 write_sheet2
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
save_plot2(p13.1, cazy_compositionpath, "cazy_flower_plot", width = 10, height = 8)
save_plot2(p13.1, cazy_compositionpath, "cazy_flower_plot", width = 10, height = 8)

# addWorksheet(cazy_composition_wb, "flower_data")  # 已合并到 write_sheet2
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
save_plot2(p14, cazy_compositionpath, "cazy_venn_network", width = 10, height = 8)
save_plot2(p14, cazy_compositionpath, "cazy_venn_network", width = 10, height = 8)

# addWorksheet(cazy_composition_wb, "venn_network_data")  # 已合并到 write_sheet2
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
save_plot2(p15.1, cazy_diffpath, "cazy_deseq2_plot1", width = 10, height = 8)
save_plot2(p15.1, cazy_diffpath, "cazy_deseq2_plot1", width = 10, height = 8)
save_plot2(p15.2, cazy_diffpath, "cazy_deseq2_plot2", width = 10, height = 8)
save_plot2(p15.2, cazy_diffpath, "cazy_deseq2_plot2", width = 10, height = 8)
save_plot2(p15.3, cazy_diffpath, "cazy_deseq2_plot3", width = 10, height = 8)
save_plot2(p15.3, cazy_diffpath, "cazy_deseq2_plot3", width = 10, height = 8)

# addWorksheet(cazy_diff_wb, "DESep2_results")  # 已合并到 write_sheet2
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
save_plot2(p16.1, cazy_diffpath, "cazy_edger_plot1", width = 10, height = 8)
save_plot2(p16.1, cazy_diffpath, "cazy_edger_plot1", width = 10, height = 8)
save_plot2(p16.2, cazy_diffpath, "cazy_edger_plot2", width = 10, height = 8)
save_plot2(p16.2, cazy_diffpath, "cazy_edger_plot2", width = 10, height = 8)
save_plot2(p16.3, cazy_diffpath, "cazy_edger_plot3", width = 10, height = 8)
save_plot2(p16.3, cazy_diffpath, "cazy_edger_plot3", width = 10, height = 8)

# addWorksheet(cazy_diff_wb, "EdgeR_results")  # 已合并到 write_sheet2
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
save_plot2(p17.1, cazy_diffpath, "cazy_ttest_plot1", width = 10, height = 8)
save_plot2(p17.1, cazy_diffpath, "cazy_ttest_plot1", width = 10, height = 8)
save_plot2(p17.2, cazy_diffpath, "cazy_ttest_plot2", width = 10, height = 8)
save_plot2(p17.2, cazy_diffpath, "cazy_ttest_plot2", width = 10, height = 8)
save_plot2(p17.3, cazy_diffpath, "cazy_ttest_plot3", width = 10, height = 8)
save_plot2(p17.3, cazy_diffpath, "cazy_ttest_plot3", width = 10, height = 8)

# addWorksheet(cazy_diff_wb, "ttest_results")  # 已合并到 write_sheet2
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
save_plot2(p18.1, cazy_diffpath, "cazy_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.1, cazy_diffpath, "cazy_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.2, cazy_diffpath, "cazy_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.2, cazy_diffpath, "cazy_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.3, cazy_diffpath, "cazy_wilcox_plot3", width = 10, height = 8)
save_plot2(p18.3, cazy_diffpath, "cazy_wilcox_plot3", width = 10, height = 8)

# addWorksheet(cazy_diff_wb, "wilcox_results")  # 已合并到 write_sheet2
writeData(cazy_diff_wb, "wilcox_results", dat, rowNames = TRUE)

# 18 stamp.metf:stamp差异分析 ----
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
save_plot2(p19, cazy_diffpath, "cazy_stamp_plot", width = 12, height = 8)
save_plot2(p19, cazy_diffpath, "cazy_stamp_plot", width = 12, height = 8)

# addWorksheet(cazy_diff_wb, "STAMP_results1")  # 已合并到 write_sheet2
# addWorksheet(cazy_diff_wb, "STAMP_results2")  # 已合并到 write_sheet2
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
save_plot2(prfcv, cazy_biomarkerpath, "cazy_rfcv_plot", width = 10, height = 8)
save_plot2(prfcv, cazy_biomarkerpath, "cazy_rfcv_plot", width = 10, height = 8)

# addWorksheet(cazy_biomarker_wb, "RFCV_results")  # 已合并到 write_sheet2
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
save_plot2(p33.1, cazy_biomarkerpath, "cazy_roc_plot1", width = 10, height = 8)
save_plot2(p33.1, cazy_biomarkerpath, "cazy_roc_plot1", width = 10, height = 8)
save_plot2(p33.2, cazy_biomarkerpath, "cazy_roc_plot2", width = 10, height = 8)
save_plot2(p33.2, cazy_biomarkerpath, "cazy_roc_plot2", width = 10, height = 8)

# addWorksheet(cazy_biomarker_wb, "ROC_results")  # 已合并到 write_sheet2
writeData(cazy_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps.cazy,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
save_plot2(p34.1, cazy_biomarkerpath, "cazy_pca_loading", width = 10, height = 8)
save_plot2(p34.1, cazy_biomarkerpath, "cazy_pca_loading", width = 10, height = 8)

# addWorksheet(cazy_biomarker_wb, "PCA_loading")  # 已合并到 write_sheet2
writeData(cazy_biomarker_wb, "PCA_loading", dat, rowNames = TRUE)

#22 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = ps.cazy %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
# addWorksheet(cazy_biomarker_wb, "SVM_AUC")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "SVM_importance")  # 已合并到 write_sheet2
writeData(cazy_biomarker_wb, "SVM_AUC", AUC, rowNames = TRUE)
writeData(cazy_biomarker_wb, "SVM_importance", importance, rowNames = TRUE)

#23 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = ps.cazy %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
# addWorksheet(cazy_biomarker_wb, "GLM_AUC")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "GLM_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cazy_biomarker_wb, "XGBoost_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "XGBoost_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cazy_biomarker_wb, "LASSO_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "LASSO_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cazy_biomarker_wb, "DecisionTree_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "DecisionTree_importance")  # 已合并到 write_sheet2
writeData(cazy_biomarker_wb, "DecisionTree_accuracy", accuracy, rowNames = TRUE)
writeData(cazy_biomarker_wb, "DecisionTree_importance", importance, rowNames = TRUE)

#27 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps = ps.cazy, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
# addWorksheet(cazy_biomarker_wb, "NaiveBayes_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "NaiveBayes_importance")  # 已合并到 write_sheet2
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
save_plot2(p35, cazy_biomarkerpath, "cazy_lda_plot", width = 12, height = 8)
save_plot2(p35, cazy_biomarkerpath, "cazy_lda_plot", width = 12, height = 8)

# addWorksheet(cazy_biomarker_wb, "LDA_results")  # 已合并到 write_sheet2
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
save_plot2(p42.1, cazy_biomarkerpath, "cazy_rf_plot1", width = 10, height = 8)
save_plot2(p42.1, cazy_biomarkerpath, "cazy_rf_plot1", width = 10, height = 8)
save_plot2(p42.2, cazy_biomarkerpath, "cazy_rf_plot2", width = 10, height = 8)
save_plot2(p42.2, cazy_biomarkerpath, "cazy_rf_plot2", width = 10, height = 8)
save_plot2(p42.4, cazy_biomarkerpath, "cazy_rf_plot4", width = 10, height = 8)
save_plot2(p42.4, cazy_biomarkerpath, "cazy_rf_plot4", width = 10, height = 8)

# addWorksheet(cazy_biomarker_wb, "RandomForest_results")  # 已合并到 write_sheet2
writeData(cazy_biomarker_wb, "RandomForest_results", dat, rowNames = TRUE)

#30 nnet.metf: 神经网络筛选特征功能  ------
res =nnet.metf(ps = ps.cazy, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
# addWorksheet(cazy_biomarker_wb, "NeuralNet_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "NeuralNet_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cazy_biomarker_wb, "Bagging_accuracy")  # 已合并到 write_sheet2
# addWorksheet(cazy_biomarker_wb, "Bagging_importance")  # 已合并到 write_sheet2
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
save_plot2(p0, cazy_networkpath, "cazy_network_main", width = 12, height = 10)
save_plot2(p0, cazy_networkpath, "cazy_network_main", width = 12, height = 10)

# 33 net_properties.4:网络属性计算 ----
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
# addWorksheet(cazy_network_wb, "network_properties")  # 已合并到 write_sheet2
writeData(cazy_network_wb, "network_properties", dat2, rowNames = TRUE)

# 34 netproperties.sample:单个样本的网络属性 ----
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
# addWorksheet(cazy_network_wb, "sample_network_properties")  # 已合并到 write_sheet2
writeData(cazy_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

# 35 node_properties:计算节点属性 ----
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
# addWorksheet(cazy_network_wb, "node_properties")  # 已合并到 write_sheet2
writeData(cazy_network_wb, "node_properties", nodepro2, rowNames = TRUE)

# 37 module.compare.net.pip:网络显著性比较 ----
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
# addWorksheet(cazy_network_wb, "network_comparison")  # 已合并到 write_sheet2
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

# 1 alpha.metf:6种alpha多样性计算 ----
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
save_plot2(p1_1, cog_diversitypath, "cog_alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_1, cog_diversitypath, "cog_alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_2, cog_diversitypath, "cog_alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_2, cog_diversitypath, "cog_alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_3, cog_diversitypath, "cog_alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_3, cog_diversitypath, "cog_alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_0, cog_diversitypath, "cog_alpha_diversity_violin", width = 12, height = 8)
save_plot2(p1_0, cog_diversitypath, "cog_alpha_diversity_violin", width = 12, height = 8)
# addWorksheet(cog_diversity_wb, "alpha_diversity_data")  # 已合并到 write_sheet2
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
save_plot2(p2_1, cog_diversitypath, "cog_rarefaction_individual", width = 10, height = 8)
save_plot2(p2_1, cog_diversitypath, "cog_rarefaction_individual", width = 10, height = 8)
save_plot2(p2_2, cog_diversitypath, "cog_rarefaction_group", width = 10, height = 8)
save_plot2(p2_2, cog_diversitypath, "cog_rarefaction_group", width = 10, height = 8)
save_plot2(p2_3, cog_diversitypath, "cog_rarefaction_group_sd", width = 10, height = 8)
save_plot2(p2_3, cog_diversitypath, "cog_rarefaction_group_sd", width = 10, height = 8)
# addWorksheet(cog_diversity_wb, "rarefaction_data")  # 已合并到 write_sheet2
writeData(cog_diversity_wb, "rarefaction_data", raretab, rowNames = TRUE)

# function beta diversity -----
# 创建COG beta多样性分析目录
cog_betapath <- file.path(metapath, "cog", "beta")
dir.create(cog_betapath, recursive = TRUE)

# 创建COG beta多样性分析总表
cog_beta_wb <- createWorkbook()

# 3 ordinate.metf:排序分析 ----
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
save_plot2(p3_1, cog_betapath, "cog_ordination_basic", width = 10, height = 8)
save_plot2(p3_1, cog_betapath, "cog_ordination_basic", width = 10, height = 8)
save_plot2(p3_2, cog_betapath, "cog_ordination_labeled", width = 10, height = 8)
save_plot2(p3_2, cog_betapath, "cog_ordination_labeled", width = 10, height = 8)
save_plot2(p3_3, cog_betapath, "cog_ordination_refined", width = 10, height = 8)
save_plot2(p3_3, cog_betapath, "cog_ordination_refined", width = 10, height = 8)
# addWorksheet(cog_beta_wb, "ordination_data")  # 已合并到 write_sheet2
writeData(cog_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

# 4 MetaTest.metf:群落功能差异检测 ----
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

# 保存群落功能差异检测结果
# addWorksheet(cog_beta_wb, "metatest_results")  # 已合并到 write_sheet2
writeData(cog_beta_wb, "metatest_results", dat1, rowNames = TRUE)

# 5 pairMetaTest.metf:两两分组群落功能差异检测 ----
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两分组群落功能差异检测结果
# addWorksheet(cog_beta_wb, "pair_metatest_results")  # 已合并到 write_sheet2
writeData(cog_beta_wb, "pair_metatest_results", dat2, rowNames = TRUE)

# 6 mantal.metf：群落功能差异检测普鲁士分析#------
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
save_plot2(p3_7, cog_betapath, "cog_mantel_test", width = 10, height = 8)
save_plot2(p3_7, cog_betapath, "cog_mantel_test", width = 10, height = 8)
# addWorksheet(cog_beta_wb, "mantel_results")  # 已合并到 write_sheet2
writeData(cog_beta_wb, "mantel_results", data, rowNames = TRUE)

#7 cluster.metf:样品聚类 ----
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
save_plot2(p4, cog_betapath, "cog_cluster_heatmap", width = 12, height = 10)
save_plot2(p4, cog_betapath, "cog_cluster_heatmap", width = 12, height = 10)
save_plot2(p4_1, cog_betapath, "cog_cluster_dendrogram1", width = 10, height = 8)
save_plot2(p4_1, cog_betapath, "cog_cluster_dendrogram1", width = 10, height = 8)
save_plot2(p4_2, cog_betapath, "cog_cluster_dendrogram2", width = 10, height = 8)
save_plot2(p4_2, cog_betapath, "cog_cluster_dendrogram2", width = 10, height = 8)
# addWorksheet(cog_beta_wb, "cluster_results")  # 已合并到 write_sheet2
writeData(cog_beta_wb, "cluster_results", dat, rowNames = TRUE)

#8 Micro_tern.metf: 三元图展示功能----
ps1 = ps %>% filter_OTU_ps(500)
res = Micro_tern.metf(ps=ps1,color = "Function"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图结果
save_plot2(p15, cog_betapath, "cog_ternary_plot", width = 10, height = 8)
save_plot2(p15, cog_betapath, "cog_ternary_plot", width = 10, height = 8)
# addWorksheet(cog_beta_wb, "ternary_data")  # 已合并到 write_sheet2
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
save_plot2(p4_1, cog_classificationpath, "cog_barplot_main", width = 12, height = 8)
save_plot2(p4_1, cog_classificationpath, "cog_barplot_main", width = 12, height = 8)
save_plot2(p4_2, cog_classificationpath, "cog_barplot_secondary", width = 12, height = 8)
save_plot2(p4_2, cog_classificationpath, "cog_barplot_secondary", width = 12, height = 8)
# addWorksheet(cog_classification_wb, "barplot_data")  # 已合并到 write_sheet2
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
save_plot2(p10.1, cog_classificationpath, "cog_venn_diagram", width = 10, height = 8)
save_plot2(p10.1, cog_classificationpath, "cog_venn_diagram", width = 10, height = 8)
save_plot2(p10.2, cog_classificationpath, "cog_upset_plot", width = 12, height = 8)
save_plot2(p10.2, cog_classificationpath, "cog_upset_plot", width = 12, height = 8)
# addWorksheet(cog_classification_wb, "venn_upset_data")  # 已合并到 write_sheet2
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
save_plot2(p7_1, cog_classificationpath, "cog_venn_detail_bar", width = 12, height = 8)
save_plot2(p7_1, cog_classificationpath, "cog_venn_detail_bar", width = 12, height = 8)
save_plot2(p7_2, cog_classificationpath, "cog_venn_detail_alluvial", width = 12, height = 8)
save_plot2(p7_2, cog_classificationpath, "cog_venn_detail_alluvial", width = 12, height = 8)
# addWorksheet(cog_classification_wb, "venn_detail_data")  # 已合并到 write_sheet2
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
save_plot2(p13.1, cog_classificationpath, "cog_flower_plot", width = 10, height = 10)
save_plot2(p13.1, cog_classificationpath, "cog_flower_plot", width = 10, height = 10)
# addWorksheet(cog_classification_wb, "flower_data")  # 已合并到 write_sheet2
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
save_plot2(p14, cog_classificationpath, "cog_venn_network", width = 10, height = 8)
save_plot2(p14, cog_classificationpath, "cog_venn_network", width = 10, height = 8)
# addWorksheet(cog_classification_wb, "venn_network_data")  # 已合并到 write_sheet2
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
save_plot2(p15.1, cog_diffpath, "cog_DESep2_plot1", width = 10, height = 8)
save_plot2(p15.1, cog_diffpath, "cog_DESep2_plot1", width = 10, height = 8)
save_plot2(p15.2, cog_diffpath, "cog_DESep2_plot2", width = 10, height = 8)
save_plot2(p15.2, cog_diffpath, "cog_DESep2_plot2", width = 10, height = 8)
save_plot2(p15.3, cog_diffpath, "cog_DESep2_plot3", width = 10, height = 8)
save_plot2(p15.3, cog_diffpath, "cog_DESep2_plot3", width = 10, height = 8)
# addWorksheet(cog_diff_wb, "DESep2_results")  # 已合并到 write_sheet2
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
save_plot2(p16.1, cog_diffpath, "cog_EdgeR_plot1", width = 10, height = 8)
save_plot2(p16.1, cog_diffpath, "cog_EdgeR_plot1", width = 10, height = 8)
save_plot2(p16.2, cog_diffpath, "cog_EdgeR_plot2", width = 10, height = 8)
save_plot2(p16.2, cog_diffpath, "cog_EdgeR_plot2", width = 10, height = 8)
save_plot2(p16.3, cog_diffpath, "cog_EdgeR_plot3", width = 10, height = 8)
save_plot2(p16.3, cog_diffpath, "cog_EdgeR_plot3", width = 10, height = 8)
# addWorksheet(cog_diff_wb, "EdgeR_results")  # 已合并到 write_sheet2
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
save_plot2(p17.1, cog_diffpath, "cog_ttest_plot1", width = 10, height = 8)
save_plot2(p17.1, cog_diffpath, "cog_ttest_plot1", width = 10, height = 8)
save_plot2(p17.2, cog_diffpath, "cog_ttest_plot2", width = 10, height = 8)
save_plot2(p17.2, cog_diffpath, "cog_ttest_plot2", width = 10, height = 8)
save_plot2(p17.3, cog_diffpath, "cog_ttest_plot3", width = 10, height = 8)
save_plot2(p17.3, cog_diffpath, "cog_ttest_plot3", width = 10, height = 8)
# addWorksheet(cog_diff_wb, "ttest_results")  # 已合并到 write_sheet2
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
save_plot2(p18.1, cog_diffpath, "cog_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.1, cog_diffpath, "cog_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.2, cog_diffpath, "cog_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.2, cog_diffpath, "cog_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.3, cog_diffpath, "cog_wilcox_plot3", width = 10, height = 8)
save_plot2(p18.3, cog_diffpath, "cog_wilcox_plot3", width = 10, height = 8)
# addWorksheet(cog_diff_wb, "wilcox_results")  # 已合并到 write_sheet2
writeData(cog_diff_wb, "wilcox_results", dat, rowNames = TRUE)

# 18 stamp.metf:stamp差异分析 ----
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
save_plot2(p19, cog_diffpath, "cog_stamp", width = 10, height = 8)
save_plot2(p19, cog_diffpath, "cog_stamp", width = 10, height = 8)
# addWorksheet(cog_diff_wb, "stamp_plot_data")  # 已合并到 write_sheet2
writeData(cog_diff_wb, "stamp_plot_data", dat1, rowNames = TRUE)
# addWorksheet(cog_diff_wb, "stamp_results")  # 已合并到 write_sheet2
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
save_plot2(prfcv, cog_biomarkerpath, "cog_rfcv", width = 10, height = 8)
save_plot2(prfcv, cog_biomarkerpath, "cog_rfcv", width = 10, height = 8)
# addWorksheet(cog_biomarker_wb, "rfcv_results")  # 已合并到 write_sheet2
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
save_plot2(p33.1, cog_biomarkerpath, "cog_ROC_plot1", width = 10, height = 8)
save_plot2(p33.1, cog_biomarkerpath, "cog_ROC_plot1", width = 10, height = 8)
save_plot2(p33.2, cog_biomarkerpath, "cog_ROC_plot2", width = 10, height = 8)
save_plot2(p33.2, cog_biomarkerpath, "cog_ROC_plot2", width = 10, height = 8)
# addWorksheet(cog_biomarker_wb, "ROC_results")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
save_plot2(p34.1, cog_biomarkerpath, "cog_loadingPCA", width = 10, height = 8)
save_plot2(p34.1, cog_biomarkerpath, "cog_loadingPCA", width = 10, height = 8)
# addWorksheet(cog_biomarker_wb, "loadingPCA_results")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "loadingPCA_results", dat, rowNames = TRUE)

#22 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = pst %>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
# addWorksheet(cog_biomarker_wb, "svm_AUC")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "svm_AUC", AUC, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "svm_importance")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#23 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
# addWorksheet(cog_biomarker_wb, "glm_AUC")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "glm_AUC", AUC, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "glm_importance")  # 已合并到 write_sheet2
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
# addWorksheet(cog_biomarker_wb, "xgboost_accuracy")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "xgboost_importance")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#25 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
# addWorksheet(cog_biomarker_wb, "lasso_accuracy")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "lasso_importance")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#26 decisiontree.metf: ----
library(rpart)
res =decisiontree.metf(ps=ps.card, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
# addWorksheet(cog_biomarker_wb, "decisiontree_accuracy")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "decisiontree_importance")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "decisiontree_importance", importance, rowNames = TRUE)

#27 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=ps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
# addWorksheet(cog_biomarker_wb, "naivebayes_accuracy")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "naivebayes_importance")  # 已合并到 write_sheet2
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
save_plot2(p35, cog_biomarkerpath, "cog_LDA", width = 10, height = 8)
save_plot2(p35, cog_biomarkerpath, "cog_LDA", width = 10, height = 8)
# addWorksheet(cog_biomarker_wb, "LDA_results")  # 已合并到 write_sheet2
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
save_plot2(p42.1, cog_biomarkerpath, "cog_randomforest_plot1", width = 10, height = 8)
save_plot2(p42.1, cog_biomarkerpath, "cog_randomforest_plot1", width = 10, height = 8)
save_plot2(p42.2, cog_biomarkerpath, "cog_randomforest_plot2", width = 10, height = 8)
save_plot2(p42.2, cog_biomarkerpath, "cog_randomforest_plot2", width = 10, height = 8)
save_plot2(p42.4, cog_biomarkerpath, "cog_randomforest_plot4", width = 10, height = 8)
save_plot2(p42.4, cog_biomarkerpath, "cog_randomforest_plot4", width = 10, height = 8)
# addWorksheet(cog_biomarker_wb, "randomforest_results")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "randomforest_results", dat, rowNames = TRUE)

#30 nnet.metf: 神经网络筛选特征功能  ------
library(nnet)
res =nnet.metf(ps=pst, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
# addWorksheet(cog_biomarker_wb, "nnet_accuracy")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "nnet_importance")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "nnet_importance", importance, rowNames = TRUE)

#31 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps =  pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
# addWorksheet(cog_biomarker_wb, "bagging_accuracy")  # 已合并到 write_sheet2
writeData(cog_biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(cog_biomarker_wb, "bagging_importance")  # 已合并到 write_sheet2
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
save_plot2(p0, cog_networkpath, "cog_network_main", width = 12, height = 10)
save_plot2(p0, cog_networkpath, "cog_network_main", width = 12, height = 10)

# 33 net_properties.4:网络属性计算 ----
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
# addWorksheet(cog_network_wb, "network_properties")  # 已合并到 write_sheet2
writeData(cog_network_wb, "network_properties", dat2, rowNames = TRUE)

# 34 netproperties.sample:单个样本的网络属性 ----
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
# addWorksheet(cog_network_wb, "sample_network_properties")  # 已合并到 write_sheet2
writeData(cog_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

# 35 node_properties:计算节点属性 ----
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
# addWorksheet(cog_network_wb, "node_properties")  # 已合并到 write_sheet2
writeData(cog_network_wb, "node_properties", nodepro2, rowNames = TRUE)

# 36 module.compare.net.pip:网络显著性比较 ----
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
# addWorksheet(cog_network_wb, "network_comparison")  # 已合并到 write_sheet2
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

# 1 alpha.metf:6种alpha多样性计算 ----
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
save_plot2(p1_1, vfdb_diversitypath, "vfdb_alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_1, vfdb_diversitypath, "vfdb_alpha_diversity_box", width = 12, height = 8)
save_plot2(p1_2, vfdb_diversitypath, "vfdb_alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_2, vfdb_diversitypath, "vfdb_alpha_diversity_bar", width = 12, height = 8)
save_plot2(p1_3, vfdb_diversitypath, "vfdb_alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_3, vfdb_diversitypath, "vfdb_alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_0, vfdb_diversitypath, "vfdb_alpha_diversity_violin", width = 12, height = 8)
save_plot2(p1_0, vfdb_diversitypath, "vfdb_alpha_diversity_violin", width = 12, height = 8)
# addWorksheet(vfdb_diversity_wb, "alpha_diversity_data")  # 已合并到 write_sheet2
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
save_plot2(p2_1, vfdb_diversitypath, "vfdb_rarefaction_individual", width = 10, height = 8)
save_plot2(p2_1, vfdb_diversitypath, "vfdb_rarefaction_individual", width = 10, height = 8)
save_plot2(p2_2, vfdb_diversitypath, "vfdb_rarefaction_group", width = 10, height = 8)
save_plot2(p2_2, vfdb_diversitypath, "vfdb_rarefaction_group", width = 10, height = 8)
save_plot2(p2_3, vfdb_diversitypath, "vfdb_rarefaction_group_sd", width = 10, height = 8)
save_plot2(p2_3, vfdb_diversitypath, "vfdb_rarefaction_group_sd", width = 10, height = 8)
# addWorksheet(vfdb_diversity_wb, "rarefaction_data")  # 已合并到 write_sheet2
writeData(vfdb_diversity_wb, "rarefaction_data", raretab, rowNames = TRUE)

# function beta diversity -----
# 创建VFDB beta多样性分析目录
vfdb_betapath <- file.path(metapath, "vfdb", "beta")
dir.create(vfdb_betapath, recursive = TRUE)

# 创建VFDB beta多样性分析总表
vfdb_beta_wb <- createWorkbook()

# 3 ordinate.metf:排序分析 ----
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
save_plot2(p3_1, vfdb_betapath, "vfdb_ordination_basic", width = 10, height = 8)
save_plot2(p3_1, vfdb_betapath, "vfdb_ordination_basic", width = 10, height = 8)
save_plot2(p3_2, vfdb_betapath, "vfdb_ordination_labeled", width = 10, height = 8)
save_plot2(p3_2, vfdb_betapath, "vfdb_ordination_labeled", width = 10, height = 8)
save_plot2(p3_3, vfdb_betapath, "vfdb_ordination_refined", width = 10, height = 8)
save_plot2(p3_3, vfdb_betapath, "vfdb_ordination_refined", width = 10, height = 8)
# addWorksheet(vfdb_beta_wb, "ordination_data")  # 已合并到 write_sheet2
writeData(vfdb_beta_wb, "ordination_data", plotdata, rowNames = TRUE)

# 4 MetaTest.metf:群落功能差异检测 ----
dat1 = MetaTest.metf(ps = ps, Micromet = "adonis", dist = "bray")
dat1

# 保存群落功能差异检测结果
# addWorksheet(vfdb_beta_wb, "metatest_results")  # 已合并到 write_sheet2
writeData(vfdb_beta_wb, "metatest_results", dat1, rowNames = TRUE)

# 5 pairMetaTest.metf:两两分组群落功能差异检测 ----
dat2 = pairMetaTest.metf(ps = ps, Micromet = "MRPP", dist = "bray")
dat2

# 保存两两分组群落功能差异检测结果
# addWorksheet(vfdb_beta_wb, "pair_metatest_results")  # 已合并到 write_sheet2
writeData(vfdb_beta_wb, "pair_metatest_results", dat2, rowNames = TRUE)

# 6 mantal.metf：群落功能差异检测普鲁士分析#------
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
save_plot2(p3_7, vfdb_betapath, "vfdb_mantel_test", width = 10, height = 8)
save_plot2(p3_7, vfdb_betapath, "vfdb_mantel_test", width = 10, height = 8)
# addWorksheet(vfdb_beta_wb, "mantel_results")  # 已合并到 write_sheet2
writeData(vfdb_beta_wb, "mantel_results", data, rowNames = TRUE)

#7 cluster.metf:样品聚类 ----
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
save_plot2(p4, vfdb_betapath, "vfdb_cluster_heatmap", width = 12, height = 10)
save_plot2(p4, vfdb_betapath, "vfdb_cluster_heatmap", width = 12, height = 10)
save_plot2(p4_1, vfdb_betapath, "vfdb_cluster_dendrogram1", width = 10, height = 8)
save_plot2(p4_1, vfdb_betapath, "vfdb_cluster_dendrogram1", width = 10, height = 8)
save_plot2(p4_2, vfdb_betapath, "vfdb_cluster_dendrogram2", width = 10, height = 8)
save_plot2(p4_2, vfdb_betapath, "vfdb_cluster_dendrogram2", width = 10, height = 8)
# addWorksheet(vfdb_beta_wb, "cluster_results")  # 已合并到 write_sheet2
writeData(vfdb_beta_wb, "cluster_results", dat, rowNames = TRUE)

#8 Micro_tern.metf: 三元图展示功能----
res = Micro_tern.metf(ps=ps %>% filter_OTU_ps(500),color = "Function"  )
p15 = res[[1]]
p15
dat =  res[[2]]
dat

# 保存三元图结果
save_plot2(p15, vfdb_betapath, "vfdb_ternary_plot", width = 10, height = 8)
save_plot2(p15, vfdb_betapath, "vfdb_ternary_plot", width = 10, height = 8)
# addWorksheet(vfdb_beta_wb, "ternary_data")  # 已合并到 write_sheet2
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
save_plot2(p4_1, vfdb_classificationpath, "vfdb_barplot_main", width = 12, height = 8)
save_plot2(p4_1, vfdb_classificationpath, "vfdb_barplot_main", width = 12, height = 8)
save_plot2(p4_2, vfdb_classificationpath, "vfdb_barplot_secondary", width = 12, height = 8)
save_plot2(p4_2, vfdb_classificationpath, "vfdb_barplot_secondary", width = 12, height = 8)
# addWorksheet(vfdb_classification_wb, "barplot_data")  # 已合并到 write_sheet2
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
save_plot2(p10.1, vfdb_classificationpath, "vfdb_venn_diagram", width = 10, height = 8)
save_plot2(p10.1, vfdb_classificationpath, "vfdb_venn_diagram", width = 10, height = 8)
save_plot2(p10.2, vfdb_classificationpath, "vfdb_upset_plot", width = 12, height = 8)
save_plot2(p10.2, vfdb_classificationpath, "vfdb_upset_plot", width = 12, height = 8)
# addWorksheet(vfdb_classification_wb, "venn_upset_data")  # 已合并到 write_sheet2
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
save_plot2(p7_1, vfdb_classificationpath, "vfdb_venn_detail_bar", width = 12, height = 8)
save_plot2(p7_1, vfdb_classificationpath, "vfdb_venn_detail_bar", width = 12, height = 8)
save_plot2(p7_2, vfdb_classificationpath, "vfdb_venn_detail_alluvial", width = 12, height = 8)
save_plot2(p7_2, vfdb_classificationpath, "vfdb_venn_detail_alluvial", width = 12, height = 8)
# addWorksheet(vfdb_classification_wb, "venn_detail_data")  # 已合并到 write_sheet2
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
save_plot2(p13.1, vfdb_classificationpath, "vfdb_flower_plot", width = 10, height = 10)
save_plot2(p13.1, vfdb_classificationpath, "vfdb_flower_plot", width = 10, height = 10)
# addWorksheet(vfdb_classification_wb, "flower_data")  # 已合并到 write_sheet2
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
save_plot2(p14, vfdb_classificationpath, "vfdb_venn_network", width = 10, height = 8)
save_plot2(p14, vfdb_classificationpath, "vfdb_venn_network", width = 10, height = 8)
# addWorksheet(vfdb_classification_wb, "venn_network_data")  # 已合并到 write_sheet2
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
save_plot2(p15.1, vfdb_diffpath, "vfdb_DESep2_plot1", width = 10, height = 8)
save_plot2(p15.1, vfdb_diffpath, "vfdb_DESep2_plot1", width = 10, height = 8)
save_plot2(p15.2, vfdb_diffpath, "vfdb_DESep2_plot2", width = 10, height = 8)
save_plot2(p15.2, vfdb_diffpath, "vfdb_DESep2_plot2", width = 10, height = 8)
save_plot2(p15.3, vfdb_diffpath, "vfdb_DESep2_plot3", width = 10, height = 8)
save_plot2(p15.3, vfdb_diffpath, "vfdb_DESep2_plot3", width = 10, height = 8)
# addWorksheet(vfdb_diff_wb, "DESep2_results")  # 已合并到 write_sheet2
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
save_plot2(p16.1, vfdb_diffpath, "vfdb_EdgeR_plot1", width = 10, height = 8)
save_plot2(p16.1, vfdb_diffpath, "vfdb_EdgeR_plot1", width = 10, height = 8)
save_plot2(p16.2, vfdb_diffpath, "vfdb_EdgeR_plot2", width = 10, height = 8)
save_plot2(p16.2, vfdb_diffpath, "vfdb_EdgeR_plot2", width = 10, height = 8)
save_plot2(p16.3, vfdb_diffpath, "vfdb_EdgeR_plot3", width = 10, height = 8)
save_plot2(p16.3, vfdb_diffpath, "vfdb_EdgeR_plot3", width = 10, height = 8)
# addWorksheet(vfdb_diff_wb, "EdgeR_results")  # 已合并到 write_sheet2
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
save_plot2(p17.1, vfdb_diffpath, "vfdb_ttest_plot1", width = 10, height = 8)
save_plot2(p17.1, vfdb_diffpath, "vfdb_ttest_plot1", width = 10, height = 8)
save_plot2(p17.2, vfdb_diffpath, "vfdb_ttest_plot2", width = 10, height = 8)
save_plot2(p17.2, vfdb_diffpath, "vfdb_ttest_plot2", width = 10, height = 8)
save_plot2(p17.3, vfdb_diffpath, "vfdb_ttest_plot3", width = 10, height = 8)
save_plot2(p17.3, vfdb_diffpath, "vfdb_ttest_plot3", width = 10, height = 8)
# addWorksheet(vfdb_diff_wb, "ttest_results")  # 已合并到 write_sheet2
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
save_plot2(p18.1, vfdb_diffpath, "vfdb_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.1, vfdb_diffpath, "vfdb_wilcox_plot1", width = 10, height = 8)
save_plot2(p18.2, vfdb_diffpath, "vfdb_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.2, vfdb_diffpath, "vfdb_wilcox_plot2", width = 10, height = 8)
save_plot2(p18.3, vfdb_diffpath, "vfdb_wilcox_plot3", width = 10, height = 8)
save_plot2(p18.3, vfdb_diffpath, "vfdb_wilcox_plot3", width = 10, height = 8)
# addWorksheet(vfdb_diff_wb, "wilcox_results")  # 已合并到 write_sheet2
writeData(vfdb_diff_wb, "wilcox_results", dat, rowNames = TRUE)

# 18 stamp.metf:stamp差异分析 ----
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
save_plot2(p19, vfdb_diffpath, "vfdb_stamp", width = 10, height = 8)
save_plot2(p19, vfdb_diffpath, "vfdb_stamp", width = 10, height = 8)
# addWorksheet(vfdb_diff_wb, "stamp_plot_data")  # 已合并到 write_sheet2
writeData(vfdb_diff_wb, "stamp_plot_data", dat1, rowNames = TRUE)
# addWorksheet(vfdb_diff_wb, "stamp_results")  # 已合并到 write_sheet2
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
save_plot2(prfcv, vfdb_biomarkerpath, "vfdb_rfcv", width = 10, height = 8)
save_plot2(prfcv, vfdb_biomarkerpath, "vfdb_rfcv", width = 10, height = 8)
# addWorksheet(vfdb_biomarker_wb, "rfcv_results")  # 已合并到 write_sheet2
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
save_plot2(p33.1, vfdb_biomarkerpath, "vfdb_ROC_plot1", width = 10, height = 8)
save_plot2(p33.1, vfdb_biomarkerpath, "vfdb_ROC_plot1", width = 10, height = 8)
save_plot2(p33.2, vfdb_biomarkerpath, "vfdb_ROC_plot2", width = 10, height = 8)
save_plot2(p33.2, vfdb_biomarkerpath, "vfdb_ROC_plot2", width = 10, height = 8)
# addWorksheet(vfdb_biomarker_wb, "ROC_results")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "ROC_results", dat, rowNames = TRUE)

#21 loadingPCA.metf:载荷矩阵筛选特征功能------
res = loadingPCA.metf(ps = ps,Top = 20)
p34.1 = res[[1]]
p34.1
dat = res[[2]]
dat

# 保存PCA载荷结果
save_plot2(p34.1, vfdb_biomarkerpath, "vfdb_loadingPCA", width = 10, height = 8)
save_plot2(p34.1, vfdb_biomarkerpath, "vfdb_loadingPCA", width = 10, height = 8)
# addWorksheet(vfdb_biomarker_wb, "loadingPCA_results")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "loadingPCA_results", dat, rowNames = TRUE)

#22 svm.metf:svm筛选特征功能 ----
res <- svm.metf(ps = ps%>% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存SVM结果
# addWorksheet(vfdb_biomarker_wb, "svm_AUC")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "svm_AUC", AUC, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "svm_importance")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "svm_importance", importance, rowNames = TRUE)

#23 glm.metf:glm筛选特征功能----
res <- glm.metf(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance

# 保存GLM结果
# addWorksheet(vfdb_biomarker_wb, "glm_AUC")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "glm_AUC", AUC, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "glm_importance")  # 已合并到 write_sheet2
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
# addWorksheet(vfdb_biomarker_wb, "xgboost_accuracy")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "xgboost_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "xgboost_importance")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "xgboost_importance", importance, rowNames = TRUE)

#25 lasso.metf: lasso筛选特征功能----
library(glmnet)
res =lasso.metf(ps =  ps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Lasso结果
# addWorksheet(vfdb_biomarker_wb, "lasso_accuracy")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "lasso_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "lasso_importance")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "lasso_importance", importance, rowNames = TRUE)

#26 decisiontree.micro: ----
library(rpart)
res =decisiontree.metf(ps=ps.card, top = 50, seed = 6358, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存决策树结果
# addWorksheet(vfdb_biomarker_wb, "decisiontree_accuracy")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "decisiontree_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "decisiontree_importance")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "decisiontree_importance", importance, rowNames = TRUE)

#27 naivebayes.metf: bayes筛选特征功能----
res = naivebayes.metf(ps=pst, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存朴素贝叶斯结果
# addWorksheet(vfdb_biomarker_wb, "naivebayes_accuracy")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "naivebayes_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "naivebayes_importance")  # 已合并到 write_sheet2
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
save_plot2(p35, vfdb_biomarkerpath, "vfdb_LDA", width = 10, height = 8)
save_plot2(p35, vfdb_biomarkerpath, "vfdb_LDA", width = 10, height = 8)
# addWorksheet(vfdb_biomarker_wb, "LDA_results")  # 已合并到 write_sheet2
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
save_plot2(p42.1, vfdb_biomarkerpath, "vfdb_randomforest_plot1", width = 10, height = 8)
save_plot2(p42.1, vfdb_biomarkerpath, "vfdb_randomforest_plot1", width = 10, height = 8)
save_plot2(p42.2, vfdb_biomarkerpath, "vfdb_randomforest_plot2", width = 10, height = 8)
save_plot2(p42.2, vfdb_biomarkerpath, "vfdb_randomforest_plot2", width = 10, height = 8)
save_plot2(p42.4, vfdb_biomarkerpath, "vfdb_randomforest_plot4", width = 10, height = 8)
save_plot2(p42.4, vfdb_biomarkerpath, "vfdb_randomforest_plot4", width = 10, height = 8)
# addWorksheet(vfdb_biomarker_wb, "randomforest_results")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "randomforest_results", dat, rowNames = TRUE)

#30 nnet.metf: 神经网络筛选特征功能  ------
library(nnet)
res =nnet.metf(ps=ps, top = 100, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存神经网络结果
# addWorksheet(vfdb_biomarker_wb, "nnet_accuracy")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "nnet_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "nnet_importance")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "nnet_importance", importance, rowNames = TRUE)

#31 bagging.metf : Bootstrap Aggregating筛选特征功能 ------
library(ipred)
res =bagging.metf(ps = ps, top = 20, seed = 1010, k = 5)
accuracy = res[[1]]
accuracy
importance = res[[2]]
importance

# 保存Bagging结果
# addWorksheet(vfdb_biomarker_wb, "bagging_accuracy")  # 已合并到 write_sheet2
writeData(vfdb_biomarker_wb, "bagging_accuracy", accuracy, rowNames = TRUE)
# addWorksheet(vfdb_biomarker_wb, "bagging_importance")  # 已合并到 write_sheet2
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
save_plot2(p0, vfdb_networkpath, "vfdb_network_main", width = 12, height = 10)
save_plot2(p0, vfdb_networkpath, "vfdb_network_main", width = 12, height = 10)

# 33 net_properties.4:网络属性计算 ----
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
# addWorksheet(vfdb_network_wb, "network_properties")  # 已合并到 write_sheet2
writeData(vfdb_network_wb, "network_properties", dat2, rowNames = TRUE)

# 34 netproperties.sample:单个样本的网络属性 ----
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
# addWorksheet(vfdb_network_wb, "sample_network_properties")  # 已合并到 write_sheet2
writeData(vfdb_network_wb, "sample_network_properties", dat3, rowNames = TRUE)

# 35 node_properties:计算节点属性 ----
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
# addWorksheet(vfdb_network_wb, "node_properties")  # 已合并到 write_sheet2
writeData(vfdb_network_wb, "node_properties", nodepro2, rowNames = TRUE)

# 36 module.compare.net.pip:网络显著性比较 ----
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
# addWorksheet(vfdb_network_wb, "network_comparison")  # 已合并到 write_sheet2
writeData(vfdb_network_wb, "network_comparison", res, rowNames = TRUE)

# 保存网络分析总表
saveWorkbook(vfdb_network_wb, file.path(vfdb_networkpath, "vfdb_network_analysis.xlsx"), overwrite = TRUE)
