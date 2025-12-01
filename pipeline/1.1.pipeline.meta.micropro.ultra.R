# ===================== 宏基因组物种分析（基于 EasyMultiOmics） =====================
rm(list = ls())

## ===================== 0. 基础设置 & 加载包 =====================
# BiocManager::install("MicrobiotaProcess")
library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)
library(ggrepel)
library(ggsci)
library(openxlsx)
library(fs)

## ===================== 1. 读入 phyloseq 对象 & 基础信息 =====================

## ===================== 1. 读入 phyloseq 对象 =====================
# 方式一：交互选择
# file_path <- tcltk::tk_choose.files(
#   caption = "请选择 ps_ITS.rds 文件",
#   multi   = FALSE,
#   filters = matrix(c("RDS files", ".rds",
#                      "All files", "*"), ncol = 2, byrow = TRUE)
# )
# 方式二：直接指定路径
# file_path <- "./data/ps_16s.rds"
# ps.16s    <- readRDS(file_path)
# ps.16s

file_path <- "./data/ps.rhi.rds"
ps.micro    <- readRDS(file_path)
ps.micro

map <- sample_data(ps.16s)
head(map)
# 假设 ps.micro 是你的 phyloseq 对象
mg_path <- create_amplicon_result_dir(ps.micro, include_time =TRUE)
mg_path

map <- sample_data(ps.micro)
head(map)

# 简单检查
sample_sums(ps.micro)
phyloseq::tax_table(ps.micro) %>% head()

# 分组数量与顺序
gnum       <- phyloseq::sample_data(ps.micro)$Group %>% unique() %>% length()
axis_order <- phyloseq::sample_data(ps.micro)$Group %>% unique()

# 分组配色
col.g <- get_group_cols(axis_order, palette = "npg")
col.g
scales::show_col(col.g)

# 主题 & 颜色
package.amp()
res     <- theme_my(ps.micro)
mytheme1 <- res[[1]]
mytheme2 <- res[[2]]
colset1  <- res[[3]]
colset2  <- res[[4]]
colset3  <- res[[5]]
colset4  <- res[[6]]



## ===================== 2. Alpha 多样性分析 =====================

mg_alpha_path <- file.path(mg_path, "01_alpha_diversity")
dir.create(mg_alpha_path, recursive = TRUE, showWarnings = FALSE)

alpha_xlsx_path <- file.path(mg_alpha_path, "alpha_diversity_results.xlsx")

if (file.exists(alpha_xlsx_path)) {
  mg_alpha_wb <- openxlsx::loadWorkbook(alpha_xlsx_path)
} else {
  mg_alpha_wb <- openxlsx::createWorkbook()
}

### ---------- 2.1 alpha.metm：Shannon/Chao1 等 ----------

all.alpha <- c(
  "Shannon", "Inv_Simpson", "Pielou_evenness",
  "Simpson_evenness", "Richness", "Chao1", "ACE"
)

tab <- alpha.metm(ps = ps.micro, group = "Group")
head(tab)

data_alpha <- cbind(
  data.frame(ID = 1:length(tab$Group), group = tab$Group),
  tab[all.alpha]
)
data_alpha$ID <- as.character(data_alpha$ID)
head(data_alpha)

# ANOVA + 多重比较
result_alpha <- MuiaovMcomper2(data = data_alpha, num = 3:6)

# 盒线图
res_box <- EasyMultiOmics::FacetMuiPlotresultBox(
  data   = data_alpha,
  num    = 3:6,
  result = result_alpha,
  sig_show = "abc",
  ncol     = 4,
  width    = 0.4
)
p1_1 <- res_box[[1]] +
  scale_fill_manual(values = col.g) +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL))

# 柱状图
res_bar <- EasyMultiOmics::FacetMuiPlotresultBar(
  data    = data_alpha,
  num     = 3:6,
  result  = result_alpha,
  sig_show = "abc",
  ncol     = 4,
  mult.y   = 0.3
)
p1_2 <- res_bar[[1]] +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

# 盒线 + 柱状
res_boxbar <- EasyMultiOmics::FacetMuiPlotReBoxBar(
  data    = data_alpha,
  num     = 3:6,
  result  = result_alpha,
  sig_show = "abc",
  ncol     = 4,
  mult.y   = 0.3,
  lab.yloc = 1.1
)
p1_3 <- res_boxbar[[1]] +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

# 小提琴图
lab_alpha <- res_box[[2]] %>% dplyr::distinct(group, name, .keep_all = TRUE)
p1_0 <- res_box[[2]] %>%
  ggplot(aes(x = group, y = dd)) +
  geom_violin(alpha = 1, aes(fill = group)) +
  geom_jitter(aes(color = group),
              position = position_jitter(0.17),
              size = 3, alpha = 0.5) +
  labs(x = "", y = "") +
  facet_wrap(. ~ name, scales = "free_y", ncol = 4) +
  geom_text(data = lab_alpha, aes(x = group, y = y, label = stat), alpha = 0.6) +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g) +
  ggplot2::scale_color_manual(values = col.g)

# ---- 保存 Alpha 图 ----
save_plot2(p1_1, mg_alpha_path, "alpha_diversity_box",    width = 12, height = 6)
save_plot2(p1_2, mg_alpha_path, "alpha_diversity_bar",    width = 12, height = 6)
save_plot2(p1_3, mg_alpha_path, "alpha_diversity_boxbar", width = 12, height = 6)
save_plot2(p1_0, mg_alpha_path, "alpha_diversity_violin", width = 12, height = 6)

# ---- 保存 Alpha 表 ----
write_sheet2(mg_alpha_wb, "alpha_diversity_data", data_alpha)
write_sheet2(mg_alpha_wb, "alpha_diversity_stat", result_alpha)
openxlsx::saveWorkbook(mg_alpha_wb, alpha_xlsx_path, overwrite = TRUE)

### ---------- 2.2 alpha_rare.metm：Alpha 稀释曲线 ----------

library(microbiome)
library(vegan)

rare_step <- mean(phyloseq::sample_sums(ps.micro)) / 20

res_rare <- alpha_rare.metm(
  ps     = ps.micro,
  group  = "Group",
  method = "Richness",
  start  = 100,
  step   = rare_step
)

p2_1 <- res_rare[[1]] +
  scale_color_manual(values = col.g) +
  theme_nature()

raretab <- res_rare[[2]]

p2_2 <- res_rare[[3]] +
  scale_color_manual(values = col.g) +
  theme_nature()

p2_3 <- res_rare[[4]] +
  scale_color_manual(values = col.g) +
  theme_nature()

# ---- 保存稀释曲线图 ----
save_plot2(p2_1, mg_alpha_path, "alpha_rarefaction_individual", width = 9,  height = 7)
save_plot2(p2_2, mg_alpha_path, "alpha_rarefaction_group",      width = 10, height = 8)
save_plot2(p2_3, mg_alpha_path, "alpha_rarefaction_group_sd",   width = 10, height = 8)

# ---- 保存稀释数据 ----
write_sheet2(mg_alpha_wb, "rarefaction_data", raretab)
openxlsx::saveWorkbook(mg_alpha_wb, alpha_xlsx_path, overwrite = TRUE)


## ===================== 3. Beta 多样性分析 =====================

mg_beta_path <- file.path(mg_path, "02_beta_diversity")
dir.create(mg_beta_path, recursive = TRUE, showWarnings = FALSE)

beta_xlsx_path <- file.path(mg_beta_path, "beta_diversity_results.xlsx")

if (file.exists(beta_xlsx_path)) {
  mg_beta_wb <- openxlsx::loadWorkbook(beta_xlsx_path)
} else {
  mg_beta_wb <- openxlsx::createWorkbook()
}

### ---------- 3.1 ordinate.metm：PCoA 排序 ----------

library(ggsci)

res_ord <- ordinate.metm(
  ps      = ps.micro,
  group   = c("Group", "Group2"),
  dist    = "bray",
  method  = "PCoA",
  Micromet = "anosim",
  pvalue.cutoff = 0.05
)

p3_1 <- res_ord[[1]] +
  scale_fill_manual(values = col.g) +
  scale_color_manual(values = col.g, guide = "none") +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

p3_2 <- res_ord[[3]] +
  scale_fill_manual(values = col.g) +
  scale_color_manual(values = col.g, guide = "none") +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

plotdata <- res_ord[[2]]
head(plotdata)

# 精修图：群心 + “蜘蛛线”
cent <- aggregate(cbind(x, y) ~ Group, data = plotdata, FUN = mean)
segs <- merge(
  plotdata,
  setNames(cent, c("Group", "oNMDS1", "oNMDS2")),
  by = "Group",
  sort = FALSE
)

p3_3 <- p3_1 +
  geom_segment(
    data    = segs,
    mapping = aes(xend = oNMDS1, yend = oNMDS2, color = Group),
    show.legend = FALSE
  ) +
  geom_point(
    data = cent,
    mapping = aes(x = x, y = y),
    size = 5, pch = 24, color = "black", fill = "yellow"
  ) +
  scale_fill_manual(values = col.g) +
  scale_color_manual(values = col.g, guide = "none") +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

# ---- 保存 PCoA 图 ----
save_plot2(p3_1, mg_beta_path, "pcoa_basic",   width = 12, height = 10)
save_plot2(p3_2, mg_beta_path, "pcoa_labeled", width = 12, height = 10)
save_plot2(p3_3, mg_beta_path, "pcoa_refined", width = 12, height = 10)

# ---- 保存排序数据 ----
write_sheet2(mg_beta_wb, "ordination_data", plotdata)
write_sheet2(mg_beta_wb, "ordination_cent", cent)
write_sheet2(mg_beta_wb, "ordination_segs", segs)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)

### ---------- 3.2 MicroTest.metm：整体群落差异 (adonis) ----------

dat_adonis <- MicroTest.metm(ps = ps.micro, Micromet = "adonis", dist = "bray")
write_sheet2(mg_beta_wb, "adonis_results", dat_adonis)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)

### ---------- 3.3 pairMicroTest.metm2：两两分组差异 (MRPP) ----------

dat_pair <- pairMicroTest.metm2(ps = ps.micro, Micromet = "MRPP", dist = "bray")
write_sheet2(mg_beta_wb, "pairwise_MRPP", dat_pair)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)

### ---------- 3.4 mantal.metm：Mantel 分析 ----------
# 计算 Group 分组及两两组合数量，并按每4个组合算一个参数
calc_group_pair_param <- function(ps,
                                  group_col = "Group",
                                  max_pairs_per_param = 4) {
  # 1. 取 sample_data
  sd <- phyloseq::sample_data(ps)

  if (!group_col %in% colnames(sd)) {
    stop("列 '", group_col, "' 在 sample_data 中不存在。")
  }

  # 2. 提取 Group 列的唯一分组，去掉 NA
  groups <- unique(as.vector(sd[[group_col]]))
  groups <- groups[!is.na(groups)]

  n_group <- length(groups)

  # 3. 计算两两组合数量：n * (n - 1) / 2
  if (n_group < 2) {
    n_pairs <- 0L
  } else {
    n_pairs <- as.integer(n_group * (n_group - 1) / 2)
  }

  # 4. 每 4 个组合算 1 个参数，不足 4 个也算 1 个
  if (n_pairs == 0L) {
    n_param <- 0L
  } else {
    n_param <- as.integer(ceiling(n_pairs / max_pairs_per_param))
  }

  # 返回一个列表，方便后续使用
  list(
    groups   = groups,   # 分组名称
    n_group  = n_group,  # 分组个数
    n_pairs  = n_pairs,  # 两两组合总数
    n_param  = n_param   # 参数个数（每4个组合一个，不足4也加1）
  )
}
res <- calc_group_pair_param(ps.micro)
# res$n_group  # 有几个 Group
# res$n_pairs  # 两两组合数量
# res$n_param  # 按“每4个组合一个，不足4也算一个”得到的参数数目
# res$groups   # 分组名称向量

res_mantel <- mantal.metm(
  ps     = ps.micro,
  method = "spearman",
  group  = "Group",
  ncol   = gnum,
  nrow   = res$n_param
)
data_mantel <- res_mantel[[1]]
p3_7        <- res_mantel[[2]]

save_plot2(p3_7, mg_beta_path, "mantel_plot", width = 10, height = 7)
write_sheet2(mg_beta_wb, "mantel_results", data_mantel)
openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)

### ---------- 3.5 cluster_metm：样品聚类 ----------

res_clust <- cluster_metm(
  ps             = ps.micro,
  hcluter_method = "complete",
  dist           = "bray",
  cuttree        = 3,
  row_cluster    = TRUE,
  col_cluster    = TRUE
)

p4   <- res_clust[[1]] + scale_fill_manual(values = col.g)
p4_1 <- res_clust[[2]]
p4_2 <- res_clust[[3]]
dat_c <- res_clust[[4]]

save_plot2(p4,   mg_beta_path, "cluster_heatmap",   width = 10, height = 8)
save_plot2(p4_1, mg_beta_path, "cluster_dendro1",   width = 10, height = 8)
save_plot2(p4_2, mg_beta_path, "cluster_scatter",   width = 10, height = 8)

write_sheet2(mg_beta_wb, "cluster_matrix", dat_c)
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

### ---------- 4.1 Ven.Upset.metm：共有/特有 OTU/ASV ----------

res_venn <- Ven.Upset.metm(
  ps    = ps.micro,
  group = "Group",
  N     = 0.5,
  size  = 3
)

p10_1 <- res_venn[[1]] + theme_void()+ scale_fill_manual(values = col.g)
p10_2 <- res_venn[[2]]
p10_3 <- res_venn[[4]]
dat10 <- res_venn[[3]]

save_plot2(p10_1, mg_comp_path, "venn_diagram",  width = 10, height = 8)
save_plot2(p10_2, mg_comp_path, "upset_plot",    width = 10, height = 8)
save_plot2(p10_3, mg_comp_path, "venn_extra",    width = 10, height = 8)

write_sheet2(mg_comp_wb, "venn_upset_data", dat10)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.2 VenSuper.metm：各部分物种组成 ----------

library(ggpubr)
library(agricolae)
library(reshape2)

sample_data(ps.micro)$Group %>% table()

res_vensuper <- VenSuper.metm(
  ps    = ps.micro,
  group = "Group",
  num   = 6
)

p7_1  <- res_vensuper[[1]]
p7_2  <- res_vensuper[[3]]
dat12 <- res_vensuper[[4]]
dat120 <- dplyr::bind_rows(dat12, .id = "Part")

save_plot2(p7_1, mg_comp_path, "vensuper_bar",     width = 10, height = 8)
save_plot2(p7_2, mg_comp_path, "vensuper_alluvial", width = 10, height = 8)

write_sheet2(mg_comp_wb, "vensuper_detail", dat120)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.3 ggflower.metm：花瓣图 ----------

res_flower <- ggflower.metm(
  ps        = ps.micro,
  group     = "ID",
  start     = 1,
  m1        = 2,
  a         = 0.3,
  b         = 1,
  lab.leaf  = 1,
  col.cir   = "yellow",
  N         = 0.5
)

p13_1 <- res_flower[[1]]
dat13 <- res_flower[[2]]

save_plot2(p13_1, mg_comp_path, "flower_plot", width = 10, height = 8)
write_sheet2(mg_comp_wb, "flower_data", dat13)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.4 ven.network.metm：Venn 网络 ----------

res_vnet <- ven.network.metm(
  ps   = ps.micro,
  N    = 0.5,
  fill = "Kingdom"
)
p14   <- res_vnet[[1]]
dat14 <- res_vnet[[2]]

save_plot2(p14, mg_comp_path, "venn_network", width = 12, height = 10)
write_sheet2(mg_comp_wb, "venn_network_data", dat14)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.5 Micro_tern.metm：三元图 ----------

ps1 <- ps.micro %>% filter_OTU_ps(500)
res_tern <- Micro_tern.metm(ps1)
p15     <- res_tern[[1]][[1]]
dat15   <- res_tern[[2]]

save_plot2(p15 + theme_classic(), mg_comp_path, "ternary_plot", width = 12, height = 10)
write_sheet2(mg_comp_wb, "ternary_data", dat15)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.6 多元 polygon 图（可选保存） ----------

p_poly1 <- ps_polygon_plot(ps.micro, group = "Group", taxrank = "Phylum")


save_plot2(p_poly1, mg_comp_path, "polygon_Group2_Phylum", width = 18, height = 12)
### ---------- 4.7 barMainplot.metm：堆积柱状图 --------

library(ggalluvial)

pst <- ps.micro %>% subset_taxa.wt("Genus", "Unassigned", TRUE)

res_barMain <- barMainplot.metm(
  ps    = pst,
  j     = "Genus",
  label = FALSE,
  sd    = FALSE,
  Top   = 10
)

p4_1 <- res_barMain[[1]] +
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

p4_2 <- res_barMain[[3]] +
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

dat16 <- res_barMain[[2]] %>%
  dplyr::group_by(Group, aa) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  as.data.frame()
colnames(dat16) <- c("Group", "Genus", "Abundance(%)")

save_plot2(p4_1, mg_comp_path, "barplot_samples", width = 10, height = 8)
save_plot2(p4_2, mg_comp_path, "barplot_groups",  width = 10, height = 8)

write_sheet2(mg_comp_wb, "barplot_raw_data",  res_barMain[[2]])
write_sheet2(mg_comp_wb, "barplot_summary",   dat16)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.8 cluMicro.bar.metm：聚类堆积柱状图 ----------

res_clubar <- cluMicro.bar.metm(
  dist           = "bray",
  ps             = ps.micro,
  j              = "Genus",
  Top            = 6,
  tran           = TRUE,
  hcluter_method = "complete",
  Group          = "Group",
  cuttree        = length(unique(phyloseq::sample_data(ps.micro)$Group))
)

p5_2    <- res_clubar[[2]]
p5_3    <- res_clubar[[3]]
p5_4    <- res_clubar[[4]]

clubardata <- res_clubar[[5]]

save_plot2(p5_2, mg_comp_path, "cluster_barplot_1", width = 10, height = 8)
save_plot2(p5_3, mg_comp_path, "cluster_barplot_2", width = 10, height = 8)
save_plot2(p5_4, mg_comp_path, "cluster_barplot_3", width = 10, height = 8)

write_sheet2(mg_comp_wb, "cluster_barplot_data", clubardata)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.9 stacked_bar_by_custom_class_v：高/中/低丰度代表菌 ----------

res_col <- stacked_bar_by_custom_class_v(
  ps         = ps.micro,
  rank       = "Genus",
  group      = "Group",
  thresholds = c(1, 0.1),   # 高(>=1%) / 中[0.1%,1%) / 低(<0.1%)
  n_each     = c(6, 6, 6),
  arrange    = "row"
)
# 这里只画图不强制保存，如需要保存可用：
save_plot2(res_col$combined, mg_comp_path, "stacked_bar_custom", width = 18, height = 8)

### ---------- 4.10 cir_barplot.metm & cir_barplot.metm2：环状堆积柱状图 ----------

library(ggtree)

res_cir <- cir_barplot.metm(
  ps             = ps.micro,
  Top            = 15,
  dist           = "bray",
  cuttree        = 3,
  hcluter_method = "complete"
)

p17  <- res_cir[[1]]
dat17 <- res_cir[[2]]

save_plot2(p17, mg_comp_path, "circular_barplot", width = 16, height = 12)
write_sheet2(mg_comp_wb, "circular_barplot_data", dat17)
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

# 高/中/低丰度分圈展示
res_cir2 <- cir_barplot.metm2(
  ps          = ps.micro,
  Top         = 20,
  rank        = "Genus",
  thresholds  = c(1, 0.8),
  n_each      = c(10, 8, 8),
  xlim_high   = c(0, 100),
  xlim_medium = c(0, 20),
  xlim_low    = c(0, 5),
  breaks_high   = seq(0, 100, 10),
  breaks_medium = seq(0, 20, 5),
  breaks_low    = seq(0, 5, 1),
  ring_pwidth   = c(2, 5, 5),
  base_offset   = 2,
  ring_gap      = 0,
  axis_pad      = 0,
  show_axes     = c(FALSE, FALSE, FALSE),
  legend_position = "right"
)
save_plot2(res_cir2$plot, mg_comp_path, "circular_barplot_multiring", width = 24, height = 20)

### ---------- 4.11 cir_plot.metm：和弦图 ----------

res_cirplot <- cir_plot.metm(ps = ps.micro, Top = 12, rank = 6)
mer_otu_mean <- res_cirplot[[1]]

# 用专用函数保存（与 16S 版本保持一致）
save_circlize_plot(
  cir_plot.metm(ps = ps.micro, Top = 12, rank = 6),
  out_dir = mg_comp_path,
  prefix  = "chord_plot",
  width   = 14,
  height  = 14,
  res     = 300
)

write_sheet2(mg_comp_wb, "chord_data", as.data.frame(mer_otu_mean))
openxlsx::saveWorkbook(mg_comp_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 4.12 Microheatmap.metm：相对丰度热图 ----------

heatnum <- 30

ps_tem <- ps.micro %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus")

id_heat <- ps.micro %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus") %>%
  ggClusterNet::filter_OTU_ps(100) %>%
  ggClusterNet::vegan_otu() %>%
  t() %>%
  as.data.frame() %>%
  rowCV %>%
  sort(decreasing = TRUE) %>%
  head(heatnum) %>%
  names()

res_heat <- Microheatmap.metm(
  ps_rela     = ps_tem,
  id          = id_heat,
  col_cluster = FALSE,
  row_cluster = FALSE,
  col1        = (ggsci::pal_gsea(alpha = 1))(12)
)

p24_1 <- res_heat[[1]]
p24_2 <- res_heat[[2]]
dat24 <- res_heat[[3]]

save_plot2(p24_1, mg_comp_path, "microbiome_heatmap_1", width = 12, height = 10)
save_plot2(p24_2, mg_comp_path, "microbiome_heatmap_2", width = 12, height = 10)

write_sheet2(mg_comp_wb, "heatmap_data", dat24)
write_sheet2(mg_comp_wb, "heatmap_selected_otus", data.frame(OTU_ID = id_heat))
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

### ---------- 5.1 EdgerSuper.metm：EdgeR 差异物种（火山图列表） ----------

res_edger <- EdgerSuper.metm(
  ps       = ps.micro,
  group    = "Group",
  artGroup = NULL,
  j        = "Species"
)

p25_1 <- res_edger[[1]][[1]]
p25_2 <- res_edger[[1]][[2]]
p25_3 <- res_edger[[1]][[3]]
dat25 <- res_edger[[2]]

save_plot2(p25_1, mg_diff_path, "edger_volcano_1", width = 10, height = 8)
save_plot2(p25_2, mg_diff_path, "edger_volcano_2", width = 10, height = 8)
save_plot2(p25_3, mg_diff_path, "edger_volcano_3", width = 10, height = 8)

write_sheet2(mg_diff_wb, "edger_results", dat25)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### ---------- 5.2 EdgerSuper2.metm：所有分组组合 EdgeR 竖表 ----------

res_edger2 <- EdgerSuper2.metm(
  ps       = ps.micro,
  group    = "Group",
  artGroup = NULL,
  j        = "Species"
)

write_sheet2(mg_diff_wb, "edger2_results", res_edger2)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### ---------- 5.3 DESep2Super.metm：DESeq2 差异物种 ----------

res_deseq <- DESep2Super.metm(
  ps       = ps.micro %>% ggClusterNet::filter_OTU_ps(500),
  group    = "Group",
  artGroup = NULL,
  j        = "Species"
)

p26_1 <- res_deseq[[1]][[1]]
p26_2 <- res_deseq[[1]][[2]]
p26_3 <- res_deseq[[1]][[3]]
dat26  <- res_deseq[[2]]

save_plot2(p26_1, mg_diff_path, "deseq2_volcano_1", width = 10, height = 8)
save_plot2(p26_2, mg_diff_path, "deseq2_volcano_2", width = 10, height = 8)
save_plot2(p26_3, mg_diff_path, "deseq2_volcano_3", width = 10, height = 8)

write_sheet2(mg_diff_wb, "deseq2_results", dat26)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### ---------- 5.4 edge_Manhattan.metm：曼哈顿图 ----------

res_manh <- edge_Manhattan.metm(
  ps     = ps.micro %>% ggClusterNet::filter_OTU_ps(500),
  pvalue = 0.05,
  lfc    = 0
)

p27_1 <- res_manh[[1]]
p27_2 <- res_manh[[2]]
p27_3 <- res_manh[[3]]
# 若有 res_manh[[4]] 可写表

save_plot2(p27_1, mg_diff_path, "manhattan_plot_1", width = 12, height = 10)
save_plot2(p27_2, mg_diff_path, "manhattan_plot_2", width = 12, height = 10)
save_plot2(p27_3, mg_diff_path, "manhattan_plot_3", width = 12, height = 10)

# 如果有数据：
# write_sheet2(mg_diff_wb, "manhattan_data", res_manh[[4]])
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### ---------- 5.5 stemp_diff.metm：STAMP 风格差异丰度 ----------

map <- sample_data(ps.micro)
allgroup <- combn(unique(map$Group), 2)
plot_list_28 <- vector("list", ncol(allgroup))

for (i in seq_len(ncol(allgroup))) {
  ps_sub <- phyloseq::subset_samples(ps.micro, Group %in% allgroup[, i])
  p_tmp  <- stemp_diff.metm(ps = ps_sub, Top = 20, ranks = 6)
  plot_list_28[[i]] <- p_tmp

  save_plot2(
    p_tmp,
    mg_diff_path,
    paste0("stamp_plot_", i),
    width  = 10,
    height = 8
  )
}
# 如 stemp_diff.metm 返回数据，可同样写表

### ---------- 5.6 Mui.Group.volcano.metm：聚类火山图（所有组） ----------

res29_input <- EdgerSuper2.metm(
  ps       = ps.micro %>% filter_OTU_ps(500),
  group    = "Group",
  artGroup = NULL,
  j        = "OTU"
)

res29 <- Mui.Group.volcano.metm(res = res29_input)


p29_1 <- res29[[1]]
p29_2 <- res29[[2]]
dat29  <- res29[[3]]

save_plot2(p29_1, mg_diff_path, "cluster_volcano_1", width = 10, height = 8)
save_plot2(p29_2, mg_diff_path, "cluster_volcano_2", width = 10, height = 8)

write_sheet2(mg_diff_wb, "cluster_volcano_results", dat29)
openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)

### ---------- 5.7 diff_all_methods + tree（可选，放到 if(FALSE) 防止默认跑） ----------

if (FALSE) {
  library(ape)
  library(ggtree)
  library(Maaslin2)
  library(metagenomeSeq)
  library(glmnet)

  # 任选一对分组做“全差异方法整合”
  map_tmp <- phyloseq::sample_data(ps.micro)
  id.g    <- map_tmp$Group %>% unique() %>% as.character() %>% combn(2)
  i       <- 2

  ps.cs <- ps.micro %>% subset_samples.wt("Group", id.g[, i])

  res_all <- diff_all_methods(
    ps.cs %>% filter_OTU_ps(400),
    group = "Group",
    alpha = 0.05
  )

  # 树 + 全差异环图
  res_tree <- plot_resall_with_tree(
    res_all,
    ps.cs,
    topN       = 175,
    type       = "3/4",
    inner_gap  = -20,
    label_size = 2.8
  )

  p_tree <- res_tree[[1]]
  save_plot2(p_tree, mg_diff_path, "diff_all_methods_tree", width = 12, height = 12)

  # 摘要表：每个方法是否显著
  details_tab <- res_all$details %>%
    dplyr::mutate(sig = ifelse(!is.na(adjust.p) & adjust.p < 0.05, 1, 0)) %>%
    dplyr::select(micro, method, sig)

  details_wide <- details_tab %>%
    group_by(micro, method) %>%
    summarise(sig = max(sig, na.rm = TRUE), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from  = method,
      values_from = sig,
      values_fill = list(sig = 0)
    ) %>%
    as.data.frame()

  write_sheet2(mg_diff_wb, "diff_all_methods_summary", details_wide)
  openxlsx::saveWorkbook(mg_diff_wb, diff_xlsx_path, overwrite = TRUE)
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
id_groups <- sample_data(ps.micro)$Group %>% unique()
aaa       <- combn(id_groups, 2)
i         <- 1
group2    <- c(aaa[1, i], aaa[2, i])
group_df  <- data.frame(group = group2)

pst <- ps.micro %>%
  subset_samples.wt("Group", group2) %>%
  filter_taxa(function(x) sum(x) > 10, prune = TRUE)

### ---------- 6.1 loadingPCA.metm：PCA 载荷筛选特征 ----------

res34 <- loadingPCA.metm(ps = pst, Top = 20)
p34_1 <- res34[[1]]
dat34 <- res34[[2]]

save_plot2(p34_1 + theme_classic(), mg_biomarker_path, "pca_loading", width = 10, height = 8)
write_sheet2(mg_biomarker_wb, "pca_loading_data", dat34)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.2 LDA.metm：LEfSe 风格 LDA 特征 ----------

tablda <- LDA.metm(
  ps      = pst,
  Top     = 10,
  p.lvl   = 0.05,
  lda.lvl = 4,
  seed    = 11,
  adjust.p = FALSE
)

lefse_tab <- tablda[[2]]
p35 <- lefse_bar(taxtree = lefse_tab) + theme_classic()

save_plot2(p35, mg_biomarker_path, "lda_barplot", width = 10, height = 8)

write_sheet2(mg_biomarker_wb, "lda_results", lefse_tab)
write_sheet2(
  mg_biomarker_wb,
  "lda_parameters",
  data.frame(
    Parameter = c("Top", "p.lvl", "lda.lvl", "seed", "adjust.p"),
    Value     = c("10", "0.05", "4", "11", "FALSE")
  )
)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.3 svm_metm：SVM 筛选特征 ----------

res_svm <- svm_metm(ps = pst %>% filter_OTU_ps(100), k = 5)
AUC_svm <- res_svm[[1]]
imp_svm <- res_svm[[2]]

write_sheet2(mg_biomarker_wb, "svm_auc", data.frame(AUC = AUC_svm))
write_sheet2(mg_biomarker_wb, "svm_importance", imp_svm)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.4 glm.metm：GLM 筛选特征 ----------

res_glm <- glm.metm(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC_glm <- res_glm[[1]]
imp_glm <- res_glm[[2]]

write_sheet2(mg_biomarker_wb, "glm_auc", data.frame(AUC = AUC_glm))
write_sheet2(mg_biomarker_wb, "glm_importance", imp_glm)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.5 lasso.metm：Lasso 特征选择 ----------

library(glmnet)

res_lasso <- lasso.metm(ps = pst, top = 20, seed = 1010, k = 5)
acc_lasso <- res_lasso[[1]]
imp_lasso <- res_lasso[[2]]

write_sheet2(mg_biomarker_wb, "lasso_accuracy", data.frame(Accuracy = acc_lasso))
write_sheet2(mg_biomarker_wb, "lasso_importance", imp_lasso)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.6 xgboost.metm：XGBoost 特征选择 ----------

library(xgboost)
library(Ckmeans.1d.dp)
library(mia)

res_xgb <- xgboost.metm(ps = pst, top = 20)
acc_xgb <- res_xgb[[1]]
imp_xgb <- as.data.frame(res_xgb[[2]]$importance)

write_sheet2(mg_biomarker_wb, "xgboost_accuracy", data.frame(Accuracy = acc_xgb))
write_sheet2(mg_biomarker_wb, "xgboost_importance", imp_xgb)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.7 decisiontree.metm：决策树 ----------

library(rpart)

res_tree <- decisiontree.metm(ps = pst, top = 50, seed = 6358, k = 5)
acc_tree <- res_tree[[1]]
imp_tree <- res_tree[[2]]

write_sheet2(mg_biomarker_wb, "tree_accuracy", data.frame(Accuracy = acc_tree))
write_sheet2(mg_biomarker_wb, "tree_importance", imp_tree)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.8 randomforest.metm：随机森林 ----------

res_rf <- randomforest.metm(
  ps     = pst %>% filter_OTU_ps(20),
  group  = "Group",
  optimal = 20
)

p42_1 <- res_rf[[1]]
p42_2 <- res_rf[[2]]
dat_rf <- res_rf[[3]]
p42_4  <- res_rf[[4]]

save_plot2(p42_1 + theme_classic(), mg_biomarker_path, "rf_importance", width = 16, height = 11)
save_plot2(p42_2, mg_biomarker_path, "rf_error",      width = 10, height = 8)
save_plot2(p42_4, mg_biomarker_path, "rf_additional", width = 6, height = 4)

write_sheet2(mg_biomarker_wb, "rf_importance", dat_rf)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.9 bagging.metm：Bagging ----------

library(ipred)

res_bag <- bagging.metm(ps = pst, top = 20, seed = 1010, k = 5)
acc_bag <- res_bag[[1]]
imp_bag <- res_bag[[2]]

write_sheet2(mg_biomarker_wb, "bagging_accuracy", data.frame(Accuracy = acc_bag))
write_sheet2(mg_biomarker_wb, "bagging_importance", imp_bag)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.10 naivebayes.metm：朴素贝叶斯 ----------

res_nb <- naivebayes.metm(ps = pst, top = 20, seed = 1010, k = 5)
acc_nb <- res_nb[[1]]
imp_nb <- res_nb[[2]]

write_sheet2(mg_biomarker_wb, "naivebayes_accuracy", data.frame(Accuracy = acc_nb))
write_sheet2(mg_biomarker_wb, "naivebayes_importance", imp_nb)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.11 nnet.metm：浅层神经网络 ----------

options(expressions = 500000)

res_nnet <- nnet.metm(ps = pst %>% filter_OTU_ps(100), seed = 1010, k = 5)
acc_nnet <- res_nnet[[1]]
imp_nnet <- res_nnet[[2]]

write_sheet2(mg_biomarker_wb, "nnet_accuracy", data.frame(Accuracy = acc_nnet))
write_sheet2(mg_biomarker_wb, "nnet_importance", imp_nnet)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.12 MLP / stacking 模型（可选，不强制写表） ----------

if (FALSE) {
  # 39 多层感知机 MLP
  res_mlp <- mlp_metm_neuralnet(
    pst %>% filter_OTU_ps(1000),
    group    = "Group",
    hidden   = c(20, 10),
    seed     = 123,
    k        = 5,
    threshold = 0.01
  )
  print(res_mlp$Accuracy)

  # 40 集成学习 1：RF + XGB + SVM
  library(caret)
  library(pROC)
  library(zCompositions)
  library(compositions)
  library(e1071)

  fit1 <- stacking_RF_XGB_SVM(pst, outcome_col = "Group", k = 5)
  print(fit1$summary_metrics)

  # 41 集成学习 2：RF + EN + XGB
  fit2 <- stacking_RF_EN_XGB(pst, outcome_col = "Group", k = 5)
  print(fit2$summary_metrics)
}

### ---------- 6.13 rfcv.metm：随机森林特征数交叉验证 ----------

library(randomForest)
library(caret)
library(ROCR)
library(e1071)

res_rfcv <- rfcv.metm(
  ps        = pst %>% filter_OTU_ps(100),
  group     = "Group",
  optimal   = 20,
  nrfcvnum  = 6
)

p_rfcv     <- res_rfcv[[1]]
rfcv_plot  <- res_rfcv[[2]]
rfcv_table <- res_rfcv[[3]]

save_plot2(p_rfcv + theme_classic(), mg_biomarker_path, "rfcv_plot", width = 10, height = 8)

write_sheet2(mg_biomarker_wb, "rfcv_plot_data", rfcv_plot)
write_sheet2(mg_biomarker_wb, "rfcv_summary",   rfcv_table)
openxlsx::saveWorkbook(mg_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

### ---------- 6.14 Roc.metm：ROC 曲线 ----------

res_roc <- Roc.metm(
  ps     = pst %>% filter_OTU_ps(500),
  group  = "Group",
  repnum = 5
)

p33_1 <- res_roc[[1]]
AUC    <- res_roc[[2]]
dat33  <- res_roc[[3]]

save_plot2(p33_1 + theme_classic(), mg_biomarker_path, "roc_curve", width = 8, height = 6)

write_sheet2(mg_biomarker_wb, "roc_results", dat33)
write_sheet2(mg_biomarker_wb, "roc_auc_summary", data.frame(AUC = AUC))
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

### ---------- 7.1 network.pip：主协同网络构建 ----------

tab_net <- network.pip(
  ps              = ps.micro,
  N               = 100,
  big             = TRUE,
  select_layout   = FALSE,
  layout_net      = "model_maptree2",
  r.threshold     = 0.8,
  p.threshold     = 0.05,
  maxnode         = 2,
  label           = FALSE,
  lab             = "elements",
  group           = "Group",
  fill            = "Phylum",
  size            = "igraph.degree",
  zipi            = TRUE,
  ram.net         = TRUE,
  clu_method      = "cluster_fast_greedy",
  step            = 100,
  R               = 10,
  ncpus           = 1
)

net_plot_list <- tab_net[[1]]
net_obj       <- tab_net[[2]]

cortab <- net_obj$net.cor.matrix$cortab   # list，每个组一张相关矩阵




## 主网络图（第一张）
p_net_main <- net_plot_list[[1]]
save_plot2(p_net_main, mg_network_path, "network_main", width = 14, height = 8)

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
        data     = cor_df,
        rownames = TRUE
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
write_sheet2(mg_network_wb, "network_parameters", network_info, rownames = FALSE)

openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


### ---------- 7.2 net_properties.4：每个组整体网络属性 ----------

cor_list <- cortab
grp_id   <- names(cor_list)
net_prop_all <- NULL

for (i in seq_along(grp_id)) {
  g_i  <- cor_list[[grp_id[i]]] %>% make_igraph()
  dat_i <- net_properties.4(g_i, n.hub = FALSE)
  colnames(dat_i) <- grp_id[i]
  if (is.null(net_prop_all)) {
    net_prop_all <- dat_i
  } else {
    net_prop_all <- cbind(net_prop_all, dat_i)
  }
}
head(net_prop_all)

write_sheet2(
  wb       = mg_network_wb,
  sheet    = "network_properties_summary",
  data     = net_prop_all,
  rownames = TRUE
)
openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


### ---------- 7.3 netproperties.sample：单样本网络属性 + metadata 合并 ----------

net_prop_sample_all <- NULL
for (i in seq_along(grp_id)) {
  pst_i <- ps.micro %>%
    subset_samples.wt("Group", grp_id[i]) %>%
    remove.zero()
  dat_f <- netproperties.sample(pst = pst_i, cor = cor_list[[grp_id[i]]])
  if (is.null(net_prop_sample_all)) {
    net_prop_sample_all <- dat_f
  } else {
    net_prop_sample_all <- rbind(net_prop_sample_all, dat_f)
  }
}

## 合并 map 信息
map_df <- phyloseq::sample_data(ps.micro) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("ID")

dat3 <- net_prop_sample_all %>%
  tibble::rownames_to_column("ID") %>%
  dplyr::inner_join(map_df, by = "ID")

write_sheet2(
  wb       = mg_network_wb,
  sheet    = "sample_network_properties",
  data     = net_prop_sample_all,
  rownames = TRUE
)
write_sheet2(
  wb       = mg_network_wb,
  sheet    = "sample_netprop_with_meta",
  data     = dat3,
  rownames = FALSE
)
openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


### ---------- 7.4 node_properties：节点属性汇总 ----------

node_prop_all <- NULL
for (i in seq_along(grp_id)) {
  g_i <- cor_list[[grp_id[i]]] %>% make_igraph()
  node_i <- node_properties(g_i) %>%
    as.data.frame()
  node_i$Group <- grp_id[i]
  node_i <- node_i %>%
    tibble::rownames_to_column("ASV.name")

  ## 为避免列名冲突给该组列名加后缀
  colnames(node_i)[-1] <- paste0(colnames(node_i)[-1], ".", grp_id[i])

  if (is.null(node_prop_all)) {
    node_prop_all <- node_i
  } else {
    node_prop_all <- dplyr::full_join(node_prop_all, node_i, by = "ASV.name")
  }
}
head(node_prop_all)

write_sheet2(
  wb       = mg_network_wb,
  sheet    = "node_properties_combined",
  data     = node_prop_all,
  rownames = FALSE
)
openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


### ---------- 7.5 module.compare.net.pip：网络模块比较 ----------

mod_comp <- module.compare.net.pip(
  ps         = NULL,
  corg       = cor_list,
  degree     = TRUE,
  zipi       = FALSE,
  r.threshold = 0.8,
  p.threshold = 0.05,
  method      = "spearman",
  padj        = FALSE,
  n           = 3
)

res_mod <- mod_comp[[1]]
head(res_mod)

write_sheet2(
  wb       = mg_network_wb,
  sheet    = "network_comparison_results",
  data     = res_mod,
  rownames = FALSE
)
openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


## ===================== 8. 机器学习网络构建 =====================

### ---------- 8.1 rf_network_from_phyloseq：随机森林网络 ----------
library(ranger)
rf_net <- rf_network_from_phyloseq(
  ps           = ps.micro %>% filter_OTU_ps(50),
  ntree        = 500,
  normalize    = TRUE,
  scale.method = "l1",
  n_perm       = 100,
  p.threshold  = 0.001,
  seed         = 123
)

edges_all  <- rf_net$edges
edges_filt <- edges_all %>% dplyr::filter(pval < 0.05)

## 转为 igraph / 邻接矩阵
g_rf <- graph_from_data_frame(
  edges_filt[, c("target", "predictor", "importance")],
  directed = FALSE
)
cor_mat_rf <- as_adjacency_matrix(g_rf, attr = "importance", sparse = FALSE)

## 模块划分 + 坐标
mod_rf <- model_maptree2(cor = cor_mat_rf,
                         method = "cluster_fast_greedy")
node_rf <- mod_rf[[1]]
edge_rf <- edgeBuild(cor = cor_mat_rf, node = node_rf)

## 网络可视化
p_rf_net <- ggplot() +
  geom_segment(aes(x = X1, y = Y1, xend = X2, yend = Y2, color = as.factor(cor)),
               data = edge_rf, size = 1, alpha = 0.1) +
  geom_point(aes(X1, X2), pch = 21, data = node_rf) +
  ggrepel::geom_text_repel(aes(X1, X2, label = elements), data = node_rf) +
  scale_colour_brewer(palette = "Set1") +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(breaks = NULL) +
  theme(panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.background = element_rect(colour = NA),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

save_plot2(p_rf_net, mg_network_path, "rf_network_plot", width = 12, height = 10)

## 写 RF 网络表
write_sheet2(mg_network_wb, "rf_edges_all",   edges_all,  rownames = FALSE)
write_sheet2(mg_network_wb, "rf_edges_filt",  edges_filt, rownames = FALSE)
write_sheet2(mg_network_wb, "rf_node_layout", node_rf,    rownames = FALSE)
write_sheet2(mg_network_wb, "rf_cor_matrix",  as.data.frame(cor_mat_rf), rownames = TRUE)

openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)


### ---------- 8.2 ML.network：多算法机器学习网络 ----------

ml_net <- ML.network(
  ps           = ps.micro %>% filter_OTU_ps(50),
  scale.method = "l1",
  method       = c("rf","xgb","glmnet","svm","lgbm","mlp","gam")[6],  # 这里是 "mlp"
  ntree        = 100,
  p.threshold  = 0.05,
  n_perm       = 10,
  n_cores      = 10
)

## 主要返回：edges + 各模型重要性
write_sheet2(mg_network_wb, "ML_network_edges", ml_net$edges,          rownames = FALSE)
write_sheet2(mg_network_wb, "ML_network_importance", ml_net$importance_all, rownames = FALSE)

openxlsx::saveWorkbook(mg_network_wb, network_xlsx_path, overwrite = TRUE)
