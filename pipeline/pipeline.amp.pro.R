# 扩增子微生物组学分析 -----
rm(list = ls())

library(EasyMultiOmics)
library(phyloseq)
library(tidyverse)
library(ggClusterNet)

library(ggsci)
library(openxlsx)
library(ape)
library(picante)
library(fs)


## ===================== 0. 基础设置 & 通用函数 =====================

# 基础路径（可根据项目修改一次即可）
# project_root   <- "E:/Shared_Folder/师晶晶氮肥梯度样本"
# amplicon_path  <- file.path(project_root, "amplicon.16S")
# dir_create(amplicon_path)
amplicon_path = "./result/"
get_group_cols <- function(groups,
                           palette = c("npg", "nejm", "lancet")) {
  palette <- match.arg(palette)


  pal_fun <- switch(
    palette,
    npg    = ggsci::pal_npg("nrc"),
    nejm   = ggsci::pal_nejm(),
    lancet = ggsci::pal_lancet()
  )

  cols <- pal_fun(length(groups))
  setNames(cols, groups)
}


save_plot2 <- function(p,
                       out_dir,
                       prefix,
                       width       = NULL,   # 如果指定，就完全按这个宽高
                       height      = NULL,
                       dpi         = 300,
                       base_width  = 8,      # 未指定宽高时，从这个基础值自动放大
                       base_height = 6,
                       x_per_inch      = 6,  # 每 inch 放多少个 x 类别
                       facets_per_row  = 4   # 每行 facet 预期数量
) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # 估算 x 轴类别数（如果拿不到就忽略）
  n_x <- NA_integer_
  if (!is.null(p$mapping$x)) {
    xvar <- rlang::as_name(p$mapping$x)
    if (!is.null(p$data) && xvar %in% names(p$data)) {
      n_x <- dplyr::n_distinct(p$data[[xvar]])
    }
  }

  # 估算 facet 数量
  n_facets <- 1L
  gb <- tryCatch(ggplot2::ggplot_build(p),
                 error = function(e) NULL)
  if (!is.null(gb) && !is.null(gb$layout$layout$PANEL)) {
    n_facets <- length(unique(gb$layout$layout$PANEL))
  }

  # 宽度：如果用户没指定，就根据 x 类别数加宽
  if (is.null(width)) {
    add_w <- if (!is.na(n_x)) max(0, n_x / x_per_inch) else 0
    width_final <- base_width + add_w
  } else {
    width_final <- width
  }

  # 高度：如果用户没指定，就根据 facet 行数加高
  if (is.null(height)) {
    facet_rows <- ceiling(n_facets / facets_per_row)
    add_h      <- max(0, (facet_rows - 1) * 2)  # 每多一行 facet 高 +2 inch
    height_final <- base_height + add_h
  } else {
    height_final <- height
  }

  # 实际保存 PNG + PDF
  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".png")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    dpi       = dpi,
    units     = "in",
    limitsize = FALSE
  )

  ggplot2::ggsave(
    filename  = file.path(out_dir, paste0(prefix, ".pdf")),
    plot      = p,
    width     = width_final,
    height    = height_final,
    units     = "in"
  )
}

# 专门用于保存 circlize / base 图的通用保存函数
save_circlize_plot <- function(expr,
                               out_dir,
                               prefix,
                               width  = 10,
                               height = 10,
                               res    = 300) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  # 保存 PNG
  png(
    filename = file.path(out_dir, paste0(prefix, ".png")),
    width    = width,
    height   = height,
    units    = "in",
    res      = res
  )
  eval(substitute(expr), envir = parent.frame())
  dev.off()

  # 保存 PDF
  pdf(
    file   = file.path(out_dir, paste0(prefix, ".pdf")),
    width  = width,
    height = height
  )
  eval(substitute(expr), envir = parent.frame())
  dev.off()
}


# 统一写入 excel sheet 的函数
write_sheet2 <- function(wb, sheet_name, df, row_names = TRUE) {
  if (sheet_name %in% names(wb)) {
    openxlsx::removeWorksheet(wb, sheet_name)
  }
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, df, rowNames = row_names)
}


## ===================== 1. 读入 phyloseq 对象 =====================

# 方式一：交互选择
# file_path <- tcltk::tk_choose.files(
#   caption = "请选择 ps_ITS.rds 文件",
#   multi   = FALSE,
#   filters = matrix(c("RDS files", ".rds",
#                      "All files", "*"), ncol = 2, byrow = TRUE)
# )

# 方式二：直接指定路径
# file_path <- "./data/ps_its.rds"
# ps.16s    <- readRDS(file_path)
ps.16s

map <- sample_data(ps.16s)
head(map)

# 提取分组因子数量与顺序
gnum       <- phyloseq::sample_data(ps.16s)$Group %>% unique() %>% length()
axis_order <- phyloseq::sample_data(ps.16s)$Group %>% unique()

col.g = get_group_cols(axis_order)

# 主题 & 颜色
package.amp()
res      <- theme_my(ps.16s)
mytheme1 <- res[[1]]
mytheme2 <- res[[2]]
colset1  <- res[[3]]
colset2  <- res[[4]]
colset3  <- res[[5]]
colset4  <- res[[6]]


## ===================== 2. Alpha 多样性分析 =====================

# 2.1 创建 Alpha 多样性目录 & 工作簿（仿照 composition/diff 的写法）
amplicon_alpha_path <- file.path(amplicon_path, "01_alpha_diversity")
dir.create(amplicon_alpha_path, recursive = TRUE, showWarnings = FALSE)

alpha_xlsx_path <- file.path(amplicon_alpha_path, "alpha_diversity_results.xlsx")

if (file.exists(alpha_xlsx_path)) {
  amplicon_alpha_wb <- openxlsx::loadWorkbook(alpha_xlsx_path)
} else {
  amplicon_alpha_wb <- openxlsx::createWorkbook()
}



### ---------- 2.1 alpha.micro：Shannon/Chao1 等 ----------

all.alpha <- c("Shannon", "Inv_Simpson", "Pielou_evenness",
               "Simpson_evenness", "Richness", "Chao1", "ACE")

# 计算 alpha 指标
tab <- alpha.micro(ps = ps.16s, group = "Group")
head(tab)

data <- cbind(
  data.frame(ID = 1:length(tab$Group), group = tab$Group),
  tab[all.alpha]
)
data$ID <- as.character(data$ID)
head(data)

# Kruskal-Wallis + 多重比较
result  <- MuiKwWlx2(data = data, num = 3:9)

# 盒线图
result1 <- FacetMuiPlotresultBox2(
  data   = data,
  num    = 3:9,
  result = result,
  sig_show = "abc",
  ncol     = 4
)
p1_1 <- result1[[1]] +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_color_manual(values = col.g)

# 柱状图
res2 <- FacetMuiPlotresultBar(
  data   = data,
  num    = 3:9,
  result = result,
  sig_show = "abc",
  ncol     = 4
)
p1_2 <- res2[[1]] +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

# 盒线 + 柱状
res3 <- FacetMuiPlotReBoxBar(
  data   = data,
  num    = 3:9,
  result = result,
  sig_show = "abc",
  ncol     = 3
)
p1_3 <- res3[[1]] +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

# 小提琴图（使用 FacetMui 输出的数据）
p1_0 <- result1[[2]] %>%
  ggplot(aes(x = group, y = dd)) +
  geom_violin(alpha = 1, aes(fill = group)) +
  geom_jitter(aes(color = group),
              position = position_jitter(0.17),
              size = 3, alpha = 0.5) +
  labs(x = "", y = "") +
  facet_wrap(. ~ name, scales = "free_y", ncol = 4) +
  geom_text(aes(x = group, y = y, label = stat))

p1_0 <- p1_0 +
  ggplot2::scale_x_discrete(limits = axis_order) +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = col.g)

## ---- 保存 alpha.micro 图 ----
save_plot2(p1_1, amplicon_alpha_path, "alpha_diversity_box",    width = 12, height = 8)
save_plot2(p1_2, amplicon_alpha_path, "alpha_diversity_bar",    width = 12, height = 8)
save_plot2(p1_3, amplicon_alpha_path, "alpha_diversity_boxbar", width = 12, height = 8)
save_plot2(p1_0, amplicon_alpha_path, "alpha_diversity_violin", width = 12, height = 8)

## ---- 保存 alpha.micro 表 ----
# 原始 alpha 数据
write_sheet2(amplicon_alpha_wb, "alpha_diversity_data", data)

# 统计分析结果（Kruskal-Wallis & 多重比较）
write_sheet2(amplicon_alpha_wb, "alpha_diversity_stat", result)
openxlsx::saveWorkbook(amplicon_alpha_wb, alpha_xlsx_path, overwrite = TRUE)

### ---------- 2.2 alpha.pd.micro：PD 多样性 ----------
# 这里修正了原来的小 bug：用 ps.16s，而不是 ps16s
tab2    <- alpha.pd.micro(ps.16s)
head(tab2)

result_pd  <- MuiKwWlx2(data = tab2, num = 3)
result_pd1 <- FacetMuiPlotresultBox(
  data   = tab2,
  num    = 3,
  result = result_pd,
  sig_show = "abc",
  ncol     = 1
)
p_pd <- result_pd1[[1]] +
  theme_nature() +
  ggplot2::guides(fill = guide_legend(title = NULL)) +
  ggplot2::scale_fill_manual(values = colset1)

## ---- 保存 PD 图 ----
save_plot2(p_pd, amplicon_alpha_path, "pd_diversity", width = 8, height = 6)

## ---- 保存 PD 表 ----
write_sheet2(amplicon_alpha_wb, "pd_diversity_data", tab2)
write_sheet2(amplicon_alpha_wb, "pd_diversity_stat", result_pd)
openxlsx::saveWorkbook(amplicon_alpha_wb, alpha_xlsx_path, overwrite = TRUE)

### ---------- 2.3 alpha.rare.line: alpha 多样性稀释曲线 ----------

rare <- mean(phyloseq::sample_sums(ps.16s)) / 10

result_rare <- alpha.rare.line.micro(
  ps     = ps.16s,
  group  = "Group",
  method = "Richness",
  start  = 100,
  step   = rare
)

# 单样本稀释曲线
p2_1 <- result_rare[[1]] +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90)) +
  scale_color_manual(values = col.g)

## 稀释表格
raretab <- result_rare[[2]]

# 分组稀释曲线
p2_2 <- result_rare[[3]] +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90)) +
  scale_color_manual(values = col.g) +
  scale_fill_manual(values = col.g)

# 分组稀释曲线（带标准差）
p2_3 <- result_rare[[4]] +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90)) +
  scale_color_manual(values = col.g) +
  scale_fill_manual(values = col.g)

## ---- 保存稀释曲线图 ----
save_plot2(p2_1, amplicon_alpha_path, "rarefaction_individual", width = 10, height = 8)
save_plot2(p2_2, amplicon_alpha_path, "rarefaction_group",      width = 10, height = 8)
save_plot2(p2_3, amplicon_alpha_path, "rarefaction_group_sd",   width = 10, height = 8)

## ---- 保存稀释曲线数据 ----
write_sheet2(amplicon_alpha_wb, "rarefaction_data", raretab)
openxlsx::saveWorkbook(amplicon_alpha_wb, alpha_xlsx_path, overwrite = TRUE)

## ===================== 3. Beta 多样性分析 =====================

# 创建 Beta 多样性分析目录 & 工作簿（仿照 composition/diff 的写法）
amplicon_beta_path <- file.path(amplicon_path, "beta_diversity")
dir.create(amplicon_beta_path, recursive = TRUE, showWarnings = FALSE)

beta_xlsx_path <- file.path(amplicon_beta_path, "beta_diversity.xlsx")

if (file.exists(beta_xlsx_path)) {
  amplicon_beta_wb <- openxlsx::loadWorkbook(beta_xlsx_path)
} else {
  amplicon_beta_wb <- openxlsx::createWorkbook()
}



### ---------- 3.1 ordinate.micro: 排序分析 ----------

result_ord <- ordinate.micro(
  ps           = ps.16s,
  group        = "Group",
  dist         = "bray",
  method       = "PCoA",
  Micromet     = "anosim",
  pvalue.cutoff = 0.05,
  pair         = FALSE
)

# 基础排序图
p3_1 <- result_ord[[1]] +
  scale_fill_manual(values = col.g) +
  scale_color_manual(values = col.g, guide = "none") +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

# 排序点坐标数据
plotdata <- result_ord[[2]]

# 带标签的排序图
p3_2 <- result_ord[[3]] +
  scale_fill_manual(values = col.g) +
  scale_color_manual(values = col.g, guide = "none") +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

# 精修图：加群心和“蜘蛛线”
cent <- aggregate(cbind(x, y) ~ Group, data = plotdata, FUN = mean)
segs <- merge(
  plotdata,
  setNames(cent, c("Group", "oNMDS1", "oNMDS2")),
  by = "Group",
  sort = FALSE
)

p3_3 <- result_ord[[1]] +
  geom_segment(
    data    = segs,
    mapping = aes(x = x, y = y, xend = oNMDS1, yend = oNMDS2, color = Group),
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

## ---- 保存排序图 ----
save_plot2(p3_1, amplicon_beta_path, "ordination_basic",   width = 10, height = 8)
save_plot2(p3_2, amplicon_beta_path, "ordination_labeled", width = 10, height = 8)
save_plot2(p3_3, amplicon_beta_path, "ordination_refined", width = 10, height = 8)

## ---- 保存排序坐标数据 ----
write_sheet2(amplicon_beta_wb, "ordination_data", plotdata)
openxlsx::saveWorkbook(amplicon_beta_wb ,beta_xlsx_path, overwrite = TRUE)


### ---------- 3.2 MicroTest.micro: 群落水平差异检测 ----------

dat1 <- MicroTest.micro(ps = ps.16s, Micromet = "adonis", dist = "bray")
write_sheet2(amplicon_beta_wb, "microtest_results", dat1)
openxlsx::saveWorkbook(amplicon_beta_wb ,beta_xlsx_path, overwrite = TRUE)


### ---------- 3.3 pairMicroTest.micro: 两两分组差异 ----------

dat2 <- pairMicroTest.micro(ps = ps.16s, Micromet = "MRPP", dist = "bray")
write_sheet2(amplicon_beta_wb, "pair_microtest_results", dat2)
openxlsx::saveWorkbook(amplicon_beta_wb ,beta_xlsx_path, overwrite = TRUE)


### ---------- 3.4 mantal.micro: Mantel 分析 ----------

result_m <- mantal.micro(
  ps     = ps.16s,
  method = "spearman",
  group  = "Group",
  ncol   = gnum,
  nrow   = 1
)

data_m <- result_m[[1]]
p3_7   <- result_m[[2]] + theme_nature()

save_plot2(p3_7, amplicon_beta_path, "mantel_test", width = 10, height = 8)
write_sheet2(amplicon_beta_wb, "mantel_results", data_m)
openxlsx::saveWorkbook(amplicon_beta_wb ,beta_xlsx_path, overwrite = TRUE)


### ---------- 3.5 cluster.micro: 样品聚类 ----------

res_clust <- cluster_micro(
  ps            = ps.16s,
  hcluter_method = "complete",
  dist          = "bray",
  cuttree       = 3,
  row_cluster   = TRUE,
  col_cluster   = TRUE
)

p4    <- res_clust[[1]]
p4_1  <- res_clust[[2]]
p4_2  <- res_clust[[3]]
dat_c <- res_clust[[4]]

save_plot2(p4,   amplicon_beta_path, "cluster_heatmap",     width = 12, height = 10)
save_plot2(p4_1, amplicon_beta_path, "cluster_dendrogram1", width = 10, height = 8)
save_plot2(p4_2, amplicon_beta_path, "cluster_dendrogram2", width = 10, height = 8)
write_sheet2(amplicon_beta_wb, "cluster_results", dat_c)
openxlsx::saveWorkbook(amplicon_beta_wb ,beta_xlsx_path, overwrite = TRUE)



# ### ---------- 3.6 distance.micro: 分组间距离比较 ----------
#
# res_dist <- distance.micro2(ps.16s, group = "Group")
#
# p5.1 <- res_dist[[1]] + mytheme1
#
# p5.2 <- res_dist[[2]] +
#   scale_fill_manual(values = col.g) +
#   scale_color_manual(values = col.g, guide = "none") +
#   theme_nature() +
#   theme(axis.title.y = element_text(angle = 90))
#
# p5.3 <- res_dist[[3]] +
#   scale_fill_manual(values = col.g) +
#   scale_color_manual(values = col.g, guide = "none") +
#   theme_nature() +
#   theme(axis.title.y = element_text(angle = 90))
#
# dat_d <- res_dist[[4]]
#
# save_plot2(p5.1, amplicon_beta_path, "distance_plot1", width = 10, height = 8)
# save_plot2(p5.2, amplicon_beta_path, "distance_plot2", width = 10, height = 8)
# save_plot2(p5.3, amplicon_beta_path, "distance_plot3", width = 10, height = 8)
#
# write_sheet2(amplicon_beta_wb, "distance_data", dat_d)
# openxlsx::saveWorkbook(amplicon_beta_wb ,beta_xlsx_path, overwrite = TRUE)


## ===================== 4. 组成（composition）分析 =====================

amplicon_composition_path <- file.path(amplicon_path, "composition")
dir.create(amplicon_composition_path, recursive = TRUE, showWarnings = FALSE)

amplicon_composition_wb <- createWorkbook()

# 创建微生物组成分析目录
amplicon_composition_path <- file.path(amplicon_path, "composition")
dir.create(amplicon_composition_path, recursive = TRUE, showWarnings = FALSE)

# Excel 保存路径
comp_xlsx_path <- file.path(amplicon_composition_path, "composition_results2.xlsx")

# 如果已经存在，就加载接着写；否则新建
if (file.exists(comp_xlsx_path)) {
  amplicon_composition_wb <- loadWorkbook(comp_xlsx_path)
} else {
  amplicon_composition_wb <- createWorkbook()
}





### ---------- 10. Ven.Upset.micro: 共有/特有 OTU/ASV ----------
# install.packages("ggvenn")
library(ggvenn)
res10 <- Ven.Upset.metm(
  ps    = ps.16s,
  group = "Group",
  N     = 0.5,
  size  = 3
)

p10_1 <- res10[[1]]
p10_2 <- res10[[2]]
# p10_3 <- res10[[4]]
dat10 <- res10[[3]]

save_plot2(p10_1, amplicon_composition_path, "10_venn_diagram",
           base_width = 10, base_height = 8)
save_plot2(p10_2, amplicon_composition_path, "10_upset_plot",
           base_width = 12, base_height = 8)

write_sheet2(amplicon_composition_wb, "10_venn_upset_data", dat10)


## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



### ---------- 11. ggVen.Upset.micro: ggplot 风格 Upset ----------

res11 <- ggVen.Upset.micro(ps = ps.16s, group = "Group")
p11   <- res11[[1]]
dat11 <- res11[[2]]

# 如果你不想每次都弹出图，可以去掉这句 grid.draw
# grid::grid.draw(p11)

save_plot2(p11, amplicon_composition_path, "11_ggVen_upset", base_width = 12, base_height = 10)
write_sheet2(amplicon_composition_wb, "11_ggVen_data", dat11)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)

### ---------- 12. VenSuper.micro: 各部分物种组成 ----------

library(ggpubr)
library(agricolae)
library(reshape2)

res12 <- VenSuper.metm(
  ps    = ps.16s,
  group = "Group",
  num   = 6
)

p12_1 <- res12[[1]]   # 所有部分门水平柱状图
p12_2 <- res12[[3]]   # 冲积图
p12_3 <- res12[[2]]   # 比例 + 差异
dat12 <- res12[[4]]   # 详细数据
dat120 <- dplyr::bind_rows(dat12, .id = "Part")


save_plot2(p12_1, amplicon_composition_path, "12_venn_detail_bar",        base_width = 12, base_height = 8)
save_plot2(p12_2, amplicon_composition_path, "12_venn_detail_alluvial",   base_width = 12, base_height = 8)
save_plot2(p12_3, amplicon_composition_path, "12_venn_detail_proportion", base_width = 10, base_height = 8)

str(dat12)

write_sheet2(amplicon_composition_wb, "12_venn_detail_data", dat120)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)




### ---------- 13. ggflower.micro: 花瓣图（按样本/ID） ----------

res13 <- ggflower.micro(
  ps    = ps.16s,
  group = "ID",   # 这里用的是 ID，确保 sample_data 里有这列
  start = 1,
  m1    = 2,
  a     = 0.3,
  b     = 1,
  lab.leaf = 1,
  col.cir  = "yellow",
  N        = 0.5
)

p13   <- res13[[1]]
dat13 <- res13[[2]]

save_plot2(p13, amplicon_composition_path, "13_flower_plot", base_width = 10, base_height = 10)
write_sheet2(amplicon_composition_wb, "13_flower_data", dat13)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



### ---------- 14. ven.network.micro: Venn 网络 ----------

res14 <- ven.network.metm(
  ps   = ps.16s,
  N    = 0.5,
  fill = "Phylum"
)

p14   <- res14[[1]]
dat14 <- res14[[2]]

save_plot2(p14, amplicon_composition_path, "14_venn_network", base_width = 10, base_height = 8)
write_sheet2(amplicon_composition_wb, "14_venn_network_data", dat14)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



### ---------- 15. Micro_tern.micro: 三元图 ----------
group_vec <- phyloseq::sample_data(ps.16s)$Group
gnum      <- length(unique(group_vec))
if (gnum >= 3) {
  # 分组数 ≥3，运行 ps_polygon_plot 并保存图
  p_poly <- ps_polygon_plot(
    ps      = ps.16s,
    group   = "Group",
    taxrank = "Genus"
  )
  save_plot2(p_poly, amplicon_composition_path, "16_polygon_plot", base_width = 10, base_height = 8)

  res15 <- Micro_tern.micro(ps.16s %>% filter_OTU_ps(100))
  p15   <- res15[[1]]   # 通常是 list，取第一个 ggplot
  dat15 <- res15[[2]]

  save_plot2(p15[[1]], amplicon_composition_path, "15_ternary_plot", base_width = 10, base_height = 8)
  write_sheet2(amplicon_composition_wb, "15_ternary_data", dat15)
}

## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)


### ---------- 16. barMainplot.micro: 堆积柱状图 ----------

library(ggalluvial)
detach("package:mia")

pst <- ps.16s %>% subset_taxa.wt("Species", "Unassigned", TRUE)
pst <- pst %>% subset_taxa.wt("Genus",   "Unassigned", TRUE)

res16 <- barMainplot.micro(
  ps   = pst,
  j    = "Genus",
  label = FALSE,
  sd    = FALSE,
  Top   = 10
)

p16_1 <- res16[[1]] +
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

p16_2 <- res16[[3]] +
  scale_fill_manual(values = colset2) +
  scale_x_discrete(limits = axis_order) +
  theme_nature() +
  theme(axis.title.y = element_text(angle = 90))

dat16 <- res16[[2]] %>%
  dplyr::group_by(Group, aa) %>%
  dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  as.data.frame()
colnames(dat16) <- c("Group", "Genus", "Abundance(%)")

save_plot2(p16_1, amplicon_composition_path, "16_barplot_main",      base_width = 12, base_height = 8)
save_plot2(p16_2, amplicon_composition_path, "16_barplot_secondary", base_width = 12, base_height = 8)
write_sheet2(amplicon_composition_wb, "16_barplot_data", dat16)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



### ---------- 17. cluMicro.bar.micro: 聚类堆积柱状图 ----------

res17 <- cluMicro.bar.micro(
  dist           = "bray",
  ps             = ps.16s,
  j              = "Genus",
  Top            = 10,
  tran           = TRUE,
  hcluter_method = "complete",
  Group          = "Group",
  cuttree        = length(unique(phyloseq::sample_data(ps.16s)$Group))
)

p17_1 <- res17[[1]]
p17_2 <- res17[[2]]
p17_3 <- res17[[3]]
p17_4 <- res17[[4]]
dat17 <- res17[[5]]

save_plot2(p17_1, amplicon_composition_path, "17_cluster_bar1", base_width = 12, base_height = 8)
save_plot2(p17_2, amplicon_composition_path, "17_cluster_bar2", base_width = 12, base_height = 8)
save_plot2(p17_3, amplicon_composition_path, "17_cluster_bar3", base_width = 12, base_height = 8)
save_plot2(p17_4, amplicon_composition_path, "17_cluster_bar4", base_width = 12, base_height = 8)

write_sheet2(amplicon_composition_wb, "17_cluster_bar_data", dat17)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



### ---------- 18. cir_barplot.micro: 环状堆积柱状图 ----------

library(ggtree)

res18 <- cir_barplot.micro(
  ps             = ps.16s,
  Top            = 15,
  dist           = "bray",
  cuttree        = 3,
  hcluter_method = "complete"
)

p18   <- res18[[1]]
dat18 <- res18$plotdata  # 比你之前 res[[2]] 更稳妥

# 保存环状堆积柱状图
save_plot2(
  p18,
  amplicon_composition_path,
  "18_circular_barplot",
  width  = 10,
  height = 10
)

# 保存对应数据表
write_sheet2(
  amplicon_composition_wb,
  "18_circular_barplot_data",
  dat18
)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)


## ---------- 19. cir_plot.micro: 和弦图展示物种组成 ----------

# 先跑一次函数，拿到矩阵（同时也会在当前设备上画一次图；如果不想在 RStudio 界面弹图，可以先 par(mfrow=c(1,1)) 或忽略）
res19 <- cir_plot.micro(ps = ps.16s, Top = 12, rank = 6)
mer_otu_mean <- res19[[1]]

# 用专用函数，把图保存为 PNG + PDF（各再画一次）
save_circlize_plot(
  cir_plot.micro(ps = ps.16s, Top = 12, rank = 6),
  out_dir = amplicon_composition_path,
  prefix  = "19_chord_plot",
  width   = 10,
  height  = 10,
  res     = 300
)

# 把 mer_otu_mean 矩阵也存进 Excel
write_sheet2(
  amplicon_composition_wb,
  "19_chord_data",
  as.data.frame(mer_otu_mean)
)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)


## ---------- 20. maptree.micro: 圈图展示物种组成 ----------
library(ggraph)
library(data.tree)
library(igraph)
tax = ps.16s %>% vegan_tax() %>%
  as.data.frame()
head(tax)
ps.bac <- subset_main_kingdom_ps(ps.16s)
ps.bac

res20 <- maptree.micro(
  ps     = ps.bac,
  Top    = 100,
  labtab = NULL,
  seed   = 11
)

p20   <- res20[[1]]
dat20 <- res20[[2]]

save_plot2(
  p20,
  amplicon_composition_path,
  "20_maptree",
  width  = 12,
  height = 12
)

write_sheet2(amplicon_composition_wb, "20_maptree_data", dat20)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



## ---------- 21. phy_tree.micro: 系统发育树 + 相关性 ----------

library(ggstar)

res21 <- EasyMultiOmics:::phy_tree.micro(ps = ps.16s,Top = 100)

p21_0 <- res21[[1]]
p21_1 <- res21[[2]]
# 中间 p2 ~ p6 如不需要可不保存
p21_7 <- res21[[8]]
dat21 <- res21[[9]]

save_plot2(
  p21_0,
  amplicon_composition_path,
  "21_phy_tree_main",
  width  = 12,
  height = 10
)
save_plot2(
  p21_1,
  amplicon_composition_path,
  "21_phy_tree_p1",
  width  = 12,
  height = 10
)
save_plot2(
  p21_7,
  amplicon_composition_path,
  "21_phy_tree_p7",
  width  = 12,
  height = 10
)

write_sheet2(amplicon_composition_wb, "21_phy_tree_data", dat21)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



## ---------- 22. sankey.micro: 桑基图展示物种组成 ----------

res22 <- sankey.micro(
  ps   = ps.16s,
  rank = 6,
  Top  = 100
)

p22   <- res22[[1]]
dat22 <- res22[[2]]

saveNetwork(
  p22,
  file.path(amplicon_composition_path, "22_sankey.html")
)




write_sheet2(amplicon_composition_wb, "22_sankey_data", dat22)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)




## ---------- 23. sankey.m.Group.micro: 按分组桑基图 ----------

res23 <- sankey.m.Group.micro(
  ps   = ps.16s %>% subset_taxa.wt("Species", "Unassigned", TRUE),
  rank = 6,
  Top  = 50
)

p23   <- res23[[1]]
dat23 <- res23[[2]]

# 保存交互 HTML
saveNetwork(
  p23,
  file.path(amplicon_composition_path, "23_sankey_Group.html")
)

# 数据表写入 Excel
write_sheet2(amplicon_composition_wb, "23_sankey_group_data", dat23)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



## ---------- 24. Microheatmap.micro: 热图展示物种相对丰度差异 ----------

heatnum <- 30

ps_tem <- ps.16s %>%
  ggClusterNet::scale_micro(method = "TMM") %>%
  ggClusterNet::tax_glom_wt(ranks = "Genus")

id <- ps.16s %>%
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

res24 <- Microheatmap.micro(
  ps_rela     = ps_tem,
  id          = id,
  col_cluster = TRUE,
  row_cluster = TRUE
)

p24_1 <- res24[[1]]
p24_2 <- res24[[2]]
dat24 <- res24[[3]]

save_plot2(
  p24_1,
  amplicon_composition_path,
  "24_heatmap1",
  width  = 12,
  height = 10
)
save_plot2(
  p24_2,
  amplicon_composition_path,
  "24_heatmap2",
  width  = 12,
  height = 10
)

write_sheet2(amplicon_composition_wb, "24_heatmap_data", dat24)
## 保存一次 Excel
# save_comp_wb(amplicon_composition_wb, comp_xlsx_path)
openxlsx::saveWorkbook(amplicon_composition_wb, comp_xlsx_path, overwrite = TRUE)



## ===================== 4. 差异分析（Differential Analysis） =====================

# 差异分析目录
amplicon_diff_path <- file.path(amplicon_path, "differential")
dir.create(amplicon_diff_path, recursive = TRUE, showWarnings = FALSE)

# 差异分析 Excel 文件路径
diff_xlsx_path <- file.path(amplicon_diff_path, "differential_results.xlsx")

# 若已有，则在原文件基础上追加；否则新建
if (file.exists(diff_xlsx_path)) {
  amplicon_composition_wb <- openxlsx::loadWorkbook(diff_xlsx_path)
} else {
  amplicon_composition_wb <- openxlsx::createWorkbook()
}



## ---------- 25. EdgerSuper.micro: EdgeR 计算差异微生物 ----------

res25 <- EdgerSuper.metm(
  ps       = ps.16s %>% ggClusterNet::filter_OTU_ps(500),
  group    = "Group",
  artGroup = NULL,
  j        = "OTU"
)

# 注意这里用 [[ ]] 取出 ggplot，而不是 [ ]
p25_1 <- res25[[1]][[1]]
dat25 <- res25[[2]]

# 保存图
save_plot2(p25_1, amplicon_diff_path, "25_EdgeR_plot1", width = 10, height = 8)

# 保存表
write_sheet2(amplicon_composition_wb, "25_EdgeR_results", dat25)
openxlsx::saveWorkbook(amplicon_composition_wb, diff_xlsx_path, overwrite = TRUE)


## ---------- 25.1 EdgerSuper2.micro: EdgeR 计算差异微生物（矩阵形式） ----------

res25_2 <- EdgerSuper2.metm(
  ps       = ps.16s,
  group    = "Group",
  artGroup = NULL,
  j        = "OTU"
)

write_sheet2(amplicon_composition_wb, "25_1_EdgeR2_results", res25_2)
openxlsx::saveWorkbook(amplicon_composition_wb, diff_xlsx_path, overwrite = TRUE)


## ---------- 26. DESep2Super.micro: DESeq2 计算差异微生物 ----------

res26 <- DESep2Super.micro(
  ps       = ps.16s %>% ggClusterNet::filter_OTU_ps(500),
  group    = "Group",
  artGroup = NULL,
  j        = "OTU"
)

p26_1 <- res26[[1]][[1]]

dat26 <- res26[[2]]

save_plot2(p26_1, amplicon_diff_path, "26_DESeq2_plot1", width = 10, height = 8)


write_sheet2(amplicon_composition_wb, "26_DESeq2_results", dat26)
openxlsx::saveWorkbook(amplicon_composition_wb, diff_xlsx_path, overwrite = TRUE)



## ---------- 27. edge_Manhattan.micro: 曼哈顿图展示差异微生物 ----------

res27 <- edge_Manhattan.metm(
  ps     = ps.16s %>% ggClusterNet::filter_OTU_ps(500),
  pvalue = 0.05,
  lfc    = 0
)

p27_1 <- res27[[1]]

save_plot2(p27_1, amplicon_diff_path, "27_Manhattan_plot1", width = 12, height = 8)

# 如果 edge_Manhattan.micro 有数据输出（例如 res27[[4]]），可以加一行：
# write_sheet2(amplicon_composition_wb, "27_Manhattan_data", res27[[4]])


## ---------- 28. stamp_diff.micro: STAMP 风格差异丰度展示 ----------

map <- phyloseq::sample_data(ps.16s)  # 确保 map 已存在
allgroup <- combn(unique(map$Group), 2)
plot_list_28 <- list()

for (i in seq_len(ncol(allgroup))) {
  ps_sub <- phyloseq::subset_samples(ps.16s, Group %in% allgroup[, i])
  p_tmp  <- stemp_diff.micro(ps = ps_sub, Top = 20, ranks = 6)
  plot_list_28[[i]] <- p_tmp
}

# 逐个保存图，不再手动写 p28.1/p28.2/p28.3（避免分组数不固定时报错）
for (i in seq_along(plot_list_28)) {
  save_plot2(
    plot_list_28[[i]],
    amplicon_diff_path,
    paste0("28_stamp_plot", i),
    width  = 10,
    height = 8
  )
}

# 如果 stemp_diff.micro 有数据输出，可在循环里顺便收集 dat_28 再写 Excel
# 这里暂时只保存图，不写表

## ---------- 29. Mui.Group.volcano.micro: 聚类火山图（所有组） ----------

res29_input <- EdgerSuper2.metm(
  ps       = ps.16s,
  group    = "Group",
  artGroup = NULL,
  j        = "OTU"
)
library(ggrepel)
res29 <- Mui.Group.volcano.micro(res = res29_input)

p29_1 <- res29[[1]]

dat29 <- res29[[3]]

save_plot2(p29_1, amplicon_diff_path, "29_volcano_cluster1", width = 12, height = 8)

write_sheet2(amplicon_composition_wb, "29_volcano_cluster_data", dat29)
openxlsx::saveWorkbook(amplicon_composition_wb, diff_xlsx_path, overwrite = TRUE)


## ---------- 30. Mui.cluster.volcano.micro: 指定分组聚类火山图 ----------


id_groups <- sample_data(ps.16s)$Group %>% unique()
aaa       <- combn(id_groups, 2)

# 取第一对作为示例（你也可以循环所有对）
group2 <- c(aaa[1, 1], aaa[2, 1])
b <- data.frame(group2)

res30_input <- EdgerSuper2.metm(
  ps       = ps.16s,
  group    = "Group",
  artGroup =b,
  j        = "OTU"
)

res30 <- Mui.cluster.volcano.micro(res = res30_input, rs.k = 6)

p30_1 <- res30[[1]] + ggtitle(paste(group2, collapse = "-"))

dat30 <- res30[[3]]

save_plot2(p30_1, amplicon_diff_path, "30_volcano_specific1", width = 12, height = 8)

write_sheet2(amplicon_composition_wb, "30_volcano_specific_data", dat30)
openxlsx::saveWorkbook(amplicon_composition_wb, diff_xlsx_path, overwrite = TRUE)



## ===================== biomarker identification =====================

# 创建生物标志物识别目录 & 总表
amplicon_biomarker_path <- file.path(amplicon_path, "biomarker")
dir.create(amplicon_biomarker_path, recursive = TRUE, showWarnings = FALSE)

biomarker_xlsx_path <- file.path(amplicon_biomarker_path, "biomarker_results.xlsx")

if (file.exists(biomarker_xlsx_path)) {
  amplicon_biomarker_wb <- openxlsx::loadWorkbook(biomarker_xlsx_path)
} else {
  amplicon_biomarker_wb <- openxlsx::createWorkbook()
}



library(randomForest)
library(caret)
library(ROCR)
library(e1071)

## ---------- 32. rfcv.micro : 交叉验证筛选重要特征 ----------

res32 <- rfcv.Micro(
  ps      = ps.16s %>% filter_OTU_ps(100),
  group   = "Group",
  optimal = 20,
  nrfcvnum = 6
)

p32_1     <- res32[[1]] + theme_classic()
rfcvtable <- res32[[3]]

# 保存图
save_plot2(
  p32_1,
  out_dir = amplicon_biomarker_path,
  prefix  = "32_rfcv",
  width   = 10,
  height  = 8
)

# 保存表
write_sheet2(amplicon_biomarker_wb, "32_rfcv_results", rfcvtable)

openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

## ---------- 33. Roc.micro : ROC 曲线绘制 ----------

id  <- sample_data(ps.16s)$Group %>% unique()
aaa <- combn(id, 2)

i      <- 1   # 仍然取第一对分组，如需多对可改为循环
group2 <- c(aaa[1, i], aaa[2, i])
group_df <- data.frame(group = group2)

ps_roc <- ps.16s %>%
  subset_taxa.wt("Family", "Unassigned", TRUE) %>%
  subset_taxa.wt("Order",  "Unassigned", TRUE) %>%
  subset_taxa.wt("Genus",  "Unassigned", TRUE) %>%
  subset_taxa.wt("Phylum", "Unassigned", TRUE) %>%
  subset_taxa.wt("class",  "Unassigned", TRUE)    # 如果分类层级是 "Class"，这里自行改成 "Class"

pst <- ps_roc %>%
  subset_samples.wt("Group", group2) %>%
  filter_taxa(function(x) sum(x) > 10, prune = TRUE)

res33 <- Roc.micro(
  ps     = pst %>% filter_OTU_ps(1000),
  group  = "Group",
  repnum = 5
)

p33_1 <- res33[[1]] + theme_classic()
p33_2 <- res33[[2]]
dat33 <- res33[[3]]

# 保存图
save_plot2(
  p33_1,
  out_dir = amplicon_biomarker_path,
  prefix  = "33_ROC_plot1",
  width   = 10,
  height  = 8
)


# 保存表
write_sheet2(amplicon_biomarker_wb, "33_ROC_results", dat33)

openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

## ---------- 34. loadingPCA.micro : PCA 载荷筛选特征微生物 ----------

res34 <- loadingPCA.micro(ps = ps.16s, Top = 20)

p34_1 <- res34[[1]] + theme_classic()
dat34 <- res34[[2]]

# 保存图
save_plot2(
  p34_1,
  out_dir = amplicon_biomarker_path,
  prefix  = "34_loadingPCA",
  width   = 10,
  height  = 8
)

# 保存表
write_sheet2(amplicon_biomarker_wb, "34_loadingPCA_results", dat34)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

## ---------- 35. LDA.micro : LDA 筛选特征微生物（LEfSe 风格） ----------

p_base <- p_base.micro(ps.16s, Top = 20)   # 如后面未直接用，可保留用于检查基础可视化

tablda <- LDA.micro(
  ps      = ps.16s,
  Top     = 20,
  p.lvl   = 0.05,
  lda.lvl = 2,
  seed    = 11,
  adjust.p = FALSE
)

lefse_tab <- tablda[[2]]

p35_1 <- lefse_bar(taxtree = lefse_tab) + theme_classic()

# 保存图
save_plot2(
  p35_1,
  out_dir = amplicon_biomarker_path,
  prefix  = "35_LDA",
  width   = 10,
  height  = 8
)

# 保存表
write_sheet2(amplicon_biomarker_wb, "35_LDA_results", lefse_tab)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#36 svm.micro:svm筛选特征微生物 ----
res <- svm_micro(ps = ps.16s %>% filter_OTU_ps(20), k = 5)
AUC        <- res[[1]]
importance <- res[[2]]

# 保存SVM结果
write_sheet2(amplicon_biomarker_wb, "svm_AUC",        AUC)
write_sheet2(amplicon_biomarker_wb, "svm_importance", importance)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#37 glm.micro:glm筛选特征微生物----
res <- glm.micro(ps = pst %>% filter_OTU_ps(50), k = 5)
AUC        <- res[[1]]
importance <- res[[2]]

# 保存GLM结果
write_sheet2(amplicon_biomarker_wb, "glm_AUC",        AUC)
write_sheet2(amplicon_biomarker_wb, "glm_importance", importance)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#38 xgboost.micro: xgboost筛选特征微生物----
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)

res <- xgboost.micro(ps = pst, top = 20)
accuracy   <- res[[1]]
importance <- res[[2]]


# 保存XGBoost结果
write_sheet2(amplicon_biomarker_wb, "xgboost_accuracy",   accuracy)
write_sheet2(amplicon_biomarker_wb, "xgboost_importance", as.data.frame(importance$importance))
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#39 lasso.micro: lasso筛选特征微生物----
library(glmnet)
res <- lasso.micro(ps = pst, top = 20, seed = 1010, k = 5)
accuracy   <- res[[1]]
importance <- res[[2]]

# 保存Lasso结果
write_sheet2(amplicon_biomarker_wb, "lasso_accuracy",   accuracy)
write_sheet2(amplicon_biomarker_wb, "lasso_importance", importance)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#40 decisiontree.micro----
library(rpart)
res <- decisiontree.micro(ps = pst, top = 50, seed = 6358, k = 5)
accuracy   <- res[[1]]
importance <- res[[2]]

# 保存决策树结果
write_sheet2(amplicon_biomarker_wb, "decisiontree_accuracy",   accuracy)
write_sheet2(amplicon_biomarker_wb, "decisiontree_importance", importance)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#41 naivebayes.micro: bayes筛选特征微生物----
res <- naivebayes.micro(ps = pst, top = 20, seed = 1010, k = 5)
accuracy   <- res[[1]]
importance <- res[[2]]

# 保存朴素贝叶斯结果
write_sheet2(amplicon_biomarker_wb, "naivebayes_accuracy",   accuracy)
write_sheet2(amplicon_biomarker_wb, "naivebayes_importance", importance)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#42 randomforest.micro: 随机森林筛选特征微生物----
res <- randomforest.micro(ps = pst, group = "Group", optimal = 50)
p42.1 <- res[[1]]
# p42.2 <- res[[2]]
dat   <- res[[3]]
p42.4 <- res[[4]]

# 保存随机森林图
save_plot2(p42.1, amplicon_biomarker_path, "randomforest_plot1", width = 10, height = 8)
# save_plot2(p42.2, amplicon_biomarker_path, "randomforest_plot2", width = 10, height = 8)
save_plot2(p42.4, amplicon_biomarker_path, "randomforest_plot4", width = 10, height = 8)

# 保存随机森林表
write_sheet2(amplicon_biomarker_wb, "randomforest_results", dat)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)


#43 bagging.micro : Bootstrap Aggregating筛选特征微生物 ------
library(ipred)
res <- bagging.micro(ps = pst, top = 20, seed = 1010, k = 5)
accuracy   <- res[[1]]
importance <- res[[2]]

# 保存Bagging结果
write_sheet2(amplicon_biomarker_wb, "bagging_accuracy",   accuracy)
write_sheet2(amplicon_biomarker_wb, "bagging_importance", importance)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)

#44 nnet.micro 神经网络筛选特征微生物  ------
res <- nnet.micro(ps = pst, top = 100, seed = 1010, k = 5)
accuracy   <- res[[1]]
importance <- res[[2]]

# 保存神经网络结果
write_sheet2(amplicon_biomarker_wb, "nnet_accuracy",   accuracy)
write_sheet2(amplicon_biomarker_wb, "nnet_importance", importance)
openxlsx::saveWorkbook(amplicon_biomarker_wb, biomarker_xlsx_path, overwrite = TRUE)


# ===================== 5. 网络分析（Network Analysis） =====================

# 创建网络分析目录
amplicon_network_path <- file.path(amplicon_path, "network")
dir.create(amplicon_network_path, recursive = TRUE, showWarnings = FALSE)

# 网络分析 Excel 文件路径
network_xlsx_path <- file.path(amplicon_network_path, "network_analysis.xlsx")

# 若已有，则在原文件基础上追加；否则新建
if (file.exists(network_xlsx_path)) {
  amplicon_network_wb <- openxlsx::loadWorkbook(network_xlsx_path)
} else {
  amplicon_network_wb <- openxlsx::createWorkbook()
}

### ---------- 51. network.pip: 网络分析主函数 ----------

library(ggClusterNet)
library(igraph)
# 确保 mia 包不会冲突
if ("package:mia" %in% search()) {
  detach("package:mia", unload = TRUE)
}

tab.r = network.pip(
  ps = ps.16s,
  N = 200,
  big = TRUE,
  select_layout = FALSE,
  layout_net = "model_maptree2",
  r.threshold = 0.6,
  p.threshold = 0.05,
  maxnode = 2,
  label = FALSE,
  lab = "elements",
  group = "Group",
  fill = "Phylum",
  size = "igraph.degree",
  zipi = TRUE,
  ram.net = TRUE,
  clu_method = "cluster_fast_greedy",
  step = 100,
  R = 10,
  ncpus = 1
)

dat = tab.r[[2]]
cortab = dat$net.cor.matrix$cortab

# 提取全部图片的存储对象
plot = tab.r[[1]]

# 提取网络图可视化结果
p0 = plot[[1]]

# 保存网络分析主图
save_plot2(
  p0,
  amplicon_network_path,
  "network_main",
  width = 12,
  height = 10
)

# 保存网络相关矩阵数据
write_sheet2(amplicon_network_wb, "network_cor_matrix", cortab)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 52. net_properties.4: 网络属性计算 ----------

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
write_sheet2(amplicon_network_wb, "network_properties", dat2)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 53. netproperties.sample: 单个样本的网络属性 ----------

for (i in 1:length(id)) {
  pst = ps.16s %>% subset_samples.wt("Group",id[i]) %>% remove.zero()
  dat.f = netproperties.sample(pst = pst,cor = cor[[id[i]]])
  if (i == 1) {
    dat.f2 = dat.f
  } else{
    dat.f2 = rbind(dat.f2,dat.f)
  }
}

map = sample_data(ps.16s) %>% as.tibble()
dat3 = dat.f2 %>% rownames_to_column("ID") %>% inner_join(map,by = "ID")

# 保存样本网络属性数据
write_sheet2(amplicon_network_wb, "sample_network_properties", dat3)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 54. node_properties: 计算节点属性 ----------

for (i in 1:length(id)) {
  igraph= cor[[id[i]]] %>% make_igraph()
  nodepro = node_properties(igraph) %>% as.data.frame()
  nodepro$Group = id[i]
  head(nodepro)
  colnames(nodepro) = paste0(colnames(nodepro),".",id[i])
  nodepro = nodepro %>%
    as.data.frame() %>%
    rownames_to_column("ASV.name")

  if (i == 1) {
    nodepro2 = nodepro
  } else{
    nodepro2 = nodepro2 %>% full_join(nodepro,by = "ASV.name")
  }
}

# 保存节点属性数据
write_sheet2(amplicon_network_wb, "node_properties", nodepro2)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 55. negative.correlation.ratio: 网络稳定性-计算负相关的比例 ----------

res4 = negative.correlation.ratio(ps = ps.16s,
                                  corg = cortab,
                                  degree = TRUE,
                                  zipi = FALSE)

p5 = res4[[1]]
p5 + theme_classic()
dat6 = res4[[2]]

# 保存负相关比例结果
save_plot2(
  p5,
  amplicon_network_path,
  "negative_correlation_ratio",
  width = 10,
  height = 8
)

write_sheet2(amplicon_network_wb, "negative_correlation_ratio", dat6)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 56. community.stability: 网络稳定性-群落稳定性 ----------

# 设置配对样本信息
treat = sample_data(ps.16s)
# treat$pair = paste( "A",c(rep(1:10,3)),sep = "")
# sample_data(ps.16s) = treat

# 群落稳定性分析
res5 = community.stability( ps = ps.16s,
                            corg = cor,
                            time = FALSE)
p6 = res5[[1]]
dat7 = res5[[2]]

# 保存群落稳定性结果
save_plot2(
  p6,
  amplicon_network_path,
  "community_stability",
  width = 10,
  height = 8
)

write_sheet2(amplicon_network_wb, "community_stability", dat7)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 57. natural.con.microp: 网络稳定性-网络抗毁性 ----------

library("pulsar")
res6 = natural.con.microp (
  ps = ps.16s,
  corg = cor,
  norm = TRUE,
  end = 150,
  start = 0
)
p7 = res6[[1]]
dat8  = res6[[2]]

# 保存网络抗毁性结果
save_plot2(
  p7,
  amplicon_network_path,
  "network_robustness",
  width = 10,
  height = 8
)

write_sheet2(amplicon_network_wb, "network_robustness", dat8)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 58. module.compare.net.pip: 网络显著性比较 ----------

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

# 保存网络显著性比较结果
write_sheet2(amplicon_network_wb, "module_compare_results", res)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 59. module.compare.m: 网络稳定性-模块相似性 ----------

library(tidyfst)
res1 = module.compare.m(
  ps = NULL,
  corg = cor,
  zipi = FALSE,
  zoom = 0.6,
  padj = FALSE,
  n = 3)

p1 = res1[[1]]
dat1 = res1[[2]]
dat2 = res1[[3]]

dat2$m1 = dat2$module1 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$m2 = dat2$module2 %>% strsplit("model") %>%
  sapply(`[`, 1)
dat2$cross = paste(dat2$m1,dat2$m2,sep = "_Vs_")
dat2 = dat2 %>% filter(module1 != "none")

p2 = ggplot(dat2) + geom_bar(aes(x = cross,fill = cross)) +
  labs(x = "",
       y = "numbers.of.similar.modules"
  )+ theme_classic()

# 保存模块相似性结果
save_plot2(
  p1,
  amplicon_network_path,
  "module_similarity1",
  width = 10,
  height = 8
)

save_plot2(
  p2,
  amplicon_network_path,
  "module_similarity2",
  width = 10,
  height = 8
)

write_sheet2(amplicon_network_wb, "module_similarity_otu", dat1)
write_sheet2(amplicon_network_wb, "module_similarity_results", dat2)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 510. Robustness.Targeted.removal: 网络稳定性-去除关键节点 ----------

library(patchwork)
res2= Robustness.Targeted.removal(ps = ps.16s,
                                  corg = cor,
                                  degree = TRUE,
                                  zipi = FALSE
)
p3 = res2[[1]]
dat4 = res2[[2]]

# 保存目标移除鲁棒性结果
save_plot2(
  p3,
  amplicon_network_path,
  "robustness_targeted",
  width = 10,
  height = 8
)

write_sheet2(amplicon_network_wb, "robustness_targeted", dat4)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

### ---------- 511. Robustness.Random.removal: 网络稳定性-随机移除节点 ----------

res3 = Robustness.Random.removal(ps = ps.16s,
                                 corg = cortab,
                                 Top = 0
)
p4 = res3[[1]]
dat5 = res3[[2]]

# 保存随机移除鲁棒性结果
save_plot2(
  p4,
  amplicon_network_path,
  "robustness_random",
  width = 10,
  height = 8
)

write_sheet2(amplicon_network_wb, "robustness_random", dat5)
openxlsx::saveWorkbook(amplicon_network_wb, network_xlsx_path, overwrite = TRUE)

# ===================== 6. 群落组装分析（Community Assembly） =====================

# 创建群落组装分析目录
amplicon_assembly_path <- file.path(amplicon_path, "06_community_assembly")
dir.create(amplicon_assembly_path, recursive = TRUE, showWarnings = FALSE)

# 群落组装分析 Excel 文件路径
assembly_xlsx_path <- file.path(amplicon_assembly_path, "community_assembly_results.xlsx")

# 若已有，则在原文件基础上追加；否则新建
if (file.exists(assembly_xlsx_path)) {
  amplicon_assembly_wb <- openxlsx::loadWorkbook(assembly_xlsx_path)
} else {
  amplicon_assembly_wb <- openxlsx::createWorkbook()
}

library(picante)
library(ape)
library(FSA)
library(eulerr)
library(minpack.lm)
library(Hmisc)

# 使用绝对丰度
# otu_table(ps16s) <- round(otu_table(ps16s) * 1000000)
psphy = filter_taxa(ps.16s, function(x) sum(x ) > 1000, TRUE)
map = sample_data(psphy)
n = map$Group %>% unique() %>% length()

### ---------- 61. neutralModel: 中性模型 ----------

result = neutralModel(ps = psphy, group  = "Group", ncol = n)
p43 =  result[[1]]
dat = result[[3]][[1]]
dat2 = result[[4]][[1]]

# 保存中性模型结果
save_plot2(
  p43,
  amplicon_assembly_path,
  "neutral_model",
  width = 12,
  height = 8
)

write_sheet2(amplicon_assembly_wb, "neutral_model_data1", dat)
write_sheet2(amplicon_assembly_wb, "neutral_model_data2", dat2)
openxlsx::saveWorkbook(amplicon_assembly_wb, assembly_xlsx_path, overwrite = TRUE)

# ### ---------- 62. nullModel: 零模型 ----------
#
# result <- nullModel(ps = psphy,
#                     group="Group",
#                     dist.method =  "bray",
#                     gamma.method = "total",
#                     transfer = "none",
#                     null.model = "ecosphere"
# )
#
# nullModeltab <- result[[1]]
# ratiotab <- result[[2]]
# aovtab <- result[[3]]
#
# # 保存零模型结果
# write_sheet2(amplicon_assembly_wb, "null_model_results", nullModeltab)
# write_sheet2(amplicon_assembly_wb, "null_model_ratio", ratiotab)
# write_sheet2(amplicon_assembly_wb, "null_model_aov", aovtab)
# openxlsx::saveWorkbook(amplicon_assembly_wb, assembly_xlsx_path, overwrite = TRUE)

### ---------- 63. bNTICul: β最近分类单元指数计算 ----------

library(parallel)
result = bNTICul(ps = psphy,
                 group  = "Group",
                 num = 10,
                 thread = 1
)
bNTI = result[[1]]

# 保存βNTI结果
write_sheet2(amplicon_assembly_wb, "bNTI_results", bNTI)
openxlsx::saveWorkbook(amplicon_assembly_wb, assembly_xlsx_path, overwrite = TRUE)

### ---------- 64. RCbary: RCbary 计算 ----------

result = RCbary(ps = psphy ,group  = "Group", num = 10, thread = 1)
RCbary = result[[1]]

# 保存RCbray结果
write_sheet2(amplicon_assembly_wb, "RCbray_results", RCbary)
openxlsx::saveWorkbook(amplicon_assembly_wb, assembly_xlsx_path, overwrite = TRUE)

### ---------- 65. bNTIRCPlot: BetaNTI和RCbray联合出图 ----------

RCb = RCbary %>%
  dplyr::mutate(Sample_1 = Site2, Sample_2 = Site1)

result = bNTIRCPlot(ps = psphy ,RCb  =RCb, bNTI = bNTI,group  = "Group")

p3 <- result[[1]]
p4 <- result[[2]]
p5 <- result[[3]]
plotdata = result[[4]]
dat = result[[5]]

# 保存βNTI和RCbray联合分析结果
save_plot2(p3, amplicon_assembly_path, "bNTI_plot", width = 10, height = 8)
save_plot2(p4, amplicon_assembly_path, "RCbray_plot", width = 10, height = 8)
save_plot2(p5, amplicon_assembly_path, "bNTI_RCbray_combined", width = 12, height = 8)

write_sheet2(amplicon_assembly_wb, "bNTI_RCbray_plotdata", plotdata)
write_sheet2(amplicon_assembly_wb, "bNTI_RCbray_summary", dat)
openxlsx::saveWorkbook(amplicon_assembly_wb, assembly_xlsx_path, overwrite = TRUE)


# ===================== 7. 其他分析（Other Analysis） =====================

# 创建其他分析目录
amplicon_other_path <- file.path(amplicon_path, "07_other_analysis")
dir.create(amplicon_other_path, recursive = TRUE, showWarnings = FALSE)

# 其他分析 Excel 文件路径
other_xlsx_path <- file.path(amplicon_other_path, "other_analysis_results.xlsx")

# 若已有，则在原文件基础上追加；否则新建
if (file.exists(other_xlsx_path)) {
  amplicon_other_wb <- openxlsx::loadWorkbook(other_xlsx_path)
} else {
  amplicon_other_wb <- openxlsx::createWorkbook()
}

##71. FEAST.micro: 溯源分析 ----------
data("ps")


result = FEAST.micro (ps = ps,
                     group = "Group",
                     sinkG = "OE",
                     sourceG = c("WT","KO")
)

result


# 保存FEAST溯源分析数据
write_sheet2(amplicon_other_wb, "FEAST_results", result)
openxlsx::saveWorkbook(amplicon_other_wb, other_xlsx_path, overwrite = TRUE)

##72. Plot_FEAST: 溯源分析可视化分组 ----------
p <- Plot_FEAST(data = result)

# 保存FEAST分组可视化结果
save_plot2(p, amplicon_other_path, "FEAST_group", width = 10, height = 8)

##73. MuiPlot_FEAST: 溯源分析可视化样品 ----------
p2 = MuiPlot_FEAST(data = result)

# 保存FEAST样品可视化结果
save_plot2(p2, amplicon_other_path, "FEAST_sample", width = 10, height = 8)


# ===================== 8. 功能预测分析（Function Prediction） =====================

# 创建功能预测分析目录
amplicon_function_path <- file.path(amplicon_path, "08_function_prediction")
dir.create(amplicon_function_path, recursive = TRUE, showWarnings = FALSE)

# 功能预测分析 Excel 文件路径
function_xlsx_path <- file.path(amplicon_function_path, "function_prediction_results.xlsx")

# 若已有，则在原文件基础上追加；否则新建
if (file.exists(function_xlsx_path)) {
  amplicon_function_wb <- openxlsx::loadWorkbook(function_xlsx_path)
} else {
  amplicon_function_wb <- openxlsx::createWorkbook()
}

library(DOSE)
library(GO.db)
library(GSEABase)
library(ggtree)
library(aplot)
library(clusterProfiler)
library("GSVA")
library(dplyr)

ps.kegg = ps.kegg %>% filter_OTU_ps(Top = 1000)
tax = ps.kegg %>% phyloseq::tax_table()
colnames(tax)[3] = "KOid"
tax_table(ps.kegg) = tax

### ---------- 81. EdgerSuper.metf: 差异分析 ----------

res = EdgerSuper.metf (ps = ps.kegg,
                       group  = "Group",
                       artGroup = NULL)
dat = res[[2]]

# 保存差异分析结果
write_sheet2(amplicon_function_wb, "differential_analysis", dat)
openxlsx::saveWorkbook(amplicon_function_wb, function_xlsx_path, overwrite = TRUE)

### ---------- 82. KEGG_enrich.micro: taxfun2功能富集分析 ----------

res2 = KEGG_enrich.micro(ps =  ps.kegg, dif = dat)
dat1= res2$`KO-OE`
dat2= res2$`KO-WT`
dat3= res2$`OE-WT`

# 保存KEGG富集分析结果
write_sheet2(amplicon_function_wb, "KEGG_enrich_KO_OE", dat1)
write_sheet2(amplicon_function_wb, "KEGG_enrich_KO_WT", dat2)
write_sheet2(amplicon_function_wb, "KEGG_enrich_OE_WT", dat3)
openxlsx::saveWorkbook(amplicon_function_wb, function_xlsx_path, overwrite = TRUE)

### ---------- 83. buplotall.micro: taxfun2功能富集分析气泡图 ----------

Desep_group <- ps.kegg %>% sample_data() %>%
  .$Group %>%
  as.factor() %>%
  levels() %>%
  as.character()
cbtab = combn(Desep_group,2)
Desep_group = cbtab[,1]
group = paste(Desep_group[1], Desep_group[2], sep = "-")
id = paste(group,"level",sep = "")

result = buplot.micro(dt=dat1, id = id)
p1 = result[[1]]
p2 = result[[2]]

# 保存功能富集气泡图结果
save_plot2(p1, amplicon_function_path, "function_bubble_plot1", width = 10, height = 8)
save_plot2(p2, amplicon_function_path, "function_bubble_plot2", width = 10, height = 8)
