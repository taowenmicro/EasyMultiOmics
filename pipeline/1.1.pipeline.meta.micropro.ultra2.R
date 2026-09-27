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

# 方式一：交互选择
# file_path <- tcltk::tk_choose.files(
#   caption = "请选择 ps_ITS.rds 文件",
#   multi   = FALSE,
#   filters = matrix(c("RDS files", ".rds",
#                      "All files", "*"), ncol = 2, byrow = TRUE)
# )
# 方式二：直接指定路径
file_path <- "../思农生信21个宏基因组西安交大/data.ps/phyloseq_species_7level_count_grouped.rds"
ps.micro    <- readRDS(file_path)
ps.micro

res <- infer_sequencing_type(ps = ps.micro)  # ps 是你的 phyloseq 对象
res$label
sample_sums(ps.micro)
map <- sample_data(ps.micro)
head(map)



### ===================== 1.1 一次性创建全部输出目录与工作簿 =====================
mg_path <- create_omics_result_dir_auto(
  ps = ps.micro,
  base_dir     = "../思农生信21个宏基因组西安交大/result2/",
  include_time = FALSE
)
# 所有模块的目录在流程最前面一次性创建。
# 后面即使只运行某个分析模块，也不会因为目标文件夹不存在而保存失败。
mg_alpha_path     <- file.path(mg_path, "01_alpha_diversity")
mg_beta_path      <- file.path(mg_path, "02_beta_diversity")
mg_comp_path      <- file.path(mg_path, "03_composition")
mg_diff_path      <- file.path(mg_path, "04_differential")
mg_biomarker_path <- file.path(mg_path, "05_biomarker")
mg_network_path   <- file.path(mg_path, "06_network")

all_result_dirs <- c(
  mg_path, mg_alpha_path, mg_beta_path, mg_comp_path,
  mg_diff_path, mg_biomarker_path, mg_network_path
)
invisible(lapply(all_result_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

alpha_xlsx_path     <- file.path(mg_alpha_path,     "alpha_diversity_results.xlsx")
beta_xlsx_path      <- file.path(mg_beta_path,      "beta_diversity_results.xlsx")
comp_xlsx_path      <- file.path(mg_comp_path,      "composition_results.xlsx")
diff_xlsx_path      <- file.path(mg_diff_path,      "differential_results.xlsx")
biomarker_xlsx_path <- file.path(mg_biomarker_path, "biomarker_results.xlsx")
network_xlsx_path   <- file.path(mg_network_path,   "network_results.xlsx")

load_or_create_wb <- function(path) {
  if (file.exists(path)) openxlsx::loadWorkbook(path) else openxlsx::createWorkbook()
}

mg_alpha_wb     <- load_or_create_wb(alpha_xlsx_path)
mg_beta_wb      <- load_or_create_wb(beta_xlsx_path)
mg_comp_wb      <- load_or_create_wb(comp_xlsx_path)
mg_diff_wb      <- load_or_create_wb(diff_xlsx_path)
mg_biomarker_wb <- load_or_create_wb(biomarker_xlsx_path)
mg_network_wb   <- load_or_create_wb(network_xlsx_path)

# 通用安全执行：某个可选模块失败时记录错误，不让整个流程中断。
run_safe <- function(label, expr) {
  tryCatch(
    expr,
    error = function(e) {
      message("[", label, "] ERROR: ", conditionMessage(e))
      structure(list(error = conditionMessage(e)), class = "analysis_error")
    }
  )
}

is_error_result <- function(x) inherits(x, "analysis_error")

# 将任意结果转换为可写入 Excel 的 data.frame。
as_df_safe <- function(x, value_name = "Value") {
  if (is.null(x)) return(data.frame())
  if (is.data.frame(x)) return(x)
  if (is.matrix(x)) return(as.data.frame(x))
  if (is.atomic(x)) {
    z <- data.frame(x, stringsAsFactors = FALSE)
    colnames(z)[1] <- value_name
    return(z)
  }
  data.frame()
}

# 保存结果对象中所有真正的绘图对象，避免只保存前3张。
# data.frame / matrix / 普通向量不会递归，防止大型结果表被逐列扫描。
is_plot_object <- function(x) {
  inherits(x, c("ggplot", "gg", "grob", "gtable", "patchwork", "ggarrange"))
}

save_all_plot_components <- function(x, out_dir, prefix, width = 10, height = 8) {
  if (is.null(x)) return(invisible(NULL))

  if (is_plot_object(x)) {
    try(save_plot2(x, out_dir, prefix, width = width, height = height), silent = TRUE)
    return(invisible(NULL))
  }

  if (!is.list(x) || is.data.frame(x)) return(invisible(NULL))

  nm <- names(x)
  for (i in seq_along(x)) {
    obj <- x[[i]]
    tag <- if (!is.null(nm) && length(nm) >= i && !is.na(nm[i]) && nzchar(nm[i])) {
      nm[i]
    } else {
      paste0("plot", i)
    }
    tag <- gsub("[^A-Za-z0-9._-]", "_", tag)
    save_all_plot_components(
      obj, out_dir, paste0(prefix, "_", tag),
      width = width, height = height
    )
  }
  invisible(NULL)
}

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
## ===================== 2. Alpha 多样性分析 =====================

mg_alpha_path <- file.path(mg_path, "01_alpha_diversity")
dir.create(mg_alpha_path, recursive = TRUE, showWarnings = FALSE)
alpha_xlsx_path <- file.path(mg_alpha_path, "alpha_diversity_results.xlsx")
mg_alpha_wb <- if (file.exists(alpha_xlsx_path)) openxlsx::loadWorkbook(alpha_xlsx_path) else openxlsx::createWorkbook()

# ---------- Alpha diversity ----------
all.alpha <- c("Shannon", "Inv_Simpson", "Pielou_evenness",
               "Simpson_evenness", "Richness", "Chao1", "ACE")

tab <- alpha.metm(ps = ps.micro, group = "Group")
data_alpha <- cbind(
  data.frame(ID = as.character(seq_len(nrow(tab))), group = tab$Group),
  tab[all.alpha]
)


alpha_res <- alpha_auto_stat.metm(
  data = data_alpha, num = 3:ncol(data_alpha),
  group_col = "group", alpha = 0.05, fallback = TRUE
)

data_alpha_use <- alpha_res$data
result_alpha   <- alpha_res$result

alpha_res$summary
alpha_res$letters


alpha_theme <- list(
  scale_x_discrete(limits = axis_order),
  scale_fill_manual(values = col.g),
  theme_nature(),
  guides(fill = guide_legend(title = NULL))
)

res_box <- EasyMultiOmics::FacetMuiPlotresultBox(
  data = data_alpha_use, num = 3:ncol(data_alpha_use),
  result = result_alpha, sig_show = "abc", ncol = 3, width = 0.4
)
p1_1 <- res_box[[1]] + alpha_theme


res_bar <- EasyMultiOmics::FacetMuiPlotresultBar(
  data = data_alpha_use, num = 3:ncol(data_alpha_use),
  result = result_alpha, sig_show = "abc", ncol = 3, mult.y = 0.3
)
p1_2 <- res_bar[[1]] + alpha_theme


res_boxbar <- EasyMultiOmics::FacetMuiPlotReBoxBar(
  data = data_alpha_use, num = 3:ncol(data_alpha_use),
  result = result_alpha, sig_show = "abc", ncol = 3,
  mult.y = 0.3, lab.yloc = 1.1
)
p1_3 <- res_boxbar[[1]] + alpha_theme


lab_alpha <- res_box[[2]] %>% dplyr::distinct(group, name, .keep_all = TRUE)

p1_0 <- res_box[[2]] %>%
  ggplot(aes(group, dd, fill = group)) +
  geom_violin(alpha = 1) +
  geom_jitter(aes(color = group), width = 0.17, size = 3, alpha = 0.5) +
  geom_text(data = lab_alpha, aes(group, y, label = stat), alpha = 0.6) +
  facet_wrap(~name, scales = "free_y", ncol = 3) +
  labs(x = NULL, y = NULL) +
  scale_x_discrete(limits = axis_order) +
  scale_fill_manual(values = col.g) +
  scale_color_manual(values = col.g) +
  theme_nature() +
  guides(fill = guide_legend(title = NULL))

plots_alpha <- list(
  alpha_diversity_box     = p1_1,
  alpha_diversity_bar     = p1_2,
  alpha_diversity_boxbar  = p1_3,
  alpha_diversity_violin  = p1_0
)

for (nm in names(plots_alpha))
  save_plot2(plots_alpha[[nm]], mg_alpha_path, nm, width = 24, height = 8)

write_sheet2(mg_alpha_wb, "alpha_diversity_data", data_alpha_use)
write_sheet2(mg_alpha_wb, "alpha_diversity_stat", as.data.frame(result_alpha))
write_sheet2(mg_alpha_wb, "alpha_test_method", alpha_res$summary)
write_sheet2(mg_alpha_wb, "alpha_games_howell", alpha_res$games_howell)
write_sheet2(mg_alpha_wb, "alpha_final_letters", alpha_res$letters)

if (nrow(alpha_res$invalid_metrics))
  write_sheet2(mg_alpha_wb, "alpha_invalid_metrics", alpha_res$invalid_metrics)

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

# ---- 保存稀释曲线图
save_plot2(p2_1, mg_alpha_path, "alpha_rarefaction_individual", width = 9,  height = 7)
save_plot2(p2_2, mg_alpha_path, "alpha_rarefaction_group",      width = 10, height = 8)
save_plot2(p2_3, mg_alpha_path, "alpha_rarefaction_group_sd",   width = 10, height = 8)

# ---- 保存稀释数据
write_sheet2(mg_alpha_wb, "rarefaction_data", raretab)
openxlsx::saveWorkbook(mg_alpha_wb, alpha_xlsx_path, overwrite = TRUE)


## ===================== 3. Beta 多样性分析 =====================

mg_beta_path <- file.path(mg_path, "02_beta_diversity")
dir.create(mg_beta_path, recursive = TRUE, showWarnings = FALSE)

beta_xlsx_path <- file.path(mg_beta_path, "beta_diversity_results.xlsx")
mg_beta_wb <- if (file.exists(beta_xlsx_path)) {
  openxlsx::loadWorkbook(beta_xlsx_path)
} else {
  openxlsx::createWorkbook()
}


### ---------- 3.1 PCoA + 95% CI ----------

res_ord <- ordinate.metm2(
  ps = ps.micro,
  group = "Group",
  dist = "bray",
  method = "PCoA",
  Micromet = "anosim",
  pvalue.cutoff = 0.05,
  ellipse = TRUE,
  ellipse.level = 0.25,
  ellipse.type = "confidence"
)

## 统一绘图风格
ord_theme <- function(p) {
  p +
    scale_fill_manual(values = col.g) +
    scale_color_manual(values = col.g, guide = "none") +
    theme_nature() +
    theme(axis.title.y = element_text(angle = 90))
}

p3_1 <- ord_theme(res_ord[[1]])
p3_2 <- ord_theme(res_ord[[3]])

## 排序坐标、群心和蜘蛛线
plotdata <- as.data.frame(res_ord[[2]])

cent <- aggregate(cbind(x, y) ~ Group, plotdata, mean)

segs <- dplyr::left_join(
  plotdata,
  dplyr::rename(cent, center_x = x, center_y = y),
  by = "Group"
)

## 精修图：95% CI + spider + centroid
p3_3 <- p3_1 +
  geom_segment(
    data = segs,
    aes(x, y, xend = center_x, yend = center_y, color = Group),
    inherit.aes = FALSE,
    linewidth = 0.6,
    alpha = 0.65,
    show.legend = FALSE
  ) +
  geom_point(
    data = cent,
    aes(x, y),
    inherit.aes = FALSE,
    shape = 24,
    size = 5,
    color = "black",
    fill = "yellow"
  )

## 保存图片
save_plot2(p3_1, mg_beta_path, "pcoa_basic_95CI",   12, 10)
save_plot2(p3_2, mg_beta_path, "pcoa_labeled_95CI", 12, 10)
save_plot2(p3_3, mg_beta_path, "pcoa_refined_95CI", 12, 10)

## 保存数据
write_sheet2(mg_beta_wb, "ordination_data", plotdata)
write_sheet2(mg_beta_wb, "ordination_cent", cent)
write_sheet2(mg_beta_wb, "ordination_segs", segs)

if (!is.null(res_ord$ellipse))
  write_sheet2(mg_beta_wb, "ordination_95CI", res_ord$ellipse)

openxlsx::saveWorkbook(
  mg_beta_wb,
  beta_xlsx_path,
  overwrite = TRUE
)

### ---------- 3.2 MicroTest.metm：整体群落差异 (adonis) ----------

dat_adonis <- MicroTest.metm(ps = ps.micro, Micromet = "adonis", dist = "bray")
write_sheet2(mg_beta_wb, "adonis_results", as.data.frame(dat_adonis))
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

res_mantel <- EasyMultiOmics::mantal.micro(
  ps     = ps.micro,
  method = "spearman",
  group  = "Group",
  ncol   = 3,
  nrow   = res$n_param
)
data_mantel <- res_mantel[[1]]
p3_7        <- res_mantel[[2]][[1]]

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



### ---------- 3.6 Procrustes：Bray-Curtis 与 Jaccard 排序一致性 ----------
# 原脚本没有实际的 Procrustes 模块。这里比较同一批样本的
# abundance-based Bray-Curtis 与 presence/absence-based Jaccard 排序。
# 若后续要比较“物种 vs 功能”等两个不同对象，可把两个坐标矩阵替换即可。
res_procrustes <- run_safe("Procrustes", {
  ps_beta <- phyloseq::transform_sample_counts(ps.micro, function(x) x / sum(x))

  d_bray <- phyloseq::distance(ps_beta, method = "bray")
  d_jacc <- phyloseq::distance(ps_beta, method = "jaccard", binary = TRUE)

  ord_bray <- stats::cmdscale(d_bray, k = 2, eig = TRUE, add = TRUE)$points
  ord_jacc <- stats::cmdscale(d_jacc, k = 2, eig = TRUE, add = TRUE)$points

  ids_proc <- Reduce(intersect, list(rownames(ord_bray), rownames(ord_jacc), sample_names(ps_beta)))
  ord_bray <- ord_bray[ids_proc, , drop = FALSE]
  ord_jacc <- ord_jacc[ids_proc, , drop = FALSE]

  proc_fit  <- vegan::procrustes(ord_bray, ord_jacc, symmetric = TRUE)
  proc_test <- vegan::protest(ord_bray, ord_jacc, permutations = 999, symmetric = TRUE)

  map_proc <- as.data.frame(sample_data(ps_beta))[ids_proc, , drop = FALSE]
  proc_dat <- data.frame(
    ID = ids_proc,
    Group = map_proc$Group,
    Bray1 = proc_fit$X[, 1],
    Bray2 = proc_fit$X[, 2],
    Jaccard1_rot = proc_fit$Yrot[, 1],
    Jaccard2_rot = proc_fit$Yrot[, 2],
    stringsAsFactors = FALSE
  )

  proc_stat <- data.frame(
    Statistic = c("Procrustes_correlation", "Permutation_P", "Permutations"),
    Value = c(
      if (!is.null(proc_test$t0)) proc_test$t0 else NA_real_,
      if (!is.null(proc_test$signif)) proc_test$signif else NA_real_,
      999
    )
  )

  p_proc <- ggplot(proc_dat) +
    geom_segment(
      aes(x = Bray1, y = Bray2, xend = Jaccard1_rot, yend = Jaccard2_rot, colour = Group),
      linewidth = 0.6, alpha = 0.65, show.legend = FALSE
    ) +
    geom_point(aes(Bray1, Bray2, fill = Group), shape = 21, size = 4, colour = "black") +
    geom_point(aes(Jaccard1_rot, Jaccard2_rot, colour = Group), shape = 4, size = 3, stroke = 1) +
    scale_fill_manual(values = col.g) +
    scale_colour_manual(values = col.g) +
    theme_nature() +
    labs(
      x = "Procrustes axis 1", y = "Procrustes axis 2",
      title = paste0(
        "Bray vs Jaccard Procrustes: r = ",
        round(if (!is.null(proc_test$t0)) proc_test$t0 else NA_real_, 3),
        ", P = ", signif(if (!is.null(proc_test$signif)) proc_test$signif else NA_real_, 3)
      )
    )

  save_plot2(p_proc, mg_beta_path, "procrustes_bray_jaccard", width = 10, height = 8)
  write_sheet2(mg_beta_wb, "procrustes_coordinates", proc_dat)
  write_sheet2(mg_beta_wb, "procrustes_test", proc_stat)
  openxlsx::saveWorkbook(mg_beta_wb, beta_xlsx_path, overwrite = TRUE)

  list(plot = p_proc, data = proc_dat, test = proc_stat, fit = proc_fit)
})


## ===================== 4. 组成（Composition）分析 =====================
## ===================== 4. Composition =====================

mg_comp_path <- file.path(mg_path, "03_composition")
dir.create(mg_comp_path, recursive = TRUE, showWarnings = FALSE)

tax_ranks <- intersect(
  c("Phylum", "Class", "Order", "Family", "Genus", "Species"),
  phyloseq::rank_names(ps.micro)
)

#
comp_dirs <- c(
  venn       = "01_venn_upset",
  vensuper   = "02_vensuper",
  flower     = "03_flower",
  vnet       = "04_venn_network",
  ternary    = "05_ternary",
  polygon    = "06_polygon",
  bar        = "07_barplot",
  multilevel = "08_multilevel",
  cluster    = "09_cluster_barplot",
  abundance  = "10_abundance_class",
  circular   = "11_circular_barplot",
  chord      = "12_chord_plot",
  heatmap    = "13_heatmap"
)

comp_dirs <- file.path(mg_comp_path, comp_dirs)
names(comp_dirs) <- c(
  "venn","vensuper","flower","vnet","ternary","polygon",
  "bar","multilevel","cluster","abundance","circular","chord","heatmap"
)
invisible(lapply(comp_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

comp_error <- data.frame(
  Module = character(), Rank = character(), Error = character()
)

log_comp_error <- function(module, rank, msg) {
  comp_error <<- rbind(
    comp_error,
    data.frame(Module = module, Rank = rank, Error = as.character(msg))
  )
  message("[SKIP] ", module, " / ", rank, ": ", msg)
}

safe_comp <- function(module, rank, expr) {
  tryCatch(
    expr,
    error = function(e) {
      log_comp_error(module, rank, conditionMessage(e))
      NULL
    }
  )
}

safe_plot_comp <- function(p, module, rank, name, w = 10, h = 8) {
  if (is.null(p)) return(invisible(NULL))
  tryCatch(
    save_plot2(
      p, comp_dirs[[module]],
      paste0(rank, "_", name), width = w, height = h
    ),
    error = function(e) log_comp_error(module, rank, conditionMessage(e))
  )
}

safe_write_comp <- function(module, rank, sheet, dat) {
  if (is.null(dat)) return(invisible(NULL))

  tryCatch({
    if (is.list(dat) && !is.data.frame(dat))
      dat <- dplyr::bind_rows(dat, .id = "Part")
    dat <- as_df_safe(dat)
    if (!nrow(dat) && !ncol(dat)) return(invisible(NULL))

    f <- file.path(comp_dirs[[module]], "data.xlsx")
    wb <- if (file.exists(f)) openxlsx::loadWorkbook(f) else openxlsx::createWorkbook()

    nm <- substr(paste0(rank, "_", sheet), 1, 31)
    write_sheet2(wb, nm, dat)
    openxlsx::saveWorkbook(wb, f, overwrite = TRUE)
  }, error = function(e) {
    log_comp_error(module, rank, conditionMessage(e))
  })
}

make_rank_ps <- function(ps, rank) {

  tx <- as.data.frame(phyloseq::tax_table(ps), stringsAsFactors = FALSE)
  x <- trimws(as.character(tx[[rank]]))

  keep <- !is.na(x) &
    nzchar(x) &
    !tolower(x) %in% c("unassigned", "unknown", "na", "unclassified")

  ps0 <- phyloseq::prune_taxa(rownames(tx)[keep], ps)
  if (phyloseq::ntaxa(ps0) < 2)
    stop("有效分类单元不足")

  psr <- phyloseq::tax_glom(ps0, taxrank = rank, NArm = TRUE)
  psr <- phyloseq::prune_taxa(phyloseq::taxa_sums(psr) > 0, psr)

  if (phyloseq::ntaxa(psr) < 2)
    stop("聚合后有效分类单元不足")

  psr
}



for (rank in tax_ranks) {

  message("\n================ ", rank, " ================")

  psr <- safe_comp(
    "prepare", rank,
    make_rank_ps(ps.micro, rank)
  )
  if (is.null(psr)) next


  ## ---------- 4.1 Venn / UpSet ----------
  x <- safe_comp(
    "venn", rank,
    Ven.Upset.metm(ps = psr, group = "Group", N = 0.5, size = 3)
  )

  if (!is.null(x)) {
    if (length(x) >= 1)
      safe_plot_comp(x[[1]] + theme_void(), "venn", rank, "venn", 10, 8)
    if (length(x) >= 2)
      safe_plot_comp(x[[2]], "venn", rank, "upset", 10, 8)
    if (length(x) >= 4)
      safe_plot_comp(x[[4]], "venn", rank, "extra", 10, 8)
    if (length(x) >= 3)
      safe_write_comp("venn", rank, "data", x[[3]])
  }


  ## ---------- 4.2 VenSuper ----------
  x <- safe_comp(
    "vensuper", rank,
    VenSuper.metm(
      ps = psr,
      group = "Group",
      num = min(
        6,
        length(unique(as.character(
          phyloseq::sample_data(psr)$Group
        )))
      )
    )
  )

  if (!is.null(x)) {

    safe_comp(
      "vensuper", rank,
      save_all_plot_components(
        x[seq_len(min(3, length(x)))],
        comp_dirs[["vensuper"]],
        paste0(rank, "_vensuper"),
        width = 10, height = 8
      )
    )

    if (length(x) >= 4)
      safe_write_comp("vensuper", rank, "detail", x[[4]])

  } else {

    ## 自动退回 UpSet
    xf <- safe_comp(
      "vensuper", rank,
      Ven.Upset.metm(ps = psr, group = "Group", N = 0.5, size = 3)
    )

    if (!is.null(xf)) {
      safe_comp(
        "vensuper", rank,
        save_all_plot_components(
          xf[c(1, 2, 4)],
          comp_dirs[["vensuper"]],
          paste0(rank, "_fallback"),
          width = 10, height = 8
        )
      )
      if (length(xf) >= 3)
        safe_write_comp("vensuper", rank, "fallback", xf[[3]])
    }
  }


  ## ---------- 4.3 Flower ----------
  flower_group <- if (
    "ID" %in% phyloseq::sample_variables(psr)
  ) "ID" else "Group"

  x <- safe_comp(
    "flower", rank,
    ggflower.metm(
      ps = psr, group = flower_group,
      start = 1, m1 = 2, a = 0.3, b = 1,
      lab.leaf = 1, col.cir = "yellow", N = 0.5
    )
  )

  if (!is.null(x)) {
    safe_plot_comp(x[[1]], "flower", rank, "flower", 10, 8)
    if (length(x) >= 2)
      safe_write_comp("flower", rank, "data", x[[2]])
  }


  ## ---------- 4.4 Venn network ----------
  x <- safe_comp(
    "vnet", rank,
    ven.network.metm(ps = psr, N = 0.5, fill = rank)
  )

  if (!is.null(x)) {
    safe_plot_comp(x[[1]], "vnet", rank, "network", 12, 10)
    if (length(x) >= 2)
      safe_write_comp("vnet", rank, "data", x[[2]])
  }


  ## ---------- 4.5 Ternary ----------
  x <- safe_comp(
    "ternary", rank,
    Micro_tern.metm(
      ggClusterNet::filter_OTU_ps(psr, 500)
    )
  )

  if (!is.null(x)) {
    safe_plot_comp(
      x[[1]][[1]] + theme_classic(),
      "ternary", rank, "ternary", 12, 10
    )
    if (length(x) >= 2)
      safe_write_comp("ternary", rank, "data", x[[2]])
  }


  ## ---------- 4.6 Polygon ----------
  x <- safe_comp(
    "polygon", rank,
    ps_polygon_plot(
      ggClusterNet::filter_OTU_ps(psr, 500),
      group = "Group",
      taxrank = rank
    )
  )

  safe_plot_comp(x, "polygon", rank, "polygon", 18, 12)


  ## ---------- 4.7 普通堆积柱 ----------
  x <- safe_comp(
    "bar", rank,
    barMainplot.metm(
      ps = psr,
      j = rank,
      label = FALSE,
      sd = FALSE,
      Top = 10
    )
  )

  if (!is.null(x)) {

    p1 <- x[[1]] +
      scale_fill_manual(values = colset2) +
      theme_nature() +
      theme(axis.title.y = element_text(angle = 90))

    p2 <- x[[3]] +
      scale_fill_manual(values = colset2) +
      theme_nature() +
      theme(axis.title.y = element_text(angle = 90))

    safe_plot_comp(p1, "bar", rank, "samples", 12, 8)
    safe_plot_comp(p2, "bar", rank, "groups", 10, 8)

    safe_write_comp("bar", rank, "raw", x[[2]])

    dat <- safe_comp(
      "bar", rank,
      x[[2]] |>
        dplyr::group_by(Group, aa) |>
        dplyr::summarise(
          Abundance = sum(Abundance),
          .groups = "drop"
        )
    )

    if (!is.null(dat)) {
      colnames(dat)[1:3] <- c(
        "Group", rank, "Abundance_percent"
      )
      safe_write_comp("bar", rank, "summary", dat)
    }
  }


  ## ---------- 4.8 聚类堆积柱 ----------
  x <- safe_comp(
    "cluster", rank,
    cluMicro.bar.metm(
      dist = "bray",
      ps = psr,
      j = rank,
      Top = 6,
      tran = TRUE,
      hcluter_method = "complete",
      Group = "Group",
      cuttree = length(unique(
        phyloseq::sample_data(psr)$Group
      )),
      group_colors = col.g
    )
  )

  if (!is.null(x)) {
    if (length(x) >= 2)
      safe_plot_comp(x[[2]], "cluster", rank, "plot1", 10, 8)
    if (length(x) >= 3)
      safe_plot_comp(x[[3]], "cluster", rank, "plot2", 10, 8)
    if (length(x) >= 4)
      safe_plot_comp(x[[4]], "cluster", rank, "plot3", 10, 8)
    if (length(x) >= 5)
      safe_write_comp("cluster", rank, "data", x[[5]])
  }


  ## ---------- 4.9 高 / 中 / 低丰度 ----------
  x <- safe_comp(
    "abundance", rank,
    stacked_bar_by_custom_class_v(
      ps = psr,
      rank = rank,
      group = "Group",
      thresholds = c(1, 0.1),
      n_each = c(6, 6, 6),
      arrange = "row"
    )
  )

  if (!is.null(x) && !is.null(x$combined))
    safe_plot_comp(
      x$combined,
      "abundance", rank,
      "high_medium_low",
      18, 8
    )


  ## ---------- 4.10 环状柱状图 ----------
  x <- safe_comp(
    "circular", rank,
    cir_barplot.metm(
      ps = psr,
      Top = 15,
      dist = "bray",
      cuttree = 3,
      hcluter_method = "complete"
    )
  )

  if (!is.null(x)) {
    safe_plot_comp(
      x[[1]], "circular", rank,
      "circular", 16, 12
    )
    if (length(x) >= 2)
      safe_write_comp("circular", rank, "data", x[[2]])
  }


  ## 多圈高 / 中 / 低丰度
  x2 <- safe_comp(
    "circular", rank,
    cir_barplot.metm2(
      ps = psr,
      Top = 20,
      rank = rank,
      thresholds = c(1, 0.8),
      n_each = c(10, 8, 8),
      xlim_high = c(0, 100),
      xlim_medium = c(0, 20),
      xlim_low = c(0, 5),
      breaks_high = seq(0, 100, 10),
      breaks_medium = seq(0, 20, 5),
      breaks_low = seq(0, 5, 1),
      ring_pwidth = c(2, 5, 5),
      base_offset = 2,
      ring_gap = 0,
      axis_pad = 0,
      show_axes = c(FALSE, FALSE, FALSE),
      legend_position = "right"
    )
  )

  if (!is.null(x2) && !is.null(x2$plot))
    safe_plot_comp(
      x2$plot,
      "circular", rank,
      "multiring",
      24, 20
    )


  ## ---------- 4.11 Chord ----------
  rank_id <- match(rank, phyloseq::rank_names(ps.micro))

  x <- safe_comp(
    "chord", rank,
    cir_plot.metm(
      ps = ps.micro,
      Top = 12,
      rank = rank_id
    )
  )

  if (!is.null(x)) {

    if (length(x) >= 1)
      safe_write_comp(
        "chord", rank,
        "data",
        as.data.frame(x[[1]])
      )

    safe_comp(
      "chord", rank,
      save_circlize_plot(
        cir_plot.metm(
          ps = ps.micro,
          Top = 12,
          rank = rank_id
        ),
        out_dir = comp_dirs[["chord"]],
        prefix = paste0(rank, "_chord"),
        width = 14,
        height = 14,
        res = 300
      )
    )
  }


  ## ---------- 4.12 Heatmap ----------
  ## ---------- 4.12 Heatmap ----------
  x <- safe_comp(
    "heatmap", rank, {

      ps_tem <- psr |>
        ggClusterNet::scale_micro(method = "TMM")

      id_heat <- ps_tem |>
        ggClusterNet::filter_OTU_ps(100) |>
        ggClusterNet::vegan_otu() |>
        t() |>
        as.data.frame() |>
        rowCV() |>
        sort(decreasing = TRUE) |>
        head(30) |>
        names()

      if (!length(id_heat))
        stop("无可用于 heatmap 的分类单元")

      res <- Microheatmap.metm(
        ps_rela = ps_tem,
        id = id_heat,
        col_cluster = FALSE,
        row_cluster = FALSE,
        col1 = colorRampPalette(
          c("#00CED1", "#FFFFFF", "#FF4500")
        )(100),
        col.g = col.g,
        y_text_size = 8
      )

      list(res = res, ids = id_heat)
    }
  )

  if (!is.null(x)) {
    safe_plot_comp(x$res[[1]], "heatmap", rank, "heatmap1", 12, 10)
    safe_plot_comp(x$res[[2]], "heatmap", rank, "heatmap2", 12, 10)

    safe_write_comp(
      "heatmap", rank, "data",
      x$res[[3]]
    )

    safe_write_comp(
      "heatmap", rank, "selected",
      data.frame(Taxon_ID = x$ids)
    )
  }

  if (!is.null(x)) {
    safe_plot_comp(
      x$res[[1]],
      "heatmap", rank,
      "heatmap1", 12, 10
    )
    safe_plot_comp(
      x$res[[2]],
      "heatmap", rank,
      "heatmap2", 12, 10
    )

    safe_write_comp(
      "heatmap", rank,
      "data",
      x$res[[3]]
    )

    safe_write_comp(
      "heatmap", rank,
      "selected",
      data.frame(Taxon_ID = x$ids)
    )
  }
}



x <- safe_comp(
  "multilevel", "AllRanks",
  barMultiLevel.metm(
    ps = ps.micro,
    group = "Group",
    Top = 10,
    label_threshold = 10
  )
)

if (!is.null(x)) {
  if (length(x) >= 1)
    safe_plot_comp(
      x[[1]], "multilevel",
      "AllRanks", "plot1", 10, 8
    )
  if (length(x) >= 3)
    safe_plot_comp(
      x[[3]], "multilevel",
      "AllRanks", "plot2", 10, 8
    )
  if (length(x) >= 4)
    safe_plot_comp(
      x[[4]], "multilevel",
      "AllRanks", "plot3", 10, 8
    )
}


write.table(
  comp_error,
  file.path(mg_comp_path, "composition_error_log.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

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
otu = ps.micro %>% vegan_otu() %>% t() %>% as.data.frame()
otu[1:5,1:5]


map = ps.micro %>% sample_data()
head(map)


res_edger <- EdgerSuper.metm(
  ps       = ps.micro %>% remove.zero(),
  group    = "Group",
  artGroup = NULL,
  j        = "Species",
  col.g    = col.g,
  gradient = TRUE,
  top_n    = 5
)


dat25 <- res_edger[[2]]
save_all_plot_components(res_edger[[1]], mg_diff_path, "edger_volcano", width = 10, height = 8)
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
  ps        = ps.micro ,
  group     = "Group",
  artGroup  = NULL,
  j         = "Species",
  col.g     = col.g,
  gradient  = TRUE,
  top_n     = 5
)


dat26 <- res_deseq[[2]]
save_all_plot_components(res_deseq[[1]], mg_diff_path, "deseq2_volcano", width = 10, height = 8)
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
  ps       = ps.micro ,
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
  ps.cs <- clean_tax_table_symbols(ps.cs)

  res_all <- diff_all_methods(
    ps.cs %>% filter_OTU_ps(500),
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

  res_all$summary$methods
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
library(caret)
library(pROC)
library(zCompositions)
library(compositions)
library(xgboost)
library(Ckmeans.1d.dp)
library(mia)
library(rpart)
library(ipred)
library(randomForest)
library(caret)
library(ROCR)
library(e1071)

## ===================== 6. Biomarker Identification (compact) =====================
## Requires:
## loadingPCA.metm2, LDA.metm2, randomforest.metm2
## svm_metm2, glm.metm2, lasso.metm2, xgboost.metm2, decisiontree.metm2
## bagging_metm2, naivebayes.metm2, nnet.metm2, Roc.metm2
## run_safe, is_error_result, as_df_safe, save_plot2, write_sheet2

## ---------- 6.0 Settings ----------
bp <- file.path(mg_path, "05_biomarker4"); dir.create(bp, TRUE, FALSE)
paths <- list(
  PCA=file.path(bp,"01_PCA"), LDA=file.path(bp,"02_LDA"),
  RF=file.path(bp,"03_RandomForest_multigroup"),
  PERM=file.path(bp,"04_Pairwise_PERMANOVA"),
  ML=file.path(bp,"05_Pairwise_ML")
)
invisible(lapply(paths, dir.create, recursive=TRUE, showWarnings=FALSE))

xlsx <- file.path(bp,"biomarker_results.xlsx")
wb <- if (file.exists(xlsx)) openxlsx::loadWorkbook(xlsx) else openxlsx::createWorkbook()

perm_n <- 999; perm_seed <- 11; padj_cut <- 0.05; fallback_n <- 3L
ml_prev <- 0.30; ml_mean <- 0.01; ml_trans <- "hellinger"; ml_seed <- 1010

ws <- function(wb, sheet, x){
  if (is.null(x)) return(invisible())
  x <- tryCatch(as_df_safe(x), error=function(e) tryCatch(as.data.frame(x), error=function(e) data.frame()))
  if (!nrow(x) && !ncol(x)) return(invisible())
  sheet <- substr(sheet,1,31); old <- names(wb); hit <- which(tolower(old)==tolower(sheet))
  if (length(hit)) openxlsx::removeWorksheet(wb, old[hit[1]])
  write_sheet2(wb, sheet, x); invisible()
}
sp <- function(p, path, name, w=10, h=8){
  if (!is.null(p)) tryCatch(save_plot2(p,path,name,w,h),
                            error=function(e) message("[plot skipped] ",name,": ",conditionMessage(e)))
}
save_book <- function() openxlsx::saveWorkbook(wb, xlsx, overwrite=TRUE)
ok <- function(x) !is.null(x) && !is_error_result(x)

save_std_ml <- function(res, tag, path, wb, file){
  if (!ok(res)) return(invisible())
  dir.create(path,TRUE,FALSE)
  tabs <- list(
    metric=if(length(res)>=1 && is.data.frame(res[[1]])) res[[1]] else NULL,
    importance=if(length(res)>=2 && is.data.frame(res[[2]])) res[[2]] else NULL,
    filter=res$filter_table, pred=res$predictions, param=res$parameters,
    stats=res$model_stats, confusion=res$confusion, performance=res$class_performance
  )
  for(nm in names(tabs)) ws(wb,paste0(tag,"_",nm),tabs[[nm]])
  plots <- list(importance=res$importance_plot,roc=res$roc_plot,
                confusion=res$confusion_plot,probability=res$probability_plot)
  for(nm in names(plots)) sp(plots[[nm]],path,paste0(tolower(tag),"_",nm),10,8)
  openxlsx::saveWorkbook(wb,file,overwrite=TRUE); invisible(res)
}

save_rf <- function(res, path, wb, file, prefix="RF"){
  if (!ok(res)) return(invisible())
  dir.create(path,TRUE,FALSE)
  pp <- list(
    variable=res$variability_plot, importance=res$importance_plot, error=res$error_plot,
    confusion=res$confusion_plot, mds=res$mds_plot, heatmap=res$abundance_heatmap,
    abundance=res$abundance_plot, combined=res$combined_plot, cv=res$rfcv_plot
  )
  sz <- list(variable=c(11,10),importance=c(12,10),error=c(11,8),confusion=c(9,8),
             mds=c(10,8),heatmap=c(12,10),abundance=c(16,12),combined=c(20,11),cv=c(9,7))
  i <- 0L
  for(nm in names(pp)){ i <- i+1L; wh <- sz[[nm]]
  sp(pp[[nm]],path,sprintf("rf_%02d_%s",i,nm),wh[1],wh[2])
  }
  tabs <- list(
    importance_top=if(length(res)>=3) res[[3]] else NULL,
    importance_all=if(length(res)>=5) res[[5]] else NULL,
    feature_filter=res$filter_table, model_stats=res$model_stats,
    confusion=res$confusion, performance=res$class_performance,
    abundance=res$abundance_summary, parameters=res$parameters, feature_CV=res$rfcv_data
  )
  for(nm in names(tabs)) ws(wb,paste0(prefix,"_",nm),tabs[[nm]])
  openxlsx::saveWorkbook(wb,file,overwrite=TRUE); invisible(res)
}

## ---------- 6.1 PCA ----------
res34 <- run_safe("loadingPCA.metm2", loadingPCA.metm2(
  ps=ps.micro, Top=20, group="Group", taxrank="Genus",
  cumulative=.80, transform="hellinger"
))
if (ok(res34)){
  P <- list(sample=res34$pca_plot,variance=res34$scree_plot,importance=res34$importance_plot,
            abundance=res34$abundance_plot,combined=res34[[1]])
  S <- list(sample=c(10,8),variance=c(10,7),importance=c(9,10),abundance=c(10,10),combined=c(18,10))
  i <- 0L; for(nm in names(P)){i<-i+1L; sp(P[[nm]],paths$PCA,sprintf("pca_%02d_%s",i,nm),S[[nm]][1],S[[nm]][2])}
  T <- list(feature_all=res34[[2]],feature_top20=res34$top_features,variance=res34$pca_variance,
            scores=res34$pca_scores,loadings=res34$loadings,abundance=res34$abundance_summary,
            group_test=res34$group_test,pairwise=res34$pairwise_test)
  for(nm in names(T)) ws(wb,paste0("PCA_",nm),T[[nm]])
  save_book()
}

## ---------- 6.2 LDA ----------
tablda <- run_safe("LDA.metm2", LDA.metm2(
  ps=ps.micro, group="Group", rank="Genus", Top=20,
  p.lvl=.5, lda.lvl=.1, seed=11, adjust.p=FALSE,
  show_n=20, group_colors=col.g
))
if (ok(tablda)){
  if (!is.null(tablda[[2]]) && nrow(tablda[[2]])>0)
    sp(tryCatch(lefse_bar(tablda[[2]])+theme_classic(),error=function(e) NULL),
       paths$LDA,"lda_01_lefse_style",10,8)
  P <- list(importance=tablda$lda_bar,significance=tablda$significance_plot,
            heatmap=tablda$abundance_heatmap,abundance=tablda$abundance_plot,
            marker_number=tablda$marker_count_plot,combined=tablda$combined_plot)
  S <- list(importance=c(11,9),significance=c(11,8),heatmap=c(12,10),
            abundance=c(16,14),marker_number=c(9,7),combined=c(20,11))
  i<-1L; for(nm in names(P)){i<-i+1L; sp(P[[nm]],paths$LDA,sprintf("lda_%02d_%s",i,nm),S[[nm]][1],S[[nm]][2])}
  T <- list(significant=tablda$taxtree,all_features=tablda$all_results,
            abundance_raw=tablda$abundance_raw,abundance_mean=tablda$group_mean,
            marker_count=tablda$marker_count,feature_filter=tablda$lda_filter,
            parameters=tablda$parameters)
  for(nm in names(T)) ws(wb,paste0("LDA_",nm),T[[nm]])
  save_book()
}

## ---------- 6.3 Global multi-group Random Forest ----------
rf_k <- min(3L, as.integer(min(table(as.character(phyloseq::sample_data(ps.micro)$Group)))))
res_rf_global <- run_safe("randomforest.metm2_global", randomforest.metm2(
  ps=ps.micro, group="Group", rank="Genus", optimal=20, ntree=1000, seed=11,
  transform="hellinger", prevalence=.30, min_mean=.01,
  feature_filter="MAD", variable_top=100,
  fill="Phylum", lab="Genus", group_colors=col.g, rfcv=TRUE, nrfcvnum=rf_k
))
save_rf(res_rf_global, paths$RF, wb, xlsx, "RF")

## ---------- 6.4 Pairwise Bray-Curtis PERMANOVA ----------
pair_adonis <- function(g){
  p <- ps.micro %>% subset_samples.wt("Group",g)
  p <- phyloseq::prune_taxa(phyloseq::taxa_sums(p)>0,p)
  m <- data.frame(phyloseq::sample_data(p),check.names=FALSE)
  m$Group <- droplevels(factor(as.character(m$Group)))
  n <- table(m$Group)
  if(length(n)!=2 || min(n)<2) return(data.frame(Group1=g[1],Group2=g[2],N1=n[1],N2=n[2],R2=NA,F=NA,P=NA,Dispersion_P=NA))
  d <- phyloseq::distance(p,"bray"); m <- m[attr(d,"Labels"),,drop=FALSE]
  a <- tryCatch(vegan::adonis2(d~Group,data=m,permutations=perm_n),error=function(e) NULL)
  dp <- tryCatch(vegan::permutest(vegan::betadisper(d,m$Group),permutations=perm_n)$tab[1,"Pr(>F)"],
                 error=function(e) NA_real_)
  if(is.null(a)) return(data.frame(Group1=g[1],Group2=g[2],N1=n[1],N2=n[2],R2=NA,F=NA,P=NA,Dispersion_P=dp))
  data.frame(Group1=g[1],Group2=g[2],N1=as.integer(n[g[1]]),N2=as.integer(n[g[2]]),
             R2=a$R2[1],F=a$F[1],P=a$`Pr(>F)`[1],Dispersion_P=dp)
}
set.seed(perm_seed)
gr <- sort(unique(as.character(phyloseq::sample_data(ps.micro)$Group)))
pairwise_adonis <- dplyr::bind_rows(lapply(combn(gr,2,simplify=FALSE),pair_adonis))
pairwise_adonis$Padj <- p.adjust(pairwise_adonis$P,"BH")
pairwise_adonis$Dispersion_warning <- !is.na(pairwise_adonis$Dispersion_P) & pairwise_adonis$Dispersion_P<.05
pairwise_adonis <- pairwise_adonis[order(pairwise_adonis$P,-pairwise_adonis$R2,na.last=TRUE),]
rownames(pairwise_adonis) <- NULL
ws(wb,"pairwise_PERMANOVA",pairwise_adonis)
write.table(pairwise_adonis,file.path(paths$PERM,"pairwise_PERMANOVA.tsv"),sep="\t",quote=FALSE,row.names=FALSE)

## ---------- 6.5 Choose pairs for ML ----------
sig <- pairwise_adonis[is.finite(pairwise_adonis$Padj) & pairwise_adonis$Padj<padj_cut,,drop=FALSE]
if(nrow(sig)){
  ml_pairs <- sig; ml_pairs$Selection_mode <- "significant_Padj"; ml_pairs$Exploratory <- FALSE
} else {
  z <- pairwise_adonis[is.finite(pairwise_adonis$P)&is.finite(pairwise_adonis$R2),,drop=FALSE]
  z <- z[order(z$P,-z$R2),,drop=FALSE]
  ml_pairs <- head(z,min(fallback_n,nrow(z)))
  ml_pairs$Selection_mode <- paste0("fallback_top",fallback_n,"_lowest_P")
  ml_pairs$Exploratory <- TRUE
}
ml_pairs$ML_selected <- TRUE
if(nrow(ml_pairs)) ml_pairs$Selection_rank <- seq_len(nrow(ml_pairs))
ws(wb,"ML_selected_pairs",ml_pairs)
write.table(ml_pairs,file.path(paths$PERM,"ML_selected_pairs.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
save_book()
message("ML pairs: ",paste(paste0(ml_pairs$Group1,"_vs_",ml_pairs$Group2),collapse=", "))

## ---------- 6.6 Pairwise ML ----------
model_cfg <- list(
  SVM=list(fun="svm_metm2",dir="01_SVM",var=10L,seed=ml_seed),
  GLM=list(fun="glm.metm2",dir="02_GLM",var=3L,seed=ml_seed),
  LASSO=list(fun="lasso.metm2",dir="03_LASSO",var=Inf,seed=ml_seed),
  XGB=list(fun="xgboost.metm2",dir="04_XGBoost",var=Inf,seed=ml_seed),
  Tree=list(fun="decisiontree.metm2",dir="05_DecisionTree",var=Inf,seed=6358),
  Bagging=list(fun="bagging_metm2",dir="07_Bagging",var=Inf,seed=ml_seed),
  NB=list(fun="naivebayes.metm2",dir="08_NaiveBayes",var=10L,seed=ml_seed,laplace=1),
  NNet=list(fun="nnet.metm2",dir="09_NNet",var=10L,seed=ml_seed)
)

for(ii in seq_len(nrow(ml_pairs))){
  g <- c(ml_pairs$Group1[ii],ml_pairs$Group2[ii]); tag <- paste(g,collapse="_vs_")
  pst <- ps.micro %>% subset_samples.wt("Group",g) %>%
    phyloseq::filter_taxa(function(x) sum(x)>10,prune=TRUE)
  phyloseq::sample_data(pst)$Group <- droplevels(factor(as.character(phyloseq::sample_data(pst)$Group)))
  gn <- table(as.character(phyloseq::sample_data(pst)$Group))
  if(length(gn)!=2 || min(gn)<2){message("[skip] ",tag); next}
  k <- max(2L,min(3L,as.integer(min(gn))))
  mt <- min(20L,max(6L,2L*phyloseq::nsamples(pst)))
  pair_path <- file.path(paths$ML,tag); dir.create(pair_path,TRUE,FALSE)
  pair_xlsx <- file.path(pair_path,"ML_results.xlsx")
  pair_wb <- if(file.exists(pair_xlsx)) openxlsx::loadWorkbook(pair_xlsx) else openxlsx::createWorkbook()

  design <- data.frame(
    Parameter=c("Group1","Group2","N_Group1","N_Group2","CV_folds","PERMANOVA_R2","PERMANOVA_F",
                "PERMANOVA_P","PERMANOVA_Padj","Dispersion_P","Dispersion_warning","Selection_mode","Exploratory"),
    Value=c(g[1],g[2],gn[g[1]],gn[g[2]],k,ml_pairs$R2[ii],ml_pairs$F[ii],ml_pairs$P[ii],
            ml_pairs$Padj[ii],ml_pairs$Dispersion_P[ii],ml_pairs$Dispersion_warning[ii],
            ml_pairs$Selection_mode[ii],ml_pairs$Exploratory[ii])
  )
  ws(pair_wb,"design",design)

  ## standard binary models
  for(nm in names(model_cfg)){
    cf <- model_cfg[[nm]]
    vtop <- if(is.finite(cf$var)) min(cf$var,mt) else mt
    args <- list(ps=pst,group="Group",rank="Genus",k=k,seed=cf$seed,top=20,
                 prevalence=ml_prev,min_mean=ml_mean,variable_top=vtop,
                 transform=ml_trans,lab="Genus",fill="Phylum",group_colors=col.g)
    if(!is.null(cf$laplace)) args$laplace <- cf$laplace
    res <- run_safe(paste0(nm,"_",tag), do.call(get(cf$fun),args))
    save_std_ml(res,nm,file.path(pair_path,cf$dir),pair_wb,pair_xlsx)
  }

  ## pairwise RF
  res_rf <- run_safe(paste0("RF_",tag), randomforest.metm2(
    ps=pst,group="Group",rank="Genus",optimal=20,ntree=1000,seed=11,
    transform=ml_trans,prevalence=ml_prev,min_mean=ml_mean,
    feature_filter="MAD",variable_top=mt,fill="Phylum",lab="Genus",
    group_colors=col.g,rfcv=TRUE,nrfcvnum=k
  ))
  save_rf(res_rf,file.path(pair_path,"06_RandomForest"),pair_wb,pair_xlsx,"RF")

  ## single-feature ROC
  res_roc <- run_safe(paste0("ROC_",tag), Roc.metm2(
    ps=pst,group="Group",rank="Genus",top=10,prevalence=ml_prev,min_mean=ml_mean,
    boot_n=500,seed=11,group_colors=col.g
  ))
  if(ok(res_roc)){
    rp <- file.path(pair_path,"10_ROC"); dir.create(rp,TRUE,FALSE)
    sp(res_roc[[1]],rp,"roc_01_top_features",10,8)
    sp(res_roc$auc_plot,rp,"roc_02_auc_ranking",10,8)
    sp(res_roc$abundance_plot,rp,"roc_03_feature_abundance",12,10)
    ws(pair_wb,"ROC_all",res_roc[[2]]); ws(pair_wb,"ROC_top",res_roc$top_results)
    ws(pair_wb,"ROC_parameters",res_roc$parameters)
  }

  openxlsx::saveWorkbook(pair_wb,pair_xlsx,overwrite=TRUE)
  message("[DONE] ",tag)
}

save_book()
message("Biomarker finished: ",normalizePath(bp,mustWork=FALSE))


## ===================== 7. 网络分析（Network Analysis） =====================



## ================= 7. Network analysis: complete compact workflow =================
## ps.micro + mg_path are required.
## Core: network -> properties -> ZiPi/random -> modules -> stability -> pairwise turnover.

library(ggClusterNet); library(igraph); library(dplyr); library(tibble); library(ggplot2)
if ("package:mia" %in% search()) detach("package:mia", unload=TRUE, character.only=TRUE)

## ---------- 7.0 settings ----------
np <- file.path(mg_path,"06_network")
dirs <- list(
  cache=file.path(np,"00_cache"), main=file.path(np,"01_main"),
  prop=file.path(np,"02_properties"), topo=file.path(np,"03_topology_roles"),
  module=file.path(np,"04_module_similarity"), stability=file.path(np,"05_stability"),
  pair=file.path(np,"06_pairwise")
)
invisible(lapply(c(np,unlist(dirs)),dir.create,recursive=TRUE,showWarnings=FALSE))
xlsx <- file.path(np,"network_results.xlsx")
wb <- openxlsx::createWorkbook()

net_N <- 200; r_thr <- .8; p_thr <- .05; net_seed <- 11
network_pairs <- NULL                # NULL=all pairwise; or list(c("CK","LA"),c("LA","MA"))
pair_col <- NULL                     # paired design only, e.g. "Block"; otherwise NULL
run_stability <- TRUE; run_module <- TRUE

ws <- function(wb,sheet,x,row_names=FALSE){
  if(is.null(x)) return(invisible())
  x <- tryCatch(as.data.frame(x),error=function(e)data.frame())
  if(!nrow(x) && !ncol(x)) return(invisible())
  sheet <- substr(sheet,1,31); hit <- which(tolower(names(wb))==tolower(sheet))
  if(length(hit)) openxlsx::removeWorksheet(wb,names(wb)[hit[1]])
  write_sheet2(wb,sheet,x,row_names=row_names); invisible()
}
sp <- function(p,path,name,w=10,h=8){
  if(!is.null(p)) tryCatch(save_plot2(p,path,name,w,h),
                           error=function(e) message("[plot skipped] ",name,": ",conditionMessage(e)))
}
safe <- function(tag,x) tryCatch(x,error=function(e){message("[",tag," ERROR] ",conditionMessage(e)); NULL})
savewb <- function(wb,file) openxlsx::saveWorkbook(wb,file,overwrite=TRUE)
wtsv <- function(x,file) if(!is.null(x) && NROW(x)) write.table(as.data.frame(x),file,sep="\t",quote=FALSE,row.names=FALSE)
pick <- function(x,keys){for(k in keys) if(!is.null(x[[k]])) return(as.data.frame(x[[k]])); data.frame()}
edge_nodes <- function(e) if(nrow(e)) unique(c(as.character(e$var1),as.character(e$var2))) else character()
tax_sum <- function(ids,tax,rank){
  if(!length(ids)||!rank%in%names(tax)) return(data.frame())
  x <- as.character(tax[match(ids,rownames(tax)),rank])
  x <- x[!is.na(x) & nzchar(x)]
  if(!length(x)) return(data.frame())
  tt <- sort(table(x),decreasing=TRUE)
  y <- data.frame(Taxon=names(tt),Nodes=as.integer(tt),stringsAsFactors=FALSE)
  names(y)[1] <- rank
  y$Fraction <- y$Nodes/sum(y$Nodes)
  y
}
top_tax_plot <- function(e,tax,rank,title){
  z<-head(tax_sum(edge_nodes(e),tax,rank),20); if(!nrow(z)) return(NULL)
  ggplot(z,aes(x=reorder(.data[[rank]],Nodes),y=Nodes))+geom_col()+coord_flip()+theme_classic()+
    labs(x=rank,y="Network nodes",title=title)
}
edge_net_plot <- function(e,tax,title){
  if(!nrow(e)||!requireNamespace("ggraph",quietly=TRUE)) return(NULL)
  ed<-e[,c("var1","var2"),drop=FALSE]; names(ed)<-c("from","to")
  g<-igraph::graph_from_data_frame(ed,FALSE); td<-tax[match(V(g)$name,rownames(tax)),,drop=FALSE]
  V(g)$Phylum<-if("Phylum"%in%names(td)) as.character(td$Phylum) else "Unknown"; V(g)$Degree<-degree(g)
  ggraph::ggraph(g,layout="fr")+ggraph::geom_edge_link(alpha=.2)+
    ggraph::geom_node_point(aes(size=Degree,fill=Phylum),shape=21)+theme_void()+labs(title=title)
}
venn_plot <- function(A,B,AB,g1,g2){
  if(!requireNamespace("VennDiagram",quietly=TRUE)||!requireNamespace("ggplotify",quietly=TRUE)) return(NULL)
  v<-VennDiagram::draw.pairwise.venn(A,B,AB,category=c(g1,g2),alpha=.55,lwd=1.2,cex=1.1,cat.cex=1.1)
  ggplotify::as.ggplot(grid::grobTree(v))
}
save_res2 <- function(res,path,prefix,w=8,h=5){
  if(is.null(res)) return(invisible())
  if(length(res)>=1) sp(res[[1]],path,prefix,w,h)
  if(length(res)>=2){d<-as.data.frame(res[[2]]); wtsv(d,file.path(path,paste0(prefix,".tsv"))); ws(wb,prefix,d)}
}

## ---------- 7.1 main network ----------
gn <- table(as.character(phyloseq::sample_data(ps.micro)$Group))
if(any(gn<6)) warning("Some groups have <6 samples: correlation networks/stability comparisons are exploratory.")
set.seed(net_seed)
tab_net <- safe("network.pip",network.pip(
  ps=ps.micro,N=net_N,big=TRUE,select_layout=FALSE,layout_net="model_maptree2",
  r.threshold=r_thr,p.threshold=p_thr,maxnode=2,method="spearman",label=FALSE,
  lab="elements",group="Group",fill="Phylum",size="igraph.degree",
  zipi=TRUE,ram.net=TRUE,clu_method="cluster_fast_greedy",step=100,R=10,ncpus=1
))
if(is.null(tab_net)) stop("network.pip failed.")
net_plots<-tab_net[[1]]; dat<-tab_net[[2]]; netmat<-dat$net.cor.matrix; cortab<-netmat$cortab

## Normalize all correlation/network matrices once, before any downstream analysis.
cortab <- lapply(cortab,function(x){
  x <- as.matrix(x)
  storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- 0
  if(nrow(x)==ncol(x)) diag(x) <- 0
  x
})

saveRDS(tab_net,file.path(dirs$cache,"network_pip.rds"))
saveRDS(cortab,file.path(dirs$cache,"cor_matrix_all_groups.rds"))

## Save network.pip plots with explicit names.
if(length(net_plots)>=1) sp(net_plots[[1]],dirs$main,"network_main",14,8)
if(length(net_plots)>=2) sp(net_plots[[2]],dirs$topo,"ZiPi_network_pip",12,8)
if(length(net_plots)>=3) sp(net_plots[[3]],dirs$topo,"random_network_network_pip",12,6)
if(length(net_plots)>3)
  for(i in 4:length(net_plots)) sp(net_plots[[i]],dirs$main,sprintf("network_extra_%02d",i),12,8)
for(nm in names(cortab)){ws(wb,paste0("cor_",nm),cortab[[nm]],TRUE); write.table(cortab[[nm]],file.path(dirs$cache,paste0("cor_",nm,".tsv")),sep="\t",quote=FALSE)}
if(!is.null(netmat$node)){ws(wb,"network_nodes",netmat$node); wtsv(netmat$node,file.path(dirs$main,"network_nodes.tsv"))}
if(!is.null(netmat$edge)){ws(wb,"network_edges",netmat$edge); wtsv(netmat$edge,file.path(dirs$main,"network_edges.tsv"))}
ws(wb,"network_parameters",data.frame(Parameter=c("N","method","r","p","layout","cluster"),
                                      Value=c(net_N,"spearman",r_thr,p_thr,"model_maptree2","cluster_fast_greedy")))
tax<-as.data.frame(phyloseq::tax_table(ps.micro)); tax$ID<-rownames(tax)
gm<-safe("group.mean",ps.micro%>%ps.group.mean()%>%left_join(tax,by="ID")); ws(wb,"group_mean_tax",gm)
savewb(wb,xlsx)

## ---------- 7.2 network properties + statistics ----------
grp<-names(cortab)
prop<-bind_rows(lapply(grp,function(g){
  z<-as.data.frame(net_properties.4(make_igraph(cortab[[g]]),n.hub=FALSE))
  data.frame(Group=g,Metric=rownames(z),Value=suppressWarnings(as.numeric(z[[1]])))
}))
ws(wb,"network_properties_summary",prop); wtsv(prop,file.path(dirs$prop,"network_properties_summary.tsv"))
prop_plot <- prop[is.finite(prop$Value),,drop=FALSE]
if(nrow(prop_plot))
  sp(ggplot(prop_plot,aes(Group,Value))+geom_col()+facet_wrap(~Metric,scales="free_y")+theme_classic()+
       theme(axis.text.x=element_text(angle=45,hjust=1)),dirs$prop,"network_properties",14,10)

sample_prop<-bind_rows(lapply(grp,function(g){
  p<-ps.micro%>%subset_samples.wt("Group",g)%>%remove.zero()
  z<-as.data.frame(netproperties.sample(pst=p,cor=cortab[[g]])); z$ID<-rownames(z); rownames(z)<-NULL; z$NetworkGroup<-g; z
}))
meta<-as.data.frame(phyloseq::sample_data(ps.micro)); if("ID"%in%names(meta)) names(meta)[names(meta)=="ID"]<-"ID0"
meta<-rownames_to_column(meta,"ID"); sample_meta<-left_join(sample_prop,meta,by="ID")
ws(wb,"sample_network_properties",sample_prop); ws(wb,"sample_netprop_with_meta",sample_meta)
wtsv(sample_prop,file.path(dirs$prop,"sample_network_properties.tsv"))
wtsv(sample_meta,file.path(dirs$prop,"sample_network_properties_with_metadata.tsv"))

num_prop<-names(sample_prop)[vapply(sample_prop,is.numeric,logical(1))]
kw<-bind_rows(lapply(num_prop,function(v){
  z<-sample_prop[is.finite(sample_prop[[v]]),,drop=FALSE]
  if(length(unique(z$NetworkGroup))<2) return(NULL)
  k<-safe(paste0("KW_",v),kruskal.test(z[[v]]~factor(z$NetworkGroup)))
  if(is.null(k)) NULL else data.frame(Metric=v,Statistic=unname(k$statistic),P=k$p.value)
}))
if(nrow(kw)){kw$Padj<-p.adjust(kw$P,"BH"); ws(wb,"sample_property_KW",kw); wtsv(kw,file.path(dirs$prop,"sample_property_KW.tsv"))}
if(length(num_prop)){
  lp<-tidyr::pivot_longer(sample_prop,all_of(num_prop),names_to="Metric",values_to="Value")
  sp(ggplot(lp,aes(NetworkGroup,Value))+geom_boxplot(outlier.shape=NA)+geom_jitter(width=.12)+
       facet_wrap(~Metric,scales="free_y")+theme_classic()+theme(axis.text.x=element_text(angle=45,hjust=1)),
     dirs$prop,"sample_network_properties",14,10)
}

node_prop<-bind_rows(lapply(grp,function(g){
  z<-as.data.frame(node_properties(make_igraph(cortab[[g]]))); z$ASV.name<-rownames(z); rownames(z)<-NULL; z$Group<-g; z
}))
node_prop<-left_join(node_prop,tax,by=c("ASV.name"="ID"))
ws(wb,"node_properties_combined",node_prop); wtsv(node_prop,file.path(dirs$prop,"node_properties_combined.tsv"))
dc<-grep("degree",names(node_prop),ignore.case=TRUE,value=TRUE)[1]
if(!is.na(dc)){hubs<-node_prop%>%group_by(Group)%>%arrange(desc(.data[[dc]]),.by_group=TRUE)%>%slice_head(n=20)%>%ungroup()
ws(wb,"top_degree_nodes",hubs); wtsv(hubs,file.path(dirs$prop,"top_degree_nodes.tsv"))}

sign_sum<-bind_rows(lapply(grp,function(g){
  m<-as.matrix(cortab[[g]]); v<-m[upper.tri(m)]; v<-v[is.finite(v)&v!=0]
  data.frame(Group=g,Edges=length(v),Positive=sum(v>0),Negative=sum(v<0),
             Positive_ratio=if(length(v))mean(v>0) else NA,Negative_ratio=if(length(v))mean(v<0) else NA)
}))
ws(wb,"edge_sign_summary",sign_sum); wtsv(sign_sum,file.path(dirs$prop,"edge_sign_summary.tsv"))
savewb(wb,xlsx)

## ---------- 7.3 ZiPi roles + random-network diagnostic ----------

zipi <- dat$zipi.data

if(
  !is.null(zipi) &&
  all(c("p", "z") %in% names(zipi))
){

  zipi <- as.data.frame(
    zipi,
    stringsAsFactors = FALSE
  )

  ## network.pip 有时会产生 list-column
  flat_chr <- function(x){

    if(is.list(x))
      return(
        vapply(
          x,
          function(y)
            paste(
              as.character(y),
              collapse = ";"
            ),
          character(1)
        )
      )

    as.character(x)
  }


  for(
    v in intersect(
      c(
        "group",
        "module",
        "label",
        "taxa",
        "roles",
        "role_7"
      ),
      names(zipi)
    )
  ){
    zipi[[v]] <- flat_chr(
      zipi[[v]]
    )
  }


  zipi$p <- suppressWarnings(
    as.numeric(zipi$p)
  )

  zipi$z <- suppressWarnings(
    as.numeric(zipi$z)
  )


  ## 四类 ZiPi role
  zipi$Role <- ifelse(
    zipi$z >= 2.5 &
      zipi$p >= 0.62,
    "Network hub",

    ifelse(
      zipi$z >= 2.5,
      "Module hub",

      ifelse(
        zipi$p >= 0.62,
        "Connector",
        "Peripheral"
      )
    )
  )


  ## 保存全部 ZiPi 数据
  ws(
    wb,
    "ZiPi_roles",
    zipi
  )

  wtsv(
    zipi,
    file.path(
      dirs$topo,
      "ZiPi_roles.tsv"
    )
  )


  ## 每组四类角色数量
  if(
    "group" %in%
    names(zipi)
  ){

    zr <- zipi %>%
      dplyr::group_by(
        .data$group,
        .data$Role
      ) %>%
      dplyr::summarise(
        N = dplyr::n(),
        .groups = "drop"
      )


    ws(
      wb,
      "ZiPi_role_counts",
      zr
    )

    wtsv(
      zr,
      file.path(
        dirs$topo,
        "ZiPi_role_counts.tsv"
      )
    )


    ## 重新绘制统一 ZiPi 图
    pz <- ggplot(
      zipi,
      aes(
        p,
        z
      )
    ) +

      geom_vline(
        xintercept = 0.62,
        linetype = 2
      ) +

      geom_hline(
        yintercept = 2.5,
        linetype = 2
      ) +

      geom_point(
        aes(
          shape = Role
        ),
        size = 2
      ) +

      facet_wrap(
        ~group
      ) +

      theme_classic() +

      labs(
        x =
          "Participation coefficient (Pi)",
        y =
          "Within-module connectivity z-score (Zi)"
      )


  } else {

    zr <- zipi %>%
      dplyr::count(
        .data$Role,
        name = "N"
      )

    ws(
      wb,
      "ZiPi_role_counts",
      zr
    )


    pz <- ggplot(
      zipi,
      aes(
        p,
        z
      )
    ) +

      geom_vline(
        xintercept = 0.62,
        linetype = 2
      ) +

      geom_hline(
        yintercept = 2.5,
        linetype = 2
      ) +

      geom_point(
        aes(
          shape = Role
        ),
        size = 2
      ) +

      theme_classic()

  }


  sp(
    pz,
    dirs$topo,
    "ZiPi_roles_replot",
    12,
    8
  )
}

## Save random-network diagnostic data returned by network.pip.
rnd <- dat$random.net.data
if(!is.null(rnd)){
  ws(wb,"random_network_data",rnd)
  wtsv(rnd,file.path(dirs$topo,"random_network_data.tsv"))
}
savewb(wb,xlsx)

## ---------- 7.4 network/module comparison ----------
cmpnet<-safe("module.compare.net.pip",module.compare.net.pip(
  ps=NULL,corg=cortab,degree=TRUE,zipi=FALSE,r.threshold=r_thr,p.threshold=p_thr,
  method="spearman",padj=FALSE,n=3))
if(!is.null(cmpnet)){ws(wb,"network_comparison_results",cmpnet[[1]])
  wtsv(cmpnet[[1]],file.path(dirs$module,"network_comparison_results.tsv"))}

if(run_module){
  ## Avoid current module.compare.m() namespace conflict (unqualified select()).
  ## Build module memberships group-by-group, then compare them with model_compare().
  mod_member <- safe("module_membership", dplyr::bind_rows(lapply(grp,function(g){
    rr <- model_maptree2(cor=as.matrix(cortab[[g]]),method="cluster_fast_greedy")

    ## Find the returned table containing ID + group, independent of list position.
    z <- NULL
    for(obj in rr){
      if(is.data.frame(obj) || is.matrix(obj)){
        zz <- as.data.frame(obj,stringsAsFactors=FALSE)
        if(all(c("ID","group") %in% names(zz))){ z <- zz; break }
      }
    }
    if(is.null(z)) stop("model_maptree2() did not return a table containing ID and group.")

    z <- z[z$group!="mother_no",c("ID","group"),drop=FALSE]
    z$group <- paste0(g,"_",z$group)
    z$Group <- g
    z
  })))

  if(!is.null(mod_member) && nrow(mod_member)){
    ws(wb,"module_membership",mod_member)
    wtsv(mod_member,file.path(dirs$module,"module_membership.tsv"))

    mod_sim <- safe("model_compare",model_compare(node_table2=mod_member,n=3,padj=FALSE))
    if(!is.null(mod_sim)){
      mod_sim <- as.data.frame(mod_sim,stringsAsFactors=FALSE)
      ws(wb,"module_similarity_table",mod_sim)
      wtsv(mod_sim,file.path(dirs$module,"module_similarity_table.tsv"))

      if(all(c("module1","module2") %in% names(mod_sim))){
        module_to_group <- setNames(mod_member$Group,mod_member$group)
        g1m <- unname(module_to_group[as.character(mod_sim$module1)])
        g2m <- unname(module_to_group[as.character(mod_sim$module2)])
        mod_sim$cross <- paste(g1m,g2m,sep="_vs_")

        mc <- mod_sim |>
          dplyr::filter(.data$module1!="none",!is.na(.data$cross)) |>
          dplyr::count(.data$cross,name="Similar_modules")

        ws(wb,"module_similarity_counts",mc)
        wtsv(mc,file.path(dirs$module,"module_similarity_counts.tsv"))
        if(nrow(mc))
          sp(ggplot(mc,aes(reorder(cross,Similar_modules),Similar_modules))+
               geom_col()+coord_flip()+theme_classic()+
               labs(x=NULL,y="Number of similar modules"),
             dirs$module,"module_similarity_counts",8,6)
      }
    }
  }
}
savewb(wb,xlsx)

## ---------- 7.5 stability / robustness ----------
if(run_stability){
  save_res2(safe("Targeted removal",
                 Robustness.Targeted.removal(ps=ps.micro,corg=cortab,degree=TRUE,zipi=FALSE)),
            dirs$stability,"targeted_removal",8,5)

  save_res2(safe("Random removal",
                 Robustness.Random.removal(ps=ps.micro,corg=cortab,Top=0)),
            dirs$stability,"random_removal",8,5)

  save_res2(safe("Negative ratio",
                 negative.correlation.ratio(ps=ps.micro,corg=cortab,degree=TRUE,zipi=FALSE)),
            dirs$stability,"negative_correlation_ratio",6,5)

  ## Natural connectivity: call pulsar explicitly instead of the broken
  ## un-namespaced natural.connectivity() call inside natural.con.microp().
  if(requireNamespace("pulsar",quietly=TRUE)){
    end_n <- min(100L,max(0L,min(vapply(cortab,nrow,integer(1)))-2L))
    if(end_n>=2){
      set.seed(net_seed)
      nc <- safe("Natural connectivity",dplyr::bind_rows(lapply(names(cortab),function(g){
        A <- as.matrix(cortab[[g]])
        ids <- rownames(A)
        dplyr::bind_rows(lapply(0:end_n,function(k){
          vv <- replicate(10,{
            keep <- if(k==0) ids else sample(ids,length(ids)-k)
            pulsar::natural.connectivity(G=A[keep,keep,drop=FALSE],norm=TRUE)
          })
          data.frame(Group=g,Removed_nodes=k,
                     Natural_connectivity=mean(vv,na.rm=TRUE),
                     SD=stats::sd(vv,na.rm=TRUE),N=10)
        }))
      })))
      if(!is.null(nc) && nrow(nc)){
        ws(wb,"natural_connectivity",nc)
        wtsv(nc,file.path(dirs$stability,"natural_connectivity.tsv"))
        sp(ggplot(nc,aes(Removed_nodes,Natural_connectivity,color=Group))+
             geom_line()+theme_classic(),
           dirs$stability,"natural_connectivity",7,5)
      }
    }
  } else message("[Natural connectivity skipped] install.packages('pulsar')")

  if(!is.null(pair_col) && pair_col%in%names(meta)){
    pp <- ps.micro
    md <- as.data.frame(phyloseq::sample_data(pp))
    md$pair <- meta[match(rownames(md),meta$ID),pair_col]
    phyloseq::sample_data(pp) <- phyloseq::sample_data(md)
    save_res2(safe("Community stability",
                   community.stability(ps=pp,corg=cortab,time=FALSE)),
              dirs$stability,"community_stability",7,5)
  }
}
savewb(wb,xlsx)

## ---------- 7.6 pairwise shared/specific network turnover ----------
## ===== Fix: group network comparison (robust, independent of align_corr_pair) =====
## Requires existing objects: cortab, tax, dirs, wb, xlsx, ws, wtsv, sp, savewb
## Optional existing plot helpers: venn_plot, edge_net_plot, top_tax_plot

library(dplyr); library(ggplot2)

## 1. normalize network matrices
cortab <- lapply(cortab,function(x){
  x <- as.matrix(x); storage.mode(x) <- "numeric"
  x[!is.finite(x)] <- 0
  if(nrow(x)==ncol(x)) diag(x) <- 0
  x
})
grp <- names(cortab)

## fresh comparison workbook: avoids inheriting broken drawing relationships from older files
cmp_xlsx <- file.path(dirs$pair,"group_network_comparison.xlsx")
cmp_wb <- openxlsx::createWorkbook()

## 2. robust helpers
edge_tab <- function(M,thr=.5){
  M <- as.matrix(M); ids <- rownames(M)
  ij <- which(upper.tri(M) & abs(M)>=thr,arr.ind=TRUE)
  if(!nrow(ij)) return(data.frame(var1=character(),var2=character(),cor=numeric(),key=character()))
  z <- data.frame(var1=ids[ij[,1]],var2=ids[ij[,2]],cor=M[ij],stringsAsFactors=FALSE)
  z$key <- paste(pmin(z$var1,z$var2),pmax(z$var1,z$var2),sep="|||")
  z
}
tax_sum2 <- function(ids,rank){
  if(!length(ids)||!rank%in%names(tax)) return(data.frame())
  x <- as.character(tax[match(ids,rownames(tax)),rank]); x <- x[!is.na(x)&nzchar(x)]
  if(!length(x)) return(data.frame())
  tt <- sort(table(x),decreasing=TRUE)
  y <- data.frame(Taxon=names(tt),Nodes=as.integer(tt),Fraction=as.integer(tt)/sum(tt))
  names(y)[1] <- rank; y
}
cmp_pair <- function(A,B,thr=.5){
  ea<-edge_tab(A,thr); eb<-edge_tab(B,thr)
  sh<-merge(ea,eb,by="key",suffixes=c("_A","_B"))
  ao<-ea[!ea$key%in%sh$key,,drop=FALSE]; bo<-eb[!eb$key%in%sh$key,,drop=FALSE]
  sh2<-data.frame(var1=sh$var1_A,var2=sh$var2_A,cor_A=sh$cor_A,cor_B=sh$cor_B,
                  same_sign=sign(sh$cor_A)==sign(sh$cor_B),key=sh$key)
  na<-unique(c(ea$var1,ea$var2)); nb<-unique(c(eb$var1,eb$var2))
  list(A=ea,B=eb,shared=sh2,Aonly=ao,Bonly=bo,nA=na,nB=nb)
}

## current cortab is a signed adjacency matrix (-1/0/1), so any non-zero entry is an edge
vals <- unique(unlist(lapply(cortab,function(x) unique(as.numeric(x)))))
edge_thr <- if(all(vals %in% c(-1,0,1))) .5 else r_thr

## 3. all pairwise comparisons
pairs <- if(exists("network_pairs") && !is.null(network_pairs)) network_pairs else combn(grp,2,simplify=FALSE)
all_sum <- list()

for(i in seq_along(pairs)){
  g<-pairs[[i]]; g1<-g[1]; g2<-g[2]; tag<-paste(g,collapse="_vs_")
  pp<-file.path(dirs$pair,tag); dir.create(pp,recursive=TRUE,showWarnings=FALSE)
  px<-file.path(pp,"network_pair_results.xlsx")
  pw<-openxlsx::createWorkbook()

  x<-cmp_pair(cortab[[g1]],cortab[[g2]],edge_thr)
  sh<-x$shared; ao<-x$Aonly; bo<-x$Bonly
  common_nodes <- intersect(x$nA,x$nB)
  shared_edge_nodes <- unique(c(as.character(sh$var1),as.character(sh$var2)))
  na<-setdiff(x$nA,x$nB); nb<-setdiff(x$nB,x$nA)
  union_e<-nrow(sh)+nrow(ao)+nrow(bo)

  sm<-data.frame(
    Group1=g1,Group2=g2,
    A_edges=nrow(x$A),B_edges=nrow(x$B),Shared_edges=nrow(sh),
    Shared_same_sign=sum(sh$same_sign),Shared_opposite_sign=sum(!sh$same_sign),
    A_only_edges=nrow(ao),B_only_edges=nrow(bo),
    Edge_Jaccard=if(union_e)nrow(sh)/union_e else NA_real_,
    A_nodes=length(x$nA),B_nodes=length(x$nB),
    Common_nodes=length(common_nodes),Shared_edge_nodes=length(shared_edge_nodes),
    A_only_nodes=length(na),B_only_nodes=length(nb),
    Node_Jaccard=if(length(union(x$nA,x$nB)))length(common_nodes)/length(union(x$nA,x$nB)) else NA_real_
  )
  all_sum[[i]]<-sm

  ws(pw,"summary",sm); ws(pw,"shared_edges",sh); ws(pw,"A_only_edges",ao); ws(pw,"B_only_edges",bo)
  ws(pw,"common_nodes",data.frame(ID=common_nodes))
  ws(pw,"shared_edge_nodes",data.frame(ID=shared_edge_nodes))
  ws(pw,"A_only_nodes",data.frame(ID=na)); ws(pw,"B_only_nodes",data.frame(ID=nb))
  for(rank in intersect(c("Phylum","Class","Genus"),names(tax))){
    ws(pw,paste0("shared_",rank),tax_sum2(shared_edge_nodes,rank))
    ws(pw,paste0("Aonly_",rank),tax_sum2(na,rank))
    ws(pw,paste0("Bonly_",rank),tax_sum2(nb,rank))
  }

  ## plots; failures do not stop the pair
  if(exists("venn_plot")) sp(venn_plot(nrow(x$A),nrow(x$B),nrow(sh),g1,g2),pp,"01_edge_venn",6,6)
  if(exists("edge_net_plot")){
    E<-list(shared=data.frame(var1=sh$var1,var2=sh$var2),
            A_only=ao[,c("var1","var2"),drop=FALSE],B_only=bo[,c("var1","var2"),drop=FALSE])
    for(nm in names(E)){
      sp(edge_net_plot(E[[nm]],tax,paste(tag,nm)),pp,paste0("02_network_",nm),12,8)
      if(exists("top_tax_plot")){
        sp(top_tax_plot(E[[nm]],tax,"Phylum",paste(tag,nm,"Phylum")),pp,paste0("03_",nm,"_Phylum"),7,6)
        sp(top_tax_plot(E[[nm]],tax,"Genus",paste(tag,nm,"Genus")),pp,paste0("04_",nm,"_Genus"),7,7)
      }
    }
  }
  openxlsx::saveWorkbook(pw,px,overwrite=TRUE)
  message("[DONE] ",tag)
}

## 4. global pairwise summary + matrices
pair_sum <- bind_rows(all_sum)
ws(cmp_wb,"pairwise_network_summary",pair_sum)
wtsv(pair_sum,file.path(dirs$pair,"pairwise_network_summary.tsv"))

mkmat <- function(value,diag_value=1){
  M<-matrix(NA_real_,length(grp),length(grp),dimnames=list(grp,grp)); diag(M)<-diag_value
  for(i in seq_len(nrow(pair_sum))){
    a<-pair_sum$Group1[i]; b<-pair_sum$Group2[i]; M[a,b]<-M[b,a]<-pair_sum[[value]][i]
  }
  M
}
for(v in c("Edge_Jaccard","Node_Jaccard","Shared_edges","Shared_same_sign")){
  mm<-mkmat(v,if(grepl("Jaccard",v))1 else NA_real_)
  ws(cmp_wb,paste0("pair_",v),as.data.frame(mm),TRUE)
}
openxlsx::saveWorkbook(cmp_wb,cmp_xlsx,overwrite=TRUE)

## 5. overview plot
sp(ggplot(pair_sum,aes(reorder(paste(Group1,Group2,sep="_vs_"),Edge_Jaccard),Edge_Jaccard))+
     geom_col()+coord_flip()+theme_classic()+labs(x=NULL,y="Edge Jaccard"),
   dirs$pair,"pairwise_edge_jaccard",9,8)
sp(ggplot(pair_sum,aes(reorder(paste(Group1,Group2,sep="_vs_"),Node_Jaccard),Node_Jaccard))+
     geom_col()+coord_flip()+theme_classic()+labs(x=NULL,y="Node Jaccard"),
   dirs$pair,"pairwise_node_jaccard",9,8)

message("Pairwise network comparison finished: ",normalizePath(dirs$pair,mustWork=FALSE))



#  查看全部结果tree#---------

# 显示完整路径
library(fs)
dir_tree(mg_path, recurse = 2)



