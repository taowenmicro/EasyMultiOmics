#' Plot Phyloseq Metabolite Feature Metrics
#'
#' 可视化代谢物的多指标特征，包括 Weighted score、丰度、显著性、网络特征等
#'
#' @param ps 一个 phyloseq 对象
#' @param dat_ms 代谢物特征数据框，必须包含列：
#'   "Rank","Metabolite","Weighted_score","Mean_intensity",
#'   "Prevalence","padj","log2FoldChange","Importance","Degree","Betweenness"
#' @param top_n 显示的前 N 个特征（默认 50）
#' @return patchwork 组合图
#' @export
plot_physeq_metabolites <- function(ps, dat_ms, top_n = 50) {
  if (!inherits(ps, "phyloseq")) stop("ps 必须是 phyloseq 对象")

  # ---- 0. 必须的列 ----
  required_cols <- c("Rank","Metabolite","Weighted_score","Mean_intensity",
                     "Prevalence","padj","log2FoldChange",
                     "Importance","Degree","Betweenness")
  missing_cols <- setdiff(required_cols, colnames(dat_ms))
  if (length(missing_cols) > 0) stop("dat_ms 缺少列: ", paste(missing_cols, collapse = ", "))

  # ---- 1. 清理代谢物名字 ----
  clean_names <- function(x) {
    x <- gsub("`", "", x)          # 去掉反引号
    x <- trimws(x)                 # 去掉前后空格
    x <- gsub("\\s+", " ", x)      # 连续空格压缩
    return(x)
  }
  dat_ms$Metabolite <- clean_names(dat_ms$Metabolite)

  # ---- 2. 计算分组均值 ----
  tem <- ps %>%
    ps.group.mean() %>%
    dplyr::arrange(dplyr::desc(mean))
  tem$ID <- clean_names(tem$ID)

  # ---- 3. 找每个代谢物的最大组 ----
  sample_cols <- setdiff(colnames(tem), c("mean","ID","group.class"))
  tem$group.class <- apply(tem[, sample_cols], 1, function(row) {
    if (all(is.na(row))) {
      return("Unknown")
    } else {
      return(names(row)[which.max(as.numeric(row))])
    }
  })

  # ---- 4. 合并 ----

  head(tem)
  head(dat_ms)
  dat2 <- dplyr::left_join(dat_ms, tem[,c("ID","group.class")],
                           by = c("Metabolite" = "ID"))
  dat2$group.class[is.na(dat2$group.class)] <- "Unknown"

  # ---- 5. 选 Top N ----
  df_top <- dat2 %>%
    dplyr::arrange(Rank) %>%
    dplyr::slice_head(n = top_n)

  # ---- 6. 固定 Group 顺序 + 排序 ----
  if ("CT" %in% df_top$group.class & "NT" %in% df_top$group.class) {
    df_top$group.class <- factor(df_top$group.class, levels = c("NT","CT"))
  }

  df_ord <- df_top %>%
    dplyr::mutate(sign_key = ifelse(log2FoldChange >= 0, 1L, 0L)) %>%
    dplyr::arrange(group.class, dplyr::desc(sign_key), dplyr::desc(Weighted_score))

  otu_levels <- df_ord$Metabolite
  df_ord <- dplyr::mutate(df_ord, Metabolite = factor(Metabolite, levels = rev(otu_levels)))

  # ---- 7. 转长表 ----
  df_long <- df_ord %>%
    dplyr::transmute(
      Metabolite,
      log2FoldChange,
      Weighted       = Weighted_score,
      Mean_intensity = Mean_intensity,
      Prevalence     = Prevalence,
      Significance   = -log10(padj),
      log2FC         = log2FoldChange,
      Importance     = Importance,
      Degree         = Degree,
      Betweenness    = Betweenness
    ) %>%
    tidyr::pivot_longer(cols = -c(Metabolite, log2FoldChange),
                        names_to = "metric", values_to = "score")

  metric_levels <- c("Weighted","Mean_intensity","Prevalence",
                     "Significance","log2FC","Importance","Degree","Betweenness")
  df_long$metric <- factor(df_long$metric, levels = metric_levels)

  # ---- 8. 格式化标签 ----
  fmt_sci <- scales::label_scientific(digits = 2)
  df_long <- dplyr::mutate(df_long,
                           label = dplyr::case_when(
                             metric == "Weighted"        ~ scales::number(score, accuracy = 0.001),
                             metric == "Mean_intensity"  ~ fmt_sci(score),
                             metric == "Prevalence"      ~ scales::percent(score, accuracy = 1),
                             metric == "Significance"    ~ scales::number(score, accuracy = 0.1),
                             metric == "log2FC"          ~ scales::number(score, accuracy = 0.01),
                             metric == "Importance"      ~ scales::number(score, accuracy = 0.1),
                             metric == "Degree"          ~ scales::number(score, accuracy = 1),
                             metric == "Betweenness"     ~ scales::number(score, accuracy = 0.1),
                             TRUE                        ~ as.character(round(score, 3))
                           ))

  fc_max <- max(abs(df_ord$log2FoldChange), na.rm = TRUE)

  # ---- 9. 左图 Weighted ----
  p_left <- ggplot2::ggplot(dplyr::filter(df_long, metric == "Weighted"),
                            ggplot2::aes(x = score, y = Metabolite, fill = log2FoldChange)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.1, size = 2.7) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-fc_max, fc_max),
                                  name = "log2FC\n(blue negative, red positive)") +
    ggplot2::labs(x = "Weighted score", y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank(),
                   legend.position    = "right")

  # ---- 10. 右图多指标 ----
  p_right <- ggplot2::ggplot(dplyr::filter(df_long, metric != "Weighted"),
                             ggplot2::aes(x = score, y = Metabolite, fill = log2FoldChange)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.1, size = 2.5) +
    ggplot2::facet_grid(cols = ggplot2::vars(metric), scales = "free_x", switch = "x") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-fc_max, fc_max),
                                  name = "log2FC\n(blue negative, red positive)") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank(),
                   strip.placement    = "outside",
                   strip.text.x       = ggplot2::element_text(face = "bold"),
                   legend.position    = "right")

  # ---- 11. Group panel ----
  df_group <- dplyr::transmute(df_ord, Metabolite, group.class, score = 1)
  p_group <- ggplot2::ggplot(df_group, ggplot2::aes(x = score, y = Metabolite, fill = group.class)) +
    ggplot2::geom_col(width = 0.8, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = group.class), hjust = -0.1, size = 2.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1.2),
                                expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::labs(x = "Group", y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank())

  # ---- 12. 合并 ----
  p_combined <- patchwork::wrap_plots(
    p_left, p_right, p_group,
    ncol = 3,
    widths = c(1, 6, 1),
    guides = "collect"
  )

  return(p_combined)
}
