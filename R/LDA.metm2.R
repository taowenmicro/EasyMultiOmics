#' LDA-based microbiome biomarker analysis
#'
#' Identifies group-associated microbial features using a Kruskal-Wallis
#' screening step followed by linear discriminant analysis. The function can
#' work at one taxonomic rank or combine multiple taxonomic ranks.
#'
#' @param ps A phyloseq object.
#' @param group Character. Sample metadata column containing group labels.
#' @param rank Character. Taxonomic rank to analyze, or `"all"` for multiple
#'   taxonomic levels plus the ASV/OTU level.
#' @param Top Integer. Initial feature filtering parameter passed to
#'   `ggClusterNet::filter_OTU_ps()`.
#' @param p.lvl Numeric. P-value or FDR threshold used for screening.
#' @param lda.lvl Numeric. Minimum LDA score retained as a biomarker.
#' @param seed Integer. Random seed.
#' @param adjust.p Logical. If TRUE, use BH-adjusted P values for screening.
#' @param show_n Integer. Maximum number of biomarkers emphasized in plots.
#' @param group_colors Optional named vector of group colors.
#'
#' @return A list containing LEfSe-style lists, significant biomarkers,
#'   all tested features, LDA filtering information, plots, abundance tables,
#'   marker counts, the fitted LDA model, and analysis parameters.
#'
#' @family biomarker analysis
#' @export
LDA.metm2 <- function(
    ps,
    group = "Group",
    rank = "all",
    Top = 10,
    p.lvl = 0.05,
    lda.lvl = 2,
    seed = 11,
    adjust.p = FALSE,
    show_n = 20,
    group_colors = NULL
) {

  ## =========================================================
  ## 1. 基础数据
  ## =========================================================
  ps <- phyloseq::filter_taxa(ps, function(x) sum(x) > 0, TRUE)

  if (!is.null(Top) && Top > 0)
    ps <- ggClusterNet::filter_OTU_ps(ps, Top)

  meta <- data.frame(
    phyloseq::sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!group %in% names(meta))
    stop("sample_data 中不存在分组列: ", group)

  samples <- phyloseq::sample_names(ps)
  meta <- meta[samples, , drop = FALSE]

  grp <- factor(as.character(meta[[group]]))
  names(grp) <- samples

  otu <- as(phyloseq::otu_table(ps), "matrix")
  if (!phyloseq::taxa_are_rows(ps))
    otu <- t(otu)

  otu <- otu[, samples, drop = FALSE]

  tax <- as.data.frame(
    phyloseq::tax_table(ps),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  tax <- tax[rownames(otu), , drop = FALSE]

  ranks <- intersect(
    c("Kingdom", "Phylum", "Class", "Order",
      "Family", "Genus", "Species"),
    colnames(tax)
  )

  prefix <- c(
    Kingdom = "k", Phylum = "p", Class = "c",
    Order = "o", Family = "f", Genus = "g",
    Species = "s"
  )

  clean_tax <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x == "" |
        tolower(x) %in% c(
          "unassigned", "unknown", "unclassified", "na"
        )] <- "Unassigned"
    x
  }

  tax[ranks] <- lapply(tax[ranks], clean_tax)


  ## =========================================================
  ## 2. 构建分析特征
  ## =========================================================
  if (tolower(rank) == "all") {

    mats <- list()
    infos <- list()

    for (i in seq_along(ranks)) {

      r <- ranks[i]

      ## 去除当前层级未注释分类单元
      ok <- tax[[r]] != "Unassigned"
      if (!any(ok)) next

      ## 用完整lineage区分同名taxon
      lineage <- apply(
        tax[ok, ranks[seq_len(i)], drop = FALSE],
        1,
        paste,
        collapse = "_Rank_"
      )

      z <- rowsum(
        otu[ok, , drop = FALSE],
        group = lineage,
        reorder = FALSE
      )

      ids <- paste0(prefix[r], "__", rownames(z))
      rownames(z) <- ids

      info <- data.frame(
        FeatureID = ids,
        Rank = r,
        Taxon = sub(".*_Rank_", "", lineage),
        Lineage = lineage,
        stringsAsFactors = FALSE
      )

      mats[[r]] <- z
      infos[[r]] <- info
    }

    ## ASV/OTU 层级
    otu_asv <- otu
    asv_ids <- paste0("st__", rownames(otu_asv))
    rownames(otu_asv) <- asv_ids

    mats[["ASV"]] <- otu_asv
    infos[["ASV"]] <- data.frame(
      FeatureID = asv_ids,
      Rank = "ASV",
      Taxon = sub("^st__", "", asv_ids),
      Lineage = sub("^st__", "", asv_ids),
      stringsAsFactors = FALSE
    )

    feature_mat <- do.call(rbind, mats)
    feature_info <- do.call(rbind, infos)

  } else {

    rank_match <- ranks[tolower(ranks) == tolower(rank)]

    if (!length(rank_match))
      stop(
        "rank 必须为 all 或现有分类层级: ",
        paste(ranks, collapse = ", ")
      )

    rank <- rank_match[1]

    ok <- tax[[rank]] != "Unassigned"
    if (!any(ok))
      stop(rank, " 层级没有有效分类信息")

    z <- rowsum(
      otu[ok, , drop = FALSE],
      group = tax[[rank]][ok],
      reorder = FALSE
    )

    ids <- paste0(prefix[rank], "__", rownames(z))
    rownames(z) <- ids

    feature_mat <- z

    feature_info <- data.frame(
      FeatureID = ids,
      Rank = rank,
      Taxon = sub("^[^_]+__", "", ids),
      Lineage = sub("^[^_]+__", "", ids),
      stringsAsFactors = FALSE
    )
  }


  ## =========================================================
  ## 3. 样本 × feature
  ## =========================================================
  count <- t(feature_mat)
  count <- count[samples, , drop = FALSE]

  ## 用原始测序深度做标准化
  lib <- phyloseq::sample_sums(ps)[samples]

  rela <- sweep(count, 1, lib, "/") * 100
  rela[!is.finite(rela)] <- 0

  ## LDA使用CPM量级，既标准化测序深度，
  ## 又保留原函数LDA score阈值的量级
  lda_norm <- rela * 10000


  ## =========================================================
  ## 4. Kruskal-Wallis
  ## =========================================================
  rawp <- vapply(
    seq_len(ncol(rela)),
    function(i) {
      x <- rela[, i]
      tryCatch(
        suppressWarnings(
          stats::kruskal.test(x ~ grp)$p.value
        ),
        error = function(e) NA_real_
      )
    },
    numeric(1)
  )

  names(rawp) <- colnames(rela)

  fdr <- stats::p.adjust(rawp, method = "BH")
  names(fdr) <- names(rawp)

  pass <- if (adjust.p) {
    names(fdr)[!is.na(fdr) & fdr <= p.lvl]
  } else {
    names(rawp)[!is.na(rawp) & rawp <= p.lvl]
  }

  ## 强制变成character，彻底避免list下标问题
  pass <- unique(as.character(unlist(pass, use.names = FALSE)))
  pass <- intersect(pass, colnames(lda_norm))


  ## =========================================================
  ## 5. 全部特征结果
  ## =========================================================
  all_results <- feature_info

  all_results$Pvalues <- rawp[
    match(all_results$FeatureID, names(rawp))
  ]

  all_results$FDR <- fdr[
    match(all_results$FeatureID, names(fdr))
  ]

  all_results$MeanAbundance <- colMeans(rela)[
    match(all_results$FeatureID, colnames(rela))
  ]

  all_results$Prevalence <- colMeans(count > 0)[
    match(all_results$FeatureID, colnames(count))
  ] * 100


  ## =========================================================
  ## 6. 无候选特征安全退出
  ## =========================================================
  empty_return <- function(msg) {
    message(msg)

    list(
      data.frame(),
      data.frame(),
      all_results = all_results,
      lda_filter = data.frame(),
      message = msg
    )
  }

  if (!length(pass))
    return(
      empty_return("No features passed P/FDR threshold.")
    )


  ## =========================================================
  ## 7. LDA前处理：去零变异 + 共线变量
  ## =========================================================
  pass_before <- pass

  X <- as.matrix(
    lda_norm[, pass, drop = FALSE]
  )

  ## 零值设极小值，避免完全零结构
  X[X == 0] <- 1

  ## 组内中心化
  Xw <- X

  for (g in levels(grp)) {
    ii <- grp == g

    Xw[ii, ] <- sweep(
      X[ii, , drop = FALSE],
      2,
      colMeans(X[ii, , drop = FALSE]),
      "-"
    )
  }

  ## 去除组内完全无变化变量
  within_sd <- apply(Xw, 2, stats::sd, na.rm = TRUE)
  keep_var <- is.finite(within_sd) & within_sd > 1e-12

  X <- X[, keep_var, drop = FALSE]
  Xw <- Xw[, keep_var, drop = FALSE]

  if (!ncol(X))
    return(
      empty_return(
        "No features remained after within-group variance filtering."
      )
    )

  ## QR去共线
  if (ncol(X) > 1) {

    qx <- qr(Xw, tol = 1e-7)

    if (qx$rank > 0) {
      keep <- qx$pivot[seq_len(qx$rank)]
      X <- X[, keep, drop = FALSE]
    }
  }

  if (!ncol(X))
    return(
      empty_return(
        "No non-collinear features remained for LDA."
      )
    )

  lda_used <- colnames(X)

  lda_filter <- data.frame(
    FeatureID = pass_before,
    Used_in_LDA = pass_before %in% lda_used,
    stringsAsFactors = FALSE
  )

  message(
    "LDA input: ",
    length(pass_before),
    " candidates -> ",
    length(lda_used),
    " non-collinear features"
  )


  ## =========================================================
  ## 8. LDA
  ## =========================================================
  set.seed(seed)

  ldares <- suppressWarnings(
    MASS::lda(
      x = as.data.frame(X),
      grouping = grp
    )
  )

  means <- as.data.frame(
    t(ldares$means),
    check.names = FALSE
  )

  gcols <- colnames(means)

  means$max <- apply(
    means[, gcols, drop = FALSE],
    1,
    max
  )

  means$min <- apply(
    means[, gcols, drop = FALSE],
    1,
    min
  )

  ## 延续原LDA.metm score定义
  means$LDAscore <- signif(
    log10(
      1 + abs(means$max - means$min) / 2
    ),
    4
  )

  means$class <- gcols[
    max.col(
      as.matrix(means[, gcols, drop = FALSE]),
      ties.method = "first"
    )
  ]

  fid <- rownames(means)

  means$Pvalues <- rawp[fid]
  means$FDR <- fdr[fid]

  m <- match(fid, feature_info$FeatureID)

  means$Rank <- feature_info$Rank[m]
  means$Taxon <- feature_info$Taxon[m]
  means$Lineage <- feature_info$Lineage[m]

  means$MeanAbundance <- colMeans(rela)[fid]
  means$Prevalence <- colMeans(count[, fid, drop = FALSE] > 0) * 100


  ## =========================================================
  ## 9. P/FDR + LDA 双阈值
  ## =========================================================
  keep <- if (adjust.p) {
    means$FDR <= p.lvl &
      means$LDAscore >= lda.lvl
  } else {
    means$Pvalues <= p.lvl &
      means$LDAscore >= lda.lvl
  }

  keep[is.na(keep)] <- FALSE

  taxtree <- means[
    keep,
    ,
    drop = FALSE
  ]

  taxtree <- taxtree[
    order(-taxtree$LDAscore),
    ,
    drop = FALSE
  ]

  message(
    "A total of ",
    nrow(taxtree),
    " significant features: ",
    if (adjust.p) "FDR" else "P",
    " <= ", p.lvl,
    " and LDA >= ", lda.lvl
  )


  ## 把LDA结果补充到all_results
  lda_result <- data.frame(
    FeatureID = rownames(means),
    LDAscore = means$LDAscore,
    EnrichedGroup = means$class,
    stringsAsFactors = FALSE
  )

  all_results <- dplyr::left_join(
    all_results,
    lda_result,
    by = "FeatureID"
  )


  if (!nrow(taxtree)) {

    return(
      list(
        data.frame(),
        taxtree,
        all_results = all_results,
        lda_filter = lda_filter,
        lda_model = ldares
      )
    )
  }


  ## =========================================================
  ## 10. 分组颜色
  ## =========================================================
  lev <- levels(grp)

  if (is.null(group_colors)) {

    cols <- setNames(
      scales::hue_pal()(length(lev)),
      lev
    )

  } else {

    cols <- group_colors

    if (is.null(names(cols)))
      names(cols) <- lev[seq_len(min(length(cols), length(lev)))]

    miss <- setdiff(lev, names(cols))

    if (length(miss)) {
      cols <- c(
        cols,
        setNames(
          scales::hue_pal()(length(miss)),
          miss
        )
      )
    }
  }


  ## =========================================================
  ## 11. lefse兼容结果
  ## =========================================================
  taxtree$color <- unname(
    cols[as.character(taxtree$class)]
  )

  lefse_lists <- data.frame(
    node = rownames(taxtree),
    color = taxtree$color,
    Group = taxtree$class,
    stringsAsFactors = FALSE
  )


  ## =========================================================
  ## 12. Top biomarkers
  ## =========================================================
  showtab <- head(
    taxtree,
    min(show_n, nrow(taxtree))
  )

  showtab$FeatureID <- rownames(showtab)

  showtab$Label <- make.unique(
    paste0(
      showtab$Taxon,
      " [",
      showtab$Rank,
      "]"
    )
  )

  feature_label <- setNames(
    showtab$Label,
    showtab$FeatureID
  )

  label_order <- rev(showtab$Label)

  showtab$Label <- factor(
    showtab$Label,
    levels = label_order
  )


  ## =========================================================
  ## 13. Top taxa真实相对丰度
  ## =========================================================
  ab <- as.data.frame(
    rela[, showtab$FeatureID, drop = FALSE],
    check.names = FALSE
  )

  ab$SampleID <- rownames(ab)
  ab$Group <- unname(
    as.character(grp[ab$SampleID])
  )

  ab_long <- tidyr::pivot_longer(
    ab,
    cols = -c(SampleID, Group),
    names_to = "FeatureID",
    values_to = "Abundance"
  )

  ab_long$Label <- feature_label[
    ab_long$FeatureID
  ]

  ab_long$Label <- factor(
    ab_long$Label,
    levels = label_order
  )

  ab_sum <- ab_long |>
    dplyr::group_by(
      FeatureID,
      Label,
      Group
    ) |>
    dplyr::summarise(
      Mean = mean(Abundance),
      SD = stats::sd(Abundance),
      SE = SD / sqrt(dplyr::n()),
      .groups = "drop"
    )

  group_mean <- ab_sum |>
    dplyr::select(
      Label,
      Group,
      Mean
    ) |>
    tidyr::pivot_wider(
      names_from = Group,
      values_from = Mean
    )


  ## =========================================================
  ## 14. 图1 LDA score
  ## =========================================================
  p_bar <- ggplot2::ggplot(
    showtab,
    ggplot2::aes(
      x = Label,
      y = LDAscore,
      fill = class
    )
  ) +
    ggplot2::geom_col(width = 0.72) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = sprintf("%.2f", LDAscore)
      ),
      hjust = -0.12,
      size = 3
    ) +
    ggplot2::coord_flip(clip = "off") +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(
      x = NULL,
      y = "LDA score",
      fill = group
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
        5.5, 30, 5.5, 5.5
      )
    )


  ## =========================================================
  ## 15. 图2 LDA × significance × abundance
  ## =========================================================
  showtab$mlogP <- -log10(
    pmax(
      showtab$Pvalues,
      .Machine$double.xmin
    )
  )

  p_dot <- ggplot2::ggplot(
    showtab,
    ggplot2::aes(
      x = LDAscore,
      y = mlogP,
      color = class,
      size = MeanAbundance
    )
  ) +
    ggplot2::geom_point(alpha = 0.85) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = Taxon),
      size = 3,
      show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      x = "LDA score",
      y = expression(-log[10](P)),
      size = "Mean abundance (%)",
      color = group
    ) +
    ggplot2::theme_classic()


  ## =========================================================
  ## 16. 图3 平均丰度heatmap
  ## =========================================================
  p_heat <- ggplot2::ggplot(
    ab_sum,
    ggplot2::aes(
      x = Group,
      y = Label,
      fill = Mean
    )
  ) +
    ggplot2::geom_tile(
      color = "white",
      linewidth = 0.3
    ) +
    ggplot2::scale_fill_viridis_c(
      trans = "sqrt"
    ) +
    ggplot2::labs(
      x = NULL,
      y = NULL,
      fill = "Mean relative\nabundance (%)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )


  ## =========================================================
  ## 17. 图4 样本真实丰度
  ## =========================================================
  p_abund <- ggplot2::ggplot(
    ab_long,
    ggplot2::aes(
      x = Group,
      y = Abundance,
      color = Group
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.6,
      linewidth = 0.4
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      size = 1.6,
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(
      ~ Label,
      scales = "free_y",
      ncol = 4
    ) +
    ggplot2::scale_color_manual(values = cols) +
    ggplot2::labs(
      x = NULL,
      y = "Relative abundance (%)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      )
    )


  ## =========================================================
  ## 18. 图5 每组marker数量
  ## =========================================================
  marker_count <- as.data.frame(
    table(taxtree$class)
  )

  colnames(marker_count) <- c(
    "Group",
    "N"
  )

  marker_count$Group <- as.character(
    marker_count$Group
  )

  p_count <- ggplot2::ggplot(
    marker_count,
    ggplot2::aes(
      x = Group,
      y = N,
      fill = Group
    )
  ) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(
      ggplot2::aes(label = N),
      vjust = -0.3
    ) +
    ggplot2::scale_fill_manual(values = cols) +
    ggplot2::labs(
      x = NULL,
      y = "Number of biomarkers"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none"
    )


  ## =========================================================
  ## 19. 组合图
  ## =========================================================
  p_combined <- NULL

  if (requireNamespace("ggpubr", quietly = TRUE)) {
    p_combined <- ggpubr::ggarrange(
      p_bar,
      p_heat,
      ncol = 2,
      widths = c(1.05, 1.35),
      align = "hv"
    )
  }


  ## =========================================================
  ## 20. 返回
  ## =========================================================
  list(
    lefse_lists,                 # [[1]]
    taxtree,                     # [[2]]

    all_results = all_results,
    lda_filter = lda_filter,

    lda_bar = p_bar,
    significance_plot = p_dot,
    abundance_heatmap = p_heat,
    abundance_plot = p_abund,
    marker_count_plot = p_count,
    combined_plot = p_combined,

    abundance_raw = ab_long,
    abundance_summary = ab_sum,
    group_mean = group_mean,
    marker_count = marker_count,

    lda_model = ldares,

    parameters = data.frame(
      Parameter = c(
        "rank", "Top", "p.lvl",
        "lda.lvl", "seed",
        "adjust.p", "show_n"
      ),
      Value = c(
        rank, Top, p.lvl,
        lda.lvl, seed,
        adjust.p, show_n
      )
    )
  )
}
