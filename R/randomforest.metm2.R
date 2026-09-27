#' Random-forest biomarker analysis
#'
#' Fits a random-forest classifier to microbiome features with optional
#' taxonomic aggregation, abundance/prevalence filtering, unsupervised
#' variability filtering, feature importance, proximity ordination, abundance
#' summaries, and random-forest feature-number cross-validation.
#'
#' @param ps A phyloseq object.
#' @param group Character. Sample metadata column containing group labels.
#' @param rank Character or NULL. Taxonomic rank used for aggregation.
#' @param optimal Integer. Number of top important features retained for display.
#' @param ntree Integer. Number of trees.
#' @param mtry Integer or NULL. Number of variables sampled at each split.
#' @param seed Integer. Random seed.
#' @param transform Character. Feature transformation method.
#' @param prevalence Numeric. Minimum feature prevalence as a proportion.
#' @param min_mean Numeric. Minimum mean relative abundance in percent.
#' @param feature_filter Character. Unsupervised variability filter.
#' @param variable_top Integer. Maximum number of variable features retained.
#' @param fill Character or NULL. Taxonomic rank used for feature colors.
#' @param lab Character or NULL. Taxonomic rank used as feature labels.
#' @param group_colors Optional named vector of group colors.
#' @param rfcv Logical. Whether to run feature-number cross-validation.
#' @param nrfcvnum Integer. Number of folds used by RF feature-number CV.
#' @param pseudocount Numeric. Pseudocount used for CLR transformation.
#'
#' @return A list containing the fitted random forest, model statistics,
#'   confusion matrices and class performance, feature-importance tables,
#'   diagnostic and abundance plots, filtering results, MDS results,
#'   RF-CV results, and analysis parameters.
#'
#' @family biomarker machine learning
#' @export
randomforest.metm2 <- function(
    ps,
    group = "Group",
    rank = NULL,
    optimal = 20,
    ntree = 1000,
    mtry = NULL,
    seed = 11,
    transform = c("hellinger", "relative", "clr", "none"),
    prevalence = 0.30,
    min_mean = 0.01,
    feature_filter = c("MAD", "IQR", "variance", "none"),
    variable_top = 100,
    fill = "Phylum",
    lab = NULL,
    group_colors = NULL,
    rfcv = TRUE,
    nrfcvnum = 3,
    pseudocount = 0.5
) {

  transform <- match.arg(transform)
  feature_filter <- match.arg(feature_filter)

  ## =========================================================
  ## 1. 分类层级聚合
  ## =========================================================
  ps <- phyloseq::filter_taxa(ps, function(x) sum(x) > 0, TRUE)

  if (!is.null(rank)) {

    tt <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
    if (is.null(tt))
      stop("phyloseq对象没有tax_table，不能按 ", rank, " 聚合")

    tax0 <- as.data.frame(tt, stringsAsFactors = FALSE)

    if (!rank %in% colnames(tax0))
      stop("tax_table中不存在分类层级: ", rank)

    x <- trimws(as.character(tax0[[rank]]))

    ok <- !is.na(x) & nzchar(x) &
      !tolower(x) %in% c(
        "unassigned", "unknown",
        "unclassified", "na"
      )

    ps <- phyloseq::prune_taxa(
      rownames(tax0)[ok], ps
    )

    ps <- phyloseq::tax_glom(
      ps,
      taxrank = rank,
      NArm = TRUE
    )
  }

  if (is.null(lab))
    lab <- if (is.null(rank)) "FeatureID" else rank


  ## =========================================================
  ## 2. OTU + metadata
  ## =========================================================
  otu <- as(
    phyloseq::otu_table(ps),
    "matrix"
  )

  if (phyloseq::taxa_are_rows(ps))
    otu <- t(otu)

  meta <- data.frame(
    phyloseq::sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!group %in% colnames(meta))
    stop("sample_data中不存在分组列: ", group)

  samples <- rownames(otu)
  meta <- meta[samples, , drop = FALSE]

  grp <- factor(
    as.character(meta[[group]])
  )

  names(grp) <- samples

  if (anyNA(grp))
    stop("分组变量中存在NA")

  if (nlevels(grp) < 2)
    stop("Random Forest至少需要2个分组")

  lib <- rowSums(otu)

  if (any(lib <= 0))
    stop("存在总reads为0的样本")


  ## =========================================================
  ## 3. 相对丰度 + 基础过滤
  ## =========================================================
  rel <- sweep(
    otu,
    1,
    lib,
    "/"
  )

  rel[!is.finite(rel)] <- 0

  prev <- colMeans(otu > 0)
  mean_abund <- colMeans(rel) * 100

  eligible <-
    colSums(otu) > 0 &
    prev >= prevalence &
    mean_abund >= min_mean

  if (!any(eligible))
    stop("prevalence/min_mean过滤后没有有效特征")


  ## =========================================================
  ## 4. 转换
  ## =========================================================
  X0 <- switch(
    transform,

    hellinger = sqrt(rel),

    relative = rel,

    clr = {
      z <- otu + pseudocount
      lz <- log(z)
      sweep(lz, 1, rowMeans(lz), "-")
    },

    none = otu
  )


  ## =========================================================
  ## 5. 无监督波动筛选
  ## =========================================================
  Xe <- X0[, eligible, drop = FALSE]

  variability <- switch(
    feature_filter,

    MAD = apply(
      Xe, 2,
      stats::mad,
      na.rm = TRUE
    ),

    IQR = apply(
      Xe, 2,
      stats::IQR,
      na.rm = TRUE
    ),

    variance = apply(
      Xe, 2,
      stats::var,
      na.rm = TRUE
    ),

    none = setNames(
      rep(NA_real_, ncol(Xe)),
      colnames(Xe)
    )
  )

  selected <- colnames(Xe)

  if (
    feature_filter != "none" &&
    !is.null(variable_top) &&
    variable_top > 0
  ) {

    ord <- order(
      variability,
      decreasing = TRUE,
      na.last = NA
    )

    selected <- names(variability)[
      head(
        ord,
        min(variable_top, length(ord))
      )
    ]
  }

  X <- X0[, selected, drop = FALSE]

  ## 最终去除零变异变量
  sdx <- apply(
    X, 2,
    stats::sd,
    na.rm = TRUE
  )

  X <- X[
    ,
    is.finite(sdx) & sdx > 1e-12,
    drop = FALSE
  ]

  if (ncol(X) < 2)
    stop("筛选后可用于Random Forest的特征少于2个")


  ## =========================================================
  ## 6. 筛选信息
  ## =========================================================
  filter_table <- data.frame(
    FeatureID = colnames(otu),
    Prevalence = prev * 100,
    MeanAbundance = mean_abund,
    Variability = NA_real_,
    Eligible = eligible,
    SelectedByVariability = colnames(otu) %in% selected,
    UsedInRF = colnames(otu) %in% colnames(X),
    stringsAsFactors = FALSE
  )

  ii <- match(
    names(variability),
    filter_table$FeatureID
  )

  filter_table$Variability[ii] <-
    variability


  ## =========================================================
  ## 7. Random Forest
  ## =========================================================
  if (is.null(mtry))
    mtry <- max(
      1,
      floor(sqrt(ncol(X)))
    )

  mtry <- min(
    mtry,
    ncol(X)
  )

  set.seed(seed)

  model <- randomForest::randomForest(
    x = X,
    y = grp,
    ntree = ntree,
    mtry = mtry,
    importance = TRUE,
    proximity = TRUE
  )


  ## =========================================================
  ## 8. 混淆矩阵 + 分类性能
  ## =========================================================
  lev <- levels(grp)

  pred <- factor(
    as.character(model$predicted),
    levels = lev
  )

  cm <- table(
    Actual = factor(grp, levels = lev),
    Predicted = pred
  )

  correct <- diag(cm)
  actual_n <- rowSums(cm)
  predicted_n <- colSums(cm)

  recall <- correct / actual_n
  precision <- correct / predicted_n
  f1 <- 2 * precision * recall /
    (precision + recall)

  recall[!is.finite(recall)] <- NA_real_
  precision[!is.finite(precision)] <- NA_real_
  f1[!is.finite(f1)] <- NA_real_

  class_performance <- data.frame(
    Group = lev,
    N = as.numeric(actual_n),
    Correct = as.numeric(correct),
    Recall = as.numeric(recall),
    Precision = as.numeric(precision),
    F1 = as.numeric(f1),
    stringsAsFactors = FALSE
  )

  oob_error <- tail(
    model$err.rate[, "OOB"],
    1
  )

  model_stats <- data.frame(
    Samples = nrow(X),
    Groups = nlevels(grp),
    Features = ncol(X),
    Trees = ntree,
    Mtry = mtry,
    OOB_error = oob_error,
    OOB_accuracy = 1 - oob_error,
    Balanced_accuracy = mean(
      recall,
      na.rm = TRUE
    ),
    Macro_F1 = mean(
      f1,
      na.rm = TRUE
    )
  )

  confusion <- data.frame(
    Actual = rownames(cm),
    as.data.frame.matrix(cm),
    check.names = FALSE
  )


  ## =========================================================
  ## 9. Importance
  ## =========================================================
  imp <- as.data.frame(
    randomForest::importance(
      model,
      scale = FALSE
    ),
    check.names = FALSE
  )

  imp$FeatureID <- rownames(imp)

  if (!"MeanDecreaseAccuracy" %in% colnames(imp))
    stop("未获得MeanDecreaseAccuracy")

  ## taxonomy
  tt <- phyloseq::tax_table(
    ps,
    errorIfNULL = FALSE
  )

  if (!is.null(tt)) {

    tax <- as.data.frame(
      tt,
      stringsAsFactors = FALSE
    )

    tax$FeatureID <- rownames(tax)

    imp <- dplyr::left_join(
      imp,
      tax,
      by = "FeatureID"
    )
  }

  ## 加入过滤信息
  imp <- dplyr::left_join(
    imp,
    filter_table,
    by = "FeatureID"
  )

  ## 显示名称
  if (
    !is.null(lab) &&
    lab %in% colnames(imp)
  ) {
    imp$Label <- as.character(
      imp[[lab]]
    )
  } else {
    imp$Label <- imp$FeatureID
  }

  bad <- is.na(imp$Label) |
    !nzchar(imp$Label) |
    tolower(imp$Label) %in%
    c(
      "unassigned",
      "unknown",
      "unclassified",
      "na"
    )

  imp$Label[bad] <-
    imp$FeatureID[bad]

  ## 分类颜色
  if (
    !is.null(fill) &&
    fill %in% colnames(imp)
  ) {

    imp$Fill <- as.character(
      imp[[fill]]
    )

    imp$Fill[
      is.na(imp$Fill) |
        !nzchar(imp$Fill)
    ] <- "Unassigned"

  } else {

    imp$Fill <- "Taxon"
  }

  imp <- imp[
    order(
      -imp$MeanDecreaseAccuracy,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]

  top <- head(
    imp,
    min(optimal, nrow(imp))
  )

  top$LabelPlot <- make.unique(
    top$Label
  )

  top$LabelPlot <- factor(
    top$LabelPlot,
    levels = rev(top$LabelPlot)
  )


  ## =========================================================
  ## 10. 颜色
  ## =========================================================
  base_cols <- setNames(
    grDevices::hcl.colors(
      length(lev),
      "Dark 3"
    ),
    lev
  )

  if (is.null(group_colors)) {

    group_colors <- base_cols

  } else {

    group_colors <-
      as.character(group_colors)

    if (is.null(names(group_colors)))
      names(group_colors) <-
        lev[seq_len(
          min(length(group_colors), length(lev))
        )]

    miss <- setdiff(
      lev,
      names(group_colors)
    )

    group_colors <- c(
      group_colors,
      base_cols[miss]
    )
  }

  fill_levels <- unique(
    top$Fill
  )

  fill_cols <- setNames(
    grDevices::hcl.colors(
      length(fill_levels),
      "Dynamic"
    ),
    fill_levels
  )


  ## =========================================================
  ## 11. Top菌真实相对丰度
  ## =========================================================
  top_ids <- top$FeatureID

  ab <- as.data.frame(
    rel[, top_ids, drop = FALSE] * 100,
    check.names = FALSE
  )

  ab$SampleID <- rownames(ab)
  ab$Group <- as.character(grp)

  ab_long <- tidyr::pivot_longer(
    ab,
    cols = -c(SampleID, Group),
    names_to = "FeatureID",
    values_to = "Abundance"
  )

  label_map <- setNames(
    as.character(top$LabelPlot),
    top$FeatureID
  )

  ab_long$Label <-
    label_map[
      ab_long$FeatureID
    ]

  ab_long$Label <- factor(
    ab_long$Label,
    levels = levels(top$LabelPlot)
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


  ## =========================================================
  ## 12. Importance图
  ## =========================================================
  p_imp <- ggplot2::ggplot(
    top,
    ggplot2::aes(
      x = MeanDecreaseAccuracy,
      y = LabelPlot,
      color = Fill
    )
  ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = MeanDecreaseAccuracy,
        yend = LabelPlot
      ),
      linewidth = 1.3
    ) +
    ggplot2::geom_point(size = 4) +
    ggplot2::scale_color_manual(
      values = fill_cols
    ) +
    ggplot2::labs(
      x = "Mean decrease accuracy",
      y = NULL,
      color = fill
    ) +
    ggplot2::theme_classic()


  ## =========================================================
  ## 13. OOB error
  ## =========================================================
  err <- as.data.frame(
    model$err.rate,
    check.names = FALSE
  )

  err$Tree <- seq_len(
    nrow(err)
  )

  err_long <- tidyr::pivot_longer(
    err,
    cols = -Tree,
    names_to = "Class",
    values_to = "Error"
  )

  err_cols <- c(
    OOB = "black",
    group_colors
  )

  p_err <- ggplot2::ggplot(
    err_long,
    ggplot2::aes(
      Tree,
      Error,
      color = Class
    )
  ) +
    ggplot2::geom_line(
      linewidth = 0.7
    ) +
    ggplot2::scale_color_manual(
      values = err_cols
    ) +
    ggplot2::labs(
      x = "Number of trees",
      y = "OOB error rate"
    ) +
    ggplot2::theme_classic()


  ## =========================================================
  ## 14. Confusion matrix
  ## =========================================================
  cm_long <- as.data.frame(
    as.table(cm)
  )

  colnames(cm_long) <- c(
    "Actual",
    "Predicted",
    "N"
  )

  row_pct <- sweep(
    cm,
    1,
    rowSums(cm),
    "/"
  ) * 100

  pct_long <- as.data.frame(
    as.table(row_pct)
  )

  cm_long$Percent <-
    pct_long$Freq

  cm_long$Label <- paste0(
    cm_long$N,
    "\n",
    round(cm_long$Percent, 1),
    "%"
  )

  p_conf <- ggplot2::ggplot(
    cm_long,
    ggplot2::aes(
      Predicted,
      Actual,
      fill = Percent
    )
  ) +
    ggplot2::geom_tile(
      color = "white"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = Label
      ),
      size = 4
    ) +
    ggplot2::scale_fill_viridis_c(
      limits = c(0, 100)
    ) +
    ggplot2::labs(
      fill = "Row %",
      x = "Predicted",
      y = "Actual"
    ) +
    ggplot2::theme_classic()


  ## =========================================================
  ## 15. RF proximity MDS
  ## =========================================================
  mds <- tryCatch({

    z <- stats::cmdscale(
      stats::as.dist(
        1 - model$proximity
      ),
      k = 2,
      eig = TRUE,
      add = TRUE
    )

    data.frame(
      MDS1 = z$points[, 1],
      MDS2 = z$points[, 2],
      Group = grp,
      SampleID = rownames(X)
    )

  }, error = function(e) NULL)

  p_mds <- NULL

  if (!is.null(mds)) {

    p_mds <- ggplot2::ggplot(
      mds,
      ggplot2::aes(
        MDS1,
        MDS2,
        color = Group,
        fill = Group
      )
    ) +
      ggplot2::geom_point(
        size = 4,
        shape = 21
      ) +
      ggplot2::scale_color_manual(
        values = group_colors
      ) +
      ggplot2::scale_fill_manual(
        values = group_colors
      ) +
      ggplot2::theme_classic()
  }


  ## =========================================================
  ## 16. 丰度Heatmap
  ## =========================================================
  p_heat <- ggplot2::ggplot(
    ab_sum,
    ggplot2::aes(
      Group,
      Label,
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
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
    )


  ## =========================================================
  ## 17. 真实样本丰度
  ## =========================================================
  show_ids <- head(
    top_ids,
    min(12, length(top_ids))
  )

  ab_show <- ab_long[
    ab_long$FeatureID %in% show_ids,
    ,
    drop = FALSE
  ]

  p_abund <- ggplot2::ggplot(
    ab_show,
    ggplot2::aes(
      Group,
      Abundance,
      color = Group
    )
  ) +
    ggplot2::geom_boxplot(
      outlier.shape = NA,
      width = 0.6
    ) +
    ggplot2::geom_jitter(
      width = 0.12,
      size = 1.5,
      alpha = 0.85
    ) +
    ggplot2::facet_wrap(
      ~ Label,
      scales = "free_y",
      ncol = 4
    ) +
    ggplot2::scale_color_manual(
      values = group_colors
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Relative abundance (%)"
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "none",
      axis.text.x =
        ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
    )


  ## =========================================================
  ## 18. 波动筛选图
  ## =========================================================
  p_filter <- NULL

  if (feature_filter != "none") {

    ft <- filter_table[
      filter_table$Eligible &
        is.finite(filter_table$Variability),
      ,
      drop = FALSE
    ]

    ft <- head(
      ft[
        order(
          -ft$Variability
        ),
        ,
        drop = FALSE
      ],
      30
    )

    if (nrow(ft)) {

      ft$FeatureID <- factor(
        ft$FeatureID,
        levels = rev(ft$FeatureID)
      )

      p_filter <- ggplot2::ggplot(
        ft,
        ggplot2::aes(
          Variability,
          FeatureID
        )
      ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0,
            xend = Variability,
            yend = FeatureID
          ),
          linewidth = 1
        ) +
        ggplot2::geom_point(
          size = 3
        ) +
        ggplot2::labs(
          x = paste0(
            feature_filter,
            " variability"
          ),
          y = NULL
        ) +
        ggplot2::theme_classic()
    }
  }


  ## =========================================================
  ## 19. 组合图
  ## =========================================================
  p_combined <- NULL

  if (
    requireNamespace(
      "patchwork",
      quietly = TRUE
    )
  ) {

    p_combined <- patchwork::wrap_plots(
      p_imp,
      p_heat,
      ncol = 2,
      widths = c(1.1, 1.3)
    )
  }


  ## =========================================================
  ## 20. RF-CV：特征数量
  ## =========================================================
  cv_result <- NULL
  cv_data <- NULL
  p_cv <- NULL

  if (rfcv && ncol(X) >= 2) {

    folds <- min(
      nrfcvnum,
      min(table(grp))
    )

    cv_result <- tryCatch(

      randomForest::rfcv(
        trainx = as.data.frame(X),
        trainy = grp,
        cv.fold = folds,
        step = 0.5
      ),

      error = function(e) {
        message(
          "rfcv failed: ",
          conditionMessage(e)
        )
        NULL
      }
    )

    if (!is.null(cv_result)) {

      cv_data <- data.frame(
        Features = cv_result$n.var,
        CV_error = cv_result$error.cv
      )

      p_cv <- ggplot2::ggplot(
        cv_data,
        ggplot2::aes(
          Features,
          CV_error
        )
      ) +
        ggplot2::geom_line() +
        ggplot2::geom_point(
          size = 3
        ) +
        ggplot2::scale_x_log10() +
        ggplot2::labs(
          x = "Number of features",
          y = "Cross-validation error"
        ) +
        ggplot2::theme_classic()
    }
  }


  ## =========================================================
  ## 21. 参数
  ## =========================================================
  parameters <- data.frame(
    Parameter = c(
      "group", "rank", "optimal",
      "ntree", "mtry", "seed",
      "transform", "prevalence",
      "min_mean_percent",
      "feature_filter",
      "variable_top",
      "rfcv", "nrfcvnum"
    ),
    Value = c(
      group,
      ifelse(is.null(rank), "NULL", rank),
      optimal,
      ntree,
      mtry,
      seed,
      transform,
      prevalence,
      min_mean,
      feature_filter,
      variable_top,
      rfcv,
      nrfcvnum
    ),
    stringsAsFactors = FALSE
  )


  ## =========================================================
  ## 22. Return
  ## =========================================================
  list(
    p_imp,                    # [[1]]
    p_err,                    # [[2]]
    top,                      # [[3]]
    p_conf,                   # [[4]]
    imp,                      # [[5]]

    model = model,
    model_stats = model_stats,
    confusion = confusion,
    confusion_long = cm_long,
    class_performance = class_performance,

    importance_plot = p_imp,
    error_plot = p_err,
    confusion_plot = p_conf,
    mds_plot = p_mds,
    abundance_heatmap = p_heat,
    abundance_plot = p_abund,
    variability_plot = p_filter,
    combined_plot = p_combined,
    rfcv_plot = p_cv,

    filter_table = filter_table,
    abundance_raw = ab_long,
    abundance_summary = ab_sum,
    mds_data = mds,
    error_data = err_long,
    rfcv_data = cv_data,
    rfcv = cv_result,

    parameters = parameters
  )
}
