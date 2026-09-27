#' Automatic statistical analysis for alpha diversity
#'
#' Performs automatic statistical testing for multiple alpha-diversity
#' indices. The primary analysis uses \code{MuiaovMcomper2}. When the
#' primary method does not detect a significant difference, the function
#' can optionally apply Welch's ANOVA followed by Games-Howell pairwise
#' comparisons as a fallback.
#'
#' Invalid alpha-diversity metrics containing no usable variation or
#' insufficient finite observations are automatically excluded.
#'
#' @param data A data frame containing sample IDs, group information,
#'   and alpha-diversity indices.
#' @param num Numeric vector specifying the columns containing alpha-diversity
#'   indices. Default is \code{3:ncol(data)}.
#' @param group_col Character string giving the grouping column name.
#'   Default is \code{"group"}.
#' @param alpha Numeric significance threshold. Default is \code{0.05}.
#' @param fallback Logical. Whether Welch ANOVA and Games-Howell tests
#'   should be used when the primary method does not detect a difference.
#'   Default is \code{TRUE}.
#'
#' @return A list containing:
#' \describe{
#'   \item{result}{Final significance-letter table for use with
#'     \code{FacetMuiPlot} functions.}
#'   \item{data}{Cleaned alpha-diversity data used for statistical analysis.}
#'   \item{primary_result}{Original result returned by
#'     \code{MuiaovMcomper2}.}
#'   \item{summary}{Summary of the statistical method used for each metric.}
#'   \item{games_howell}{Complete Games-Howell pairwise-comparison results.}
#'   \item{letters}{Final compact-letter display results.}
#'   \item{invalid_metrics}{Metrics excluded from analysis.}
#' }
#'
#' @details
#' Non-finite values including \code{NaN}, \code{Inf}, and \code{-Inf}
#' are converted to missing values before analysis.
#'
#' The primary analysis is performed using \code{MuiaovMcomper2}.
#' If no real difference is represented by the compact-letter display,
#' and \code{fallback = TRUE}, Welch's ANOVA and Games-Howell tests are
#' evaluated. Games-Howell-derived letters replace the primary letters
#' only when both the Welch test and at least one Games-Howell comparison
#' are significant.
#'
#' @examples
#' \dontrun{
#' res <- alpha_auto_stat.metm(
#'   data = alpha_data,
#'   num = 3:ncol(alpha_data),
#'   group_col = "group"
#' )
#'
#' res$summary
#' res$letters
#' }
#'
#' @export
alpha_auto_stat.metm <- function(
    data,
    num = 3:ncol(data),
    group_col = "group",
    alpha = 0.05,
    fallback = TRUE
) {

  if (!requireNamespace("rstatix", quietly = TRUE)) {
    stop("请先安装 rstatix: install.packages('rstatix')")
  }

  if (!requireNamespace("multcompView", quietly = TRUE)) {
    stop("请先安装 multcompView: install.packages('multcompView')")
  }

  data <- as.data.frame(data)

  metric_names <- colnames(data)[num]

  # ----------------------------------------------------------
  # NaN / Inf -> NA
  # ----------------------------------------------------------

  for (m in metric_names) {

    data[[m]] <- suppressWarnings(
      as.numeric(data[[m]])
    )

    data[[m]][
      !is.finite(data[[m]])
    ] <- NA
  }

  # ----------------------------------------------------------
  # 自动筛选有效指标
  # ----------------------------------------------------------

  valid_metrics <- metric_names[
    vapply(
      metric_names,
      function(m) {

        x <- data[[m]]
        ok <- !is.na(x)

        # 至少有有效数据
        enough_data <- sum(ok) >= 2

        # 至少两个组有数据
        enough_group <-
          length(
            unique(
              data[[group_col]][ok]
            )
          ) >= 2

        # 指标本身不能完全恒定
        variable <-
          length(
            unique(x[ok])
          ) >= 2

        enough_data &&
          enough_group &&
          variable
      },
      logical(1)
    )
  ]

  invalid_metrics <- setdiff(
    metric_names,
    valid_metrics
  )

  if (length(valid_metrics) == 0) {
    stop("没有可用于 Alpha 统计的有效指标")
  }

  # ----------------------------------------------------------
  # 构造统计数据
  # ----------------------------------------------------------

  data_use <- data.frame(
    ID = data[[1]],
    group = data[[group_col]],
    data[, valid_metrics, drop = FALSE],
    check.names = FALSE
  )

  # 保持原始分组顺序
  if (is.factor(data[[group_col]])) {

    group_levels <- levels(
      droplevels(data[[group_col]])
    )

  } else {

    group_levels <- unique(
      as.character(data[[group_col]])
    )
  }

  data_use$group <- factor(
    data_use$group,
    levels = group_levels
  )

  groups <- levels(data_use$group)

  # ==========================================================
  # 1. Primary：原 MuiaovMcomper2
  # ==========================================================

  result_primary <- MuiaovMcomper2(
    data = data_use,
    num = 3:ncol(data_use)
  )

  result_final <- result_primary

  result_df <- as.data.frame(
    result_primary,
    check.names = FALSE
  )

  # ----------------------------------------------------------
  # 找 result 中的 Group 顺序
  # ----------------------------------------------------------

  if ("group" %in% colnames(result_df)) {

    result_groups <- as.character(
      result_df$group
    )

  } else if (
    !is.null(rownames(result_df)) &&
    all(rownames(result_df) %in% groups)
  ) {

    result_groups <- rownames(result_df)

  } else if (
    nrow(result_df) == length(groups)
  ) {

    result_groups <- groups

  } else {

    stop(
      "无法识别 MuiaovMcomper2 返回结果中的分组顺序"
    )
  }

  # ----------------------------------------------------------
  # 判断字母是否真正存在差异
  #
  # a vs b  -> 显著
  # a vs ab -> 不显著
  # ab vs b -> 不显著
  # ----------------------------------------------------------

  has_real_difference <- function(x) {

    x <- as.character(x)

    x <- x[
      !is.na(x) &
        nzchar(x)
    ]

    if (length(x) < 2) {
      return(FALSE)
    }

    if (length(unique(x)) == 1) {
      return(FALSE)
    }

    cmb <- combn(
      seq_along(x),
      2
    )

    for (i in seq_len(ncol(cmb))) {

      a <- strsplit(
        x[cmb[1, i]],
        ""
      )[[1]]

      b <- strsplit(
        x[cmb[2, i]],
        ""
      )[[1]]

      # 完全没有共同字母 -> 真正显著
      if (
        length(
          intersect(a, b)
        ) == 0
      ) {
        return(TRUE)
      }
    }

    FALSE
  }

  # ==========================================================
  # 2. 逐指标判断是否需要 fallback
  # ==========================================================

  stat_list <- list()
  gh_list <- list()
  letters_list <- list()

  for (m in valid_metrics) {

    cat(
      "\n========================================\n",
      "Alpha metric: ", m, "\n",
      "========================================\n",
      sep = ""
    )

    # --------------------------------------------------------
    # Primary 字母
    # --------------------------------------------------------

    primary_sig <- FALSE
    primary_letters <- NULL

    if (m %in% colnames(result_df)) {

      primary_letters <- result_df[[m]]

      primary_sig <- has_real_difference(
        primary_letters
      )
    }

    cat(
      "Primary significant:",
      primary_sig,
      "\n"
    )

    # --------------------------------------------------------
    # Welch ANOVA
    # --------------------------------------------------------

    tmp <- data.frame(
      group = data_use$group,
      value = data_use[[m]]
    )

    tmp <- tmp[
      complete.cases(tmp),
      ,
      drop = FALSE
    ]

    tmp$group <- droplevels(
      tmp$group
    )

    welch <- tryCatch(

      stats::oneway.test(
        value ~ group,
        data = tmp,
        var.equal = FALSE
      ),

      error = function(e) NULL
    )

    welch_p <- if (
      is.null(welch)
    ) {
      NA_real_
    } else {
      welch$p.value
    }

    # --------------------------------------------------------
    # Games-Howell
    # 每组至少2个有效重复才运行
    # --------------------------------------------------------

    group_n <- table(tmp$group)

    can_gh <-
      length(group_n) >= 2 &&
      min(group_n) >= 2

    if (can_gh) {

      gh <- tryCatch(

        rstatix::games_howell_test(
          tmp,
          value ~ group
        ),

        error = function(e) NULL
      )

    } else {

      gh <- NULL
    }

    if (!is.null(gh)) {

      gh$Metric <- m

      gh_list[[m]] <- gh

      gh_sig <- any(
        gh$p.adj < alpha,
        na.rm = TRUE
      )

      min_gh <- if (
        all(is.na(gh$p.adj))
      ) {
        NA_real_
      } else {
        min(
          gh$p.adj,
          na.rm = TRUE
        )
      }

    } else {

      gh_sig <- FALSE
      min_gh <- NA_real_
    }

    fallback_sig <-
      !is.na(welch_p) &&
      welch_p < alpha &&
      gh_sig

    # --------------------------------------------------------
    # 决定最终方法
    # --------------------------------------------------------

    method_used <- "MuiaovMcomper2"
    fallback_applied <- FALSE

    if (
      fallback &&
      !primary_sig &&
      fallback_sig &&
      !is.null(gh)
    ) {

      cat(
        "Fallback: Welch + Games-Howell\n"
      )

      pvec <- gh$p.adj

      names(pvec) <- paste(
        gh$group1,
        gh$group2,
        sep = "-"
      )

      letters <- multcompView::multcompLetters(
        pvec,
        threshold = alpha
      )$Letters

      letters_df <- data.frame(
        Metric = m,
        group = names(letters),
        Letter = unname(letters),
        stringsAsFactors = FALSE
      )

      # 按实验分组顺序整理
      letters_df <- letters_df[
        match(
          groups,
          letters_df$group
        ),
        ,
        drop = FALSE
      ]

      new_letters <- letters_df$Letter[
        match(
          result_groups,
          letters_df$group
        )
      ]

      # 替换最终 result 中该指标字母
      if (is.data.frame(result_final)) {

        result_final[[m]] <- new_letters

      } else if (is.matrix(result_final)) {

        result_final[, m] <- new_letters
      }

      # result_df 同步
      result_df[[m]] <- new_letters

      method_used <- "Welch + Games-Howell"
      fallback_applied <- TRUE

    } else {

      # 保留原方法字母
      if (!is.null(primary_letters)) {

        letters_df <- data.frame(
          Metric = m,
          group = result_groups,
          Letter = as.character(primary_letters),
          stringsAsFactors = FALSE
        )

      } else {

        letters_df <- data.frame(
          Metric = m,
          group = groups,
          Letter = NA_character_,
          stringsAsFactors = FALSE
        )
      }
    }

    letters_list[[m]] <- letters_df

    stat_list[[m]] <- data.frame(
      Metric = m,
      Primary_significant = primary_sig,
      Welch_P = welch_p,
      GamesHowell_min_Padj = min_gh,
      Fallback_significant = fallback_sig,
      Fallback_applied = fallback_applied,
      Method_used = method_used,
      stringsAsFactors = FALSE
    )
  }

  # ==========================================================
  # 3. 汇总
  # ==========================================================

  stat_summary <- dplyr::bind_rows(
    stat_list
  )

  # 同时保存跨 Alpha 指标的 BH-FDR，供查看
  stat_summary$Welch_FDR <- p.adjust(
    stat_summary$Welch_P,
    method = "BH"
  )

  games_howell <- dplyr::bind_rows(
    gh_list
  )

  final_letters <- dplyr::bind_rows(
    letters_list
  )

  invalid_table <- data.frame(
    Metric = invalid_metrics,
    Reason = "No finite data or no usable variation",
    stringsAsFactors = FALSE
  )

  # ==========================================================
  # 4. 返回
  # ==========================================================

  list(
    result = result_final,
    data = data_use,
    primary_result = result_primary,
    summary = stat_summary,
    games_howell = games_howell,
    letters = final_letters,
    invalid_metrics = invalid_table
  )
}
