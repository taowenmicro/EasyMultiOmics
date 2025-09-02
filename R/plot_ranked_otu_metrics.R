#' Ranked multi-metric panel plot for OTU (or family) features
#'
#' @description
#' Create a two-part panel plot that ranks features by sign (non-negative
#' first) and weighted score, shows the primary *Weighted* column on the left
#' (with value labels and OTU names), and displays additional metrics as
#' column facets on the right. Fill colors map to `log2FoldChange` with a
#' diverging palette centered at 0.
#'
#' @param data A data.frame or tibble containing the required columns.
#' @param col_otu Name of the feature ID column. Default `"OTU"`.
#' @param col_weighted Name of the weighted score column. Default `"Weighted_score"`.
#' @param col_mean Name of the mean abundance column. Default `"Mean_abundance"`.
#' @param col_prev Name of the prevalence column. Default `"Prevalence"`.
#' @param col_padj Name of the adjusted p-value (FDR) column. Default `"padj"`.
#' @param col_log2fc Name of the log2 fold change column. Default `"log2FoldChange"`.
#' @param col_importance Name of the importance score column. Default `"Importance"`.
#' @param col_degree Name of the network degree column. Default `"Degree"`.
#' @param col_betweenness Name of the betweenness column. Default `"Betweenness"`.
#' @param legend_title Title for the fill legend. Default
#'   `"log2FC(blue = negative, red = positive)"`.
#' @param palette_low,palette_mid,palette_high Colors for the diverging palette.
#'   Defaults are `#2166AC`, `"white"`, `#B2182B`.
#' @param metric_order Character vector controlling the facet order on the right;
#'   defaults to `c("Weighted","Mean_abundance","Prevalence","Significance",
#'   "log2FC","Importance","Degree","Betweenness")`.
#' @param base_size Base text size for `theme_minimal()`. Default `11`.
#' @param left_right_widths Numeric vector of length 2 for patchwork widths.
#'   Default `c(1, 6)`.
#' @param show_value_labels Logical; whether to print value labels on bars.
#'   Default `TRUE`.
#' @param min_p Smallest p used when computing `-log10(padj)` to avoid Inf.
#'   Default `1e-300`.
#'
#' @return A list with elements:
#'   \item{plot}{A patchwork object ready to print.}
#'   \item{df_ordered}{The ordered data.frame used for plotting (left panel order).}
#'   \item{df_long}{The long-format data used for plotting.}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' p_out <- plot_ranked_otu_metrics(df_top)
#' p_out$plot
#' }
plot_ranked_otu_metrics <- function(
    data,
    col_otu         = "OTU",
    col_weighted    = "Weighted_score",
    col_mean        = "Mean_abundance",
    col_prev        = "Prevalence",
    col_padj        = "padj",
    col_log2fc      = "log2FoldChange",
    col_importance  = "Importance",
    col_degree      = "Degree",
    col_betweenness = "Betweenness",
    legend_title    = "log2FC(blue = negative, red = positive)",
    palette_low     = "#2166AC",
    palette_mid     = "white",
    palette_high    = "#B2182B",
    metric_order    = c("Weighted","Mean_abundance","Prevalence",
                        "Significance","log2FC","Importance","Degree","Betweenness"),
    base_size       = 11,
    left_right_widths = c(1, 6),
    show_value_labels = TRUE,
    min_p           = 1e-300
) {
  # -- Dependencies
  requireNamespace("dplyr", quietly = TRUE)
  requireNamespace("tidyr", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("patchwork", quietly = TRUE)
  requireNamespace("scales", quietly = TRUE)

  # -- Column existence checks
  req_cols <- c(col_otu, col_weighted, col_mean, col_prev,
                col_padj, col_log2fc, col_importance, col_degree, col_betweenness)
  missing_cols <- setdiff(req_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  df <- data

  # -- Ordering: sign (>=0 first) then weighted score (desc)
  df <- dplyr::mutate(
    df,
    .sign_key = ifelse(df[[col_log2fc]] >= 0, 1L, 0L)
  )
  df_ord <- dplyr::arrange(df, dplyr::desc(.sign_key), dplyr::desc(.data[[col_weighted]]))
  otu_levels <- df_ord[[col_otu]]
  df_ord[[col_otu]] <- factor(df_ord[[col_otu]], levels = rev(otu_levels)) # top-most = most important

  # -- Build long table
  # Significance from padj with floor to avoid Inf
  padj_clean <- df_ord[[col_padj]]
  padj_clean <- ifelse(is.na(padj_clean), NA_real_, pmax(padj_clean, min_p))

  df_long <- dplyr::transmute(
    df_ord,
    OTU             = .data[[col_otu]],
    log2FoldChange  = .data[[col_log2fc]],
    Weighted        = .data[[col_weighted]],
    Mean_abundance  = .data[[col_mean]],
    Prevalence      = .data[[col_prev]],
    Significance    = -log10(padj_clean),
    log2FC          = .data[[col_log2fc]],
    Importance      = .data[[col_importance]],
    Degree          = .data[[col_degree]],
    Betweenness     = .data[[col_betweenness]]
  )

  df_long <- tidyr::pivot_longer(
    df_long,
    cols = -c(OTU, log2FoldChange),
    names_to = "metric",
    values_to = "score"
  )

  # -- Metric facet order
  df_long$metric <- factor(df_long$metric, levels = metric_order)

  # -- Value labels (metric-specific formats)
  fmt_sci <- scales::label_scientific(digits = 2)
  df_long <- dplyr::mutate(
    df_long,
    label = dplyr::case_when(
      metric == "Weighted"        ~ scales::number(score, accuracy = 0.001),
      metric == "Mean_abundance"  ~ fmt_sci(score),
      metric == "Prevalence"      ~ scales::percent(score, accuracy = 1),
      metric == "Significance"    ~ scales::number(score, accuracy = 0.1),
      metric == "log2FC"          ~ scales::number(score, accuracy = 0.01),
      metric == "Importance"      ~ scales::number(score, accuracy = 0.1),
      metric == "Degree"          ~ scales::number(score, accuracy = 1),
      metric == "Betweenness"     ~ scales::number(score, accuracy = 0.1),
      TRUE                        ~ as.character(round(score, 3))
    )
  )

  # -- Symmetric fill limits by |log2FC|
  fc_max <- max(abs(df_ord[[col_log2fc]]), na.rm = TRUE)

  # -- Left panel (Weighted with OTU names)
  p_left <- ggplot2::ggplot(
    dplyr::filter(df_long, .data$metric == "Weighted"),
    ggplot2::aes(x = .data$score, y = .data$OTU, fill = .data$log2FoldChange)
  ) +
    ggplot2::geom_col(width = 0.8) +
    {
      if (show_value_labels) ggplot2::geom_text(
        ggplot2::aes(label = .data$label), hjust = -0.1, size = 2.7
      )
    } +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::scale_fill_gradient2(
      low = palette_low, mid = palette_mid, high = palette_high,
      midpoint = 0, limits = c(-fc_max, fc_max),
      name = legend_title
    ) +
    ggplot2::labs(x = "Weighted score", y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      legend.position    = "right"
    )

  # -- Right panels (other metrics as facets)
  p_right <- ggplot2::ggplot(
    dplyr::filter(df_long, .data$metric != "Weighted"),
    ggplot2::aes(x = .data$score, y = .data$OTU, fill = .data$log2FoldChange)
  ) +
    ggplot2::geom_col(width = 0.8) +
    {
      if (show_value_labels) ggplot2::geom_text(
        ggplot2::aes(label = .data$label), hjust = -0.1, size = 2.5
      )
    } +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$metric),
                        scales = "free_x", switch = "x") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::scale_fill_gradient2(
      low = palette_low, mid = palette_mid, high = palette_high,
      midpoint = 0, limits = c(-fc_max, fc_max),
      name = legend_title
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      axis.text.y          = ggplot2::element_blank(),
      axis.ticks.y         = ggplot2::element_blank(),
      panel.grid.major.y   = ggplot2::element_blank(),
      panel.grid.minor     = ggplot2::element_blank(),
      strip.placement      = "outside",
      strip.text.x         = ggplot2::element_text(face = "bold"),
      legend.position      = "right"
    )

  p <- p_left + p_right + patchwork::plot_layout(
    widths = left_right_widths, guides = "collect"
  )

  out <- list(
    plot = p,
    df_ordered = df_ord,
    df_long = df_long
  )
  return(out)
}
