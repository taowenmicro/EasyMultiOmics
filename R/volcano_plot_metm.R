#' Volcano plot with symmetric x-axis for edgeR/DE-style results
#'
#' @description
#' Draws a publication-grade volcano plot for differential-abundance results.
#' The x-axis is made symmetric around 0 using a robust quantile, and the
#' top labels from each direction (enriched/depleted) are placed with
#' collision-avoiding repelled text.
#'
#' @param x1 A `data.frame` with at least the columns `logFC`, `p`, and `level`.
#'   If a `Genus` column is present it will be used for labels; otherwise row names are used.
#' @param p_cut Significance threshold used to draw the horizontal reference line. Default `0.05`.
#' @param lfc_cut Absolute log2 fold-change threshold used to draw the vertical reference lines. Default `1`.
#' @param top_n Maximum number of labeled points per direction (enriched/depleted). Default `8`.
#' @param shade Logical; if `TRUE`, adds light shaded panels beyond \code{±lfc_cut}.
#'   Default `FALSE` (no background shading).
#' @param sym_q Quantile used to set the symmetric x-limits robustly (e.g. `0.995` = 99.5th percentile). Default `0.995`.
#' @param base_size Base font size passed to `ggplot2::theme_minimal()`. Default `12`.
#'
#' @details
#' - Expected `level` values are `"enriched"`, `"depleted"`, and `"nosig"`.
#' - If any p-values are zero (leading to `Inf` on `-log10(p)`), those points are capped
#'   to the maximum finite y-value + 1 to keep the plot finite.
#' - If a function named `theme_cell()` exists in the search path, it will be added on top
#'   of the theme (useful for lab-specific styling). Otherwise it is ignored.
#'
#' @return
#' Prints a `ggplot` object and (invisibly) returns a list with:
#' \describe{
#'   \item{plot}{The ggplot object.}
#'   \item{data}{Processed plotting data.}
#'   \item{data_sig}{The subset of labeled significant points.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose 'res' is a data.frame with columns: logFC, p, level, Genus
#' out <- volcano_plot_metm(res, p_cut = 0.05, lfc_cut = 1, top_n = 10)
#' out$plot
#' }
#'
#' @export
volcano_plot_metm <- function(
    x1,
    group = group,
    p_cut   = 0.05,
    lfc_cut = 1,
    top_n   = 8,
    shade   = FALSE,
    sym_q   = 0.995,
    base_size = 12
){
  stopifnot(is.data.frame(x1))
  req_cols <- c("logFC","p","level")
  miss <- setdiff(req_cols, colnames(x1))
  if (length(miss)) stop("x1 缺少必要列：", paste(miss, collapse=", "))

  df <- x1 %>%
    dplyr::mutate(
      y = -log10(.data$p),
      level = factor(.data$level, levels = c("depleted","nosig","enriched")),
      Genus = if ("Genus" %in% colnames(x1)) .data$Genus else rownames(x1)
    )

  # Handle p=0 (Inf on -log10)
  if (any(!is.finite(df$y))) {
    ymax_finite <- max(df$y[is.finite(df$y)], na.rm = TRUE)
    df$y[!is.finite(df$y)] <- ymax_finite + 1
  }

  # Select labels: top_n per direction by significance × effect size
  df_sig <- df %>%
    dplyr::filter(.data$level != "nosig") %>%
    dplyr::mutate(dir = dplyr::if_else(.data$logFC >= 0, "enriched", "depleted")) %>%
    dplyr::arrange(dplyr::desc(.data$y * abs(.data$logFC))) %>%
    dplyr::group_by(.data$dir) %>% dplyr::slice_head(n = top_n) %>% dplyr::ungroup()

  n_enriched <- sum(df$level == "enriched", na.rm = TRUE)
  n_depleted <- sum(df$level == "depleted", na.rm = TRUE)

  # Colors (NPG-like)
  col_enriched <- "#E64B35FF"  # red
  col_depleted <- "#4DBBD5FF"  # blue
  col_nosig    <- "grey78"

  # Symmetric x-axis (robust) & stable y-axis
  x_max <- max(stats::quantile(abs(df$logFC), sym_q, na.rm = TRUE), lfc_cut) * 1.05
  y_max <- max(df$y, na.rm = TRUE) * 1.05

  p <- ggplot2::ggplot()

  # Optional side shading
  if (isTRUE(shade)) {
    p <- p +
      ggplot2::annotate("rect", xmin = -Inf, xmax = -lfc_cut, ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(col_depleted, 0.05), color = NA) +
      ggplot2::annotate("rect", xmin =  lfc_cut, xmax =  Inf,     ymin = -Inf, ymax = Inf,
                        fill = scales::alpha(col_enriched, 0.05), color = NA)
  }

  p <- p +
    # Non-significant points: small and light
    ggplot2::geom_point(
      data = dplyr::filter(df, .data$level == "nosig"),
      ggplot2::aes(x = .data$logFC, y = .data$y),
      color = col_nosig, size = 1.4, alpha = 0.5, stroke = 0
    ) +
    # Significant points: colored by direction
    ggplot2::geom_point(
      data = dplyr::filter(df, .data$level != "nosig"),
      ggplot2::aes(x = .data$logFC, y = .data$y, color = .data$level),
      size = 2.2, alpha = 0.9, stroke = 0.2
    ) +
    # Threshold lines
    ggplot2::geom_hline(yintercept = -log10(p_cut), linetype = 2, linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = 3, linewidth = 0.4) +
    # Smart labels
    ggrepel::geom_text_repel(
      data = df_sig,
      ggplot2::aes(x = .data$logFC, y = .data$y, label = .data$Genus, color = .data$level),
      size = 3, min.segment.length = 0.1, box.padding = 0.35, point.padding = 0.25,
      max.overlaps = Inf, show.legend = FALSE
    ) +
    # Colors & legend
    ggplot2::scale_color_manual(
      values = c(depleted = col_depleted, nosig = col_nosig, enriched = col_enriched),
      breaks = c("depleted","enriched"),
      labels = c(
        paste0("Depleted (n=", n_depleted, ")"),
        paste0("Enriched (n=", n_enriched, ")")
      ),
      name = NULL
    ) +
    ggplot2::labs(x = "log2 fold change", y = expression(-log[10](p))) +
    ggplot2::coord_cartesian(xlim = c(-x_max, x_max), ylim = c(0, y_max)) +
    ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.line = element_blank(),
      legend.position = "right"
    ) +
    ggtitle(group)


  if (exists("theme_cell", mode = "function")) {
    p <- p + theme_cell()
  }

  out <- list(plot = p, data = df, data_sig = df_sig)
  print(p)
  invisible(out)
}
