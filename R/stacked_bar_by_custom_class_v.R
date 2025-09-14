#' @title Stacked Barplots by Custom Abundance Classes (Column/Row layout, Zoomed Medium/Low)
#'
#' @description
#' Aggregate a \code{phyloseq} object to a given taxonomic \code{rank}, convert to
#' relative abundance (percent), split taxa into High/Medium/Low classes by user
#' thresholds, select representative taxa per class, and draw Nature-style stacked
#' barplots. Panels can be arranged in a single column (default) or a single row.
#' Medium/Low panels use tighter auto-zoom on the abundance axis to highlight trends.
#' Legends are collected at the bottom; only the first panel keeps y-axis tick labels.
#'
#' @param ps A \code{phyloseq} object.
#' @param rank Character. Taxonomic rank to aggregate (e.g. "Phylum","Genus").
#' @param group Character. Sample metadata column for grouping (default "Group").
#' @param thresholds Numeric length-2 vector in percent, e.g. \code{c(1, 0.1)}:
#'   High >= max(thresholds); Medium in [min, max); Low < min. Order-insensitive.
#' @param n_each Integer length-3, number of representatives to keep in \code{c(High,Medium,Low)}.
#'   Default \code{c(5, 5, 5)}.
#' @param tran Logical; if TRUE, transform counts to per-sample relative abundances first. Default TRUE.
#' @param y_max_high,y_max_medium,y_max_low Numerics or NULL; panel x upper limits（bars are horizontal）.
#'   If NULL for Medium/Low, limits are tightly auto-estimated.
#' @param zoom_exclude_other Logical; when auto-zooming Medium/Low, exclude \code{other_label}
#'   from computing the max to better magnify selected taxa. Default TRUE.
#' @param other_label Character; label for non-selected taxa (default "Other").
#' @param arrange Character, \code{"column"} (3 panels stacked vertically) or \code{"row"} (3 panels in a row).
#'   Default \code{"column"}.
#' @param panel_heights Numeric length-3; relative heights for \code{c(High,Medium,Low)} when \code{arrange="column"}.
#'   Default \code{c(1.2, 1, 1)}.
#' @param panel_widths Numeric length-3; relative widths for \code{c(High,Medium,Low)} when \code{arrange="row"}.
#'   Default \code{c(1, 1.2, 1.2)} so Medium/Low get more space.
#'
#' @return A list with \code{high}, \code{medium}, \code{low}, \code{combined}, \code{table},
#'   and \code{y_limits} (actually the x-limits used by horizontal bars).
#' @export
stacked_bar_by_custom_class_v <- function(
    ps, rank = "Phylum", group = "Group",
    thresholds = c(1, 0.1), n_each = c(5, 5, 5),
    tran = TRUE,
    y_max_high = 100, y_max_medium = NULL, y_max_low = NULL,
    zoom_exclude_other = TRUE, other_label = "Other",
    arrange = c("column","row"),
    panel_heights = c(1.2, 1, 1),
    panel_widths  = c(1, 1.2, 1.2)
){
  arrange <- match.arg(arrange)
  if (is.null(ps)) stop("`ps` must be a phyloseq object.")
  tax_table_ps <- phyloseq::tax_table(ps)
  if (!rank %in% colnames(tax_table_ps)) stop("`rank` not found in taxonomy table.")

  ## 1) Relative abundance per sample
  ps_rela <- if (isTRUE(tran)) {
    phyloseq::transform_sample_counts(ps, function(x) x / sum(x))
  } else ps

  ## 2) Collapse to rank
  ps_rank <- ggClusterNet::tax_glom_wt(ps = ps_rela, ranks = rank)

  ## 3) Mean abundance (%) per taxon across samples
  otu_tab <- phyloseq::otu_table(ps_rank)
  if (!phyloseq::taxa_are_rows(ps_rank)) otu_tab <- t(otu_tab)
  otu_mat <- as.matrix(otu_tab)
  mean_abund_pct <- rowMeans(otu_mat) * 100

  ## 4) Class by thresholds (robust numeric)
  thr <- suppressWarnings(as.numeric(thresholds))
  if (length(thr) != 2 || anyNA(thr)) {
    stop("`thresholds` must be a numeric length-2 vector in percent, e.g. c(1, 0.1).")
  }
  hi_cut <- max(thr); lo_cut <- min(thr)

  taxa_names <- names(mean_abund_pct)
  class_vec  <- rep("Low", length(taxa_names))
  class_vec[mean_abund_pct >= hi_cut] <- "High"
  class_vec[mean_abund_pct < hi_cut & mean_abund_pct >= lo_cut] <- "Medium"
  names(class_vec) <- taxa_names

  ## 5) Pick representatives per class
  if (length(n_each) != 3) stop("`n_each` must be length 3: c(High, Medium, Low)")
  names(n_each) <- c("High","Medium","Low")
  pick_top <- function(class_name, n_keep){
    idx <- which(class_vec == class_name)
    if (!length(idx)) return(character(0))
    sel <- taxa_names[idx]
    sel_sorted <- sel[order(mean_abund_pct[sel], decreasing = TRUE)]
    utils::head(sel_sorted, n_keep)
  }
  keep_high   <- pick_top("High",   n_each["High"])
  keep_medium <- pick_top("Medium", n_each["Medium"])
  keep_low    <- pick_top("Low",    n_each["Low"])
  keep_all    <- unique(c(keep_high, keep_medium, keep_low))

  ## 6) Long table (percent) and ensure group column
  Tax <- phyloseq::psmelt(ps_rank)
  Tax$Abundance <- Tax$Abundance * 100
  if (!group %in% colnames(Tax)) {
    sam <- data.frame(phyloseq::sample_data(ps_rank))
    if (!group %in% colnames(sam)) stop("`group` column not found in sample_data(ps).")
    Tax[[group]] <- sam[match(Tax$Sample, rownames(sam)), group]
  }
  ## Non-representatives -> Other
  Tax[[rank]] <- ifelse(Tax[[rank]] %in% keep_all, as.character(Tax[[rank]]), other_label)

  ## Keep group order if factor
  smeta <- data.frame(phyloseq::sample_data(ps_rank))
  if (group %in% colnames(smeta) && is.factor(smeta[[group]])) {
    group_levels <- levels(smeta[[group]])
  } else {
    group_levels <- unique(as.character(Tax[[group]]))
  }

  ## 7) Group-wise mean composition (%)
  df_tmp <- data.frame(
    grp = Tax[[group]],
    tax = Tax[[rank]],
    Abundance = Tax$Abundance,
    stringsAsFactors = FALSE
  )
  comp_group <- stats::aggregate(Abundance ~ grp + tax, data = df_tmp, FUN = mean)
  names(comp_group) <- c(group, rank, "Abundance")

  ## 8) Per-class tables (representatives + Other)
  tbl_for <- function(kept_vec){
    if (!length(kept_vec)) return(NULL)
    comp_group[comp_group[[rank]] %in% c(kept_vec, other_label), , drop = FALSE]
  }
  tbl_high   <- tbl_for(keep_high)
  tbl_medium <- tbl_for(keep_medium)
  tbl_low    <- tbl_for(keep_low)

  ## 9) Tighter auto x-limits for Medium/Low (exclude Other if needed)
  auto_xlim <- function(tbl, fixed_max = NULL, class_name = c("High","Medium","Low")){
    class_name <- match.arg(class_name)
    if (!is.null(fixed_max)) return(c(0, fixed_max))
    if (is.null(tbl) || !nrow(tbl))  return(c(0, 100))
    tmp <- tbl
    if (isTRUE(zoom_exclude_other)) tmp <- tmp[tmp[[rank]] != other_label, , drop = FALSE]
    if (!nrow(tmp)) tmp <- tbl
    sums_df <- data.frame(
      grp = tmp[[group]],
      Abundance = tmp$Abundance,
      stringsAsFactors = FALSE
    )
    sums <- stats::aggregate(Abundance ~ grp, data = sums_df, FUN = sum)
    m <- max(sums$Abundance, na.rm = TRUE); if (!is.finite(m)) m <- 100
    headroom <- c(High = 0.05, Medium = 0.03, Low = 0.02)[[class_name]]
    step     <- c(High = 5,    Medium = 1,     Low = 0.5)[[class_name]]
    floorcap <- c(High = 5,    Medium = 0.5,   Low = 0.2)[[class_name]]
    upper <- m * (1 + headroom)
    if (step > 0) upper <- ceiling(upper / step) * step
    upper <- min(100, max(floorcap, upper))
    c(0, upper)
  }
  xlim_high   <- auto_xlim(tbl_high,   fixed_max = y_max_high,   class_name = "High")
  xlim_medium <- auto_xlim(tbl_medium, fixed_max = y_max_medium, class_name = "Medium")
  xlim_low    <- auto_xlim(tbl_low,    fixed_max = y_max_low,    class_name = "Low")

  ## 10) Titles (English) with threshold ranges
  fmt <- function(x) {
    x <- suppressWarnings(as.numeric(x))
    if (!is.finite(x)) return("?")
    if (x >= 10)       sprintf("%.1f", x)
    else if (x >= 1)   sprintf("%.2f", x)
    else if (x >= 0.1) sprintf("%.3f", x)
    else               sprintf("%.4f", x)
  }
  title_high   <- sprintf("High (≥ %s%%)",     fmt(hi_cut))
  title_medium <- sprintf("Medium (%s%%–<%s%%)", fmt(lo_cut), fmt(hi_cut))
  title_low    <- sprintf("Low (< %s%%)",      fmt(lo_cut))

  ## 11) Palettes (Nature-like)
  pal_high_base   <- RColorBrewer::brewer.pal(8, "Set1")
  pal_medium_base <- RColorBrewer::brewer.pal(8, "Set2")
  pal_low_base    <- RColorBrewer::brewer.pal(9, "Pastel1")

  ## 12) Plotting helper (horizontal bars: x=Abundance, y=group)
  make_plot <- function(tbl, class_title, pal_base, xlim_vec){
    if (is.null(tbl) || !nrow(tbl)) return(NULL)
    taxa_levels <- sort(unique(tbl[[rank]]))
    n_cols <- length(taxa_levels)
    pal <- grDevices::colorRampPalette(pal_base)(n_cols)
    ## set group factor order
    tbl[[group]] <- factor(tbl[[group]], levels = group_levels)

    ggplot2::ggplot(
      tbl, ggplot2::aes_string(y = group, x = "Abundance", fill = rank)
    ) +
      ggplot2::geom_col(position = "stack", color = NA) +
      ggplot2::scale_fill_manual(values = stats::setNames(pal, taxa_levels), name = rank) +
      ggplot2::scale_x_continuous(limits = xlim_vec, breaks = scales::pretty_breaks(n = 5)) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::labs(title = class_title, x = "Relative abundance (%)", y = NULL) +
      ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0))
  }

  p_high   <- make_plot(tbl_high,   title_high,   pal_high_base,   xlim_high)
  p_medium <- make_plot(tbl_medium, title_medium, pal_medium_base, xlim_medium)
  p_low    <- make_plot(tbl_low,    title_low,    pal_low_base,    xlim_low)

  ## 13) Hide y-axis text/ticks on medium & low
  if (!is.null(p_medium)) {
    p_medium <- p_medium +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }
  if (!is.null(p_low)) {
    p_low <- p_low +
      ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                     axis.ticks.y = ggplot2::element_blank())
  }

  ## 14) Combine; legend at bottom; column/row layout
  panels <- Filter(Negate(is.null), list(p_high, p_medium, p_low))
  if (!length(panels)) stop("No panels to plot (no taxa selected in any class).")

  if (arrange == "row") {
    combined <- patchwork::wrap_plots(panels, nrow = 1) +
      patchwork::plot_layout(
        widths = panel_widths[seq_along(panels)],
        guides = "collect"
      ) &
      ggplot2::theme(legend.position = "bottom")
  } else {
    combined <- patchwork::wrap_plots(panels, ncol = 1) +
      patchwork::plot_layout(
        heights = panel_heights[seq_along(panels)],
        guides  = "collect"
      ) &
      ggplot2::theme(legend.position = "bottom")
  }

  list(
    high = p_high, medium = p_medium, low = p_low,
    combined = combined, table = comp_group,
    y_limits = list(high = xlim_high, medium = xlim_medium, low = xlim_low)
  )
}
