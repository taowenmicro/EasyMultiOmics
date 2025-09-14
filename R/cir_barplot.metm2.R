#' @title Circular dendrogram + tight three concentric stacked-bar rings (separate palettes & legends)
#'
#' @description
#' Add **three** stacked-bar rings (High / Medium / Low) outside a circular hierarchical
#' clustering dendrogram. To create a near “no-gap” appearance by default, ring axes are
#' turned off and both the inter-ring gap and axis padding are set to 0, while the
#' tree-to-inner-ring distance (`base_offset`) is kept small.
#' Each ring uses its **own** fill scale via \strong{ggnewscale::new_scale_fill} and a
#' dedicated \strong{ggplot2::scale_fill_manual}, yielding **three different palettes and
#' three separate legends** (High / Medium / Low).
#'
#' @param ps A \code{phyloseq} object.
#' @param Top Integer. Number of globally top taxa to keep for binning/plotting. Default: 15.
#' @param dist Character. Distance metric passed to \code{phyloseq::distance}. Default: "bray".
#' @param cuttree Integer. Number of clusters to cut the dendrogram into. Default: 3.
#' @param hcluter_method Character. Hierarchical clustering method. Default: "complete".
#' @param rank Character. Taxonomic rank for aggregation (e.g., "Phylum", "Genus"); must exist in \code{tax_table}.
#' @param thresholds Numeric length-2 (percent), e.g., \code{c(1, 0.1)} meaning
#'   High (≥ max), Medium [min, max), Low (< min). Order-insensitive.
#' @param n_each Integer length-3. Number of representative taxa for \code{c(High, Medium, Low)}.
#'   Default \code{c(10, 8, 8)}.
#' @param xlim_high,xlim_medium,xlim_low Numeric length-2 or \code{NULL}. Display ranges (in %) for the three rings.
#'   If \code{NULL}, ranges are automatically determined.
#' @param breaks_high,breaks_medium,breaks_low Numeric vectors or \code{NULL}. Tick breaks for the three rings.
#'   If \code{NULL}, breaks are auto-calculated from ranges.
#' @param clip Logical. Whether to truncate values outside x-limits (display only). Default: \code{TRUE}.
#' @param ring_pwidth Numeric length-3. Radial width of rings \code{c(High, Medium, Low)}. Default \code{c(1.2, 1.1, 1.1)}.
#' @param ring_bar_width Numeric length-3. Bar width in each ring. Default \code{c(0.5, 0.5, 0.5)}.
#' @param base_offset Numeric. Distance from the tree to the inner ring. Default: 0.9.
#' @param ring_gap Numeric. Base gap between adjacent rings. Default: 0.
#' @param axis_pad Numeric. Extra padding if a ring shows its axis. Default: 0.
#' @param show_axes Logical length-3. Whether to show abundance (x) axis for each ring. Default \code{c(FALSE, FALSE, FALSE)}.
#' @param legend_position Character. Legend position for \code{ggplot2}. Default: "right".
#'
#' @return A list with:
#' \describe{
#'   \item{plot}{The circular dendrogram with three stacked-bar rings (High → Medium → Low),
#'               using three distinct palettes and three separate legends.}
#'   \item{data_high,data_medium,data_low}{Data frames used for the High/Medium/Low rings.}
#'   \item{selected_taxa}{A list of representative taxa: \code{list(high, medium, low)}.}
#'   \item{cuts}{A list with thresholds: \code{list(hi_cut, lo_cut)}.}
#' }
#'
#' @examples
#' \dontrun{
#' library(phyloseq); data(GlobalPatterns)
#' ps_demo <- prune_samples(sample_sums(GlobalPatterns) > 0, GlobalPatterns)
#' set.seed(1); ps_demo <- prune_samples(sample(sample_names(ps_demo), 20), ps_demo)
#'
#' res <- cir_barplot.metm2(
#'   ps = ps_demo, Top = 20, rank = "Phylum",
#'   thresholds = c(1, 0.1), n_each = c(10, 8, 8),
#'   xlim_high = c(0,100), xlim_medium = c(0,20), xlim_low = c(0,5),
#'   breaks_high = seq(0,100,10), breaks_medium = seq(0,20,5), breaks_low = seq(0,5,1),
#'   ring_pwidth = c(1.2,1.1,1.1),
#'   base_offset = 0.9, ring_gap = 0, axis_pad = 0,
#'   show_axes = c(FALSE, FALSE, FALSE),
#'   legend_position = "right"
#' )
#' print(res$plot)
#' }
#'
#' @author Tao Wen, Peng-Hao Xie
#' @export
cir_barplot.metm2 <- function(
    ps,
    Top = 15,
    dist = "bray",
    cuttree = 3,
    hcluter_method = "complete",
    rank = "Phylum",
    thresholds = c(1, 0.1),
    n_each = c(10, 8, 8),
    xlim_high = NULL, xlim_medium = NULL, xlim_low = NULL,
    breaks_high = NULL, breaks_medium = NULL, breaks_low = NULL,
    clip = TRUE,
    ring_pwidth = c(1.2, 1.1, 1.1),
    ring_bar_width = c(0.5, 0.5, 0.5),
    base_offset = 0.9,
    ring_gap = 0,
    axis_pad = 0,
    show_axes = c(FALSE, FALSE, FALSE),
    legend_position = "right"
){
  if (missing(ps) || is.null(ps)) stop("`ps` must be a phyloseq object.")
  if (!inherits(ps, "phyloseq")) stop("`ps` must be a phyloseq object.")
  tax_tab0 <- phyloseq::tax_table(ps)
  if (!rank %in% colnames(tax_tab0)) stop("`rank` not found in taxonomy table.")
  if (length(n_each) != 3) stop("`n_each` must be length 3: c(High,Medium,Low).")
  if (length(ring_pwidth) != 3) stop("`ring_pwidth` must be length 3.")
  if (length(ring_bar_width) != 3) stop("`ring_bar_width` must be length 3.")
  if (length(show_axes) != 3) stop("`show_axes` must be length 3 (logicals).")

  ## helpers
  .pal_npg <- function(n){
    if (requireNamespace("ggsci", quietly = TRUE)) {
      base <- ggsci::pal_npg("nrc")(10)
    } else {
      base <- c("#E64B35FF","#4DBBD5FF","#00A087FF","#3C5488FF","#F39B7FFF",
                "#8491B4FF","#91D1C2FF","#DC0000FF","#7E6148FF","#B09C85FF")
    }
    if (n <= length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
  }
  .pal_alt <- function(n, which = c("lancet","nejm","jama")){
    which <- match.arg(which)
    if (requireNamespace("ggsci", quietly = TRUE)) {
      fun <- switch(which,
                    lancet = ggsci::pal_lancet(),
                    nejm   = ggsci::pal_nejm(),
                    jama   = ggsci::pal_jama())
      fun(n)
    } else {
      pal <- switch(which,
                    lancet = RColorBrewer::brewer.pal(8,"Set1"),
                    nejm   = RColorBrewer::brewer.pal(8,"Set2"),
                    jama   = RColorBrewer::brewer.pal(9,"Pastel1"))
      if (n <= length(pal)) pal[seq_len(n)] else grDevices::colorRampPalette(pal)(n)
    }
  }
  .fmt <- function(x){ sub("\\.?0+$","", formatC(x, format = "f", digits = 3)) }
  .auto_xlim <- function(dat){
    if (is.null(dat) || !nrow(dat)) return(c(0, 1))
    m <- max(dat$Abundance, na.rm = TRUE); if (!is.finite(m)) m <- 1
    hi <- min(100, max(5, ceiling(m * 1.05 / 5) * 5))
    c(0, hi)
  }
  .auto_breaks <- function(xlim){
    lo <- xlim[1]; hi <- xlim[2]; span <- hi - lo
    step <- if (span <= 10) 1 else if (span <= 20) 2 else if (span <= 50) 5 else 10
    seq(lo, hi, by = step)
  }

  ## robust sample_data
  sn <- phyloseq::sample_names(ps)
  sdat <- phyloseq::sample_data(ps)
  if (is.null(sdat) || nrow(sdat) == 0 || ncol(sdat) == 0) {
    sd <- data.frame(Group = factor(rep("All", length(sn)), levels = "All"),
                     row.names = sn, check.names = FALSE)
  } else {
    sd <- as.data.frame(sdat, stringsAsFactors = FALSE)
    if (is.null(rownames(sd)) || !all(sn %in% rownames(sd))) {
      rownames(sd) <- sn
    } else {
      sd <- sd[sn, , drop = FALSE]
    }
    if (!"Group" %in% colnames(sd)) sd$Group <- factor("All", levels = "All")
  }

  ## tree
  ps_scaled <- ggClusterNet::scale_micro(ps)
  unif <- phyloseq::distance(ps_scaled, method = dist)
  hc   <- stats::hclust(unif, method = hcluter_method)
  clus <- stats::cutree(hc, cuttree)
  d  <- data.frame(label = names(clus), member = factor(clus))
  dd <- data.frame(label = d$label, Group = sd[d$label, "Group", drop = TRUE])

  p_tree <- ggtree::ggtree(hc, layout = "circular") %<+% dd +
    ggplot2::geom_point(size = 3.2, shape = 21, ggplot2::aes(fill = .data$Group, x = x)) +
    ggtree::geom_tiplab(ggplot2::aes(color = .data$Group, x = x * 1.2), hjust = 1, offset = 0.25) +
    ggplot2::xlim(-0.5, NA)

  ## rank aggregation + relative abundance + Top prefilter
  ps_rank <- ggClusterNet::tax_glom_wt(ps = ps_scaled, ranks = rank)
  ps_rank <- phyloseq::transform_sample_counts(ps_rank, function(x) x / sum(x))

  otu <- phyloseq::otu_table(ps_rank)
  if (!phyloseq::taxa_are_rows(ps_rank)) otu <- t(otu)
  otu_mat <- as.matrix(otu)
  top_ids <- names(sort(rowSums(otu_mat), decreasing = TRUE))[seq_len(min(Top, nrow(otu_mat)))]

  tax_long <- phyloseq::psmelt(ps_rank)
  tax_long$Abundance <- tax_long$Abundance * 100
  tax_long$id <- tax_long$Sample
  tax_long <- tax_long[tax_long[[rank]] %in% top_ids, , drop = FALSE]

  ## classing by global mean (%)
  mean_pct <- tapply(tax_long$Abundance, tax_long[[rank]], mean, na.rm = TRUE)
  thr <- suppressWarnings(as.numeric(thresholds))
  if (length(thr) != 2 || anyNA(thr)) stop("`thresholds` must be numeric length-2, e.g. c(1,0.1).")
  hi_cut <- max(thr); lo_cut <- min(thr)

  taxa_all <- names(mean_pct)
  class_vec <- rep("Low", length(taxa_all))
  class_vec[mean_pct >= hi_cut] <- "High"
  class_vec[mean_pct < hi_cut & mean_pct >= lo_cut] <- "Medium"
  names(class_vec) <- taxa_all

  pick_top <- function(class_name, n_keep){
    sel <- names(which(class_vec == class_name))
    if (!length(sel)) return(character(0))
    sel_sorted <- sel[order(mean_pct[sel], decreasing = TRUE)]
    utils::head(sel_sorted, n_keep)
  }
  keep_high   <- pick_top("High",   n_each[1])
  keep_medium <- pick_top("Medium", n_each[2])
  keep_low    <- pick_top("Low",    n_each[3])

  ## ring data
  dat_high   <- tax_long[tax_long[[rank]] %in% keep_high,   c("id","Abundance", rank)]
  dat_medium <- tax_long[tax_long[[rank]] %in% keep_medium, c("id","Abundance", rank)]
  dat_low    <- tax_long[tax_long[[rank]] %in% keep_low,    c("id","Abundance", rank)]
  names(dat_high)   <- c("id","Abundance","RankName")
  names(dat_medium) <- c("id","Abundance","RankName")
  names(dat_low)    <- c("id","Abundance","RankName")

  ## ranges / breaks / clipping
  prep_ring <- function(dat, xlim, breaks){
    if (is.null(dat) || !nrow(dat)) return(list(dat=NULL, xvar=NULL, lim=NULL, br=NULL))
    if (is.null(xlim)) xlim <- .auto_xlim(dat)
    if (is.null(breaks)) breaks <- .auto_breaks(xlim)
    if (isTRUE(clip)) {
      dat$Abundance_plot <- pmin(pmax(dat$Abundance, xlim[1]), xlim[2])
      xvar <- "Abundance_plot"
    } else xvar <- "Abundance"
    list(dat=dat, xvar=xvar, lim=xlim, br=breaks)
  }
  Rhi <- prep_ring(dat_high,   xlim_high,   breaks_high)
  Rme <- prep_ring(dat_medium, xlim_medium, breaks_medium)
  Rlo <- prep_ring(dat_low,    xlim_low,    breaks_low)

  ## offsets（tight; no gaps, no axis padding by default）
  off_high <- base_offset
  # off_med  <- base_offset + ring_pwidth[1] + (if (isTRUE(show_axes[1])) axis_pad else 0) + ring_gap
  # off_low  <- off_med     + ring_pwidth[2] + (if (isTRUE(show_axes[2])) axis_pad else 0) + ring_gap
  off_med  = 0.5
  off_low  = 0.2

  ## palettes (per ring): High=NPG, Medium=Lancet, Low=NEJM (fallback to Brewer if ggsci missing)
  pal_high <- {
    n <- if (is.null(Rhi$dat)) 0 else length(unique(Rhi$dat$RankName))
    if (n) .pal_npg(n) else character(0)
  }
  pal_med <- {
    n <- if (is.null(Rme$dat)) 0 else length(unique(Rme$dat$RankName))
    if (n) .pal_alt(n, "lancet") else character(0)
  }
  pal_low <- {
    n <- if (is.null(Rlo$dat)) 0 else length(unique(Rlo$dat$RankName))
    if (n) .pal_alt(n, "nejm") else character(0)
  }

  ## legend titles (different titles so legends stay separate)
  title_high <- paste0(rank, " (High ≥ ", .fmt(hi_cut), "%)")
  title_med  <- paste0(rank, " (Medium ", .fmt(lo_cut), "–", .fmt(hi_cut), "%)")
  title_low  <- paste0(rank, " (Low < ", .fmt(lo_cut), "%)")

  ## assemble: each ring gets its own scale (new_scale_fill + scale_fill_manual)
  add_ring <- function(p, R, pwidth, bwidth, off, show_axis){
    if (is.null(R$dat) || !nrow(R$dat)) return(p)
    ap <- if (isTRUE(show_axis)) {
      list(axis = "x", breaks = R$br, text.angle = -45, hjust = 0, vjust = 0.5)
    } else {
      list(axis = "none")
    }
    p + ggtreeExtra::geom_fruit(
      inherit.aes = FALSE,
      data  = R$dat,
      geom  = geom_col,  # uses @importFrom ggplot2 geom_col
      mapping = ggplot2::aes_string(x = R$xvar, y = "id", fill = "RankName"),
      width = bwidth,
      orientation = "y",
      offset = off,
      pwidth = pwidth,
      axis.params = ap
    )
  }

  p_out <- p_tree

  ## High ring
  if (!is.null(Rhi$dat) && nrow(Rhi$dat)) {
    p_out <- p_out + ggnewscale::new_scale_fill()
    p_out <- add_ring(p_out, Rhi, ring_pwidth[1], ring_bar_width[1], off_high, show_axes[1])
    p_out <- p_out + ggplot2::scale_fill_manual(values = pal_high, name = title_high)
  }

  ## Medium ring
  if (!is.null(Rme$dat) && nrow(Rme$dat)) {
    p_out <- p_out + ggnewscale::new_scale_fill()
    p_out <- add_ring(p_out, Rme, ring_pwidth[2], ring_bar_width[2], off_med, show_axes[2])
    p_out <- p_out + ggplot2::scale_fill_manual(values = pal_med, name = title_med)
  }

  ## Low ring
  if (!is.null(Rlo$dat) && nrow(Rlo$dat)) {
    p_out <- p_out + ggnewscale::new_scale_fill()
    p_out <- add_ring(p_out, Rlo, ring_pwidth[3], ring_bar_width[3], off_low, show_axes[3])
    p_out <- p_out + ggplot2::scale_fill_manual(values = pal_low, name = title_low)
  }

  p_out <- p_out +
    ggplot2::theme(
      legend.position = legend_position,
      legend.box = "vertical"
    )

  list(
    plot = p_out,
    data_high = Rhi$dat, data_medium = Rme$dat, data_low = Rlo$dat,
    selected_taxa = list(high = keep_high, medium = keep_medium, low = keep_low),
    cuts = list(hi_cut = hi_cut, lo_cut = lo_cut)
  )
}
