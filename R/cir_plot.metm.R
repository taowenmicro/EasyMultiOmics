#' @title Generate Circular Plot for Microbial Composition
#'
#' @description
#' This function creates a circular plot to visualize microbial community composition based on OTU abundance.
#' It summarizes the microbial taxa at the specified rank, calculates the mean abundance within groups, and displays the data in a circular chord diagram.
#'
#' @param ps A `phyloseq` object containing microbiome data.
#' @param Top An integer specifying the number of top taxa to display. Taxa not in the top are grouped as `"others"`. Default is `10`.
#' @param rank A numeric or character value specifying the taxonomic rank for aggregation (e.g., `"Phylum"`, `"Genus"`, etc.). Default is `7`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{mer_otu_mean}{A matrix summarizing the mean abundance of OTUs across groups.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Transforms the data into relative abundances using `phyloseq::transform_sample_counts`.
#'   \item Aggregates microbial taxa at the specified rank using `ggClusterNet::tax_glom_wt`.
#'   \item Calculates the mean abundance of OTUs for each group and selects the top `Top` taxa.
#'   \item Groups less abundant taxa into a category labeled `"others"`.
#'   \item Visualizes the relationships between microbial taxa and groups using a circular chord diagram with the `circlize` package.
#' }
#'
#' @examples
#' \dontrun{
#' res = cir_plot.metm(ps  = ps.16s,Top = 12,rank = 6)
#'
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export


cir_plot.metm=function (ps = ps, Top = 10, rank = 7)
{
  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  ps_rela
  ps_P <- ps_rela %>% ggClusterNet::tax_glom_wt(rank = rank)
  ps_P
  otu_P = as.data.frame((ggClusterNet::vegan_otu(ps_P)))
  head(otu_P)
  tax_P = as.data.frame(ggClusterNet::vegan_tax(ps_P))
  sub_design <- as.data.frame(phyloseq::sample_data(ps_P))
  count2 = otu_P
  iris.split <- split(count2, as.factor(sub_design$Group))
  iris.apply <- lapply(iris.split, function(x) colSums(x[]))
  iris.combine <- do.call(rbind, iris.apply)
  ven2 = t(iris.combine)
  if (is.numeric(rank)) {
    lev = phyloseq::rank_names(ps)[rank]
  }
  else {
    lev = rank
  }
  Taxonomies <- ps %>% ggClusterNet::tax_glom_wt(rank = rank) %>%
    phyloseq::transform_sample_counts(function(x) {
      x/sum(x)
    }) %>% phyloseq::psmelt() %>% dplyr::arrange(!!sym(lev))
  iris_groups <- dplyr::group_by(Taxonomies, !!sym(lev))
  ps0_sum <- dplyr::summarise(iris_groups, mean(Abundance),
                              sd(Abundance))
  ps0_sum[is.na(ps0_sum)] <- 0
  head(ps0_sum)
  colnames(ps0_sum) = c("ID", "mean", "sd")
  ps0_sum <- dplyr::arrange(ps0_sum, desc(mean))
  ps0_sum$mean <- ps0_sum$mean * 100
  ps0_sum <- as.data.frame(ps0_sum)
  head(ps0_sum)
  top_P = ps0_sum$ID[1:Top]
  top_P
  otu_P = as.data.frame(t(otu_P))
  otu_tax = merge(ven2, tax_P, by = "row.names", all = F)
  dim(otu_tax)
  otu_tax[, lev] = as.character(otu_tax[, lev])
  otu_tax[, lev][is.na(otu_tax[, lev])] = "others"
  i = 1
  for (i in 1:nrow(otu_tax)) {
    if (otu_tax[, lev][i] %in% top_P) {
      otu_tax[, lev][i] = otu_tax[, lev][i]
    }
    else if (!otu_tax[, lev][i] %in% top_P) {
      otu_tax[, lev][i] = "others"
    }
  }
  otu_tax[, lev] = as.factor(otu_tax[, lev])
  head(otu_tax)
  otu_mean = otu_tax[as.character(unique(sub_design$Group))]
  head(otu_mean)
  row.names(otu_mean) = row.names(otu_tax)
  iris.split <- split(otu_mean, as.factor(otu_tax[, lev]))
  iris.apply <- lapply(iris.split, function(x) colSums(x[]))
  iris.combine <- do.call(rbind, iris.apply)
  mer_otu_mean = t(iris.combine)
  head(mer_otu_mean)
  mi_sam = RColorBrewer::brewer.pal(9, "Set1")
  mi_tax = colorRampPalette(RColorBrewer::brewer.pal(9, "Set3"))(length(row.names(mer_otu_mean)))
  grid.col = NULL
  grid.col[as.character(unique(sub_design$Group))] = mi_sam
  grid.col[row.names(mer_otu_mean)] = mi_tax
  circlize::circos.par(gap.degree = c(rep(2, nrow(mer_otu_mean) -
                                            1), 10, rep(2, ncol(mer_otu_mean) - 1), 10), start.degree = 180)
  circlize::chordDiagram(mer_otu_mean, directional = F, diffHeight = 0.06,
                         grid.col = grid.col, reduce = 0, transparency = 0.5,
                         annotationTrack = c("grid", "axis"), preAllocateTracks = 2)
  circlize::circos.track(track.index = 1, panel.fun = function(x,
                                                               y) {
    circlize::circos.text(circlize::CELL_META$xcenter, circlize::CELL_META$ylim[1],
                          circlize::CELL_META$sector.index, facing = "clockwise",
                          niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA)
  circlize::circos.clear()
  return(list(mer_otu_mean))
}
