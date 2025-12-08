#' @title Clustered Bar Plots with Hierarchical Clustering for Microbial Data
#'
#' @description
#' Perform hierarchical clustering on microbial community data at both the sample
#' level and the group-mean level, and visualize results by combining dendrograms
#' with stacked bar plots at a specified taxonomic rank.
#'
#' @param dist Character string specifying the distance method (default: "bray").
#' @param otu,tax,map,tree Optional input tables. If NULL, values are extracted
#'   from \code{ps}.
#' @param j Character string specifying the taxonomic rank to visualize
#'   (e.g., "Phylum", "Genus"). Default is "Phylum".
#' @param ps A \code{phyloseq} object containing microbiome data.
#' @param rep Integer, number of replicates (unused, kept for compatibility).
#' @param Top Integer, number of top taxa to keep (by abundance). Others are
#'   merged into "Other". Default = 10.
#' @param tran Logical, if TRUE transform abundances to relative values.
#'   Default = TRUE.
#' @param hcluter_method Character string specifying hierarchical clustering
#'   method. Default = "complete".
#' @param Group Character string specifying the grouping column in sample metadata.
#'   Default = "Group".
#' @param cuttree Integer, number of clusters to cut the dendrogram into.
#'   Default = 3.
#'
#' @details
#' Steps performed:
#' \itemize{
#'   \item Transform OTU table to relative abundances (if \code{tran = TRUE}).
#'   \item Perform hierarchical clustering at sample level and group-average level.
#'   \item Collapse taxa to the specified rank, retain top \code{Top} taxa and merge others.
#'   \item Create stacked bar plots aligned with dendrograms using \code{facet_plot()}.
#' }
#'
#' @examples
#' \dontrun{
#' res <- cluMicro.bar.metm(
#'   dist = "bray", ps = ps.16s, j = "Genus", Top = 10,
#'   tran = TRUE, hcluter_method = "complete", Group = "Group",
#'   cuttree = length(unique(phyloseq::sample_data(ps.16s)[["Group"]]))
#' )
#' res[[1]]  # sample-level dendrogram
#' res[[2]]  # sample-level dendrogram + stacked barplot
#' res[[3]]  # group-level dendrogram
#' res[[4]]  # group-level dendrogram + stacked barplot
#' head(res[[5]])
#' }
#'
#' @author Tao Wen \email{2018203048@njau.edu.cn},
#'   Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @export
cluMicro.bar.metm <- function(
    dist = "bray",
    otu = NULL,
    tax = NULL,
    map = NULL,
    tree = NULL,
    j = "Phylum",
    ps = NULL,
    rep = 6,
    Top = 10,
    tran = TRUE,
    hcluter_method = "complete",
    Group = "Group",
    cuttree = 3,
    group_colors = NULL    # 新增：Group颜色参数
) {

  stopifnot(!is.null(ps))

  ps1_rela <- phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  dmat <- phyloseq::distance(ps1_rela, method = dist)
  hc <- stats::hclust(dmat, method = hcluter_method)
  sd <- data.frame(phyloseq::sample_data(ps1_rela))

  # 如果未提供group_colors，使用默认颜色
  if (is.null(group_colors)) {
    group_levels <- unique(sd[[Group]])
    default_colors <- c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF",
                        "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF")
    group_colors <- setNames(default_colors[seq_along(group_levels)], group_levels)
  }

  # 样本树图 - 添加 Group 颜色设置
  p <- ggtree::`%<+%`(ggtree::ggtree(hc), sd) +
    ggtree::geom_tippoint(size = 3.5, shape = 21,
                          ggplot2::aes(fill = .data[[Group]], x = x)) +
    ggtree::geom_tiplab(ggplot2::aes(color = .data[[Group]], x = x * 1.2),
                        hjust = 1, size = 3) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::theme_minimal()

  # 物种处理部分保持不变
  psdata <- ggClusterNet::tax_glom_wt(ps = ps1_rela, ranks = j)

  if (isTRUE(tran)) {
    psdata <- phyloseq::transform_sample_counts(psdata, function(x) x/sum(x))
  }

  otu_tab <- phyloseq::otu_table(psdata)
  if (!phyloseq::taxa_are_rows(psdata))
    otu_tab <- t(otu_tab)
  otu_mat <- as.matrix(otu_tab)

  top_taxa <- names(sort(rowSums(otu_mat), decreasing = TRUE))[seq_len(min(Top, nrow(otu_mat)))]

  tax_tab <- as.data.frame(phyloseq::tax_table(psdata))
  tax_tab[[j]][!(rownames(tax_tab) %in% top_taxa)] <- "Other"
  phyloseq::tax_table(psdata) <- as.matrix(tax_tab)

  Taxonomies <- phyloseq::psmelt(psdata)
  Taxonomies$Abundance <- Taxonomies$Abundance * 100
  Taxonomies$OTU <- NULL
  colnames(Taxonomies)[1] <- "id"

  if (!Group %in% colnames(Taxonomies)) {
    Taxonomies[[Group]] <- phyloseq::sample_data(psdata)[[Group]]
  }

  p <- p + ggnewscale::new_scale_fill()

  # 堆叠条形图 - 保持原始 Set2 调色板
  p1 <- ggtree::facet_plot(p, panel = "Stacked Barplot", data = Taxonomies,
                           geom = ggstance::geom_barh,
                           mapping = ggplot2::aes(x = Abundance, fill = .data[[j]]),
                           color = NA, stat = "identity") +
    ggplot2::scale_fill_brewer(palette = "Set2") +
    ggplot2::theme_minimal()

  # 组水平聚类
  otu_rel <- phyloseq::otu_table(ps1_rela)
  if (!phyloseq::taxa_are_rows(ps1_rela))
    otu_rel <- t(otu_rel)
  otu_rel <- as.matrix(otu_rel)

  smeta <- data.frame(phyloseq::sample_data(ps1_rela))
  grp_vec <- smeta[[Group]]

  otuG <- EasyMultiOmics:::avg_by_group_matrix(otu_rel, grp_vec)

  ps_group <- phyloseq::phyloseq(
    phyloseq::otu_table(otuG, taxa_are_rows = TRUE),
    phyloseq::tax_table(phyloseq::tax_table(ps1_rela))
  )

  hc_g <- stats::hclust(phyloseq::distance(ps_group, method = dist),
                        method = hcluter_method)
  clus_g <- stats::cutree(hc_g, k = cuttree)
  ddf <- data.frame(label = names(clus_g), member = factor(clus_g))
  rownames(ddf) <- ddf$label

  # 组树图 - 按 label (即Group名称) 设置颜色
  p3 <- ggtree::`%<+%`(ggtree::ggtree(hc_g), ddf) +
    ggtree::geom_tippoint(size = 3.5, shape = 21,
                          ggplot2::aes(fill = label, x = x)) +
    ggtree::geom_tiplab(ggplot2::aes(color = label, x = x * 1.2),
                        hjust = 1, size = 3) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::theme_minimal()

  grotax <- dplyr::summarise(
    dplyr::group_by(Taxonomies, .data[[Group]], .data[[j]]),
    Abundance = mean(Abundance),
    .groups = "drop"
  )

  p3 <- p3 + ggnewscale::new_scale_fill()

  # 组堆叠条形图 - 保持原始 Set1 调色板
  p4 <- ggtree::facet_plot(p3, panel = "Stacked Barplot", data = grotax,
                           geom = ggstance::geom_barh,
                           mapping = ggplot2::aes(x = Abundance, fill = .data[[j]]),
                           color = NA, stat = "identity") +
    ggplot2::scale_fill_brewer(palette = "Set1") +
    ggplot2::theme_minimal()

  return(list(p, p1, p3, p4, Taxonomies))
}


#' @keywords internal
#' @noRd
avg_by_group_matrix <- function(otu_mat, group_vector) {
  # otu_mat: taxa x samples (relative abundance)
  # 返回 taxa x groups (按组均值)
  stopifnot(is.matrix(otu_mat), length(group_vector) == ncol(otu_mat))
  groups <- unique(group_vector)
  res <- sapply(groups, function(g) {
    cols <- which(group_vector == g)
    if (length(cols) == 1) {
      otu_mat[, cols, drop = FALSE]
    } else {
      rowMeans(otu_mat[, cols, drop = FALSE])
    }
  })
  if (is.vector(res)) res <- matrix(res, ncol = length(groups))
  colnames(res) <- groups
  res
}

