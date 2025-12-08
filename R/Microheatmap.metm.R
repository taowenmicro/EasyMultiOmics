#' @title Generate a Heatmap for Microbial Abundance Data
#'
#' @description
#' This function generates a heatmap to visualize microbial abundance across samples. The heatmap includes options for clustering, scaling, and annotations of rows and columns. The data is sourced from a `phyloseq` object.
#'
#' @param ps_rela A `phyloseq` object containing relative abundance data.
#' @param id A vector specifying the taxa or features to include in the heatmap.
#' @param label Logical. If `TRUE`, displays sample labels as a color bar. Default is `TRUE`.
#' @param col_cluster Logical. If `TRUE`, clusters columns (samples). Default is `TRUE`.
#' @param row_cluster Logical. If `TRUE`, clusters rows (taxa or features). Default is `TRUE`.
#' @param ord.col Logical. If `TRUE`, orders columns based on `axis_order.s`. Default is `TRUE`.
#' @param scale Logical. If `TRUE`, scales abundance data row-wise to mitigate the influence of high-abundance taxa. Default is `TRUE`.
#' @param axis_order.s A vector specifying the order of sample IDs for column arrangement. Default is `axis_order.s`.
#' @param row.lab A character string specifying a taxonomic rank (e.g., `"Phylum"`, `"Genus"`) to annotate rows.
#' @param col1 A vector of colors for the heatmap gradient. Default is `ggsci::pal_gsea(alpha = 1)(12)`.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{Heatmap Tile Plot}{A heatmap displaying abundance data as a tile plot.}
#'   \item{Bubble Heatmap}{A heatmap displaying abundance data as a bubble plot.}
#'   \item{Processed Data}{A data frame containing the processed data used for the heatmap.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts and processes abundance and taxonomic data for the specified taxa (`id`).
#'   \item Scales the data row-wise if `scale = TRUE` to normalize abundances.
#'   \item Constructs hierarchical clusters for rows and/or columns based on user input.
#'   \item Optionally annotates rows with taxonomic information and columns with sample group information.
#'   \item Generates two types of heatmaps:
#'     \itemize{
#'       \item A tile heatmap (p1).
#'       \item A bubble heatmap (p2).
#'     }
#'   \item Adds clustering dendrograms and side plots (e.g., row mean bar plots) as optional annotations.
#' }
#'
#' This function provides a highly customizable way to visualize microbial abundance data with integrated clustering and annotations.
#'
#' @examples
#' \dontrun{
#' result <- Microheatmap.metm(ps_rela = ps_tem, id = id, col_cluster = FALSE, row_cluster = FALSE)
#' p24.1 <- result[[1]]
#' p24.2 <- result[[2]]
#' dat <- result[[3]]
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#'
#' @export

Microheatmap.metm <- function(
    ps_rela,
    id,
    label = TRUE,
    col_cluster = TRUE,
    row_cluster = TRUE,
    ord.col = TRUE,
    scale = TRUE,
    axis_order.s = NULL,
    row.lab = NULL,
    col1 = (ggsci::pal_gsea(alpha = 1))(12),
    y_text_size = 8
) {

  map = phyloseq::sample_data(ps_rela)
  map$ID = row.names(map)
  phyloseq::sample_data(ps_rela) = map

  otu = as.data.frame(t(ggClusterNet::vegan_otu(ps_rela)))
  otu = as.matrix(otu[id, ])

  ps_heatm = phyloseq::phyloseq(
    phyloseq::otu_table(otu, taxa_are_rows = TRUE),
    phyloseq::tax_table(ps_rela),
    phyloseq::sample_data(ps_rela)
  )

  datah <- as.data.frame(t(ggClusterNet::vegan_otu(ps_heatm)))
  tax = as.data.frame(ggClusterNet::vegan_tax(ps_heatm))
  otutaxh = cbind(datah, tax)
  otutaxh$id = paste(row.names(otutaxh))
  row.names(otutaxh) = otutaxh$id

  data <- otutaxh[, c("id", phyloseq::sample_names(ps_rela))]

  rig <- data[, phyloseq::sample_names(ps_rela)] %>%
    rowMeans() %>%
    as.data.frame()
  colnames(rig) = "MeanAbundance"
  rig$id = row.names(rig)
  rig = rig %>% dplyr::arrange(MeanAbundance)
  rig$id = factor(rig$id, levels = rig$id)

  p_rig = ggplot(rig) +
    geom_bar(aes(y = id, x = MeanAbundance), fill = "#A54657", stat = "identity") +
    theme_void()

  tem = data[, phyloseq::sample_names(ps_rela)] %>% as.matrix()
  if (scale == TRUE) {
    tem = scale(t(tem)) %>% t() %>% as.data.frame()
  } else {
    tem = tem
  }
  data[, phyloseq::sample_names(ps_rela)] = tem

  mat <- data[, -1]

  if (row_cluster == TRUE) {
    clust <- hclust(dist(mat %>% as.matrix()))
    ggtree_plot <- ggtree::ggtree(clust)
  }

  if (col_cluster == TRUE) {
    v_clust <- hclust(dist(mat %>% as.matrix() %>% t()))
    ggtree_plot_col <- ggtree::ggtree(v_clust) + ggtree::layout_dendrogram()
  }

  if (label == TRUE) {
    map_label = data.frame(phyloseq::sample_data(ps_rela), check.names = FALSE)
    map_label$ID = row.names(map_label)
    labels = ggplot(map_label, aes(x = ID, y = 1, fill = Group)) +
      geom_tile() +
      scale_fill_brewer(palette = "Set1", name = "Cell Type") +
      theme_void()
  }

  if (!is.null(row.lab)) {
    tax = ps_rela %>% ggClusterNet::vegan_tax() %>% as.data.frame()
    tax$ID = row.names(tax)
    p.row.lab = ggplot(tax, aes(y = ID, x = 1, fill = !!sym(row.lab))) +
      geom_tile() +
      scale_fill_brewer(palette = "Set3", name = "Cell Type") +
      theme_void()
  }

  # 修复：正确转换 sample_data 为 data.frame
  map_order = data.frame(phyloseq::sample_data(ps_rela), check.names = FALSE)
  map_order$ID = row.names(map_order)
  map_order = map_order %>% dplyr::arrange(Group)

  pcm = reshape2::melt(data, id = c("id"))
  pcm$variable = factor(pcm$variable, levels = map_order$ID)

  if (ord.col == TRUE && !is.null(axis_order.s)) {
    pcm$variable = factor(pcm$variable, levels = axis_order.s)
  }

  pcm$id = factor(pcm$id, levels = rig$id)

  # p1 - 修改 y 轴字号
  p1 = ggplot(pcm, aes(y = id, x = variable)) +
    geom_tile(aes(size = value, fill = value)) +
    labs(y = "", x = "", size = "Relative Abundance (%)", fill = "") +
    scale_x_discrete(limits = rev(levels(pcm$variable))) +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours = col1) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = y_text_size),
      axis.text.x = element_text(colour = "black", angle = 90)
    )

  # p2 - 添加 y 轴字号设置
  p2 = ggplot(pcm, aes(y = id, x = variable)) +
    geom_point(aes(size = value, fill = value), alpha = 0.75, shape = 21) +
    scale_size_continuous(limits = c(0.001, 1000), range = c(2, 25), breaks = c(0.1, 0.5, 1)) +
    labs(y = "", x = "", size = "Relative Abundance (%)", fill = "") +
    scale_x_discrete(limits = rev(levels(pcm$variable))) +
    scale_y_discrete(position = "right") +
    scale_fill_gradientn(colours = col1) +
    theme(
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.text.y = element_text(size = y_text_size),
      axis.text.x = element_text(colour = "black", angle = 90)
    )

  p1 <- p1 %>% aplot::insert_right(p_rig, width = 0.2)
  p2 <- p2 %>% aplot::insert_right(p_rig, width = 0.2)

  if (!is.null(row.lab)) {
    p1 <- p1 %>% aplot::insert_left(p.row.lab, width = 0.02)
    p2 <- p2 %>% aplot::insert_left(p.row.lab, width = 0.02)
  }

  if (row_cluster == TRUE) {
    p1 <- p1 %>% aplot::insert_left(ggtree_plot, width = 0.2)
    p2 <- p2 %>% aplot::insert_left(ggtree_plot, width = 0.2)
  }

  if (label == TRUE) {
    p1 <- p1 %>% aplot::insert_top(labels, height = 0.02)
    p2 <- p2 %>% aplot::insert_top(labels, height = 0.02)
  }

  if (col_cluster == TRUE) {
    p1 <- p1 %>% aplot::insert_top(ggtree_plot_col, height = 0.1)
    p2 <- p2 %>% aplot::insert_top(ggtree_plot_col, height = 0.1)
  }

  return(list(p1, p2, plotdata = pcm))
}
