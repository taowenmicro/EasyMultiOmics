# library(ggplot2)
# library(ggtree)
# library(ape)
# library(dplyr)
# library(tidyr)
#
#
# set.seed(123)
# mat <- matrix(rnorm(100*8), nrow=100, ncol=8)
# rownames(mat) <- paste0("ASV_",1:100)
# colnames(mat) <- paste0("Sample",1:8)
#
#
#
# p <- circular_heatmap_with_tree0(mat, cluster = TRUE, open = 90,
#                                 fill_colors = c("white","orange","red"),
#                                 label_offset = 2
# )
# print(p)
#
#
#
# library(dplyr)
# library(tidyr)

#' @title Circular heatmap with taxonomy tree and abundance
#' @description
#' Draw a circular (or fan) tree with an outer ring heatmap showing microbial abundance across samples.
#'
#' @param mat numeric matrix (rows = taxa/OTUs/ASVs, cols = samples).
#'            rownames must be taxa IDs, colnames = sample names.
#' @param cluster logical, whether to perform hierarchical clustering on taxa.
#' @param open numeric, open angle (0 = full circle, 180 = half circle, etc.)
#' @param fill_colors vector of length 3, low-mid-high colors for abundance scale.
#' @param label_size numeric, size of tip labels.
#' @param label_offset numeric, distance of labels from heatmap.
#'
#' @return ggplot object
#' @export



circular_heatmap_with_tree0 <- function(mat,
                                       cluster = TRUE,
                                       open = 0,
                                       fill_colors = c("blue","white","red"),
                                       label_size = 2.5,
                                       label_offset = 1) {
  # === 1. clustering on taxa ===
  if (cluster) {
    row_hc <- hclust(dist(mat))        # 聚类基于丰度
    phylo_tree <- as.phylo(row_hc)
  } else {
    # 如果不聚类，直接构造一棵星型树
    phylo_tree <- as.phylo(hclust(dist(diag(nrow(mat)))))
    phylo_tree$tip.label <- rownames(mat)
  }

  # === 2. base tree ===
  p_tree <- ggtree(phylo_tree, layout = "fan", open.angle = open)

  # === 3. reorder abundance matrix by tip order ===
  tip_order <- p_tree$data$label[p_tree$data$isTip]
  mat2 <- mat[match(tip_order, rownames(mat)), , drop = FALSE]

  # === 4. attach heatmap (abundance) ===
  p_out <- gheatmap(p_tree, mat2,
                    width = 0.5,              # 热图厚度
                    colnames = TRUE,
                    colnames_angle = 90,
                    colnames_offset_y = 0.3,
                    font.size = label_size) +
    scale_fill_gradient2(low = fill_colors[1],
                         mid = fill_colors[2],
                         high = fill_colors[3],
                         midpoint = mean(mat2, na.rm = TRUE),
                         name = "Abundance") +
    theme(legend.position = "right")

  # === 5. add tip labels ===
  p_out <- p_out +
    geom_tiplab(size = label_size,
                offset = label_offset)

  return(p_out)
}


