#' @title Polygonal composition plot from phyloseq
#' @description
#' Generate polygon plots (k-gon) from a phyloseq object, where k is the number of groups.
#' Each point represents a taxon (OTU/ASV/Genus...), positioned inside the polygon
#' according to its relative abundance distribution across groups.
#'
#' @param ps A \code{phyloseq} object containing microbiome data.
#' @param group A character string specifying the grouping variable in sample metadata.
#'        The number of unique groups determines the polygon edges (e.g. 4 groups → quadrilateral).
#' @param taxrank A taxonomy rank (e.g. "Phylum", "Class", "Order", "Family", "Genus")
#'        used for coloring the taxa points. Default = "Phylum".
#' @param rotate Rotation angle in radians for the polygon. Default = \code{-pi/2}
#'        (first vertex points upwards).
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' library(phyloseq)
#' # Example: polygon plot at Phylum level
#' ps_polygon_plot(ps = ps.16s, group = "Group", taxrank = "Phylum")
#' }
#' @author Tao Wen
#' @export
ps_polygon_plot <- function(ps, group = "Group", taxrank = "Phylum", rotate = -pi/2) {
  # 提取 metadata
  map <- as.data.frame(phyloseq::sample_data(ps))
  stopifnot(group %in% colnames(map))
  groups <- factor(map[[group]])
  gnames <- levels(groups)

  # 提取 OTU 表
  otu <- as.data.frame(phyloseq::otu_table(ps))
  if (!phyloseq::taxa_are_rows(ps)) otu <- t(otu)  # taxa 在行

  # 聚合样本 → 分组
  otu_group <- rowsum(t(otu), group = groups)   # (group × taxa)
  otu_group <- t(otu_group)                     # (taxa × group)

  # 去掉全0行
  otu_group <- otu_group[rowSums(otu_group) > 0, , drop = FALSE]

  # 归一化（每行和=1）
  otu_prop <- otu_group / rowSums(otu_group)
  otu_prop[is.na(otu_prop)] <- 0

  # 构造数据框
  df <- as.data.frame(otu_prop)
  df$TaxaID <- rownames(df)

  # taxonomy 表
  tax <- as.data.frame(phyloseq::tax_table(ps))
  if (!taxrank %in% colnames(tax)) stop("taxrank not found in taxonomy table")
  df$TaxRank <- tax[df$TaxaID, taxrank]

  # 自动生成调色板
  n_col <- length(unique(df$TaxRank))
  base_cols <- c(
    "#1f77b4","#ff7f0e","#2ca02c",
    "#d62728","#9467bd","#8c564b",
    "#17becf","#bcbd22","#7f7f7f"
  )
  palette <- grDevices::colorRampPalette(base_cols)(n_col)

  # 调用绘图函数
  polygon_composition_plot(
    df, comp_cols = gnames, id_col = "TaxRank", rotate = rotate, palette = palette
  )
}


#' @title General polygon composition plot
#' @description
#' Plot points inside a regular k-gon (triangle, quadrilateral, pentagon, ...),
#' where each point is positioned according to its composition across k groups
#' (barycentric coordinates).
#'
#' @param df A data.frame containing composition data and optional grouping variable.
#' @param comp_cols Character vector, names of columns representing the compositions
#'        (values must sum to 1 per row after normalization).
#' @param id_col Character string, optional column name used for coloring points.
#' @param rotate Rotation angle in radians. Default = \code{-pi/2}.
#' @param palette Vector of colors or \code{NULL}. If \code{NULL}, colors are auto-generated.
#' @param point_size Point size. Default = 2.5.
#' @param alpha_point Transparency of points. Default = 0.85.
#'
#' @return A \code{ggplot} object.
#' @examples
#' \dontrun{
#' # Example: 4-group polygon
#' df <- data.frame(G1 = c(0.7,0.1,0.2), G2 = c(0.2,0.6,0.1),
#'                  G3 = c(0.1,0.2,0.6), G4 = c(0,0.1,0.1),
#'                  Phylum = c("Firmicutes","Proteobacteria","Actinobacteria"))
#' polygon_composition_plot(df, comp_cols = c("G1","G2","G3","G4"), id_col = "Phylum")
#' }
#' @author Tao Wen
#' @export
polygon_composition_plot <- function(df, comp_cols, id_col = NULL,
                                     rotate = -pi/2, palette = NULL,
                                     point_size = 2.5, alpha_point = 0.85) {
  stopifnot(all(comp_cols %in% names(df)))
  k <- length(comp_cols)
  if (k < 3) stop("需要 >=3 个组分")

  # 1) 行归一化
  W <- as.matrix(df[, comp_cols, drop = FALSE])
  row_sums <- rowSums(W)
  W <- W[row_sums > 0, , drop = FALSE]  # 删除全0行
  row_sums <- rowSums(W)
  W <- W / row_sums

  # 2) 正 k 边形顶点坐标
  angles <- rotate + 2*pi*(0:(k-1))/k
  vx <- cos(angles)
  vy <- sin(angles)
  V <- cbind(vx, vy)

  # 3) 投影到二维
  XY <- W %*% V
  colnames(XY) <- c("x","y")

  # 4) 绘图数据
  plot_df <- data.frame(XY, df[row.names(W), comp_cols, drop = FALSE])
  if (!is.null(id_col) && id_col %in% names(df)) {
    plot_df[[id_col]] <- df[row.names(W), id_col]
  }

  # 5) 顶点坐标
  poly_df <- data.frame(
    x = vx, y = vy, label = comp_cols, idx = factor(seq_len(k), levels = seq_len(k))
  )
  poly_df_closed <- rbind(poly_df, poly_df[1,])

  # 6) 自动配色
  if (is.null(palette)) {
    base_cols <- c(
      "#1f77b4", "#ff7f0e", "#2ca02c",
      "#d62728", "#9467bd", "#8c564b",
      "#17becf", "#bcbd22", "#7f7f7f"
    )
    n_colors <- if (!is.null(id_col)) length(unique(df[[id_col]])) else 1
    palette <- grDevices::colorRampPalette(base_cols)(n_colors)
  }

  # 7) 绘制
  p <- ggplot2::ggplot() +
    ggplot2::geom_path(
      data = poly_df_closed,
      ggplot2::aes(x = x, y = y),
      linewidth = 0.6, color = "grey30"
    ) +
    ggplot2::geom_segment(
      data = poly_df,
      ggplot2::aes(x = 0, y = 0, xend = x, yend = y),
      linewidth = 0.3, color = "grey70", linetype = 2
    ) +
    ggplot2::geom_point(
      data = poly_df,
      ggplot2::aes(x = x, y = y, color = idx),
      size = 3, stroke = 0.2
    ) +
    ggplot2::geom_label(
      data = poly_df,
      ggplot2::aes(x = 1.08*x, y = 1.08*y, label = label, color = idx),
      size = 3.5, label.padding = ggplot2::unit(0.15, "lines"),
      label.size = 0.1, label.r = ggplot2::unit(0.1, "lines"),
      fill = ggplot2::alpha("white", 0.8)
    ) +
    {
      if (!is.null(id_col) && id_col %in% names(plot_df)) {
        ggplot2::geom_point(
          data = plot_df,
          ggplot2::aes(x = x, y = y, fill = .data[[id_col]]),
          shape = 21, size = point_size, color = "grey25", alpha = alpha_point
        )
      } else {
        ggplot2::geom_point(
          data = plot_df,
          ggplot2::aes(x = x, y = y),
          shape = 21, size = point_size, fill = "grey40", color = "grey15", alpha = alpha_point
        )
      }
    } +
    ggplot2::coord_equal(xlim = c(-1.25, 1.25), ylim = c(-1.25, 1.25), expand = TRUE) +
    ggplot2::scale_color_manual(values = palette[seq_len(k)], guide = "none") +
    ggplot2::scale_fill_manual(values = palette, guide = ggplot2::guide_legend(title = id_col)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank()
    )

  return(p)
}
