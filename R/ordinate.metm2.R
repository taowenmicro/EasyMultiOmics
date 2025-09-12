#' Enhanced Beta Diversity Ordination Analysis
#'
#' @description
#' Performs beta diversity ordination analysis with multiple group support and enhanced visualization.
#' This function extends traditional ordination methods by supporting up to three grouping variables
#' with distinct aesthetic mappings (color, shape, fill) for comprehensive multivariate visualization.
#'
#' @param otu OTU/ASV abundance table with samples as rows and taxa as columns. Can be NULL if \code{ps} is provided.
#' @param tax Taxonomic table with taxa as rows and taxonomic ranks as columns. Can be NULL if \code{ps} is provided.
#' @param map Sample metadata table with samples as rows and variables as columns. Can be NULL if \code{ps} is provided.
#' @param ps A phyloseq object containing OTU table, taxonomic table, and sample metadata. Alternative input to \code{otu}, \code{tax}, and \code{map}.
#' @param group Character vector specifying grouping variable name(s) from sample metadata. Supports 1-3 variables:
#'   \itemize{
#'     \item 1 variable: mapped to point color
#'     \item 2 variables: mapped to color and shape
#'     \item 3 variables: mapped to color, shape, and fill
#'   }
#' @param dist Distance metric for beta diversity calculation. Options include:
#'   \itemize{
#'     \item \strong{Phylogenetic}: "unifrac", "wunifrac"
#'     \item \strong{Ecological}: "bray", "jaccard", "kulczynski", "gower", "morisita", "horn"
#'     \item \strong{Geometric}: "euclidean", "manhattan", "canberra"
#'     \item \strong{Other}: "jsd", "cao", "binomial", "chao" and 30+ additional metrics
#'   }
#'   Default: "bray"
#' @param method Ordination method. Options:
#'   \itemize{
#'     \item \strong{"PCoA"}: Principal Coordinates Analysis (default)
#'     \item \strong{"NMDS"}: Non-metric Multidimensional Scaling
#'     \item \strong{"PCA"}: Principal Component Analysis
#'     \item \strong{"DCA"}: Detrended Correspondence Analysis
#'     \item \strong{"CCA"}: Canonical Correspondence Analysis
#'     \item \strong{"RDA"}: Redundancy Analysis
#'     \item \strong{"MDS"}: Multidimensional Scaling
#'     \item \strong{"LDA"}: Linear Discriminant Analysis
#'     \item \strong{"t-sne"}: t-Distributed Stochastic Neighbor Embedding
#'   }
#' @param Micromet Statistical test method for group comparison: "adonis", "anosim", or "MRPP". Default: "adonis"
#' @param pvalue.cutoff P-value threshold for statistical significance. Default: 0.05
#' @param pair Logical value indicating whether to perform pairwise group comparisons. Default: FALSE
#'
#' @details
#' The function performs several key steps:
#' \enumerate{
#'   \item Converts abundance data to relative abundances
#'   \item Calculates distance matrix using specified metric
#'   \item Performs ordination analysis using specified method
#'   \item Creates publication-ready ggplot2 visualization with intelligent aesthetic mapping
#'   \item Supports multiple grouping variables with distinct visual encodings
#' }
#'
#' \strong{Aesthetic Mapping Strategy:}
#' \itemize{
#'   \item \strong{1 group}: Color differentiation with solid shapes
#'   \item \strong{2 groups}: Color (group 1) + Shape (group 2) with expanded shape palette
#'   \item \strong{3 groups}: Color (group 1) + Shape (group 2) + Fill (group 3) with fillable shapes
#' }
#'
#' \strong{Shape Selection:}
#' \itemize{
#'   \item 2 groups: Uses extended shape palette including solid and hollow shapes (16,17,15,18,8,11,13,14,0,1,2,5,6)
#'   \item 3 groups: Automatically uses fillable shapes (21-25) to support the fill aesthetic
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{plot}}{Main ggplot2 object with ordination visualization}
#'   \item{\code{points}}{Data frame with ordination coordinates and metadata}
#'   \item{\code{plot_labeled}}{ggplot2 object with sample labels added}
#'   \item{\code{eigenvalues/stress}}{Eigenvalues (for most methods) or stress value (for NMDS)}
#' }
#'
#' @note
#' \itemize{
#'   \item Requires \code{phyloseq}, \code{ggplot2}, \code{ggrepel} packages
#'   \item For LDA method, automatically reduces feature dimensionality if > 10 features
#'   \item NMDS method returns stress value instead of eigenvalues
#'   \item Function automatically handles sample-wise relative abundance transformation
#' }
#'
#' @seealso
#' \code{\link[phyloseq]{distance}}, \code{\link[phyloseq]{ordinate}}, \code{\link[vegan]{vegdist}}
#'
#' @examples
#' \dontrun{
#' # Load required packages and data
#' library(phyloseq)
#' library(ggplot2)
#'
#' # Single grouping variable
#' result1 <- ordinate.metm2(ps = ps.16s,
#'                         group = "Treatment",
#'                         dist = "bray",
#'                         method = "PCoA")
#' print(result1[[1]])  # Main plot
#'
#' # Two grouping variables
#' result2 <- ordinate.metm2(ps = ps.16s,
#'                         group = c("Treatment", "Timepoint"),
#'                         dist = "unifrac",
#'                         method = "PCoA")
#' print(result2[[1]])
#'
#' # Three grouping variables
#' result3 <- ordinate.metm2(ps = ps.16s,
#'                         group = c("Treatment", "Timepoint", "Subject"),
#'                         dist = "bray",
#'                         method = "PCoA")
#' print(result3[[1]])
#'
#' # NMDS analysis
#' nmds_result <- ordinate.metm2(ps = ps.16s,
#'                             group = "Treatment",
#'                             dist = "bray",
#'                             method = "NMDS")
#' print(nmds_result[[1]])  # Plot includes stress value
#' }
#'
#' @author Contact: taowen \email{taowen@@njau.edu.cn}
#' @export

ordinate.metm2 <- function(otu = NULL, tax = NULL, map = NULL,
                          ps = NULL, group = "Group",
                          dist = "bray", method = "PCoA", Micromet = "adonis",
                          pvalue.cutoff = 0.05, pair = FALSE,
                          width = 0.005, height = 0.005,stroke = 2
                          ) {
  jit <- position_jitter(width = width, height = height)
  # Input validation
  if (is.null(ps) && (is.null(otu) || is.null(map))) {
    stop("Either 'ps' must be provided, or both 'otu' and 'map' must be provided.")
  }

  if (!method %in% c("DCA", "CCA", "RDA", "MDS", "PCoA", "PCA", "LDA", "NMDS", "t-sne")) {
    stop("Method must be one of: DCA, CCA, RDA, MDS, PCoA, PCA, LDA, NMDS, t-sne")
  }

  if (length(group) > 3) {
    warning("Only the first 3 grouping variables will be used.")
    group <- group[1:3]
  }

  # Transform to relative abundance
  ps1_rela <- phyloseq::transform_sample_counts(ps, function(x) x/sum(x))

  # Calculate ordination based on method
  if (method == "DCA") {
    ordi <- phyloseq::ordinate(ps1_rela, method = method, distance = dist)
    points <- ordi$rproj[, 1:2]
    colnames(points) <- c("x", "y")
    eig <- ordi$evals^2
  }

  if (method == "CCA") {
    ordi <- phyloseq::ordinate(ps1_rela, method = method, distance = dist)
    points <- ordi$CA$v[, 1:2]
    colnames(points) <- c("x", "y")
    eig <- ordi$CA$eig^2
  }

  if (method == "RDA") {
    ordi <- phyloseq::ordinate(ps1_rela, method = method, distance = dist)
    points <- ordi$CA$v[, 1:2]
    colnames(points) <- c("x", "y")
    eig <- ordi$CA$eig
  }

  if (method == "MDS") {
    ordi <- phyloseq::ordinate(ps1_rela, method = method, distance = dist)
    points <- ordi$vectors[, 1:2]
    colnames(points) <- c("x", "y")
    eig <- ordi$values[, 1]
  }

  if (method == "PCoA") {
    unif <- phyloseq::distance(ps1_rela, method = dist, type = "samples")
    pcoa <- stats::cmdscale(unif, k = 2, eig = TRUE)
    points <- as.data.frame(pcoa$points)
    colnames(points) <- c("x", "y")
    eig <- pcoa$eig
  }

  # Get OTU table for PCA and LDA
  otu_table <- as.data.frame(t(ggClusterNet::vegan_otu(ps1_rela)))

  if (method == "PCA") {
    otu.pca <- stats::prcomp(t(otu_table), scale.default = TRUE)
    points <- otu.pca$x[, 1:2]
    colnames(points) <- c("x", "y")
    eig <- otu.pca$sdev
    eig <- eig * eig
  }

  if (method == "LDA") {
    data <- t(otu_table)
    data <- as.data.frame(data)
    data <- scale(data, center = TRUE, scale = TRUE)
    # Reduce dimensionality if too many features
    if(ncol(data) > 10) {
      data <- data[, 1:10]
    }
    map_temp <- as.data.frame(phyloseq::sample_data(ps1_rela))
    # Use first grouping variable for LDA
    model <- MASS::lda(data, map_temp[[group[1]]])
    ord_in <- model
    axes <- c(1:2)
    points <- data.frame(predict(ord_in)$x[, axes])
    colnames(points) <- c("x", "y")
    eig <- ord_in$svd^2
  }

  if (method == "NMDS") {
    ordi <- phyloseq::ordinate(ps1_rela, method = method, distance = dist)
    points <- ordi$points[, 1:2]
    colnames(points) <- c("x", "y")
    stress <- ordi$stress
  }

  if (method == "t-sne") {
    data <- t(otu_table)
    data <- as.data.frame(data)
    data <- scale(data, center = TRUE, scale = TRUE)
    map_temp <- as.data.frame(phyloseq::sample_data(ps1_rela))
    tsne <- Rtsne::Rtsne(data, perplexity = 3)
    points <- as.data.frame(tsne$Y)
    row.names(points) <- row.names(map_temp)
    colnames(points) <- c("x", "y")
  }

  # Prepare metadata
  map <- as.data.frame(phyloseq::sample_data(ps1_rela))
  group <- unique(group)  # Remove duplicates
  for (g in group) {
    map[[g]] <- as.factor(map[[g]])
  }
  points <- cbind(points, map[match(rownames(points), rownames(map)), ])
  points$ID <- row.names(points)

  # Create visualization based on number of grouping variables
  if (length(group) == 1) {
    # Single group: color mapping with nature-style dark colors
    n_colors <- length(unique(points[[group[1]]]))
    # Nature-style dark color palette
    color_values <- c("#1f4e79", "#2d5016", "#8b1538", "#5d4e75", "#b45f04", "#1c4a3a", "#8b4513", "#2f4f4f")[1:n_colors]


    p2 <- ggplot(points, aes(x = x, y = y, color = !!sym(group[1]))) +
      geom_point(alpha = 0.7, size = 5, stroke = 1.2, shape = 19, position = jit) +
      scale_color_manual(values = color_values)

  } else if (length(group) == 2) {
    # Two groups: color + shape mapping
    n_shapes <- length(unique(points[[group[2]]]))
    n_colors <- length(unique(points[[group[1]]]))
    # Extended shape palette for non-fill usage
    shape_values <- c(16, 17, 15, 18, 8, 11, 13, 14, 0, 1, 2, 5, 6)[1:n_shapes]
    # Nature-style dark colors for first group
    color_values <- c("#1f4e79", "#2d5016", "#8b1538", "#5d4e75", "#b45f04", "#1c4a3a", "#8b4513", "#2f4f4f")[1:n_colors]

    p2 <- ggplot(points, aes(x = x, y = y,
                             color = !!sym(group[1]),
                             shape = !!sym(group[2]))) +
      geom_point(alpha = 0.7, size = 5, stroke = 1.2, position = jit) +
      scale_shape_manual(values = shape_values) +
      scale_color_manual(values = color_values)

  } else if (length(group) >= 3) {
    # Three groups: color + shape + fill mapping
    n_shapes <- length(unique(points[[group[2]]]))
    n_colors <- length(unique(points[[group[1]]]))
    n_fills <- length(unique(points[[group[3]]]))

    # Must use fillable shapes when using fill aesthetic
    shape_values <- c(21, 22, 23, 24, 25)[1:n_shapes]

    # Nature-style dark colors for first group (edge colors)
    color_values <- c("#1f4e79", "#2d5016", "#8b1538", "#5d4e75", "#b45f04", "#1c4a3a", "#8b4513", "#2f4f4f")[1:n_colors]

    # Bright, soft fill colors for third group that contrast with dark edges
    fill_values <- c("#FFE5B4", "#E5F3FF", "#E8F5E8", "#FFF0E6", "#F0E5FF", "#E5F9F6",
                     "#FFE5E5", "#F0F8E5", "#E5E5FF", "#FFE5F0")[1:n_fills]

    p2 <- ggplot(points, aes(x = x, y = y,
                             color = !!sym(group[1]),
                             shape = !!sym(group[2]),
                             fill = !!sym(group[3]))) +
      geom_point(alpha = 0.7, size = 5, stroke = stroke, position = jit) +
      scale_shape_manual(values = shape_values) +
      scale_color_manual(values = color_values) +
      scale_fill_manual(values = fill_values) +
      guides(
        color = guide_legend(
          override.aes = list(shape = 19, fill = NA, stroke = 1.5, alpha = 1)
        ),
        shape = guide_legend(
          override.aes = list(color = "black", fill = "grey", stroke = 1, alpha = 1)  # 图例中也增加边框厚度
        ),
        fill = guide_legend(
          override.aes = list(shape = 21, color = "black", stroke = 1, alpha = 1)  # 图例中也增加边框厚度
        )
      )
  }

  # Add axis labels based on method
  if (method %in% c("DCA", "CCA", "RDA", "MDS", "PCoA", "PCA", "LDA")) {
    p2 <- p2 +
      labs(x = paste0(method, " 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)"),
           y = paste0(method, " 2 (", format(100 * eig[2]/sum(eig), digits = 4), "%)"))
  } else if (method == "NMDS") {
    p2 <- p2 +
      labs(x = paste0(method, "1"),
           y = paste0(method, "2"),
           title = paste("Stress:", round(stress, 3)))
  } else if (method == "t-sne") {
    p2 <- p2 +
      labs(x = paste0(method, "1"),
           y = paste0(method, "2"))
  }

  # Apply theme
  p2 <- p2 + theme_bw()

  # Create labeled version
  p3 <- p2 + ggrepel::geom_text_repel(aes(label = ID), size = 4)

  # Return results
  if (method == "NMDS") {
    return(list(plot = p2, points = points, plot_labeled = p3, stress = stress))
  } else {
    return(list(plot = p2, points = points, plot_labeled = p3, eigenvalues = eig))
  }
}

ordinate.metm0 = function (otu = NULL, tax = NULL, map = NULL,
                          ps = NULL, group = "Group",
                          dist = "bray", method = "PCoA", Micromet = "adonis",
                          pvalue.cutoff = 0.05, pair = FALSE)
{
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))

  # --- 计算 ordination ---
  if (method == "PCoA") {
    unif = phyloseq::distance(ps1_rela, method = dist, type = "samples")
    pcoa = stats::cmdscale(unif, k = 2, eig = TRUE)
    points = as.data.frame(pcoa$points)
    colnames(points) = c("x", "y")
    eig = pcoa$eig
  }
  # 省略其他 method（DCA, RDA, PCA, NMDS ...）

  # --- 准备 metadata ---
  map = as.data.frame(phyloseq::sample_data(ps1_rela))
  group = unique(group)  # 防止重复
  for (g in group) {
    map[[g]] = as.factor(map[[g]])
  }
  points = cbind(points, map[match(rownames(points), rownames(map)), ])
  points$ID = row.names(points)

  # --- 作图 ---
  if (length(group) == 1) {
    # 单分组：第一个分组用颜色
    p2 = ggplot(points, aes(x = x, y = y, color = !!sym(group[1]))) +
      geom_point(alpha = 0.7, size = 5, stroke = 1.2, shape = 19)

  } else if (length(group) == 2) {
    # 双分组：第一个分组用颜色 + 第二个分组用形状
    n_shapes <- length(unique(points[[group[2]]]))
    # 不使用fill时，可以用所有形状
    shape_values <- c(16, 17, 15, 18, 8, 11, 13, 14, 0, 1, 2, 5, 6)[1:n_shapes]

    p2 = ggplot(points, aes(x = x, y = y,
                            color = !!sym(group[1]),
                            shape = !!sym(group[2]))) +
      geom_point(alpha = 0.7, size = 5, stroke = 1.2) +
      scale_shape_manual(values = shape_values)

  } else if (length(group) >= 3) {
    # 三分组：第一个分组用颜色 + 第二个分组用形状 + 第三个分组用填充
    n_shapes <- length(unique(points[[group[2]]]))
    # 使用fill时，必须用有填充属性的形状（21-25）
    shape_values <- c(21, 22, 23, 24, 25)[1:n_shapes]

    # 为第三个分组设置填充颜色调色板
    n_fills <- length(unique(points[[group[3]]]))
    # 使用柔和的填充色，与边框色区分
    fill_values <- c("#FFB6C1", "#87CEEB", "#98FB98", "#F0E68C", "#DDA0DD", "#F0F8FF")[1:n_fills]

    p2 = ggplot(points, aes(x = x, y = y,
                            color = !!sym(group[1]),
                            shape = !!sym(group[2]),
                            fill = !!sym(group[3]))) +
      geom_point(alpha = 0.7, size = 5, stroke = 1.2) +
      scale_shape_manual(values = shape_values) +
      scale_fill_manual(values = fill_values) +
      guides(
        color = guide_legend(
          override.aes = list(shape = 19, fill = NA, stroke = 1.2, alpha = 1)
        ),
        shape = guide_legend(
          override.aes = list(color = "black", fill = "grey", stroke = 1.2, alpha = 1)
        ),
        fill = guide_legend(
          override.aes = list(shape = 21, color = "black", stroke = 1.2, alpha = 1)
        )
      )
  }

  # 添加坐标轴标签和主题
  p2 = p2 +
    labs(x = paste0(method, " 1 (", format(100 * eig[1]/sum(eig), digits = 4), "%)"),
         y = paste0(method, " 2 (", format(100 * eig[2]/sum(eig), digits = 4), "%)")) +
    theme_bw()

  # 添加标签
  p3 = p2 + ggrepel::geom_text_repel(aes(label = ID), size = 4)

  return(list(p2, points, p3, eig))
}
