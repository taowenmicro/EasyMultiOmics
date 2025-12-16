#' @title Differential Abundance Analysis Using DESeq2
#'
#' @description
#' The \code{DESep2Super.metm} function performs differential abundance analysis between
#' microbial groups using the \code{DESeq2} package. It identifies significantly
#' enriched or depleted taxa across multiple group comparisons by calculating log2
#' fold changes, raw p-values, and adjusted p-values (FDR).
#'
#' The input data can be provided as a \code{phyloseq} object or as separate OTU,
#' taxonomic, and metadata tables.
#'
#' @param otu A data frame containing OTU counts. Optional if \code{ps} is provided.
#' @param tax A data frame containing taxonomic annotations. Optional if \code{ps} is provided.
#' @param map A data frame containing sample metadata. Optional if \code{ps} is provided.
#' @param tree A phylogenetic tree. Optional if \code{ps} is provided.
#' @param ps A \code{phyloseq} object containing microbiome data. If provided,
#'   overrides \code{otu}, \code{tax}, \code{map}, and \code{tree}.
#' @param j A string or integer specifying the taxonomic rank to perform the analysis on.
#'   Can be a numeric rank (1-7), a taxonomic name (e.g., \code{"Phylum"}), or \code{"OTU"}.
#'   Default is \code{"Genus"}.
#' @param group A string specifying the grouping variable in the sample metadata.
#'   Default is \code{"Group"}.
#' @param pvalue A numeric value specifying the significance threshold for adjusted
#'   p-values (FDR). Default is \code{0.05}.
#' @param artGroup A custom matrix specifying pairwise group comparisons. Optional.
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{Plots}{A list of volcano plots for each pairwise comparison.}
#'   \item{Results}{A data frame containing differential abundance results, including
#'     log2 fold changes, p-values, adjusted p-values, and group-specific normalized
#'     abundances.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Prepares OTU, taxonomic, and metadata tables from the \code{phyloseq} object
#'     or individual inputs.
#'   \item Aggregates data to the specified taxonomic rank (if applicable).
#'   \item Constructs a \code{DESeq2} dataset object and normalizes OTU counts.
#'   \item Performs pairwise group comparisons using contrasts in the \code{DESeq2} framework.
#'   \item Calculates log2 fold changes, raw p-values, and FDR-adjusted p-values for each taxon.
#'   \item Labels taxa as \code{"enriched"}, \code{"depleted"}, or \code{"nosig"} based on thresholds.
#'   \item Outputs a combined data frame with differential abundance results, normalized
#'     abundances, and taxonomic annotations.
#' }
#'
#' The results include volcano plots that visualize the differential abundance results
#' for each group comparison. Significant taxa are highlighted, and their log2 fold
#' changes and p-values are displayed.
#'
#' @examples
#' \dontrun{
#' res <- DESep2Super.metm(
#'   ps = ps.16s %>% ggClusterNet::filter_OTU_ps(500),
#'   group = "Group",
#'   artGroup = NULL,
#'   j = "OTU"
#' )
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @export


DESep2Super.metm <- function(otu = NULL, tax = NULL, map = NULL, tree = NULL, ps = NULL,
                             j = "Genus", group = "Group", pvalue = 0.05, artGroup = NULL,
                             col.g = NULL, gradient = TRUE, mid_color = "auto", top_n = 5) {

  get_gradient_colors <- function(color1, color2) {
    rgb1 <- col2rgb(color1)
    rgb2 <- col2rgb(color2)

    mid_rgb <- (rgb1 + rgb2) / 2
    mid_col <- rgb(mid_rgb[1], mid_rgb[2], mid_rgb[3], maxColorValue = 255)

    q1_rgb <- (rgb1 + mid_rgb) / 2
    q1_col <- rgb(q1_rgb[1], q1_rgb[2], q1_rgb[3], maxColorValue = 255)

    q3_rgb <- (rgb2 + mid_rgb) / 2
    q3_col <- rgb(q3_rgb[1], q3_rgb[2], q3_rgb[3], maxColorValue = 255)

    return(c(color2, q3_col, mid_col, q1_col, color1))
  }

  ps = ggClusterNet::inputMicro(otu, tax, map, tree, ps, group = group)

  if (j %in% c("OTU", "gene", "meta")) {
    ps = ps
  } else if (j %in% c(1:7)) {
    ps = ps %>% ggClusterNet::tax_glom_wt(ranks = j)
  } else if (j %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) {
    ps = ps %>% ggClusterNet::tax_glom_wt(ranks = j)
  } else {
    ps = ps
    print("unknown j, checked please")
  }

  Desep_group <- ps %>%
    phyloseq::sample_data() %>%
    .$Group %>%
    as.factor() %>%
    levels() %>%
    as.character()

  if (is.null(artGroup)) {
    aaa = combn(Desep_group, 2)
  } else {
    aaa = as.matrix(artGroup)
  }

  if (is.null(col.g)) {
    library(ggsci)
    col.g <- pal_npg()(length(Desep_group))
    names(col.g) <- Desep_group
  }

  count <- ps %>%
    ggClusterNet::vegan_otu() %>%
    round(0) %>%
    t()

  map = ps %>%
    phyloseq::sample_data() %>%
    data.frame(check.names = FALSE)

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count,
    colData = map,
    design = ~Group
  )
  dds2 <- DESeq2::DESeq(dds)

  plot_list <- list()

  for (i in 1:dim(aaa)[2]) {
    Desep_group = aaa[, i]
    print(Desep_group)

    group_name = paste(Desep_group[1], Desep_group[2], sep = "-")
    group1 = Desep_group[1]
    group2 = Desep_group[2]

    res <- DESeq2::results(dds2, contrast = c("Group", Desep_group), alpha = 0.05)
    x = res

    x$level = as.factor(ifelse(
      as.vector(x$padj) < pvalue & x$log2FoldChange > 0, "enriched",
      ifelse(as.vector(x$padj) < pvalue & x$log2FoldChange < 0, "depleted", "nosig")
    ))

    x = data.frame(
      row.names = row.names(x),
      logFC = x$log2FoldChange,
      level = x$level,
      p = x$pvalue
    )

    x$taxa = row.names(x)
    x$p[is.na(x$p)] = 1
    x$log10p = -log10(x$p)

    df_nosig <- x %>% dplyr::filter(level == "nosig")
    df_sig <- x %>% dplyr::filter(level != "nosig")

    df_top <- df_sig %>%
      dplyr::mutate(ord = logFC^2) %>%
      dplyr::arrange(desc(ord)) %>%
      head(n = top_n)

    df_sig_no_border <- df_sig %>%
      dplyr::filter(!taxa %in% df_top$taxa)

    color1 <- as.character(col.g[group1])
    color2 <- as.character(col.g[group2])

    grad_colors <- get_gradient_colors(color1, color2)
    cat("  渐变色:", paste(grad_colors, collapse = " -> "), "\n")

    if (gradient) {
      max_abs_fc <- max(abs(x$logFC), na.rm = TRUE)

      p <- ggplot2::ggplot() +
        ggplot2::annotate("rect", xmin = -Inf, xmax = -1, ymin = -Inf, ymax = Inf,
                          fill = color2, alpha = 0.05) +
        ggplot2::annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
                          fill = color1, alpha = 0.05) +
        ggplot2::geom_point(data = df_nosig,
                            ggplot2::aes(x = logFC, y = log10p),
                            color = "grey80", size = 2, alpha = 0.5) +
        ggplot2::geom_point(data = df_sig_no_border,
                            ggplot2::aes(x = logFC, y = log10p, color = logFC),
                            size = 2.5, alpha = 0.8) +
        ggplot2::geom_point(data = df_top,
                            ggplot2::aes(x = logFC, y = log10p, fill = logFC),
                            shape = 21,
                            color = "black",
                            size = 3.5,
                            stroke = 0.8,
                            alpha = 0.9) +
        ggplot2::scale_color_gradientn(
          colours = grad_colors,
          values = seq(0, 1, 0.25),
          limits = c(-max_abs_fc, max_abs_fc),
          name = "log2FC"
        ) +
        ggplot2::scale_fill_gradientn(
          colours = grad_colors,
          values = seq(0, 1, 0.25),
          limits = c(-max_abs_fc, max_abs_fc),
          guide = "none"
        ) +
        ggplot2::geom_hline(yintercept = -log10(pvalue), linetype = "dashed", color = "grey40") +
        ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
        ggrepel::geom_text_repel(
          data = df_top,
          ggplot2::aes(x = logFC, y = log10p, label = taxa),
          size = 3, max.overlaps = 20, segment.color = "grey50"
        ) +
        ggplot2::labs(
          title = paste(group1, "vs", group2),
          x = "log2 Fold Change",
          y = "-log10(P-value)"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right",
          panel.grid.minor = ggplot2::element_blank()
        )

    } else {
      my_colors <- c(
        "enriched" = color1,
        "depleted" = color2
      )

      p <- ggplot2::ggplot() +
        ggplot2::annotate("rect", xmin = -Inf, xmax = -1, ymin = -Inf, ymax = Inf,
                          fill = color2, alpha = 0.05) +
        ggplot2::annotate("rect", xmin = 1, xmax = Inf, ymin = -Inf, ymax = Inf,
                          fill = color1, alpha = 0.05) +
        ggplot2::geom_point(data = df_nosig,
                            ggplot2::aes(x = logFC, y = log10p),
                            color = "grey70", size = 2, alpha = 0.7) +
        ggplot2::geom_point(data = df_sig_no_border,
                            ggplot2::aes(x = logFC, y = log10p, color = level),
                            size = 2, alpha = 0.7) +
        ggplot2::geom_point(data = df_top,
                            ggplot2::aes(x = logFC, y = log10p, fill = level),
                            shape = 21,
                            color = "black",
                            size = 3.5,
                            stroke = 0.8,
                            alpha = 0.9) +
        ggplot2::scale_color_manual(
          values = my_colors,
          labels = c(
            "enriched" = paste0(group1, " enriched"),
            "depleted" = paste0(group2, " enriched")
          )
        ) +
        ggplot2::scale_fill_manual(
          values = my_colors,
          guide = "none"
        ) +
        ggplot2::geom_hline(yintercept = -log10(pvalue), linetype = "dashed", color = "grey40") +
        ggplot2::geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
        ggrepel::geom_text_repel(
          data = df_top,
          ggplot2::aes(x = logFC, y = log10p, label = taxa),
          size = 3, max.overlaps = 20, segment.color = "grey50"
        ) +
        ggplot2::labs(
          title = paste(group1, "vs", group2),
          x = "log2 Fold Change",
          y = "-log10(P-value)",
          color = "Regulation"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
          plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
          legend.position = "right",
          panel.grid.minor = ggplot2::element_blank()
        )
    }

    plot_list[[i]] <- p

    x_out <- x[, c("logFC", "level", "p")]
    colnames(x_out) = paste(group_name, colnames(x_out), sep = "")

    if (i == 1) {
      table = x_out
    } else {
      table = cbind(table, x_out)
    }
  }

  count = as.matrix(count)
  norm = t(t(count) / colSums(count))
  norm1 = norm %>% t() %>% as.data.frame()

  iris.split <- split(norm1, as.factor(map$Group))
  iris.apply <- lapply(iris.split, function(x) colMeans(x))
  iris.combine <- do.call(rbind, iris.apply)
  norm2 = t(iris.combine) %>% as.data.frame()

  x = cbind(table, norm2)

  if (!is.null(ps@tax_table)) {
    taxonomy = as.data.frame(ggClusterNet::vegan_tax(ps))

    if (length(colnames(taxonomy)) == 6) {
      colnames(taxonomy) = c("kingdom", "phylum", "class", "order", "family", "genus")
    } else if (length(colnames(taxonomy)) == 7) {
      colnames(taxonomy) = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
    } else if (length(colnames(taxonomy)) == 8) {
      colnames(taxonomy) = c("kingdom", "phylum", "class", "order", "family", "genus", "species", "rep")
    }

    taxonomy$id = rownames(taxonomy)
    tax = taxonomy[row.names(x), ]
    x = x[rownames(tax), ]

    if (length(colnames(taxonomy)) >= 7) {
      x$phylum = tax$phylum
      x$class = tax$class
      x$order = tax$order
      x$family = tax$family
      x$genus = tax$genus
    }
    if (length(colnames(taxonomy)) >= 8) {
      x$species = tax$species
    }
  }

  return(list(plot_list, x))
}
