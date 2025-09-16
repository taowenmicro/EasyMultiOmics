
# # === 把 p 值转成显著性 ===
#
# details_tab <- res_all$details %>%
#   mutate(sig = ifelse(!is.na(adjust.p) & adjust.p < 0.05, 1, 0)) %>%
#   dplyr::select(micro, method, sig)
#
# details_wide <- details_tab %>%
#   group_by(micro, method) %>%
#   summarise(sig = max(sig, na.rm = TRUE), .groups = "drop") %>%
#   tidyr::pivot_wider(
#     names_from = method,
#     values_from = sig,
#     values_fill = list(sig = 0)
#   ) %>%
#   as.data.frame()
#
# rownames(details_wide) <- details_wide$micro
# mat <- details_wide[, -1, drop = FALSE]
#
#
#
# p <- circular_heatmap_with_tree(
#   ps.cs,
#   mat,
#   open = 90,
#   label_size = 2.5,
#   label_offset = 60
# )
#
# print(p)


#' Circular heatmap with taxonomy tree
#'
#' This function plots a circular taxonomy tree from a \code{phyloseq} object
#' and attaches a binary (0/1) heatmap around the tree to visualize feature
#' detection results (e.g., significance across methods).
#'
#' @param ps A \code{phyloseq} object containing taxonomy information.
#' @param mat A matrix or data.frame with rows corresponding to OTUs/ASVs
#'   (rownames must match taxonomy IDs in \code{ps}), and columns representing samples or methods.
#'   Values should be binary (0/1) or factors convertible to 0/1.
#' @param open Numeric, the opening angle of the fan tree (default = 0 for full circle).
#' @param label_size Numeric, size of tip labels (default = 2.5).
#' @param label_offset Numeric, distance of tip labels from the heatmap (default = 20).
#' @param fill_colors Named vector of colors for 0/1 values, default:
#'   \code{c("0" = "grey90", "1" = "steelblue")}.
#'
#' @return A \code{ggplot} object with the circular tree + heatmap + tip labels.
#' @export
#' @examples
#' \dontrun{
#' library(phyloseq)
#' library(ggtree)
#' # ps: your phyloseq object
#' # mat: a binary matrix (OTUs x methods)
#' p <- circular_heatmap_with_tree(ps, mat)
#' print(p)
#' }
circular_heatmap_with_tree <- function(ps,
                                       mat,
                                       open = 0,
                                       label_size = 2.5,
                                       label_offset = 20,
                                       fill_colors = c("0" = "grey90", "1" = "steelblue")) {

  # --- 1. Subset phyloseq by taxa present in matrix ---
  ps_sub <- phyloseq::subset_taxa(ps, taxa_names(ps) %in% rownames(mat))

  # --- 2. Extract taxonomy table ---
  tax <- as.data.frame(phyloseq::tax_table(ps_sub))
  tax$ID <- rownames(tax)

  # --- 3. Construct taxonomy path ---
  tax$pathString <- apply(tax, 1, function(x) {
    paste(c("Root", na.omit(x), x["ID"]), collapse = "/")
  })

  # --- 4. Build phylo tree ---
  tax_tree <- data.tree::as.Node(tax)
  phylo_tree <- ape::as.phylo(tax_tree)

  # --- 5. Plot taxonomy tree ---
  p_tree <- ggtree::ggtree(phylo_tree, layout = "fan", open.angle = open)

  # --- 6. Reorder matrix by tip order ---
  tip_order <- p_tree$data$label[p_tree$data$isTip]
  mat2 <- mat[match(tip_order, rownames(mat)), , drop = FALSE]

  # --- 7. Convert values to factors (0/1) ---
  mat2[] <- lapply(as.data.frame(mat2), function(x) factor(x, levels = c(0, 1)))

  # --- 8. Heatmap + discrete fill ---
  p_out <- ggtree::gheatmap(p_tree, mat2,
                            width = 0.5,
                            colnames = TRUE,
                            colnames_angle = 90,
                            colnames_offset_y = 0.4,
                            hjust = 1,
                            font.size = label_size) +
    ggplot2::scale_fill_manual(values = fill_colors, na.value = "white") +
    ggplot2::theme(legend.position = "right")

  # --- 9. Add tip labels ---
  p_out <- p_out +
    ggtree::geom_tiplab(size = label_size, offset = label_offset)

  return(p_out)
}
