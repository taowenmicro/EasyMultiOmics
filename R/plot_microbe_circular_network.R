
#' Plot Circular Correlation Network for a Microbe
#'
#' This function visualizes correlations between a target microbe and metabolites
#' in a circular dendrogram-style network. Positive and negative correlations
#' are colored differently, and the microbe and group nodes are visually distinct.
#'
#' @param dat A data.frame containing at least columns:
#'   metabolite, spearman_cor, integrated_score
#' @param microbe_name Character, name of the target microbe (default: "Target_microbe")
#' @param thr Numeric, threshold for absolute correlation to include metabolite (default: 0.0)
#' @param label_offset Numeric, radial offset for metabolite labels (default: 0.05)
#' @return A ggraph ggplot object
#' @export
#' @examples
#' \dontrun{
#' p <- plot_microbe_circular_network(dat, "Microbe_X", thr = 0.1)
#' print(p)
#' }




plot_microbe_circular_network <- function(dat, microbe_name = "Target_microbe",
                                          thr = 0.0, label_offset = 0.05) {
  # ---- 0. Prepare metabolite correlation data ----
  id <- ifelse(dat$spearman_cor >= 0, 1, -1)
  dat1 <- tibble::tibble(
    metabolite = dat$metabolite,
    group      = ifelse(dat$spearman_cor >= 0, "pos", "neg"),
    r          = dat$integrated_score * id
  )

  dat_use <- dat1[is.finite(dat1$r) & abs(dat1$r) >= thr, , drop = FALSE]
  if (nrow(dat_use) == 0) stop("No metabolites pass the threshold.")

  # ---- 1. Define nodes (microbe, groups, metabolites) ----
  groups <- unique(dat_use$group)
  nodes <- tibble::tibble(
    name = c(microbe_name, groups, dat_use$metabolite),
    type = c("microbe", rep("group", length(groups)), rep("metabolite", nrow(dat_use))),
    group = c(NA, groups, dat_use$group)
  )

  # ---- 2. Define edges ----
  edges <- rbind(
    data.frame(from = microbe_name, to = groups, r = NA_real_),
    data.frame(from = dat_use$group, to = dat_use$metabolite, r = dat_use$r)
  )

  # ---- 3. Construct tbl_graph ----
  g <- tidygraph::tbl_graph(nodes = nodes, edges = edges, directed = TRUE) %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(
      sign  = dplyr::case_when(is.na(r) ~ "struct", r >= 0 ~ "pos", TRUE ~ "neg"),
      alpha = ifelse(sign == "struct", 0.35, 0.65)
    )

  # ---- 4. Define colors ----
  col_pos <- "#C43E2F"
  col_neg <- "#2B6CB0"
  col_struct <- "grey80"

  # ---- 5. Build ggraph plot ----
  p <- ggraph::ggraph(g, layout = "dendrogram", circular = TRUE) +
    # edges
    ggraph::geom_edge_diagonal(aes(edge_colour = sign, edge_alpha = alpha),
                               edge_width = 0.28, lineend = "round") +
    ggraph::scale_edge_colour_manual(
      values = c(pos = col_pos, neg = col_neg, struct = col_struct),
      breaks = c("pos", "neg"),
      labels = c("Positive", "Negative"),
      name   = "Correlation"
    ) +
    ggplot2::guides(edge_alpha = "none") +

    # nodes
    ggraph::geom_node_point(aes(shape = type, fill = type),
                            size = 2, colour = "white", stroke = 0.5) +
    ggplot2::scale_shape_manual(values = c(microbe = 21, group = 21, metabolite = 21),
                                guide = "none") +
    ggplot2::scale_fill_manual(values = c(microbe = "#8C564B", group = "#EADBC8", metabolite = "#C43E2F"),
                               guide = "none") +

    # metabolite labels with radial offset
    ggraph::geom_node_text(
      data = function(n) {
        n %>%
          dplyr::filter(type == "metabolite") %>%
          dplyr::mutate(
            ang    = atan2(y, x) * 180 / pi,
            flip   = ang > 90 & ang < 270,
            lab_ang = ifelse(flip, ang - 180, ang),
            lab_h   = ifelse(flip, 1, 0),
            r       = sqrt(x^2 + y^2),
            ux      = ifelse(r > 0, x / r, 0),
            uy      = ifelse(r > 0, y / r, 0),
            x_lab   = x + ux * label_offset,
            y_lab   = y + uy * label_offset
          )
      },
      aes(x = x_lab, y = y_lab, label = name, angle = lab_ang, hjust = lab_h),
      size = 2.2, family = "sans", check_overlap = TRUE
    ) +

    # group labels (center, invisible)
    ggraph::geom_node_label(data = function(n) dplyr::filter(n, type == "group"),
                            aes(label = name),
                            label.size = 0, size = 3, fill = "white", alpha = 0) +

    # microbe label (center)
    ggraph::geom_node_label(data = function(n) dplyr::filter(n, type == "microbe"),
                            aes(label = name),
                            label.size = 0, size = 3.6, fill = "white", alpha = 0) +

    # theme
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.margin = grid::unit(c(60, 60, 60, 60), "pt")
    )

  return(p)
}
