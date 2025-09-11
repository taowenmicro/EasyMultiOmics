#' Plot Circular Correlation Network for a Microbe with Enhanced Visualization
#'
#' This function visualizes correlations between a target microbe and metabolites
#' in a circular dendrogram-style network. Node colors represent positive/negative
#' correlations, node sizes are scaled by integrated_score, and an outer ring of
#' circular bar plots shows the integrated_score values.
#'
#' @param dat A data.frame containing at least columns:
#'   metabolite, spearman_cor, integrated_score
#' @param microbe_name Character, name of the target microbe (default: "Target_microbe")
#' @param thr Numeric, threshold for absolute correlation to include metabolite (default: 0.0)
#' @param label_offset Numeric, radial offset for metabolite labels (default: 0.05)
#' @param bar_offset Numeric, radial offset for bar plots (default: 0.15)
#' @param bar_width Numeric, width of bars in radians (default: 0.08)
#' @param size_range Vector of length 2, range for node sizes (default: c(1, 6))
#' @return A ggplot object with combined network and bar plots
#' @export
#' @examples
#' \dontrun{
#' p <- plot_microbe_circular_network2(dat, "Microbe_X", thr = 0.1)
#' print(p)
#' }

plot_microbe_circular_network2 <- function(dat, microbe_name = "Target_microbe",
                                          thr = 0.0, label_offset = 0.05,
                                          bar_offset = 0.15, bar_width = 0.08,
                                          size_range = c(1, 6)) {

  # Load required libraries
  if (!requireNamespace("ggraph", quietly = TRUE)) stop("Package 'ggraph' is required")
  if (!requireNamespace("tidygraph", quietly = TRUE)) stop("Package 'tidygraph' is required")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required")
  if (!requireNamespace("tibble", quietly = TRUE)) stop("Package 'tibble' is required")

  # ---- 0. Prepare metabolite correlation data ----
  id <- ifelse(dat$spearman_cor >= 0, 1, -1)
  dat1 <- tibble::tibble(
    metabolite = dat$metabolite,
    group      = ifelse(dat$spearman_cor >= 0, "pos", "neg"),
    r          = dat$integrated_score * id,
    integrated_score = dat$integrated_score,
    spearman_cor = dat$spearman_cor
  )

  dat_use <- dat1[is.finite(dat1$r) & abs(dat1$r) >= thr, , drop = FALSE]
  if (nrow(dat_use) == 0) stop("No metabolites pass the threshold.")

  # Sort by integrated_score for consistent ordering
  dat_use <- dat_use[order(-dat_use$integrated_score), ]

  # ---- 1. Define nodes (microbe, groups, metabolites) ----
  groups <- unique(dat_use$group)
  nodes <- tibble::tibble(
    name = c(microbe_name, groups, dat_use$metabolite),
    type = c("microbe", rep("group", length(groups)), rep("metabolite", nrow(dat_use))),
    group = c("microbe", groups, dat_use$group),
    integrated_score = c(max(dat_use$integrated_score, na.rm = TRUE),
                         rep(mean(dat_use$integrated_score, na.rm = TRUE), length(groups)),
                         dat_use$integrated_score),
    spearman_cor = c(0, rep(0, length(groups)), dat_use$spearman_cor)
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
  col_microbe <- "#8C564B"
  col_group <- "#EADBC8"

  # ---- 5. Create layout and extract coordinates ----
  layout_coords <- ggraph::create_layout(g, layout = "dendrogram", circular = TRUE)

  # Extract metabolite positions for bar plots
  metabolite_coords <- layout_coords[layout_coords$type == "metabolite", ]

  # Calculate angles and positions for bars
  metabolite_coords$angle <- atan2(metabolite_coords$y, metabolite_coords$x)
  metabolite_coords$radius <- sqrt(metabolite_coords$x^2 + metabolite_coords$y^2)

  # Normalize integrated_score for bar height
  max_score <- max(metabolite_coords$integrated_score, na.rm = TRUE)
  min_score <- min(metabolite_coords$integrated_score, na.rm = TRUE)
  score_range <- max_score - min_score
  if (score_range == 0) score_range <- 1  # Handle case where all scores are equal

  metabolite_coords$bar_height <- (metabolite_coords$integrated_score - min_score) /
    score_range * 0.12 + 0.02  # Scale to 0.02-0.14

  # Create circular bar plot data
  bar_data <- do.call(rbind, lapply(1:nrow(metabolite_coords), function(i) {
    coord <- metabolite_coords[i, ]
    base_radius <- max(layout_coords$radius, na.rm = TRUE) + bar_offset

    # Create sector for bar (circular bar)
    n_points <- 30
    angles <- seq(coord$angle - bar_width/2, coord$angle + bar_width/2, length.out = n_points)

    # Inner edge (base of bar)
    inner_x <- base_radius * cos(angles)
    inner_y <- base_radius * sin(angles)

    # Outer edge (top of bar)
    outer_radius <- base_radius + coord$bar_height
    outer_x <- outer_radius * cos(angles)
    outer_y <- outer_radius * sin(angles)

    # Create polygon data for the bar
    x_coords <- c(inner_x, rev(outer_x))
    y_coords <- c(inner_y, rev(outer_y))

    data.frame(
      x = x_coords,
      y = y_coords,
      metabolite = coord$name,
      group = coord$group,
      integrated_score = coord$integrated_score,
      polygon_id = i
    )
  }))

  # ---- 6. Build ggraph plot ----
  p <- ggraph::ggraph(layout_coords) +
    # Add circular bar plots first (background layer)
    ggplot2::geom_polygon(data = bar_data,
                          aes(x = x, y = y, group = polygon_id, fill = group),
                          alpha = 0.7, color = "white", linewidth = 0.3) +

    # edges
    ggraph::geom_edge_diagonal(aes(edge_colour = sign, edge_alpha = alpha),
                               edge_width = 0.28, lineend = "round") +
    ggraph::scale_edge_colour_manual(
      values = c(pos = col_pos, neg = col_neg, struct = col_struct),
      breaks = c("pos", "neg"),
      labels = c("Positive", "Negative"),
      name   = "Edge Type"
    ) +
    ggplot2::guides(edge_alpha = "none") +

    # nodes with size mapped to integrated_score and color to correlation type
    ggraph::geom_node_point(aes(size = integrated_score, fill = group),
                            shape = 21, colour = "white", stroke = 0.6) +
    ggplot2::scale_size_continuous(range = size_range,
                                   name = "Integrated\nScore",
                                   guide = guide_legend(
                                     override.aes = list(shape = 21, stroke = 0.6, colour = "white"),
                                     order = 2)) +

    # Color scale for both nodes and bars
    ggplot2::scale_fill_manual(
      values = c(microbe = col_microbe,
                 group = col_group,
                 pos = col_pos,
                 neg = col_neg),
      breaks = c("pos", "neg"),
      labels = c("Positive correlation", "Negative correlation"),
      name = "Correlation Type",
      guide = guide_legend(
        override.aes = list(shape = 21, size = 4, stroke = 0.6, colour = "white"),
        order = 1)
    ) +

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
      size = 2.2, family = "sans", check_overlap = TRUE, fontface = "bold"
    ) +

    # group labels (center, invisible for structure)
    ggraph::geom_node_label(data = function(n) dplyr::filter(n, type == "group"),
                            aes(label = name),
                            label.size = 0, size = 3, fill = "white", alpha = 0) +

    # microbe label (center)
    ggraph::geom_node_label(data = function(n) dplyr::filter(n, type == "microbe"),
                            aes(label = name),
                            label.size = 0, size = 3.6, fill = "white", alpha = 0,
                            fontface = "bold") +

    # Add annotation for bar plots
    ggplot2::annotate("text",
                      x = 0,
                      y = -(max(layout_coords$radius, na.rm = TRUE) + bar_offset + 0.08),
                      label = "Integrated Score",
                      size = 2.8, angle = 0, hjust = 0.5, vjust = 1,
                      color = "grey40", fontface = "italic") +

    # theme
    ggplot2::coord_equal(clip = "off") +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      plot.margin = grid::unit(c(90, 90, 90, 90), "pt"),
      legend.box = "vertical",
      legend.spacing.y = unit(0.5, "cm"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey40")
    )

  return(p)
}
