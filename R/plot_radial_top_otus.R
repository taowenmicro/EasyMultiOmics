#' Plot Radial Bar Chart of Top OTUs by Weighted Score
#'
#' This function creates a circular bar plot (radial) for the top OTUs
#' based on their Weighted_score, with OTU names cleaned from f__, g__, s__ prefixes.
#'
#' @param dat A data.frame containing at least columns: Rank, OTU, Weighted_score
#' @param top_n Integer, number of top OTUs to display (default 50)
#' @param inner_ratio Numeric, inner radius proportion (0 = solid, 0.7 = thin ring)
#' @param bar_width Numeric, width of bars (default 0.96)
#' @param text_size Numeric, size of OTU text labels (default 2)
#' @return A ggplot2 object of the radial bar chart
#' @export
#' @examples
#' \dontrun{
#' p <- plot_radial_top_otus(dat, top_n = 50, inner_ratio = 0.45)
#' print(p)
#' }




plot_radial_top_otus <- function(dat, top_n = 50, inner_ratio = 0.45,
                                 bar_width = 0.96, text_size = 2) {

  # ---- 0. Parameter checks ----
  if (!all(c("Rank", "OTU", "Weighted_score") %in% colnames(dat))) {
    stop("dat must contain columns: Rank, OTU, Weighted_score")
  }

  if (!is.numeric(top_n) || top_n <= 0) top_n <- 50
  if (!is.numeric(inner_ratio) || inner_ratio < 0 || inner_ratio > 1) inner_ratio <- 0.45

  # ---- 1. Select top N features ----
  df_top <- dat[order(dat$Rank), , drop = FALSE][seq_len(min(top_n, nrow(dat))), ]

  # ---- 2. Clean OTU names ----
  df_top$OTU <- gsub("^[fgs]__+", "", df_top$OTU)    # remove prefixes
  df_top$OTU[df_top$OTU == ""] <- "Unknown"          # replace empty names

  # ---- 3. Order by Weighted_score descending ----
  df_ord <- df_top[order(-df_top$Weighted_score), ]
  df_ord$id <- seq_len(nrow(df_ord))

  # ---- 4. Compute polar angles ----
  n <- nrow(df_ord)
  df_ord$angle_base <- 360 * (df_ord$id - 0.5) / n
  df_ord$angle_rad  <- 90 - df_ord$angle_base
  df_ord$angle_show <- ifelse(df_ord$angle_rad < -90 | df_ord$angle_rad > 90,
                              df_ord$angle_rad + 180,
                              df_ord$angle_rad)

  # ---- 5. Define inner radius and label properties ----
  y_max <- max(df_ord$Weighted_score, na.rm = TRUE)
  inner <- y_max * inner_ratio
  df_ord$score_norm <- scales::rescale(df_ord$Weighted_score, to = c(0,1))
  df_ord$label_y    <- pmax(df_ord$Weighted_score * 0.65,
                            df_ord$Weighted_score - y_max * 0.08)
  df_ord$label_col  <- ifelse(df_ord$score_norm > 0.65, "white", "#2B2B2B")
  df_ord$label_txt  <- paste0(df_ord$OTU, "\n", sprintf("%.3f", df_ord$Weighted_score))

  # ---- 6. Define color palette ----
  warm_pal <- c("#FDE6D8", "#F9C7A9", "#F28C5B", "#D95F0E", "#7F2704")

  # ---- 7. Create ggplot ----
  p <- ggplot2::ggplot(df_ord, ggplot2::aes(x = id, y = Weighted_score, fill = Weighted_score)) +
    ggplot2::geom_col(width = bar_width, color = "white", linewidth = 0.2) +
    ggplot2::coord_polar(theta = "x") +
    ggplot2::ylim(-inner, y_max * 1.05) +
    ggplot2::scale_fill_gradientn(colors = warm_pal, name = "Weighted score") +
    ggplot2::geom_text(
      ggplot2::aes(x = id, y = label_y * 0.8, label = label_txt, angle = angle_show),
      color = df_ord$label_col, lineheight = 0.9,
      hjust = 0.5, vjust = 0.5,
      size = text_size, family = "sans", inherit.aes = FALSE
    ) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      legend.position = "right",
      legend.title    = ggplot2::element_text(size = 10),
      legend.text     = ggplot2::element_text(size = 9),
      plot.margin     = ggplot2::margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
    )

  return(p)
}
