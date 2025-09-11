#' Plot Phyloseq Feature Metrics
#'
#' This function takes a phyloseq object and a feature data frame and
#' generates a multi-panel plot showing Weighted score, other metrics,
#' and group classification.
#'
#' @param ps A phyloseq object.
#' @param dat A data frame with feature information. Required columns:
#'   "OTU", "Rank", "Weighted_score", "Mean_abundance", "Prevalence",
#'   "padj", "log2FoldChange", "Importance", "Degree", "Betweenness".
#' @param top_n Integer. Number of top features to display (default 50).
#' @return A patchwork plot object.
#' @export
plot_physeq_features <- function(ps, dat, top_n = 50) {

  # ---- 0. Input validation ----
  if (!inherits(ps, "phyloseq")) stop("ps must be a phyloseq object")
  required_cols <- c("OTU","Rank","Weighted_score","Mean_abundance",
                     "Prevalence","padj","log2FoldChange",
                     "Importance","Degree","Betweenness")
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) stop("dat is missing columns: ", paste(missing_cols, collapse = ", "))

  # ---- 1. Clean OTU names in dat ----
  dat <- dplyr::mutate(dat,
                       OTU = gsub("^[fgs]__+", "", OTU),
                       OTU = ifelse(OTU == "", "Unknown.plot", OTU))

  # ---- 2. Compute group mean from phyloseq object ----
  tem <- ps %>% ps.group.mean() %>% dplyr::arrange(dplyr::desc(mean))

  # ---- 3. Clean OTU names in tem ----
  tem <- dplyr::mutate(tem,
                       ID = gsub("^[fgs]__+", "", ID),
                       ID = ifelse(ID == "", "Unknown", ID))
  head(tem)

  # Initialize vector
  A <- character(nrow(tem))

  # Identify the sample columns (exclude 'mean' and 'ID')
  sample_cols <- setdiff(colnames(tem), c("mean", "ID", "group.class"))

  # For each row, assign the group with the highest abundance
  for (j in seq_len(nrow(tem))) {
    row_values <- tem[j, sample_cols]
    max_col <- sample_cols[which.max(as.numeric(row_values))]
    # If all values are NA, assign "Unknown"
    if (all(is.na(row_values))) {
      A[j] <- "Unknown"
    } else {
      A[j] <- max_col
    }
  }

  tem$group.class <- A

  # ---- 5. Join with original dat ----
  tem$ID2 = row.names(tem)
  dat2 <- dplyr::left_join(dat, tem, by = c("Taxonomy" = "ID2"))

  # ---- 6. Select top N features ----
  df_top <- dplyr::arrange(dat2, Rank) %>% dplyr::slice_head(n = top_n)

  # ---- 7. Ordering for plotting ----
  df_ord <- dplyr::mutate(df_top,
                          sign_key = ifelse(log2FoldChange >= 0, 1L, 0L)) %>%
    dplyr::arrange(dplyr::desc(sign_key), dplyr::desc(Weighted_score))

  otu_levels <- df_ord$OTU
  df_ord <- dplyr::mutate(df_ord, OTU = factor(OTU, levels = rev(otu_levels)))

  # ---- 8. Convert to long format ----
  df_long <- df_ord %>%
    dplyr::transmute(
      OTU,
      log2FoldChange,
      Weighted       = Weighted_score,
      Mean_abundance = Mean_abundance,
      Prevalence     = Prevalence,
      Significance   = -log10(padj),
      log2FC         = log2FoldChange,
      Importance     = Importance,
      Degree         = Degree,
      Betweenness    = Betweenness
    ) %>%
    tidyr::pivot_longer(cols = -c(OTU, log2FoldChange),
                        names_to = "metric", values_to = "score")

  metric_levels <- c("Weighted","Mean_abundance","Prevalence",
                     "Significance","log2FC","Importance","Degree","Betweenness")
  df_long$metric <- factor(df_long$metric, levels = metric_levels)

  # ---- 9. Format labels ----
  fmt_sci <- scales::label_scientific(digits = 2)
  df_long <- dplyr::mutate(df_long,
                           label = dplyr::case_when(
                             metric == "Weighted"        ~ scales::number(score, accuracy = 0.001),
                             metric == "Mean_abundance"  ~ fmt_sci(score),
                             metric == "Prevalence"      ~ scales::percent(score, accuracy = 1),
                             metric == "Significance"    ~ scales::number(score, accuracy = 0.1),
                             metric == "log2FC"          ~ scales::number(score, accuracy = 0.01),
                             metric == "Importance"      ~ scales::number(score, accuracy = 0.1),
                             metric == "Degree"          ~ scales::number(score, accuracy = 1),
                             metric == "Betweenness"     ~ scales::number(score, accuracy = 0.1),
                             TRUE                        ~ as.character(round(score, 3))
                           ))

  fc_max <- max(abs(df_ord$log2FoldChange), na.rm = TRUE)

  # ---- 10. Left panel ----
  p_left <- ggplot2::ggplot(dplyr::filter(df_long, metric == "Weighted"),
                            ggplot2::aes(x = score, y = OTU, fill = log2FoldChange)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.1, size = 2.7) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-fc_max, fc_max),
                                  name = "log2FC\n(blue negative, red positive)") +
    ggplot2::labs(x = "Weighted score", y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank(),
                   legend.position    = "right")

  # ---- 11. Right panel ----
  p_right <- ggplot2::ggplot(dplyr::filter(df_long, metric != "Weighted"),
                             ggplot2::aes(x = score, y = OTU, fill = log2FoldChange)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = -0.1, size = 2.5) +
    ggplot2::facet_grid(cols = ggplot2::vars(metric), scales = "free_x", switch = "x") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                                  midpoint = 0, limits = c(-fc_max, fc_max),
                                  name = "log2FC\n(blue negative, red positive)") +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank(),
                   strip.placement    = "outside",
                   strip.text.x       = ggplot2::element_text(face = "bold"),
                   legend.position    = "right")

  # ---- 12. Group panel ----
  df_group <- dplyr::transmute(df_ord, OTU, group.class, score = 1)
  p_group <- ggplot2::ggplot(df_group, ggplot2::aes(x = score, y = OTU, fill = group.class)) +
    ggplot2::geom_col(width = 0.8, show.legend = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = group.class), hjust = -0.1, size = 2.5) +
    ggplot2::scale_x_continuous(limits = c(0, 1.2), expand = ggplot2::expansion(mult = c(0.01, 0.12))) +
    ggplot2::labs(x = "Group", y = NULL) +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor   = ggplot2::element_blank())

  # ---- 13. Combine panels with patchwork ----
  p_combined <- patchwork::wrap_plots(
    p_left, p_right, p_group,
    ncol = 3,
    widths = c(1, 6, 1),
    guides = "collect"
  )

  return(p_combined)
}

