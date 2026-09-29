#' Plot multi-level metabolite taxonomy alluvial diagrams
#'
#' Generate group-wise Sankey/alluvial plots for HMDB-style metabolite
#' taxonomy levels (Kingdom, Super_Class, Class and Sub_Class) from a
#' phyloseq object.
#'
#' @param ps A phyloseq object.
#' @param group Character. Column name in sample_data(ps) used for grouping.
#' @param Top Integer. Number of top taxa/metabolite classes retained at each rank.
#' @param col.g Optional named character vector of colors.
#' @param width Numeric. Width of alluvial strata/flows.
#' @param label_threshold Numeric. Minimum mean relative abundance (%) for labels.
#'
#' @return A list containing:
#' \describe{
#'   \item{plot_list}{A named list of ggplot objects, one for each group.}
#'   \item{summary}{Overall abundance summary.}
#'   \item{summary_by_group}{Group-wise abundance summary.}
#'   \item{colors}{Named color vector used in the plots.}
#'   \item{raw_data}{Long-format abundance table.}
#'   \item{flow_data_list}{Alluvial plotting data for each group.}
#' }
#'
#' @export
barMultiLevel <- function(ps = NULL,
                          group = "Group",
                          Top = 10,
                          col.g = NULL,
                          width = 0.6,
                          label_threshold = 1) {

  # -----------------------------
  # 0. Basic checks
  # -----------------------------
  if (is.null(ps)) {
    stop("`ps` must be a phyloseq object.", call. = FALSE)
  }

  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must inherit from class 'phyloseq'.", call. = FALSE)
  }

  if (!is.character(group) || length(group) != 1L) {
    stop("`group` must be a single character string.", call. = FALSE)
  }

  if (!is.numeric(Top) || length(Top) != 1L || Top < 1) {
    stop("`Top` must be a positive integer.", call. = FALSE)
  }
  Top <- as.integer(Top)

  # -----------------------------
  # 1. Taxonomy cleaning
  # -----------------------------
  tax_mat <- as.matrix(phyloseq::tax_table(ps))
  tax_mat[is.na(tax_mat)] <- "Unknown"
  tax_mat[tax_mat == ""] <- "Unknown"
  phyloseq::tax_table(ps) <- phyloseq::tax_table(tax_mat)

  tax <- as.data.frame(
    phyloseq::tax_table(ps),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # HMDB taxonomy columns accepted by the function
  required_cols <- c("Kingdom", "Super_Class", "Class", "Sub_Class")

  col_mapping <- list(
    Kingdom     = c("Kingdom", "kingdom"),
    Super_Class = c("Super_Class", "Super_class", "super_class"),
    Class       = c("Class", "class"),
    Sub_Class   = c("Sub_Class", "Sub_class", "sub_class")
  )

  # Normalize taxonomy column names
  for (target_name in names(col_mapping)) {
    possible_names <- col_mapping[[target_name]]
    found_name <- possible_names[possible_names %in% colnames(tax)]

    if (length(found_name) > 0L) {
      found_name <- found_name[[1L]]
      if (!identical(found_name, target_name)) {
        colnames(tax)[colnames(tax) == found_name] <- target_name
      }
    }
  }

  missing_cols <- setdiff(required_cols, colnames(tax))

  if (length(missing_cols) > 0L) {
    stop(
      "Missing required taxonomy columns: ",
      paste(missing_cols, collapse = ", "),
      "\nAvailable columns: ",
      paste(colnames(tax), collapse = ", "),
      "\n\nPlease ensure your data has been annotated with HMDB taxonomy.",
      call. = FALSE
    )
  }

  phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tax))

  # -----------------------------
  # 2. Check grouping variable
  # -----------------------------
  map <- data.frame(
    phyloseq::sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  if (!group %in% colnames(map)) {
    stop(
      "Grouping variable `", group, "` was not found in sample_data(ps).\n",
      "Available columns: ", paste(colnames(map), collapse = ", "),
      call. = FALSE
    )
  }

  # -----------------------------
  # 3. Transform to relative abundance (%)
  # -----------------------------
  ps <- phyloseq::transform_sample_counts(
    ps,
    function(x) {
      s <- sum(x, na.rm = TRUE)
      if (s == 0) {
        return(rep(0, length(x)))
      }
      x / s * 100
    }
  )

  # -----------------------------
  # 4. Extract abundance/taxonomy/metadata
  # -----------------------------
  otu <- as.data.frame(
    phyloseq::otu_table(ps),
    check.names = FALSE
  )

  if (!phyloseq::taxa_are_rows(ps)) {
    otu <- as.data.frame(t(otu), check.names = FALSE)
  }

  tax <- as.data.frame(
    phyloseq::tax_table(ps),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  map <- data.frame(
    phyloseq::sample_data(ps),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  map$Sample <- rownames(map)

  otu$OTU <- rownames(otu)
  tax$OTU <- rownames(tax)

  # -----------------------------
  # 5. Build long table
  # -----------------------------
  otu_long <- tidyr::pivot_longer(
    otu,
    cols = -dplyr::all_of("OTU"),
    names_to = "Sample",
    values_to = "Abundance"
  )

  otu_long <- dplyr::left_join(otu_long, tax, by = "OTU")
  otu_long <- dplyr::left_join(
    otu_long,
    map[, c("Sample", group), drop = FALSE],
    by = "Sample"
  )

  colnames(otu_long)[colnames(otu_long) == group] <- "Group"

  ranks_list <- c("Kingdom", "Super_Class", "Class", "Sub_Class")

  # -----------------------------
  # 6. Select Top N at each rank
  # -----------------------------
  get_top_taxa <- function(data, rank, Top) {
    tmp <- dplyr::group_by(data, rlang::.data[[rank]])
    tmp <- dplyr::summarise(
      tmp,
      Total = sum(rlang::.data$Abundance, na.rm = TRUE),
      .groups = "drop"
    )
    tmp <- dplyr::arrange(tmp, dplyr::desc(rlang::.data$Total))
    tmp <- dplyr::filter(
      tmp,
      !is.na(rlang::.data[[rank]]),
      rlang::.data[[rank]] != "Unknown"
    )
    tmp <- dplyr::slice_head(tmp, n = Top)

    tmp[[rank]]
  }

  for (rank in ranks_list) {
    top_taxa <- get_top_taxa(otu_long, rank, Top)
    new_col <- paste0(rank, "_plot")

    otu_long[[new_col]] <- ifelse(
      otu_long[[rank]] %in% top_taxa,
      otu_long[[rank]],
      "Others"
    )
  }

  groups <- unique(otu_long$Group)
  groups <- groups[!is.na(groups)]

  if (length(groups) == 0L) {
    stop("No valid groups were found in sample_data(ps).", call. = FALSE)
  }

  # -----------------------------
  # 7. Global color assignment
  # -----------------------------
  all_taxa <- character(0)

  for (rank in ranks_list) {
    col_name <- paste0(rank, "_plot")
    all_taxa <- c(all_taxa, unique(otu_long[[col_name]]))
  }

  all_taxa <- unique(all_taxa)
  all_taxa <- c(setdiff(all_taxa, "Others"), "Others")

  if (is.null(col.g)) {
    kingdom_taxa <- setdiff(unique(otu_long$Kingdom_plot), "Others")
    n_kingdom <- length(kingdom_taxa)

    super_taxa <- setdiff(
      unique(otu_long$Super_Class_plot),
      c("Others", kingdom_taxa)
    )

    class_taxa <- setdiff(
      unique(otu_long$Class_plot),
      c("Others", kingdom_taxa, super_taxa)
    )

    sub_taxa <- setdiff(
      unique(otu_long$Sub_Class_plot),
      c("Others", kingdom_taxa, super_taxa, class_taxa)
    )

    if (n_kingdom > 0L) {
      if (n_kingdom <= 10L) {
        kingdom_colors <- ggsci::pal_jco()(n_kingdom)
      } else {
        kingdom_colors <- c(
          ggsci::pal_jco()(10),
          grDevices::colorRampPalette(
            c("#0073C2FF", "#EFC000FF", "#868686FF")
          )(n_kingdom - 10L)
        )
      }
      names(kingdom_colors) <- kingdom_taxa
    } else {
      kingdom_colors <- character(0)
    }

    if (length(super_taxa) > 0L) {
      super_colors <- grDevices::colorRampPalette(
        c("#E64B35", "#4DBBD5", "#00A087")
      )(length(super_taxa))
      names(super_colors) <- super_taxa
    } else {
      super_colors <- character(0)
    }

    if (length(class_taxa) > 0L) {
      class_colors <- grDevices::colorRampPalette(
        c("#F39B7F", "#8491B4", "#91D1C2")
      )(length(class_taxa))
      names(class_colors) <- class_taxa
    } else {
      class_colors <- character(0)
    }

    if (length(sub_taxa) > 0L) {
      sub_colors <- grDevices::colorRampPalette(
        c("#DC0000", "#7E6148", "#B09C85")
      )(length(sub_taxa))
      names(sub_colors) <- sub_taxa
    } else {
      sub_colors <- character(0)
    }

    col.g <- c(
      kingdom_colors,
      super_colors,
      class_colors,
      sub_colors,
      Others = "grey80",
      Unknown = "grey90"
    )
  } else {
    if (is.null(names(col.g))) {
      stop("`col.g` must be a named character vector.", call. = FALSE)
    }

    missing_colors <- setdiff(all_taxa, names(col.g))
    if (length(missing_colors) > 0L) {
      warning(
        "No color supplied for: ",
        paste(missing_colors, collapse = ", "),
        ". These values will use grey50."
      )
    }
  }

  # -----------------------------
  # 8. Generate plot for each group
  # -----------------------------
  plot_list <- list()
  flow_data_list <- list()

  for (g in groups) {

    group_otu_long <- dplyr::filter(
      otu_long,
      rlang::.data$Group == g
    )

    group_data <- dplyr::group_by(
      group_otu_long,
      rlang::.data$Kingdom_plot,
      rlang::.data$Super_Class_plot,
      rlang::.data$Class_plot,
      rlang::.data$Sub_Class_plot
    )

    group_data <- dplyr::summarise(
      group_data,
      Abundance = sum(rlang::.data$Abundance, na.rm = TRUE),
      .groups = "drop"
    )

    n_samples <- length(unique(group_otu_long$Sample))

    if (n_samples == 0L) {
      next
    }

    group_data$Abundance <- group_data$Abundance / n_samples
    group_data$flow_id <- seq_len(nrow(group_data))
    group_data$Kingdom_for_color <- group_data$Kingdom_plot

    flow_data <- tidyr::pivot_longer(
      group_data,
      cols = dplyr::all_of(
        c("Kingdom_plot", "Super_Class_plot", "Class_plot", "Sub_Class_plot")
      ),
      names_to = "Rank",
      values_to = "Taxon"
    )

    flow_data$Rank <- sub("_plot$", "", flow_data$Rank)
    flow_data$Rank <- factor(flow_data$Rank, levels = ranks_list)

    kingdom_levels <- intersect(
      names(col.g),
      unique(group_data$Kingdom_plot)
    )

    flow_data$Kingdom_color <- factor(
      flow_data$Kingdom_for_color,
      levels = kingdom_levels
    )

    # Preserve named color order where possible, but keep unexpected levels too
    taxon_levels <- unique(c(names(col.g), as.character(flow_data$Taxon)))
    flow_data$Taxon <- factor(flow_data$Taxon, levels = taxon_levels)

    stratum_abundance <- dplyr::group_by(
      flow_data,
      rlang::.data$Rank,
      rlang::.data$Taxon
    )

    stratum_abundance <- dplyr::summarise(
      stratum_abundance,
      Total_Abundance = sum(rlang::.data$Abundance, na.rm = TRUE),
      .groups = "drop"
    )

    flow_data <- dplyr::left_join(
      flow_data,
      stratum_abundance,
      by = c("Rank", "Taxon")
    )

    flow_data_list[[as.character(g)]] <- flow_data

    p <- ggplot2::ggplot(
      flow_data,
      ggplot2::aes(
        x = rlang::.data$Rank,
        y = rlang::.data$Abundance,
        alluvium = rlang::.data$flow_id,
        stratum = rlang::.data$Taxon
      )
    ) +
      ggalluvial::geom_flow(
        ggplot2::aes(fill = rlang::.data$Kingdom_color),
        width = width,
        alpha = 0.6,
        curve_type = "cubic"
      ) +
      ggalluvial::geom_stratum(
        ggplot2::aes(fill = rlang::.data$Taxon),
        width = width,
        color = "white",
        linewidth = 0.3
      ) +
      ggplot2::geom_text(
        stat = "stratum",
        ggplot2::aes(
          label = ifelse(
            rlang::.data$Total_Abundance > label_threshold,
            as.character(rlang::.data$Taxon),
            ""
          )
        ),
        size = 3,
        fontface = "bold"
      ) +
      ggplot2::scale_fill_manual(
        values = col.g,
        na.value = "grey50"
      ) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::scale_x_discrete(
        expand = c(0.05, 0.05),
        labels = c(
          Kingdom = "Kingdom",
          Super_Class = "Super Class",
          Class = "Class",
          Sub_Class = "Sub Class"
        )
      ) +
      ggplot2::labs(
        x = "",
        y = "Mean Relative Abundance (%)",
        title = paste0("Metabolite Hierarchy - ", g)
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          face = "bold",
          size = 16,
          color = "#2C3E50"
        ),
        axis.text.x = ggplot2::element_text(
          size = 12,
          face = "bold",
          color = "black",
          angle = 0
        ),
        axis.text.y = ggplot2::element_text(size = 10),
        axis.title.y = ggplot2::element_text(size = 12, face = "bold"),
        legend.position = "none",
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.background = ggplot2::element_rect(fill = "white", color = NA),
        plot.margin = ggplot2::margin(t = 20, r = 10, b = 10, l = 10)
      )

    plot_list[[as.character(g)]] <- p
  }

  # -----------------------------
  # 9. Summary tables
  # -----------------------------
  summary_by_group <- dplyr::select(
    otu_long,
    dplyr::all_of(
      c(
        "Sample", "Group",
        "Kingdom_plot", "Super_Class_plot",
        "Class_plot", "Sub_Class_plot",
        "Abundance"
      )
    )
  )

  summary_by_group <- tidyr::pivot_longer(
    summary_by_group,
    cols = dplyr::all_of(
      c("Kingdom_plot", "Super_Class_plot", "Class_plot", "Sub_Class_plot")
    ),
    names_to = "Rank",
    values_to = "Taxon"
  )

  summary_by_group$Rank <- sub("_plot$", "", summary_by_group$Rank)

  summary_by_group <- dplyr::group_by(
    summary_by_group,
    rlang::.data$Group,
    rlang::.data$Rank,
    rlang::.data$Taxon
  )

  summary_by_group <- dplyr::summarise(
    summary_by_group,
    Mean_Abundance = mean(rlang::.data$Abundance, na.rm = TRUE),
    .groups = "drop"
  )

  summary_by_group <- dplyr::arrange(
    summary_by_group,
    rlang::.data$Group,
    rlang::.data$Rank,
    dplyr::desc(rlang::.data$Mean_Abundance)
  )

  summary_total <- dplyr::select(
    otu_long,
    dplyr::all_of(
      c(
        "Sample",
        "Kingdom_plot", "Super_Class_plot",
        "Class_plot", "Sub_Class_plot",
        "Abundance"
      )
    )
  )

  summary_total <- tidyr::pivot_longer(
    summary_total,
    cols = dplyr::all_of(
      c("Kingdom_plot", "Super_Class_plot", "Class_plot", "Sub_Class_plot")
    ),
    names_to = "Rank",
    values_to = "Taxon"
  )

  summary_total$Rank <- sub("_plot$", "", summary_total$Rank)

  summary_total <- dplyr::group_by(
    summary_total,
    rlang::.data$Rank,
    rlang::.data$Taxon
  )

  summary_total <- dplyr::summarise(
    summary_total,
    Total_Abundance = sum(rlang::.data$Abundance, na.rm = TRUE),
    Mean_Abundance = mean(rlang::.data$Abundance, na.rm = TRUE),
    .groups = "drop"
  )

  summary_total <- dplyr::arrange(
    summary_total,
    rlang::.data$Rank,
    dplyr::desc(rlang::.data$Total_Abundance)
  )

  # -----------------------------
  # 10. Return
  # -----------------------------
  list(
    plot_list = plot_list,
    summary = summary_total,
    summary_by_group = summary_by_group,
    colors = col.g,
    raw_data = otu_long,
    flow_data_list = flow_data_list
  )
}
