#' @export
barMultiLevel.metm <- function(ps = NULL, group = "Group", Top = 10,
                               col.g = NULL, width = 0.6,
                               label_threshold = 5) {

  library(ggplot2)
  library(ggalluvial)
  library(dplyr)
  library(tidyr)
  library(patchwork)

  ps <- ps %>% phyloseq::transform_sample_counts(function(x) {
    x / sum(x, na.rm = TRUE)
  })

  otu <- as.data.frame(phyloseq::otu_table(ps))
  if (!phyloseq::taxa_are_rows(ps)) {
    otu <- t(otu)
  }
  otu <- as.data.frame(otu)

  tax <- as.data.frame(phyloseq::tax_table(ps))
  map <- data.frame(phyloseq::sample_data(ps), check.names = FALSE)
  map$Sample <- rownames(map)

  otu$OTU <- rownames(otu)
  tax$OTU <- rownames(tax)

  otu_long <- otu %>%
    tidyr::pivot_longer(cols = -OTU, names_to = "Sample", values_to = "Abundance") %>%
    dplyr::left_join(tax, by = "OTU") %>%
    dplyr::left_join(map[, c("Sample", group)], by = "Sample")

  colnames(otu_long)[colnames(otu_long) == group] <- "Group"

  otu_long$Abundance <- otu_long$Abundance * 100

  ranks_list <- c("Genus", "Family", "Class", "Phylum")

  get_top_taxa <- function(data, rank, Top) {
    top <- data %>%
      dplyr::group_by(.data[[rank]]) %>%
      dplyr::summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(desc(Total)) %>%
      dplyr::slice(1:Top) %>%
      dplyr::pull(.data[[rank]])
    return(top)
  }

  for (rank in ranks_list) {
    top_taxa <- get_top_taxa(otu_long, rank, Top)
    new_col <- paste0(rank, "_plot")
    otu_long[[new_col]] <- ifelse(otu_long[[rank]] %in% top_taxa,
                                  otu_long[[rank]], "Others")
  }

  groups <- unique(otu_long$Group)
  plot_list <- list()


  all_taxa <- c()
  for (rank in ranks_list) {
    col_name <- paste0(rank, "_plot")
    all_taxa <- c(all_taxa, unique(otu_long[[col_name]]))
  }
  all_taxa <- unique(all_taxa)
  all_taxa <- c(setdiff(all_taxa, "Others"), "Others")


  phylum_taxa <- unique(otu_long$Phylum_plot)
  phylum_taxa <- phylum_taxa[phylum_taxa != "Others"]
  n_phylum <- length(phylum_taxa)

  if (is.null(col.g)) {
    library(ggsci)


    if (n_phylum <= 10) {
      phylum_base_colors <- c(
        "#E64B35", "#4DBBD5", "#00A087", "#F39B7F", "#8491B4",
        "#91D1C2", "#DC0000", "#FFD700", "#7E6148", "#FF69B4"
      )[1:n_phylum]
    } else {
      colors_npg <- pal_npg()(10)
      colors_d3 <- pal_d3()(10)
      colors_jco <- pal_jco()(10)
      all_colors <- c(colors_npg, colors_d3, colors_jco)
      indices <- round(seq(1, length(all_colors), length.out = n_phylum))
      phylum_base_colors <- all_colors[indices]
    }

    names(phylum_base_colors) <- phylum_taxa


    other_taxa <- setdiff(all_taxa, c(phylum_taxa, "Others"))
    n_other <- length(other_taxa)

    if (n_other > 0) {
      other_colors <- colorRampPalette(c(pal_npg()(10), pal_d3()(10)))(n_other)
      names(other_colors) <- other_taxa
    } else {
      other_colors <- c()
    }


    col.g <- c(phylum_base_colors, other_colors)
    col.g["Others"] <- "grey70"

  } else {

    if (!"Others" %in% names(col.g)) {
      col.g["Others"] <- "grey70"
    }


    missing_taxa <- setdiff(all_taxa, names(col.g))
    if (length(missing_taxa) > 0) {

      library(ggsci)
      extra_colors <- colorRampPalette(c(pal_npg()(10), pal_d3()(10)))(length(missing_taxa))
      names(extra_colors) <- missing_taxa
      col.g <- c(col.g, extra_colors)
    }
  }


  phylum_colors <- col.g[c(phylum_taxa, "Others")]

  for (g in groups) {

    group_data <- otu_long %>%
      dplyr::filter(Group == g) %>%
      dplyr::group_by(Phylum_plot, Class_plot, Family_plot, Genus_plot) %>%
      dplyr::summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")

    n_samples <- otu_long %>%
      dplyr::filter(Group == g) %>%
      dplyr::pull(Sample) %>%
      unique() %>%
      length()

    group_data$Abundance <- group_data$Abundance / n_samples
    group_data$flow_id <- 1:nrow(group_data)


    group_data$Phylum_for_color <- group_data$Phylum_plot

    flow_data <- group_data %>%
      tidyr::pivot_longer(
        cols = c(Genus_plot, Family_plot, Class_plot, Phylum_plot),
        names_to = "Rank",
        values_to = "Taxon"
      ) %>%
      dplyr::mutate(
        Rank = gsub("_plot", "", Rank),
        Rank = factor(Rank, levels = ranks_list)
      )

    flow_data$Phylum_color <- factor(flow_data$Phylum_for_color, levels = names(phylum_colors))
    flow_data$Taxon <- factor(flow_data$Taxon, levels = names(col.g))


    stratum_abundance <- flow_data %>%
      group_by(Rank, Taxon) %>%
      summarise(Total_Abundance = sum(Abundance), .groups = "drop")

    flow_data <- flow_data %>%
      left_join(stratum_abundance, by = c("Rank", "Taxon"))

    p <- ggplot(flow_data, aes(x = Rank, y = Abundance,
                               alluvium = flow_id,
                               stratum = Taxon)) +

      ggalluvial::geom_flow(aes(fill = Phylum_color),
                            width = width,
                            alpha = 0.7,
                            curve_type = "cubic",
                            lode.guidance = "frontback") +

      ggalluvial::geom_stratum(aes(fill = Taxon),
                               width = width,
                               color = "black",
                               size = 0.3) +

      geom_text(stat = "stratum",
                aes(label = ifelse(Total_Abundance > label_threshold,
                                   as.character(Taxon), "")),
                size = 3, fontface = "bold") +

      scale_fill_manual(
        values = c(col.g, phylum_colors),
        na.value = "grey50"
      ) +
      scale_y_continuous(expand = c(0, 0)) +
      scale_x_discrete(expand = c(0.15, 0.15)) +
      labs(x = "", y = "Relative Abundance (%)", title = g) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.title.y = element_text(size = 12, face = "bold"),
        legend.position = "none",
        panel.grid = element_blank()
      )

    plot_list[[g]] <- p
  }


  summary_data <- otu_long %>%
    dplyr::select(Sample, Group, Phylum_plot, Class_plot, Family_plot, Genus_plot, Abundance) %>%
    tidyr::pivot_longer(
      cols = c(Phylum_plot, Class_plot, Family_plot, Genus_plot),
      names_to = "Rank",
      values_to = "Taxon"
    ) %>%
    dplyr::mutate(Rank = gsub("_plot", "", Rank))


  summary_by_group <- summary_data %>%
    dplyr::group_by(Group, Rank, Taxon) %>%
    dplyr::summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    dplyr::arrange(Group, Rank, desc(Abundance))


  return(list(
    plot_list[[1]],
    summary_data,
    plot_list[[2]],
    plot_list[[3]],
    raw_data = otu_long,
    summary = summary_by_group,
    colors = col.g,
    plot_list = plot_list
  ))
}
