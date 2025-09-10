# results <- identify_minimal_change_units(
#   physeq = ps,
#   group_col = "Group",
#   consist_thresh = 0.8,
#   p_thresh = 0.05
# )
# p1 <- plot_change_units(results, physeq = ps)
#
# p3 <- plot_change_units(results, plot_type = "treemap")





#' Identify Minimal Change Units in Microbiome Data
#'
#' This function identifies taxonomic units that show consistent and significant
#' changes between two groups in microbiome data. It performs centered log-ratio (CLR)
#' transformation and identifies the most specific taxonomic level where changes occur.
#'
#' @param physeq A \code{phyloseq} object containing OTU table, taxonomy table, and sample data.
#' @param group_col Character string specifying the column name in sample data that contains
#'   the grouping variable. Must be a binary factor with exactly 2 levels.
#' @param levels Character vector of taxonomic levels to analyze, ordered from most specific
#'   to least specific. Default is \code{c("Species", "Genus", "Family", "Order", "Class", "Phylum")}.
#' @param consist_thresh Numeric value between 0 and 1 specifying the minimum proportion
#'   of sub-units that must change in the same direction for consistency. Default is 0.8.
#' @param p_thresh Numeric value specifying the p-value threshold for statistical significance.
#'   Default is 0.05.
#' @param abund_thresh Numeric value specifying the minimum mean relative abundance
#'   threshold for filtering low-abundance taxa. Default is 0.001.
#' @param use_fdr Logical indicating whether to apply False Discovery Rate (FDR) correction.
#'   Default is \code{TRUE}.
#' @param min_sub_units Integer specifying the minimum number of sub-units required
#'   for analysis at each taxonomic level. Default is 2.
#' @param verbose Logical indicating whether to print progress messages. Default is \code{TRUE}.
#'
#' @return A data.frame containing the identified minimal change units with the following columns:
#' \describe{
#'   \item{Level}{Taxonomic level of the identified unit}
#'   \item{Taxon}{Name of the taxonomic unit}
#'   \item{Taxonomy_Path}{Full taxonomic path separated by "|"}
#'   \item{N_Sub_Units}{Number of sub-units within this taxonomic unit}
#'   \item{Consist_Prop}{Proportion of sub-units changing in the same direction}
#'   \item{Direction}{Direction of change ("Up_in_Group1" or "Down_in_Group1")}
#'   \item{P_Value}{P-value from Wilcoxon rank-sum test}
#'   \item{P_Value_FDR}{FDR-adjusted p-value (if use_fdr = TRUE)}
#'   \item{Upper_Consistency}{Logical indicating if upper taxonomic level also shows consistency}
#'   \item{Mean_LFC}{Mean log fold change across sub-units}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Applies CLR transformation to the abundance data
#'   \item Filters low-abundance taxa based on \code{abund_thresh}
#'   \item For each taxonomic level, aggregates taxa and identifies sub-units
#'   \item Calculates consistency of change direction among sub-units
#'   \item Tests for statistical significance using Wilcoxon rank-sum test
#'   \item Identifies minimal change units by selecting the most specific taxonomic level
#'   \item Optionally applies FDR correction for multiple testing
#' }
#'
#' @examples
#' \dontrun{
#' # Load example data
#' data("GlobalPatterns", package = "phyloseq")
#'
#' # Add binary grouping variable
#' sample_data(GlobalPatterns)$Treatment <- factor(
#'   sample(c("Control", "Treatment"), nsamples(GlobalPatterns), replace = TRUE)
#' )
#'
#' # Run analysis
#' results <- identify_minimal_change_units(
#'   physeq = GlobalPatterns,
#'   group_col = "Treatment",
#'   consist_thresh = 0.8,
#'   p_thresh = 0.05
#' )
#'
#' # View results
#' print(results)
#' }
#'
#' @seealso \code{\link{plot_change_units}} for visualization of results
#'
#' @export
identify_minimal_change_units <- function(physeq,
                                          group_col,
                                          levels = c("Species", "Genus", "Family", "Order", "Class", "Phylum"),
                                          consist_thresh = 0.8,
                                          p_thresh = 0.05,
                                          abund_thresh = 0.001,
                                          use_fdr = TRUE,
                                          min_sub_units = 2,
                                          verbose = TRUE) {

  # Input validation
  .validate_inputs(physeq, group_col, levels, consist_thresh, p_thresh,
                   abund_thresh, min_sub_units)

  # Extract and validate grouping variable
  meta <- phyloseq::sample_data(physeq)
  groups <- factor(meta[[group_col]])
  group_levels <- levels(groups)

  if (length(group_levels) != 2) {
    stop("Group variable must have exactly 2 levels, found: ",
         paste(group_levels, collapse = ", "))
  }

  # Validate taxonomic levels
  available_ranks <- colnames(phyloseq::tax_table(physeq))
  missing_levels <- setdiff(levels, available_ranks)
  if (length(missing_levels) > 0) {
    warning("Missing taxonomic levels: ", paste(missing_levels, collapse = ", "))
    levels <- intersect(levels, available_ranks)
  }

  if (length(levels) == 0) {
    stop("No valid taxonomic levels found in the data")
  }

  # Data preprocessing
  if (verbose) cat("Preprocessing data...\n")

  # CLR transformation
  physeq_clr <- microbiome::transform(physeq, 'clr')

  # Filter low-abundance taxa
  mean_abund <- phyloseq::taxa_sums(physeq) / sum(phyloseq::taxa_sums(physeq))
  physeq_filtered <- phyloseq::prune_taxa(mean_abund > abund_thresh, physeq)
  physeq_clr_filtered <- phyloseq::prune_taxa(mean_abund > abund_thresh, physeq_clr)

  if (verbose) {
    cat("Filtered from", phyloseq::ntaxa(physeq), "to",
        phyloseq::ntaxa(physeq_filtered), "taxa\n")
  }

  # Main analysis
  results <- .analyze_taxonomic_levels(physeq_clr_filtered, groups, group_levels,
                                       levels, consist_thresh, p_thresh,
                                       min_sub_units, verbose)

  # Process results
  if (length(results) == 0) {
    warning("No significant consistent change units found")
    return(.create_empty_result())
  }

  results_df <- purrr::map_dfr(results, ~ data.frame(.x, stringsAsFactors = FALSE))

  # Apply FDR correction
  if (use_fdr && nrow(results_df) > 1) {
    results_df$P_Value_FDR <- stats::p.adjust(results_df$P_Value, method = "fdr")
    results_df <- results_df[results_df$P_Value_FDR < p_thresh, ]
  }

  if (nrow(results_df) == 0) {
    warning("No significant results after FDR correction")
    return(.create_empty_result())
  }

  # Select minimal change units
  final_results <- .select_minimal_units(results_df, levels, verbose)

  return(final_results)
}

#' Plot Minimal Change Units as Phylogenetic Tree
#'
#' Creates a tree-like visualization of the identified minimal change units,
#' displaying the taxonomic hierarchy and detailed statistics for each unit.
#'
#' @param results A data.frame output from \code{identify_minimal_change_units}.
#' @param physeq The original phyloseq object used in the analysis (optional, for extracting full taxonomy).
#' @param plot_type Character string specifying plot type: "tree" (default), "sunburst", or "treemap".
#' @param show_all Logical indicating whether to show all results or filter by significance. Default is \code{TRUE}.
#' @param p_threshold P-value threshold for filtering when show_all is FALSE. Default is 0.05.
#' @param title Character string for the plot title. If NULL, a default title is used.
#' @param interactive Logical indicating whether to create an interactive plot using plotly. Default is \code{FALSE}.
#'
#' @return A ggplot2 object (or plotly object if interactive=TRUE) showing the taxonomic tree with change units.
#'
#' @details
#' The plot creates a tree-like structure showing:
#' \itemize{
#'   \item Taxonomic hierarchy from Kingdom to Species
#'   \item Node sizes proportional to consistency percentage
#'   \item Colors indicating direction of change
#'   \item Labels showing number of consistent sub-units and percentages
#'   \item Statistical significance information
#' }
#'
#' @examples
#' \dontrun{
#' # Basic tree plot
#' p1 <- plot_change_units(results, physeq = your_phyloseq)
#'
#' # Interactive sunburst plot
#' p2 <- plot_change_units(results, plot_type = "sunburst", interactive = TRUE)
#'
#' # Treemap visualization
#' p3 <- plot_change_units(results, plot_type = "treemap")
#' }
#'
#' @export

plot_change_units <- function(results,
                              physeq = NULL,
                              plot_type = "tree",
                              show_all = TRUE,
                              p_threshold = 0.05,
                              title = NULL,
                              interactive = FALSE) {

  # Check required packages
  required_packages <- c("ggplot2", "dplyr")
  if (plot_type == "sunburst" && interactive) required_packages <- c(required_packages, "plotly", "sunburstR")
  if (plot_type == "treemap") required_packages <- c(required_packages, "treemap")
  if (interactive && plot_type != "sunburst") required_packages <- c(required_packages, "plotly")

  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }

  # Input validation
  if (!is.data.frame(results) || nrow(results) == 0) {
    warning("No results to plot")
    return(NULL)
  }

  # Filter results if needed
  if (!show_all) {
    p_col <- ifelse("P_Value_FDR" %in% colnames(results), "P_Value_FDR", "P_Value")
    results <- results[results[[p_col]] < p_threshold, ]
    if (nrow(results) == 0) {
      warning("No significant results to plot after filtering")
      return(NULL)
    }
  }

  # Route to appropriate plotting function
  switch(plot_type,
         "tree" = .plot_phylogenetic_tree(results, physeq, title, interactive),
         "sunburst" = .plot_sunburst(results, physeq, title, interactive),
         "treemap" = .plot_treemap(results, physeq, title),
         stop("Invalid plot_type. Choose 'tree', 'sunburst', or 'treemap'")
  )
}

#' Create Summary Statistics for Minimal Change Units
#'
#' Generates summary statistics and diagnostic information for the identified
#' minimal change units.
#'
#' @param results A data.frame output from \code{identify_minimal_change_units}.
#' @param verbose Logical indicating whether to print summary to console. Default is \code{TRUE}.
#'
#' @return A list containing summary statistics:
#' \describe{
#'   \item{n_total}{Total number of minimal change units identified}
#'   \item{level_distribution}{Table of units by taxonomic level}
#'   \item{direction_distribution}{Table of units by direction of change}
#'   \item{consistency_summary}{Summary statistics of consistency proportions}
#'   \item{significance_summary}{Summary statistics of p-values}
#' }
#'
#' @examples
#' \dontrun{
#' # Assuming 'results' is output from identify_minimal_change_units
#' summary_stats <- summarize_change_units(results)
#' }
#'
#' @export
summarize_change_units <- function(results, verbose = TRUE) {

  if (!is.data.frame(results) || nrow(results) == 0) {
    if (verbose) cat("No results to summarize\n")
    return(list(n_total = 0))
  }

  # Calculate summary statistics
  summary_list <- list(
    n_total = nrow(results),
    level_distribution = table(results$Level),
    direction_distribution = table(results$Direction),
    consistency_summary = summary(results$Consist_Prop),
    significance_summary = summary(results$P_Value)
  )

  # Add FDR summary if available
  if ("P_Value_FDR" %in% colnames(results)) {
    summary_list$fdr_summary <- summary(results$P_Value_FDR)
  }

  # Print summary if requested
  if (verbose) {
    cat("\n=== Minimal Change Units Summary ===\n")
    cat("Total units identified:", summary_list$n_total, "\n\n")

    cat("Distribution by taxonomic level:\n")
    print(summary_list$level_distribution)
    cat("\n")

    cat("Distribution by direction:\n")
    print(summary_list$direction_distribution)
    cat("\n")

    cat("Consistency proportion summary:\n")
    print(summary_list$consistency_summary)
    cat("\n")

    cat("P-value summary:\n")
    print(summary_list$significance_summary)

    if ("fdr_summary" %in% names(summary_list)) {
      cat("\nFDR-adjusted p-value summary:\n")
      print(summary_list$fdr_summary)
    }
    cat("\n")
  }

  return(summary_list)
}

# Plotting helper functions ----

#' Create Phylogenetic Tree Plot
#' @keywords internal
.plot_phylogenetic_tree <- function(results, physeq, title, interactive) {

  # Prepare tree data
  tree_data <- .prepare_tree_data(results, physeq)

  if (is.null(title)) {
    title <- paste("Taxonomic Tree of Minimal Change Units (n =", nrow(results), ")")
  }

  # Calculate node positions for tree layout
  tree_layout <- .calculate_tree_layout(tree_data)

  # Create base plot
  p <- ggplot2::ggplot(tree_layout, ggplot2::aes(x = x, y = y)) +
    # Draw tree branches
    ggplot2::geom_segment(ggplot2::aes(xend = parent_x, yend = parent_y),
                          color = "gray60", size = 0.5, alpha = 0.7) +
    # Draw nodes
    ggplot2::geom_point(ggplot2::aes(color = Direction, size = Consist_Prop_Numeric),
                        alpha = 0.8) +
    # Add text labels
    ggplot2::geom_text(ggplot2::aes(label = Node_Label),
                       hjust = -0.1, vjust = 0.5, size = 3.5, check_overlap = TRUE) +
    # Styling
    ggplot2::scale_color_manual(
      values = c("Up_in_Group1" = "#E31A1C", "Down_in_Group1" = "#1F78B4", "Root" = "#999999"),
      na.value = "#999999"
    ) +
    ggplot2::scale_size_continuous(
      name = "Consistency %",
      range = c(3, 12),
      breaks = c(0.8, 0.85, 0.9, 0.95, 1.0),
      labels = c("80%", "85%", "90%", "95%", "100%")
    ) +
    ggplot2::labs(
      title = title,
      subtitle = paste("Node size = consistency percentage; Color = direction of change"),
      color = "Direction of Change"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 14, hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = 11, hjust = 0.5, color = "gray60"),
      legend.position = "bottom"
    )

  if (interactive) {
    if (!requireNamespace("plotly", quietly = TRUE)) {
      warning("plotly package required for interactive plots. Returning static plot.")
      return(p)
    }

    # Add hover information
    tree_layout$hover_text <- paste(
      "Taxon:", tree_layout$Taxon,
      "<br>Level:", tree_layout$Level,
      "<br>Consistency:", paste0(round(tree_layout$Consist_Prop_Numeric * 100, 1), "%"),
      "<br>N Sub-units:", tree_layout$N_Sub_Units,
      "<br>P-value:", formatC(tree_layout$P_Value, format = "e", digits = 2),
      "<br>Mean LFC:", round(tree_layout$Mean_LFC, 3)
    )

    p <- p + ggplot2::aes(text = hover_text)
    p <- plotly::ggplotly(p, tooltip = "text")
  }

  return(p)
}

#' Create Sunburst Plot
#' @keywords internal
.plot_sunburst <- function(results, physeq, title, interactive) {

  if (!requireNamespace("sunburstR", quietly = TRUE)) {
    stop("sunburstR package required for sunburst plots")
  }

  # Prepare sunburst data
  sunburst_data <- .prepare_sunburst_data(results, physeq)

  if (is.null(title)) {
    title <- "Taxonomic Sunburst of Minimal Change Units"
  }

  # Create color palette
  colors <- c("#E31A1C", "#1F78B4", "#999999")
  names(colors) <- c("Up_in_Group1", "Down_in_Group1", "Root")

  # Create sunburst
  sb <- sunburstR::sunburst(
    data = sunburst_data,
    legend = list(w = 150, h = 25, s = 5, t = 25),
    colors = colors,
    explanation = paste0("<span style='font-size:1.2em;'>", title, "</span>")
  )

  return(sb)
}

#' Create Treemap Plot
#' @keywords internal
.plot_treemap <- function(results, physeq, title) {

  if (!requireNamespace("treemap", quietly = TRUE)) {
    stop("treemap package required for treemap plots")
  }

  # Prepare treemap data
  treemap_data <- .prepare_treemap_data(results, physeq)

  if (is.null(title)) {
    title <- "Treemap of Minimal Change Units"
  }

  # Create treemap
  treemap::treemap(
    treemap_data,
    index = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    vSize = "N_Sub_Units",
    vColor = "Consist_Prop_Numeric",
    type = "value",
    palette = "RdYlBu",
    title = title,
    fontsize.labels = c(15, 12, 10, 8, 6, 4, 2),
    fontcolor.labels = c("white", "white", "black", "black", "black", "black", "black"),
    border.col = "white",
    border.lwds = c(4, 3, 2, 1, 0.5, 0.2, 0.1)
  )
}

#' Prepare Tree Data Structure
#' @keywords internal
.prepare_tree_data <- function(results, physeq) {

  # Standard taxonomic levels (7 levels as specified)
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Extract or parse taxonomy
  if (!is.null(physeq)) {
    full_tax <- .extract_full_taxonomy(results, physeq)
  } else {
    full_tax <- .parse_taxonomy_path(results)
  }

  # Create unique taxonomy paths and build tree structure
  tree_nodes <- list()
  node_counter <- 1
  node_map <- list()  # To track existing nodes

  # Add root node
  root_key <- "Root"
  tree_nodes[[node_counter]] <- list(
    id = node_counter,
    name = "Microbiome",
    level = "Root",
    parent_id = NA,
    taxon = "Root",
    Direction = "Root",
    Consist_Prop_Numeric = 0.5,
    N_Sub_Units = 0,
    P_Value = 1,
    Mean_LFC = 0
  )
  node_map[[root_key]] <- node_counter
  node_counter <- node_counter + 1

  # Process each result to build the tree
  for (i in 1:nrow(results)) {
    result_row <- results[i, ]
    tax_path <- full_tax[i, ]

    current_parent_id <- 1  # Start from root
    current_path <- "Root"

    # Build path through taxonomic levels
    for (level_idx in seq_along(tax_levels)) {
      level <- tax_levels[level_idx]
      taxon_name <- tax_path[[level]]

      # Skip if taxon name is missing
      if (is.na(taxon_name) || taxon_name == "" || taxon_name == "NA") next

      # Create unique key for this node
      node_key <- paste(current_path, level, taxon_name, sep = "::")

      # Check if this node already exists
      if (!node_key %in% names(node_map)) {
        # Create new node
        is_target_node <- (level == result_row$Level && taxon_name == result_row$Taxon)

        tree_nodes[[node_counter]] <- list(
          id = node_counter,
          name = taxon_name,
          level = level,
          parent_id = current_parent_id,
          taxon = taxon_name,
          Direction = if (is_target_node) result_row$Direction else NA,
          Consist_Prop_Numeric = if (is_target_node) result_row$Consist_Prop else 0.5,
          N_Sub_Units = if (is_target_node) result_row$N_Sub_Units else 0,
          P_Value = if (is_target_node) result_row$P_Value else 1,
          Mean_LFC = if (is_target_node) result_row$Mean_LFC else 0
        )

        node_map[[node_key]] <- node_counter
        current_parent_id <- node_counter
        node_counter <- node_counter + 1
      } else {
        # Use existing node
        existing_node_id <- node_map[[node_key]]
        current_parent_id <- existing_node_id

        # Update if this is a target node (in case of multiple matches)
        is_target_node <- (level == result_row$Level && taxon_name == result_row$Taxon)
        if (is_target_node) {
          tree_nodes[[existing_node_id]]$Direction <- result_row$Direction
          tree_nodes[[existing_node_id]]$Consist_Prop_Numeric <- result_row$Consist_Prop
          tree_nodes[[existing_node_id]]$N_Sub_Units <- result_row$N_Sub_Units
          tree_nodes[[existing_node_id]]$P_Value <- result_row$P_Value
          tree_nodes[[existing_node_id]]$Mean_LFC <- result_row$Mean_LFC
        }
      }

      current_path <- node_key
    }
  }

  return(tree_nodes)
}

#' Calculate Tree Layout Positions
#' @keywords internal
.calculate_tree_layout <- function(tree_nodes) {

  # Convert to data frame
  tree_df <- dplyr::bind_rows(tree_nodes)

  # Initialize positions
  tree_df$depth <- NA
  tree_df$x <- NA
  tree_df$y <- NA

  # Set root position
  root_idx <- which(tree_df$id == 1)
  tree_df$depth[root_idx] <- 0
  tree_df$x[root_idx] <- 0
  tree_df$y[root_idx] <- 0

  # Calculate positions using breadth-first approach
  max_depth <- 8  # Accommodate all taxonomic levels plus root

  for (depth in 0:(max_depth-1)) {
    # Get all nodes at current depth
    current_level_indices <- which(tree_df$depth == depth)

    if (length(current_level_indices) == 0) next

    # Process each parent node
    for (parent_idx in current_level_indices) {
      parent_node <- tree_df[parent_idx, ]

      # Find children of this parent
      children_indices <- which(tree_df$parent_id == parent_node$id & is.na(tree_df$depth))

      if (length(children_indices) > 0) {
        # Set children positions
        n_children <- length(children_indices)
        child_depth <- depth + 1
        child_x <- parent_node$x + 1

        if (n_children == 1) {
          child_y_positions <- parent_node$y
        } else {
          # Spread children vertically around parent
          spread <- (n_children - 1) * 0.8  # Adjust spacing as needed
          child_y_positions <- seq(
            from = parent_node$y - spread/2,
            to = parent_node$y + spread/2,
            length.out = n_children
          )
        }

        # Update children positions
        for (i in seq_along(children_indices)) {
          child_idx <- children_indices[i]
          tree_df$depth[child_idx] <- child_depth
          tree_df$x[child_idx] <- child_x
          tree_df$y[child_idx] <- child_y_positions[i]
        }
      }
    }
  }

  # Add parent coordinates for drawing lines
  tree_df$parent_x <- NA
  tree_df$parent_y <- NA

  for (i in 1:nrow(tree_df)) {
    if (!is.na(tree_df$parent_id[i])) {
      parent_idx <- which(tree_df$id == tree_df$parent_id[i])
      if (length(parent_idx) > 0) {
        tree_df$parent_x[i] <- tree_df$x[parent_idx[1]]
        tree_df$parent_y[i] <- tree_df$y[parent_idx[1]]
      }
    }
  }

  # Create node labels with consistency information
  tree_df$Node_Label <- ifelse(
    is.na(tree_df$Direction) | tree_df$Direction == "Root",
    tree_df$name,
    paste0(
      tree_df$name, "\n(",
      round(tree_df$Consist_Prop_Numeric * 100, 1), "%, n=",
      tree_df$N_Sub_Units, ")"
    )
  )

  return(tree_df)
}

#' Extract Full Taxonomy from Phyloseq Object
#' @keywords internal
.extract_full_taxonomy <- function(results, physeq) {

  tax_table_full <- as.data.frame(phyloseq::tax_table(physeq))
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Initialize result matrix
  full_taxonomy <- matrix(NA, nrow = nrow(results), ncol = length(tax_levels))
  colnames(full_taxonomy) <- tax_levels

  # For each result, find matching taxa and extract full taxonomy
  for (i in 1:nrow(results)) {
    result_row <- results[i, ]
    target_level <- result_row$Level
    target_taxon <- result_row$Taxon

    # Find matching rows in tax table
    if (target_level %in% colnames(tax_table_full)) {
      matching_rows <- which(tax_table_full[[target_level]] == target_taxon)

      if (length(matching_rows) > 0) {
        # Take first match and extract full taxonomy
        match_row <- tax_table_full[matching_rows[1], ]

        for (level in tax_levels) {
          if (level %in% colnames(tax_table_full)) {
            full_taxonomy[i, level] <- match_row[[level]]
          }
        }
      }
    }
  }

  return(as.data.frame(full_taxonomy))
}

#' Parse Taxonomy Path from Results
#' @keywords internal
.parse_taxonomy_path <- function(results) {

  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  full_taxonomy <- matrix(NA, nrow = nrow(results), ncol = length(tax_levels))
  colnames(full_taxonomy) <- tax_levels

  for (i in 1:nrow(results)) {
    if ("Taxonomy_Path" %in% colnames(results)) {
      # Parse the taxonomy path: Species|Genus|Family|Order|Class|Phylum|Kingdom
      path_parts <- strsplit(results$Taxonomy_Path[i], "\\|")[[1]]
      path_parts <- path_parts[path_parts != "NA" & path_parts != "" & !is.na(path_parts)]

      # The taxonomy path seems to be in reverse order (most specific to least specific)
      # Let's try to map them correctly
      if (length(path_parts) > 0) {
        # Reverse the order to go from Kingdom to Species
        path_parts <- rev(path_parts)

        # Map to appropriate levels (starting from the most general available)
        for (j in seq_along(path_parts)) {
          if (j <= length(tax_levels)) {
            # Start from Kingdom level and work down
            level_idx <- length(tax_levels) - length(path_parts) + j
            if (level_idx > 0 && level_idx <= length(tax_levels)) {
              full_taxonomy[i, level_idx] <- path_parts[j]
            }
          }
        }
      }
    }

    # Ensure the target taxon is correctly placed
    target_level <- results$Level[i]
    target_taxon <- results$Taxon[i]

    if (target_level %in% tax_levels && !is.na(target_taxon) && target_taxon != "") {
      level_idx <- which(tax_levels == target_level)
      full_taxonomy[i, level_idx] <- target_taxon
    }
  }

  # Fill in missing higher-level taxonomy with generic names
  for (i in 1:nrow(full_taxonomy)) {
    for (j in 1:ncol(full_taxonomy)) {
      if (is.na(full_taxonomy[i, j]) || full_taxonomy[i, j] == "" || full_taxonomy[i, j] == "NA") {
        # Use a generic name based on the lowest known taxonomy
        lowest_known <- NULL
        for (k in ncol(full_taxonomy):1) {
          if (!is.na(full_taxonomy[i, k]) && full_taxonomy[i, k] != "" && full_taxonomy[i, k] != "NA") {
            lowest_known <- full_taxonomy[i, k]
            break
          }
        }
        if (!is.null(lowest_known)) {
          full_taxonomy[i, j] <- paste0("Unclassified_", tax_levels[j], "_of_", lowest_known)
        } else {
          full_taxonomy[i, j] <- paste0("Unknown_", tax_levels[j])
        }
      }
    }
  }

  return(as.data.frame(full_taxonomy, stringsAsFactors = FALSE))
}

#' Prepare Sunburst Data
#' @keywords internal
.prepare_sunburst_data <- function(results, physeq) {

  full_tax <- if (!is.null(physeq)) .extract_full_taxonomy(results, physeq) else .parse_taxonomy_path(results)
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Create sunburst paths
  sunburst_paths <- c()

  for (i in 1:nrow(results)) {
    path_parts <- c()

    for (level in tax_levels) {
      taxon <- full_tax[i, level]
      if (!is.na(taxon) && taxon != "") {
        # Add consistency info for target level
        if (level == results$Level[i] && taxon == results$Taxon[i]) {
          consistency_info <- paste0("(", round(results$Consist_Prop[i] * 100, 1), "%)")
          path_parts <- c(path_parts, paste0(taxon, " ", consistency_info))
        } else {
          path_parts <- c(path_parts, taxon)
        }
      }
    }

    if (length(path_parts) > 0) {
      sunburst_paths <- c(sunburst_paths, paste(path_parts, collapse = "-"))
    }
  }

  # Create data frame for sunburst
  sunburst_df <- data.frame(
    path = sunburst_paths,
    value = results$N_Sub_Units,
    stringsAsFactors = FALSE
  )

  return(sunburst_df)
}

#' Prepare Treemap Data
#' @keywords internal
.prepare_treemap_data <- function(results, physeq) {

  full_tax <- if (!is.null(physeq)) .extract_full_taxonomy(results, physeq) else .parse_taxonomy_path(results)
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

  # Combine taxonomy with results
  treemap_df <- cbind(full_tax, results)

  # Fill missing taxonomic levels
  for (level in tax_levels) {
    if (level %in% colnames(treemap_df)) {
      treemap_df[[level]][is.na(treemap_df[[level]]) | treemap_df[[level]] == ""] <- paste0("Unknown_", level)
    } else {
      treemap_df[[level]] <- paste0("Unknown_", level)
    }
  }

  return(treemap_df)
}

#' Test Function for Plot Change Units
#'
#' A simplified test function to debug the plotting without requiring phyloseq object.
#'
#' @param results A data.frame with the results structure shown in the example
#' @return A ggplot2 object or error message
#'
#' @examples
#' \dontrun{
#' # Test with the provided results
#' p <- test_plot_change_units(results)
#' print(p)
#' }
#'
#' @export
test_plot_change_units <- function(results) {

  # Validate input
  if (!is.data.frame(results) || nrow(results) == 0) {
    return("Error: Invalid or empty results data frame")
  }

  required_cols <- c("Level", "Taxon", "Taxonomy_Path", "N_Sub_Units",
                     "Consist_Prop", "Direction", "P_Value", "Mean_LFC")
  missing_cols <- setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0) {
    return(paste("Error: Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Try the tree plot
  tryCatch({
    plot_change_units(results, physeq = NULL, plot_type = "tree")
  }, error = function(e) {
    paste("Error in tree plot:", e$message)
  })
}

# Test data creation function for debugging
create_test_data <- function() {
  data.frame(
    Level = c("Species", "Genus", "Family", "Order", "Class"),
    Taxon = c("Test_species", "Test_genus", "Test_family", "Test_order", "Test_class"),
    Taxonomy_Path = c("Test_species|Test_genus|Test_family|Test_order|Test_class|Test_phylum",
                      "Test_genus|Test_family|Test_order|Test_class|Test_phylum",
                      "Test_family|Test_order|Test_class|Test_phylum",
                      "Test_order|Test_class|Test_phylum",
                      "Test_class|Test_phylum"),
    N_Sub_Units = c(2, 3, 4, 5, 6),
    Consist_Prop = c(1.0, 0.9, 0.8, 0.85, 0.95),
    Direction = c("Up_in_Group1", "Down_in_Group1", "Up_in_Group1", "Down_in_Group1", "Up_in_Group1"),
    P_Value = c(0.01, 0.005, 0.02, 0.001, 0.03),
    Upper_Consistency = c(FALSE, TRUE, FALSE, TRUE, FALSE),
    Mean_LFC = c(1.2, -0.8, 0.9, -1.1, 0.7),
    stringsAsFactors = FALSE
  )
}

# Internal helper functions ----

#' Validate Input Parameters
#' @keywords internal
.validate_inputs <- function(physeq, group_col, levels, consist_thresh,
                             p_thresh, abund_thresh, min_sub_units) {

  if (!inherits(physeq, "phyloseq")) {
    stop("'physeq' must be a phyloseq object")
  }

  if (!is.character(group_col) || length(group_col) != 1) {
    stop("'group_col' must be a single character string")
  }

  if (!group_col %in% colnames(phyloseq::sample_data(physeq))) {
    stop("'group_col' not found in sample data")
  }

  if (!is.numeric(consist_thresh) || consist_thresh < 0 || consist_thresh > 1) {
    stop("'consist_thresh' must be between 0 and 1")
  }

  if (!is.numeric(p_thresh) || p_thresh <= 0 || p_thresh >= 1) {
    stop("'p_thresh' must be between 0 and 1")
  }

  if (!is.numeric(abund_thresh) || abund_thresh < 0) {
    stop("'abund_thresh' must be non-negative")
  }

  if (!is.numeric(min_sub_units) || min_sub_units < 2) {
    stop("'min_sub_units' must be at least 2")
  }
}

#' Calculate Consistency and Direction
#' @keywords internal
.calculate_consistency <- function(otu_matrix, groups, group_levels) {

  # Calculate log fold changes
  lfc <- apply(otu_matrix, 2, function(x) {
    mean1 <- mean(x[groups == group_levels[1]], na.rm = TRUE)
    mean2 <- mean(x[groups == group_levels[2]], na.rm = TRUE)
    mean1 - mean2
  })

  # Calculate consistency
  signs <- sign(lfc[!is.na(lfc)])
  if (length(signs) == 0) {
    return(list(consist_prop = 0, direction = "None", lfc = lfc))
  }

  pos_prop <- sum(signs > 0) / length(signs)
  neg_prop <- sum(signs < 0) / length(signs)
  consist_prop <- max(pos_prop, neg_prop)
  direction <- ifelse(pos_prop > neg_prop, "Up_in_Group1", "Down_in_Group1")

  return(list(consist_prop = consist_prop, direction = direction, lfc = lfc))
}

#' Calculate Statistical Significance
#' @keywords internal
.calculate_significance <- function(otu_matrix, groups, group_levels) {

  if (ncol(otu_matrix) == 1) {
    group1_vals <- otu_matrix[groups == group_levels[1], 1]
    group2_vals <- otu_matrix[groups == group_levels[2], 1]
  } else {
    group1_vals <- rowMeans(otu_matrix[groups == group_levels[1], ], na.rm = TRUE)
    group2_vals <- rowMeans(otu_matrix[groups == group_levels[2], ], na.rm = TRUE)
  }

  tryCatch({
    test_result <- stats::wilcox.test(group1_vals, group2_vals)
    return(test_result$p.value)
  }, error = function(e) {
    return(1.0)
  })
}

#' Analyze Taxonomic Levels
#' @keywords internal
.analyze_taxonomic_levels <- function(physeq_clr_filtered, groups, group_levels,
                                      levels, consist_thresh, p_thresh,
                                      min_sub_units, verbose) {

  results <- list()
  otu_clr <- as(phyloseq::otu_table(physeq_clr_filtered), "matrix")
  if (phyloseq::taxa_are_rows(physeq_clr_filtered)) otu_clr <- t(otu_clr)

  tax_table_mat <- as(phyloseq::tax_table(physeq_clr_filtered), "matrix")

  for (level_idx in seq_along(levels)) {
    current_level <- levels[level_idx]
    if (verbose) cat("Processing level:", current_level, "\n")

    # Aggregate to current taxonomic level
    physeq_agg <- tryCatch({
      phyloseq::tax_glom(physeq_clr_filtered, taxrank = current_level)
    }, error = function(e) {
      if (verbose) warning("Failed to aggregate at level ", current_level, ": ", e$message)
      return(NULL)
    })

    if (is.null(physeq_agg)) next

    tax_agg <- as(phyloseq::tax_table(physeq_agg), "matrix")

    # Analyze each aggregated unit
    for (i in 1:nrow(tax_agg)) {
      taxon_name <- tax_agg[i, current_level]
      if (is.na(taxon_name) || taxon_name == "") next

      # Find sub-units
      sub_indices <- which(tax_table_mat[, current_level] == taxon_name)
      if (length(sub_indices) < min_sub_units) next

      sub_otu <- otu_clr[, sub_indices, drop = FALSE]

      # Calculate consistency and significance
      consistency_result <- .calculate_consistency(sub_otu, groups, group_levels)
      p_value <- .calculate_significance(sub_otu, groups, group_levels)

      if (consistency_result$consist_prop >= consist_thresh && p_value < p_thresh) {

        # Check upper level consistency
        upper_consistency <- .check_upper_consistency(
          tax_agg, i, tax_table_mat, otu_clr, groups, group_levels,
          levels, level_idx, consist_thresh, p_thresh, min_sub_units
        )

        # Build result entry
        taxonomy_path <- paste(tax_agg[i, levels[levels %in% colnames(tax_agg)]],
                               collapse = "|")

        result_entry <- list(
          Level = current_level,
          Taxon = taxon_name,
          Taxonomy_Path = taxonomy_path,
          N_Sub_Units = length(sub_indices),
          Consist_Prop = round(consistency_result$consist_prop, 3),
          Direction = consistency_result$direction,
          P_Value = p_value,
          Upper_Consistency = upper_consistency,
          Mean_LFC = round(mean(consistency_result$lfc, na.rm = TRUE), 3)
        )

        result_key <- paste(current_level, taxon_name, sep = "_")
        results[[result_key]] <- result_entry
      }
    }
  }

  return(results)
}

#' Check Upper Level Consistency
#' @keywords internal
.check_upper_consistency <- function(tax_agg, i, tax_table_mat, otu_clr,
                                     groups, group_levels, levels, level_idx,
                                     consist_thresh, p_thresh, min_sub_units) {

  if (level_idx >= length(levels)) return(FALSE)

  upper_level <- levels[level_idx + 1]
  if (!upper_level %in% colnames(tax_table_mat)) return(FALSE)

  upper_taxon <- tax_agg[i, upper_level]
  if (is.na(upper_taxon) || upper_taxon == "") return(FALSE)

  upper_indices <- which(tax_table_mat[, upper_level] == upper_taxon)
  if (length(upper_indices) < min_sub_units) return(FALSE)

  upper_otu <- otu_clr[, upper_indices, drop = FALSE]
  upper_consist <- .calculate_consistency(upper_otu, groups, group_levels)
  upper_p <- .calculate_significance(upper_otu, groups, group_levels)

  return(upper_consist$consist_prop >= consist_thresh && upper_p < p_thresh)
}

#' Select Minimal Units
#' @keywords internal
.select_minimal_units <- function(results_df, levels, verbose) {

  results_df$Level_Order <- match(results_df$Level, levels)

  final_results <- results_df %>%
    dplyr::group_by(Taxonomy_Path) %>%
    dplyr::filter(Level_Order == min(Level_Order)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(Level_Order, P_Value) %>%
    dplyr::select(-Level_Order)

  if (verbose) {
    cat("\n=== Analysis Summary ===\n")
    cat("Total candidates found:", nrow(results_df), "\n")
    cat("Final minimal change units:", nrow(final_results), "\n")
    cat("Distribution by taxonomic level:\n")
    print(table(final_results$Level))
    cat("\n")
  }

  return(as.data.frame(final_results))
}

#' Create Empty Result Data Frame
#' @keywords internal
.create_empty_result <- function() {
  data.frame(
    Level = character(0),
    Taxon = character(0),
    Taxonomy_Path = character(0),
    N_Sub_Units = integer(0),
    Consist_Prop = numeric(0),
    Direction = character(0),
    P_Value = numeric(0),
    Upper_Consistency = logical(0),
    Mean_LFC = numeric(0),
    stringsAsFactors = FALSE
  )
}

#' Shorten Taxon Names for Plotting
#' @keywords internal
.shorten_taxon_names <- function(taxon_names, max_length = 30) {
  ifelse(nchar(taxon_names) > max_length,
         paste0(substr(taxon_names, 1, max_length - 3), "..."),
         taxon_names)
}
