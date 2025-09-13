# RMT-based Microbiome Network Construction Algorithm
# Author: Tao Wen
# Date: 2025
#
# This file implements functions for constructing microbiome correlation
# networks based on Random Matrix Theory (RMT) and scale-free properties.
# Functions are designed for integration into an R package.
#
# Required packages: phyloseq, igraph, ggplot2, RColorBrewer, poweRlaw

#' Random Matrix Theory (RMT) Network Construction for Microbiome Data
#'
#' Applies RMT theory to determine an optimal correlation threshold and build
#' a microbiome association network from a phyloseq object.
#'
#' @param phyloseq_obj A \code{phyloseq} object containing OTU/ASV abundance data.
#' @param method Correlation method: "pearson", "spearman", or "kendall".
#' @param min_abundance Minimum relative abundance filter (default = 0.0001).
#' @param min_prevalence Minimum prevalence filter (default = 0.1).
#' @param transformation Data transformation: "clr", "log", "sqrt", or "none".
#' @param plot_results Logical, whether to generate diagnostic plots.
#'
#' @return A list containing:
#' \describe{
#'   \item{correlation_matrix}{Computed correlation matrix.}
#'   \item{processed_data}{Filtered and transformed abundance table.}
#'   \item{rmt_results}{Results of RMT threshold analysis.}
#'   \item{network}{igraph object representing the network.}
#'   \item{network_statistics}{Basic network statistics.}
#'   \item{parameters}{Input parameter settings.}
#' }
#'
#' @examples
#' \dontrun{
#'   data("GlobalPatterns")
#'   res <- rmt_microbiome_network(GlobalPatterns, method = "spearman")
#'   res$optimal_threshold
#'   plot(res$network)
#' }
#' @export
rmt_microbiome_network <- function(phyloseq_obj,
                                   method = "spearman",
                                   min_abundance = 0.0001,
                                   min_prevalence = 0.1,
                                   transformation = "clr",
                                   plot_results = TRUE) {

  cat("Starting RMT-based microbiome network analysis...\n")

  # Step 1: Data preprocessing
  cat("Step 1: Preprocessing phyloseq data...\n")
  processed_data <- preprocess_phyloseq(phyloseq_obj, min_abundance, min_prevalence, transformation)

  # Step 2: Calculate correlation matrix
  cat("Step 2: Computing correlation matrix...\n")
  cor_matrix <- compute_correlation_matrix(processed_data, method)

  # Step 3: Apply RMT analysis
  cat("Step 3: Applying RMT analysis...\n")
  rmt_results <- apply_rmt_analysis(cor_matrix, plot_results)

  # Step 4: Construct network with optimal threshold
  cat("Step 4: Constructing network with RMT-derived threshold...\n")
  network <- construct_network(cor_matrix, rmt_results$optimal_threshold)

  # Step 5: Network analysis
  cat("Step 5: Analyzing network properties...\n")
  network_stats <- analyze_network_properties(network)

  # Compile results
  results <- list(
    correlation_matrix = cor_matrix,
    processed_data = processed_data,
    rmt_results = rmt_results,
    network = network,
    network_statistics = network_stats,
    parameters = list(
      method = method,
      min_abundance = min_abundance,
      min_prevalence = min_prevalence,
      transformation = transformation
    )
  )

  cat("RMT analysis completed successfully!\n")
  cat("Optimal threshold:", round(rmt_results$optimal_threshold, 4), "\n")
  cat("Network nodes:", vcount(network), "\n")
  cat("Network edges:", ecount(network), "\n")

  return(results)
}

#' Preprocess phyloseq data
preprocess_phyloseq <- function(phyloseq_obj, min_abundance, min_prevalence, transformation) {

  # Extract OTU table
  otu_table <- as.data.frame(otu_table(phyloseq_obj))

  # Convert to relative abundance
  rel_abundance <- sweep(otu_table, 2, colSums(otu_table), FUN = "/")

  # Filter by abundance and prevalence
  abundance_filter <- rowMeans(rel_abundance) >= min_abundance
  prevalence_filter <- rowSums(rel_abundance > 0) / ncol(rel_abundance) >= min_prevalence

  filtered_data <- rel_abundance[abundance_filter & prevalence_filter, ]

  cat("  Filtered from", nrow(rel_abundance), "to", nrow(filtered_data), "taxa\n")

  # Data transformation
  if (transformation == "clr") {
    # Centered log-ratio transformation
    # Add pseudocount to handle zeros
    pseudocount <- min(filtered_data[filtered_data > 0]) / 2
    filtered_data <- filtered_data + pseudocount

    # CLR transformation
    log_data <- log(filtered_data)
    clr_data <- log_data - rowMeans(log_data)
    transformed_data <- clr_data

  } else if (transformation == "log") {
    pseudocount <- min(filtered_data[filtered_data > 0]) / 2
    transformed_data <- log(filtered_data + pseudocount)

  } else if (transformation == "sqrt") {
    transformed_data <- sqrt(filtered_data)

  } else {
    transformed_data <- filtered_data
  }

  cat("  Applied", transformation, "transformation\n")

  return(transformed_data)
}

#' Compute correlation matrix
compute_correlation_matrix <- function(data, method) {

  # Transpose so that taxa are columns (variables)
  data_t <- t(data)

  # Calculate correlation matrix
  cor_matrix <- cor(data_t, method = method, use = "complete.obs")

  # Remove self-correlations (diagonal)
  diag(cor_matrix) <- 0

  cat("  Correlation matrix dimensions:", dim(cor_matrix)[1], "x", dim(cor_matrix)[2], "\n")

  return(cor_matrix)
}
#' Apply RMT analysis to find optimal threshold
apply_rmt_analysis <- function(cor_matrix, plot_results = TRUE) {

  # Convert correlation matrix to absolute values
  abs_cor_matrix <- abs(cor_matrix)

  # Thresholds to scan
  thresholds <- seq(0.1, 0.9, by = 0.02)
  rmt_metrics <- data.frame(
    threshold = thresholds,
    nnsd = numeric(length(thresholds)),
    eigenvalue_ratio = numeric(length(thresholds)),
    spectral_density = numeric(length(thresholds)),
    edge_density = numeric(length(thresholds))
  )

  cat("  Analyzing", length(thresholds), "threshold values...\n")

  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]

    # Adjacency matrix
    adj_matrix <- ifelse(abs_cor_matrix >= thresh, 1, 0)
    diag(adj_matrix) <- 0

    # Edge density
    edge_density <- sum(adj_matrix) / (nrow(adj_matrix) * (nrow(adj_matrix) - 1))
    rmt_metrics$edge_density[i] <- edge_density

    # Skip if too sparse/dense
    if (edge_density < 0.01 || edge_density > 0.8) {
      rmt_metrics$nnsd[i] <- NA
      rmt_metrics$eigenvalue_ratio[i] <- NA
      rmt_metrics$spectral_density[i] <- NA
      next
    }

    # Eigenvalues
    eigenvals <- tryCatch({
      eigen(adj_matrix, symmetric = TRUE, only.values = TRUE)$values
    }, error = function(e) { NA })

    if (any(is.na(eigenvals))) {
      rmt_metrics$nnsd[i] <- NA
      rmt_metrics$eigenvalue_ratio[i] <- NA
      rmt_metrics$spectral_density[i] <- NA
      next
    }

    # Filter small eigenvalues
    eigenvals <- eigenvals[abs(eigenvals) > 1e-10]

    if (length(eigenvals) < 10) {
      rmt_metrics$nnsd[i] <- NA
      rmt_metrics$eigenvalue_ratio[i] <- NA
      rmt_metrics$spectral_density[i] <- NA
      next
    }

    # RMT metrics
    rmt_metrics$nnsd[i] <- calculate_nnsd_simple(eigenvals)
    rmt_metrics$eigenvalue_ratio[i] <- max(abs(eigenvals)) / mean(abs(eigenvals))
    rmt_metrics$spectral_density[i] <- calculate_spectral_density(eigenvals)

    if (i %% 10 == 0) cat("    Processed", i, "/", length(thresholds), "thresholds\n")
  }

  # Valid thresholds
  valid_indices <- which(!is.na(rmt_metrics$nnsd) &
                           rmt_metrics$edge_density >= 0.05 &
                           rmt_metrics$edge_density <= 0.5)

  if (length(valid_indices) < 3) {
    warning("Insufficient valid thresholds. Using edge density fallback.")
    target_density <- 0.15
    density_diff <- abs(rmt_metrics$edge_density - target_density)
    optimal_idx <- which.min(density_diff)
    optimal_threshold <- rmt_metrics$threshold[optimal_idx]
    rmt_method <- "edge_density_fallback"
  } else {
    valid_nnsd <- rmt_metrics$nnsd[valid_indices]
    valid_thresholds <- rmt_metrics$threshold[valid_indices]

    # Smooth NNSD and find maximum (Poisson peak)
    if (length(valid_nnsd) >= 5) {
      smooth_result <- tryCatch({
        smooth.spline(valid_thresholds, valid_nnsd, spar = 0.3)
      }, error = function(e) NULL)

      if (!is.null(smooth_result)) {
        predicted_nnsd <- predict(smooth_result, valid_thresholds)$y
        optimal_idx_in_valid <- which.max(predicted_nnsd)
        optimal_threshold <- valid_thresholds[optimal_idx_in_valid]
        rmt_method <- "NNSD_maximization"
      } else {
        optimal_idx_in_valid <- which.max(valid_nnsd)
        optimal_threshold <- valid_thresholds[optimal_idx_in_valid]
        rmt_method <- "NNSD_direct_max"
      }
    } else {
      optimal_idx_in_valid <- which.max(valid_nnsd)
      optimal_threshold <- valid_thresholds[optimal_idx_in_valid]
      rmt_method <- "NNSD_direct_max"
    }
  }

  cat("  RMT analysis completed. Method:", rmt_method, "\n")
  cat("  Valid thresholds:", length(valid_indices), "\n")

  if (plot_results) {
    plot_rmt_analysis(rmt_metrics, optimal_threshold)
  }

  return(list(
    threshold_analysis = rmt_metrics,
    optimal_threshold = optimal_threshold,
    rmt_method = rmt_method,
    valid_range = range(rmt_metrics$threshold[valid_indices])
  ))
}



#' Simplified NNSD calculation - Alternative method
calculate_nnsd_simple <- function(eigenvals) {
  # Basic filtering
  eigenvals <- eigenvals[!is.na(eigenvals) & abs(eigenvals) > 1e-10]

  # Remove duplicates (critical for network adjacency matrices)
  eigenvals <- unique(eigenvals)

  if (length(eigenvals) < 5) return(NA)  # 降低最小要求

  # Sort eigenvalues
  sorted_eigs <- sort(eigenvals)

  # Calculate spacings
  spacings <- diff(sorted_eigs)
  spacings <- spacings[spacings > 1e-12]
  if (length(spacings) < 3) return(NA)  # 降低最小要求

  # Normalize spacings
  mean_spacing <- mean(spacings)
  if (mean_spacing <= 0) return(NA)
  normalized_spacings <- spacings / mean_spacing

  # Remove extreme outliers using IQR method
  q75 <- quantile(normalized_spacings, 0.75, na.rm = TRUE)
  q25 <- quantile(normalized_spacings, 0.25, na.rm = TRUE)
  iqr <- q75 - q25
  upper_bound <- q75 + 1.5 * iqr
  lower_bound <- max(0, q25 - 1.5 * iqr)  # spacing不能为负

  normalized_spacings <- normalized_spacings[normalized_spacings >= lower_bound &
                                               normalized_spacings <= upper_bound]
  if (length(normalized_spacings) < 3) return(NA)

  # Simple statistic: variance of normalized spacings
  # For Poisson: Var(s) = 1, for GOE: Var(s) ≈ 0.18
  # For more ordered systems: Var(s) < 0.18
  spacing_var <- var(normalized_spacings)

  # Return variance as measure of "randomness"
  # Higher values = more random (Poisson-like)
  # Lower values = more structured (GOE-like)
  return(spacing_var)
}



#' Apply RMT + Scale-free threshold analysis
#'
#' This function performs Random Matrix Theory (RMT) analysis combined with
#' scale-free property evaluation to select an optimal correlation threshold
#' for network construction. It scans a series of thresholds, calculates RMT
#' metrics (nearest-neighbor spacing distribution, eigenvalue ratios, spectral
#' density) and fits power-law to degree distributions.
#'
#' @param cor_matrix A numeric correlation matrix (square, symmetric).
#' @param plot_results Logical, whether to plot analysis results (default = TRUE).
#'
#' @details
#' The function works in several steps:
#' \enumerate{
#'   \item \strong{RMT analysis:} For each threshold, build adjacency matrix
#'         and compute eigenvalue-based RMT metrics (NNSD, eigenvalue ratio,
#'         spectral density).
#'   \item \strong{Scale-free fitting:} Fit a power-law distribution to the
#'         degree distribution of the network and estimate gamma exponent.
#'   \item \strong{Threshold selection:} Integrate RMT and scale-free metrics
#'         to select an optimal threshold, maximizing
#'         \deqn{NNSD / |gamma - 2.5|}.
#'   \item \strong{Visualization:} If \code{plot_results = TRUE}, plot NNSD,
#'         gamma, and edge density across thresholds, with vertical line
#'         indicating optimal threshold.
#' }
#'
#' @return A list containing:
#' \describe{
#'   \item{threshold_analysis}{A data.frame with threshold and computed metrics.}
#'   \item{optimal_threshold}{The selected optimal threshold value.}
#'   \item{p}{ggplot object of the analysis visualization.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Example with random correlation matrix
#'   mat <- matrix(rnorm(100*100), 100, 100)
#'   cor_mat <- cor(mat)
#'   res <- apply_rmt_scale_free2(cor_mat)
#'   res$optimal_threshold
#'   res$p
#' }
#'
#' @import poweRlaw
#' @export
apply_rmt_scale_free2 <- function(cor_matrix, plot_results = TRUE) {


  # Step 1: RMT分析
  rmt_res <- apply_rmt_analysis(cor_matrix, plot_results = FALSE)
  rmt_metrics <- rmt_res$threshold_analysis

  # 扩展阈值范围 (0.1 ~ 0.99)
  thresholds <- seq(0.1, 0.99, by = 0.02)
  rmt_metrics <- data.frame(
    threshold = thresholds,
    nnsd = NA,
    eigenvalue_ratio = NA,
    spectral_density = NA,
    edge_density = NA,
    gamma = NA
  )

  cat("  Analyzing", length(thresholds), "threshold values...\n")

  for (i in seq_along(thresholds)) {
    thresh <- thresholds[i]

    # Adjacency matrix
    adj_matrix <- ifelse(abs(cor_matrix) >= thresh, 1, 0)
    diag(adj_matrix) <- 0

    # Edge density
    edge_density <- sum(adj_matrix) / (nrow(adj_matrix) * (nrow(adj_matrix) - 1))
    rmt_metrics$edge_density[i] <- edge_density

    # Skip too sparse/dense cases
    if (edge_density < 0.01 || edge_density > 0.8) next

    # Eigenvalues
    eigenvals <- tryCatch({
      eigen(adj_matrix, symmetric = TRUE, only.values = TRUE)$values
    }, error = function(e) NA)

    if (any(is.na(eigenvals))) next

    # Filter small eigenvalues
    eigenvals <- eigenvals[abs(eigenvals) > 1e-10]
    if (length(eigenvals) < 10) next

    # RMT metrics
    rmt_metrics$nnsd[i] <- calculate_nnsd_simple(eigenvals)
    rmt_metrics$eigenvalue_ratio[i] <- max(abs(eigenvals)) / mean(abs(eigenvals))
    rmt_metrics$spectral_density[i] <- calculate_spectral_density(eigenvals)

    # Scale-free fitting (power-law degree distribution)
    g <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected", diag = FALSE)
    deg <- igraph::degree(g)
    deg <- deg[deg > 0]

    if (length(deg) >= 10) {
      pl <- displ$new(deg)
      est <- tryCatch(estimate_xmin(pl), error = function(e) NULL)
      if (!is.null(est)) {
        pl$setXmin(est)
        rmt_metrics$gamma[i] <- estimate_pars(pl)$pars
      }
    }

    if (i %% 10 == 0) cat("    Processed", i, "/", length(thresholds), "thresholds\n")
  }

  # Step 2: Select optimal threshold
  valid_indices <- which(!is.na(rmt_metrics$nnsd) & !is.na(rmt_metrics$gamma))
  if (length(valid_indices) > 0) {
    candidate_thresholds <- rmt_metrics$threshold[valid_indices]
    gamma_values <- rmt_metrics$gamma[valid_indices]
    nnsd_values <- rmt_metrics$nnsd[valid_indices]

    combined_score <- nnsd_values / abs(gamma_values - 2.5)
    optimal_idx <- which.max(combined_score)
    optimal_threshold <- candidate_thresholds[optimal_idx]
  } else {
    optimal_threshold <- NA
    warning("No valid thresholds with both NNSD and gamma available.")
  }

  # Step 3: Visualization
  if (plot_results) {
    df_plot <- rmt_metrics %>%
      select(threshold, nnsd, gamma, edge_density) %>%
      pivot_longer(cols = c(nnsd, gamma, edge_density),
                   names_to = "metric", values_to = "value")

    p <- ggplot(df_plot, aes(x = threshold, y = value, color = metric)) +
      geom_line() + geom_point() +
      scale_y_continuous(
        name = "NNSD / Edge density",
        sec.axis = sec_axis(~ ., name = "Gamma")
      ) +
      geom_vline(xintercept = optimal_threshold, linetype = 2) +
      labs(title = "RMT + Scale-free Analysis", x = "Threshold") +
      theme_minimal() +
      scale_color_manual(values = c("blue", "red", "green"))
    print(p)
  }

  return(list(
    threshold_analysis = rmt_metrics,
    optimal_threshold = optimal_threshold,
    p = if (exists("p")) p else NULL
  ))
}

#' Spectral density deviation from semicircle law
#'
#' Computes deviation between observed eigenvalue density and Wigner semicircle law.
#'
#' @param eigenvals Numeric vector of eigenvalues.
#'
#' @return Deviation score (smaller = closer to semicircle law).
#' @keywords internal
#' @export
calculate_spectral_density <- function(eigenvals) {
  if (length(eigenvals) < 5) return(NA)

  # Normalize eigenvalues
  normalized_eigs <- eigenvals / max(eigenvals)

  # Calculate density using kernel density estimation
  density_est <- density(normalized_eigs, n = 50)

  # Compare with semicircle law (Wigner distribution)
  x <- density_est$x
  theoretical_density <- ifelse(abs(x) <= 1, (2/pi) * sqrt(1 - x^2), 0)

  # Calculate deviation
  deviation <- sum((density_est$y - theoretical_density)^2)

  return(deviation)
}
