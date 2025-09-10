#' MultiOmicsIntegration: Tools for Multi-Omics Data Integration and Analysis
#'
#' This package provides comprehensive tools for integrating and analyzing
#' microbiome and metabolome data, including cross-omics alignment,
#' association analysis, and intelligent feature selection.
#'
#' _PACKAGE
#' @name MultiOmicsIntegration
NULL

#' Multi-Omics Data Alignment
#'
#' Aligns microbiome and metabolome data into a common latent space using
#' graph-regularized joint matrix factorization or CCA-based methods.
#'
#' @param microbiome_data Matrix or data.frame with samples as rows and microbes as columns
#' @param metabolome_data Matrix or data.frame with samples as rows and metabolites as columns
#' @param n_factors Integer, number of latent factors (default: 10)
#' @param method Character, alignment method: "auto", "fast", or "full" (default: "auto")
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A list of class "multi_omics_alignment" containing:
#' \itemize{
#'   \item aligned_microbiome: Aligned microbiome data
#'   \item aligned_metabolome: Aligned metabolome data
#'   \item latent_factors: Shared latent factor matrix
#'   \item quality_score: Overall alignment quality score (0-1)
#'   \item quality_details: Detailed quality metrics
#' }
#'
#' @export
#' @examples
#' \dontrun{
#' alignment <- multi_omics_alignment(microbiome, metabolome, n_factors = 10)
#' print(alignment$quality_score)
#' }
multi_omics_alignment <- function(microbiome_data,
                                  metabolome_data,
                                  n_factors = 10,
                                  method = c("auto", "fast", "full"),
                                  verbose = TRUE) {

  method <- match.arg(method)

  if(verbose) cat("Starting multi-omics alignment...\n")

  # Data preprocessing
  processed_data <- preprocess_omics(microbiome_data, metabolome_data)

  # Build sample relationship graph
  sample_graph <- build_sample_graph(processed_data$microbiome,
                                     processed_data$metabolome)

  # Joint embedding
  if(method == "auto") {
    method <- ifelse(ncol(microbiome_data) > 500, "fast", "full")
  }

  if(method == "fast") {
    aligned_data <- fast_alignment(processed_data, sample_graph, n_factors)
  } else {
    aligned_data <- full_alignment(processed_data, sample_graph, n_factors)
  }

  # Post-processing and quality assessment
  results <- postprocess_alignment(aligned_data, processed_data)

  if(verbose) {
    cat("\nAlignment complete!\n")
    cat(sprintf("Alignment quality score: %.3f\n", results$quality_score))
  }

  return(results)
}

#' Find Associated Metabolites
#'
#' Identifies metabolites associated with a target microbe using multiple
#' analytical approaches including correlation, network analysis, differential
#' analysis, and machine learning.
#'
#' @param phyloseq_obj A phyloseq object containing microbiome data
#' @param metabolome_mat Matrix or data.frame of metabolome data
#' @param target_microbe_name Character, name of the target microbe
#' @param top_n Integer, number of top metabolites to return (default: 20)
#'
#' @return A list containing:
#' \itemize{
#'   \item top_metabolites: Data.frame of top-ranked metabolites
#'   \item all_scores: Complete scoring results for all metabolites
#'   \item correlation_results: Correlation analysis results
#'   \item network_metrics: Network topology metrics
#'   \item diff_results: Differential analysis results
#'   \item ml_results: Machine learning feature importance
#' }
#'
#' @export

find_associated_metabolites <- function(phyloseq_obj,
                                        metabolome_mat,
                                        target_microbe_name,
                                        top_n = 20) {

  cat("Step 1: Preparing data...\n")
  data_list <- prepare_data(phyloseq_obj, metabolome_mat, target_microbe_name)

  cat("Step 2: Cross-omics normalization...\n")
  norm_data <- cross_omics_normalization(data_list$target_abundance,
                                         data_list$metabolome)

  cat("Step 3: Correlation analysis...\n")
  correlation_results <- comprehensive_correlation(norm_data$microbe_norm,
                                                   norm_data$metabolome_norm)

  cat("Step 4: Network analysis...\n")
  network_metrics <- metabolite_network_analysis(norm_data$metabolome_norm)

  cat("Step 5: Differential analysis...\n")
  diff_results <- differential_analysis(data_list$target_abundance,
                                        data_list$metabolome)

  cat("Step 6: Machine learning feature selection...\n")
  ml_results <- ml_feature_selection(data_list$target_abundance,
                                     data_list$metabolome)

  cat("Step 7: Integrated scoring...\n")
  final_scores <- integrated_scoring(correlation_results, network_metrics,
                                     diff_results, ml_results,
                                     data_list$metabolome)

  # Select top metabolites
  top_metabolites <- head(final_scores, top_n)

  # Create result object
  results <- list(
    top_metabolites = top_metabolites,
    all_scores = final_scores,
    correlation_results = correlation_results,
    network_metrics = network_metrics,
    diff_results = diff_results,
    ml_results = ml_results,
    normalized_data = norm_data
  )

  cat(paste("\nAnalysis complete! Top", top_n, "metabolites identified.\n"))

  return(results)
}

#' Smart Metabolite Selection
#'
#' Intelligently determines the optimal number of metabolites to select
#' based on multiple criteria and data characteristics.
#'
#' @param results Output from find_associated_metabolites function
#' @param sample_size Integer, number of samples (optional)
#' @param validation_available Logical, whether validation cohort exists
#' @param resource_level Character, resource availability: "low", "medium", or "high"
#'
#' @return A list containing:
#' \itemize{
#'   \item final_recommendation: Recommended number of metabolites
#'   \item confidence: Confidence level of recommendation
#'   \item tiered_recommendations: Recommendations for different scenarios
#'   \item metabolites: Lists of selected metabolites (core, standard, extended)
#' }
#'
#' @export
smart_metabolite_selection <- function(results,
                                       sample_size = NULL,
                                       validation_available = FALSE,
                                       resource_level = c("medium", "low", "high")) {

  resource_level <- match.arg(resource_level)

  # Extract scoring data
  scores_df <- results$all_scores
  scores <- sort(scores_df$integrated_score, decreasing = TRUE)
  n_total <- length(scores)

  # Initialize recommendation system
  recommendations <- list()
  weights <- list()

  # Multiple recommendation methods
  recommendations$natural_breaks <- find_natural_breaks(scores)
  weights$natural_breaks <- 0.15

  recommendations$info_gain <- information_gain_analysis(scores)
  weights$info_gain <- 0.20

  recommendations$sample_constraint <- sample_size_constraint(scores, sample_size, results)
  weights$sample_constraint <- 0.25

  recommendations$correlation <- independent_components(results)
  weights$correlation <- 0.15

  recommendations$network <- network_modularity(results)
  weights$network <- 0.10

  recommendations$smoothness <- smoothness_analysis(scores)
  weights$smoothness <- 0.15

  # Weighted average
  weighted_recommendation <- sum(unlist(recommendations) * unlist(weights)) / sum(unlist(weights))
  base_recommendation <- round(weighted_recommendation)

  # Apply resource adjustment
  final_recommendation <- resource_adjustment(base_recommendation, resource_level, validation_available)

  # Calculate confidence
  cv_methods <- sd(unlist(recommendations)) / mean(unlist(recommendations))
  confidence <- ifelse(cv_methods < 0.2, "High", ifelse(cv_methods < 0.4, "Medium", "Low"))

  # Generate report
  report <- list(
    final_recommendation = final_recommendation,
    confidence = confidence,
    method_recommendations = recommendations,
    method_weights = weights,
    base_recommendation = base_recommendation,
    weighted_average = weighted_recommendation,
    tiered_recommendations = list(
      core_set = min(5, final_recommendation),
      standard_set = final_recommendation,
      extended_set = min(final_recommendation * 1.5, 50),
      maximum_set = min(final_recommendation * 2, 100)
    ),
    metabolites = list(
      core = head(scores_df, min(5, final_recommendation)),
      standard = head(scores_df, final_recommendation),
      extended = head(scores_df, min(final_recommendation * 1.5, 50))
    )
  )

  print_selection_report(report, n_total)

  return(report)
}

# ==================== Internal Helper Functions ====================

#' @keywords internal
preprocess_omics <- function(microbiome, metabolome) {
  common_samples <- intersect(rownames(microbiome), rownames(metabolome))
  microbiome <- microbiome[common_samples, ]
  metabolome <- metabolome[common_samples, ]

  # CLR transform for microbiome
  microbiome_clr <- apply(microbiome + 1, 1, function(x) {
    log(x) - mean(log(x))
  })
  microbiome_clr <- t(microbiome_clr)

  # Log transform and scale metabolome
  metabolome_log <- log2(metabolome + 1)
  metabolome_scaled <- scale(metabolome_log)

  # Remove low quality features
  microbiome_clr <- microbiome_clr[, apply(microbiome_clr, 2, function(x) {
    sd(x, na.rm = TRUE) > 0.1 & sum(is.finite(x)) == length(x)
  })]

  metabolome_scaled <- metabolome_scaled[, apply(metabolome_scaled, 2, function(x) {
    sd(x, na.rm = TRUE) > 0.1 & sum(is.finite(x)) == length(x)
  })]

  # Handle NA/Inf values
  microbiome_clr[!is.finite(microbiome_clr)] <- 0
  metabolome_scaled[!is.finite(metabolome_scaled)] <- 0

  return(list(
    microbiome = microbiome_clr,
    metabolome = metabolome_scaled,
    samples = common_samples
  ))
}

#' @keywords internal
build_sample_graph <- function(microbiome, metabolome, k = 5) {
  n_samples <- nrow(microbiome)

  # Calculate within-omics similarity
  micro_dist <- as.matrix(dist(microbiome, method = "euclidean"))
  metab_dist <- as.matrix(dist(metabolome, method = "euclidean"))

  # KNN graph construction
  get_knn_matrix <- function(dist_mat, k) {
    n <- nrow(dist_mat)
    knn_mat <- matrix(0, n, n)

    for(i in 1:n) {
      neighbors <- order(dist_mat[i, ])[2:(k+1)]
      knn_mat[i, neighbors] <- 1
      knn_mat[neighbors, i] <- 1
    }
    return(knn_mat)
  }

  micro_knn <- get_knn_matrix(micro_dist, k)
  metab_knn <- get_knn_matrix(metab_dist, k)

  # Combined graph
  combined_graph <- (micro_knn + metab_knn) / 2

  # Laplacian matrix
  degree <- rowSums(combined_graph)
  laplacian <- diag(degree) - combined_graph

  # Normalized Laplacian
  degree_sqrt_inv <- diag(1/sqrt(degree + 1e-10))
  laplacian_norm <- degree_sqrt_inv %*% laplacian %*% degree_sqrt_inv

  return(list(
    adjacency = combined_graph,
    laplacian = laplacian,
    laplacian_norm = laplacian_norm,
    micro_knn = micro_knn,
    metab_knn = metab_knn
  ))
}

#' @keywords internal
fast_alignment <- function(processed_data, sample_graph, n_factors) {
  microbiome <- processed_data$microbiome
  metabolome <- processed_data$metabolome

  # Clean data
  microbiome <- microbiome[, colSums(is.na(microbiome)) == 0]
  metabolome <- metabolome[, colSums(is.na(metabolome)) == 0]

  # Remove zero variance columns
  microbiome <- microbiome[, apply(microbiome, 2, var, na.rm = TRUE) > 1e-10]
  metabolome <- metabolome[, apply(metabolome, 2, var, na.rm = TRUE) > 1e-10]

  n_comp <- min(n_factors, ncol(microbiome), ncol(metabolome), nrow(microbiome) - 1)

  # Try different CCA methods
  cca_result <- tryCatch({
    if(requireNamespace("mixOmics", quietly = TRUE)) {
      mixOmics::rcc(microbiome, metabolome, ncomp = n_comp, lambda1 = 0.1, lambda2 = 0.1)
    } else {
      cancor(microbiome, metabolome)
    }
  }, error = function(e) {
    list(
      xcoef = prcomp(microbiome, scale = TRUE)$rotation[, 1:n_comp],
      ycoef = prcomp(metabolome, scale = TRUE)$rotation[, 1:n_comp]
    )
  })

  # Extract common components
  if(inherits(cca_result, "rcc")) {
    micro_scores <- cca_result$variates$X
    metab_scores <- cca_result$variates$Y
  } else {
    micro_scores <- microbiome %*% cca_result$xcoef[, 1:n_comp]
    metab_scores <- metabolome %*% cca_result$ycoef[, 1:n_comp]
  }

  # Graph regularization
  graph_reg_scores <- (micro_scores + metab_scores) / 2
  smoothed_scores <- graph_smoothing(graph_reg_scores, sample_graph$laplacian_norm)

  # Reconstruct aligned data
  micro_loadings <- t(microbiome) %*% smoothed_scores / nrow(microbiome)
  metab_loadings <- t(metabolome) %*% smoothed_scores / nrow(metabolome)

  aligned_microbiome <- smoothed_scores %*% t(micro_loadings)
  aligned_metabolome <- smoothed_scores %*% t(metab_loadings)

  return(list(
    microbiome = aligned_microbiome,
    metabolome = aligned_metabolome,
    latent_factors = smoothed_scores,
    micro_loadings = micro_loadings,
    metab_loadings = metab_loadings
  ))
}

#' @keywords internal
full_alignment <- function(processed_data, sample_graph, n_factors) {
  microbiome <- processed_data$microbiome
  metabolome <- processed_data$metabolome
  n_samples <- nrow(microbiome)

  # Initialize matrices
  set.seed(123)
  U <- matrix(rnorm(n_samples * n_factors), n_samples, n_factors)
  V1 <- matrix(rnorm(ncol(microbiome) * n_factors), ncol(microbiome), n_factors)
  V2 <- matrix(rnorm(ncol(metabolome) * n_factors), ncol(metabolome), n_factors)

  # Iterative optimization
  max_iter <- 100
  tolerance <- 1e-4
  alpha <- 0.1

  for(iter in 1:max_iter) {
    U_old <- U

    # Update U
    micro_term <- microbiome %*% V1
    metab_term <- metabolome %*% V2
    graph_term <- alpha * sample_graph$laplacian_norm %*% U

    U <- (micro_term + metab_term - graph_term) /
      (V1 %*% t(V1) + V2 %*% t(V2) + alpha * diag(n_factors) + 1e-10)

    # Update V1 and V2
    V1 <- t(microbiome) %*% U %*% solve(t(U) %*% U + 1e-10 * diag(n_factors))
    V2 <- t(metabolome) %*% U %*% solve(t(U) %*% U + 1e-10 * diag(n_factors))

    # Check convergence
    if(max(abs(U - U_old)) < tolerance) {
      break
    }
  }

  # Reconstruct aligned data
  aligned_microbiome <- U %*% t(V1)
  aligned_metabolome <- U %*% t(V2)

  return(list(
    microbiome = aligned_microbiome,
    metabolome = aligned_metabolome,
    latent_factors = U,
    micro_loadings = V1,
    metab_loadings = V2,
    n_iterations = iter
  ))
}

#' @keywords internal
graph_smoothing <- function(scores, laplacian_norm, alpha = 0.5) {
  n <- nrow(scores)
  smoothed <- solve(diag(n) + alpha * laplacian_norm) %*% scores
  return(smoothed)
}

#' @keywords internal
postprocess_alignment <- function(aligned_data, original_data) {
  # Preserve dimension names
  rownames(aligned_data$microbiome) <- rownames(original_data$microbiome)
  colnames(aligned_data$microbiome) <- colnames(original_data$microbiome)
  rownames(aligned_data$metabolome) <- rownames(original_data$metabolome)
  colnames(aligned_data$metabolome) <- colnames(original_data$metabolome)

  # Calculate alignment quality
  quality_metrics <- assess_alignment_quality(
    aligned_data$microbiome,
    aligned_data$metabolome,
    original_data$microbiome,
    original_data$metabolome
  )

  # Compile results
  results <- list(
    aligned_microbiome = aligned_data$microbiome,
    aligned_metabolome = aligned_data$metabolome,
    latent_factors = aligned_data$latent_factors,
    micro_loadings = aligned_data$micro_loadings,
    metab_loadings = aligned_data$metab_loadings,
    quality_score = quality_metrics$overall_score,
    quality_details = quality_metrics,
    original_data = original_data
  )

  class(results) <- "multi_omics_alignment"
  return(results)
}

#' @keywords internal
assess_alignment_quality <- function(aligned_micro, aligned_metab, orig_micro, orig_metab) {
  # Information retention
  micro_reconstruction <- 1 - mean((aligned_micro - orig_micro)^2) / var(as.vector(orig_micro))
  metab_reconstruction <- 1 - mean((aligned_metab - orig_metab)^2) / var(as.vector(orig_metab))

  # Cross-omics correlation enhancement
  n_features <- min(100, ncol(aligned_micro), ncol(aligned_metab))

  orig_cors <- sapply(1:n_features, function(i) {
    cor(orig_micro[,i], orig_metab[,i], use = "pairwise.complete.obs")
  })

  aligned_cors <- sapply(1:n_features, function(i) {
    cor(aligned_micro[,i], aligned_metab[,i], use = "pairwise.complete.obs")
  })

  correlation_improvement <- mean(abs(aligned_cors), na.rm = TRUE) -
    mean(abs(orig_cors), na.rm = TRUE)

  # Overall score
  overall_score <- (micro_reconstruction + metab_reconstruction) / 2 * 0.6 +
    max(0, correlation_improvement) * 0.3 +
    0.5 * 0.1  # Default batch removal score

  return(list(
    overall_score = overall_score,
    micro_reconstruction = micro_reconstruction,
    metab_reconstruction = metab_reconstruction,
    correlation_improvement = correlation_improvement,
    mean_correlation_after = mean(abs(aligned_cors), na.rm = TRUE)
  ))
}

#' @keywords internal
print_selection_report <- function(report, n_total) {
  cat("\n=====================================\n")
  cat("    METABOLITE SELECTION REPORT\n")
  cat("=====================================\n\n")

  cat("📊 Method-specific Recommendations:\n")
  for(method in names(report$method_recommendations)) {
    cat(sprintf("   %-20s: %3d metabolites (weight: %.0f%%)\n",
                method,
                report$method_recommendations[[method]],
                report$method_weights[[method]] * 100))
  }

  cat("\n🎯 Tiered Recommendations:\n")
  cat(sprintf("   Core Set (must-have)      : %d metabolites\n",
              report$tiered_recommendations$core_set))
  cat(sprintf("   Standard Set (recommended): %d metabolites\n",
              report$tiered_recommendations$standard_set))
  cat(sprintf("   Extended Set (if feasible): %d metabolites\n",
              report$tiered_recommendations$extended_set))

  cat("\n🏆 FINAL RECOMMENDATION:\n")
  cat(sprintf("   >>> %d metabolites <<<\n", report$final_recommendation))
  cat(sprintf("   Confidence Level: %s\n", report$confidence))

  if(report$confidence == "Low") {
    cat("\n⚠️  Warning: Low confidence due to high variability between methods.\n")
    cat("   Consider manual review of the results.\n")
  }

  cat("\n=====================================\n")
}

# Additional internal helper functions would continue here...
# Note: Some helper functions from the original code are omitted for brevity
# but would be included in the full package file
