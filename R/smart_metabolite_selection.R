#' Smart Metabolite Selection
#'
#' Intelligently determines the optimal number of metabolites to select based on
#' multiple criteria including score distribution, information gain, sample size
#' constraints, correlation patterns, and network topology.
#'
#' @param results Results object from find_associated_metabolites function
#' @param sample_size Integer, number of samples in the study (optional)
#' @param validation_available Logical, whether a validation cohort is available (default: FALSE)
#' @param resource_level Character, resource availability level: "low", "medium", or "high" (default: "medium")
#'
#' @return A list of class "metabolite_selection_report" containing:
#' \describe{
#'   \item{final_recommendation}{Integer, recommended number of metabolites}
#'   \item{confidence}{Character, confidence level: "High", "Medium", or "Low"}
#'   \item{method_recommendations}{List of recommendations from each method}
#'   \item{method_weights}{List of weights for each method}
#'   \item{tiered_recommendations}{List with core, standard, extended, and maximum set sizes}
#'   \item{metabolites}{List containing actual metabolite selections for each tier}
#'   \item{statistics}{Summary statistics of the analysis}
#'   \item{plots}{Visualization plots (if generated)}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # After running find_associated_metabolites
#' selection_report <- smart_metabolite_selection(
#'   results = results,
#'   sample_size = 30,
#'   validation_available = FALSE,
#'   resource_level = "medium"
#' )
#'
#' # Get recommended metabolites
#' selected <- selection_report$metabolites$standard
#' print(selection_report$final_recommendation)
#' }
smart_metabolite_selection <- function(results,
                                       sample_size = NULL,
                                       validation_available = FALSE,
                                       resource_level = c("medium", "low", "high")) {

  resource_level <- match.arg(resource_level)

  # Extract and validate scoring data
  if(is.null(results$all_scores)) {
    stop("Results object must contain 'all_scores' component")
  }

  scores_df <- results$all_scores
  scores <- sort(scores_df$integrated_score, decreasing = TRUE)
  n_total <- length(scores)

  # Initialize recommendation system
  recommendations <- list()
  weights <- list()

  # 1. Natural breaks analysis
  natural_break_n <- find_natural_breaks(scores)
  recommendations$natural_breaks <- natural_break_n
  weights$natural_breaks <- 0.15

  # 2. Information gain analysis
  info_gain_n <- information_gain_analysis(scores)
  recommendations$info_gain <- info_gain_n
  weights$info_gain <- 0.20

  # 3. Sample size constraint
  sample_based_n <- sample_size_constraint(scores, sample_size, results)
  recommendations$sample_constraint <- sample_based_n
  weights$sample_constraint <- 0.25

  # 4. Correlation-based independent components
  correlation_based_n <- independent_components(results)
  recommendations$correlation <- correlation_based_n
  weights$correlation <- 0.15

  # 5. Network modularity
  network_based_n <- network_modularity(results)
  recommendations$network <- network_based_n
  weights$network <- 0.10

  # 6. Smoothness analysis
  smooth_n <- smoothness_analysis(scores)
  recommendations$smoothness <- smooth_n
  weights$smoothness <- 0.15

  # Calculate weighted recommendation
  weighted_recommendation <- sum(unlist(recommendations) * unlist(weights)) / sum(unlist(weights))
  base_recommendation <- round(weighted_recommendation)

  # Apply resource adjustment
  final_recommendation <- resource_adjustment(
    base_recommendation,
    resource_level,
    validation_available
  )

  # Assess confidence
  cv_methods <- sd(unlist(recommendations)) / mean(unlist(recommendations))
  confidence <- ifelse(cv_methods < 0.2, "High",
                       ifelse(cv_methods < 0.4, "Medium", "Low"))

  # Generate comprehensive report
  report <- compile_selection_report(
    final_recommendation = final_recommendation,
    confidence = confidence,
    recommendations = recommendations,
    weights = weights,
    base_recommendation = base_recommendation,
    weighted_recommendation = weighted_recommendation,
    scores_df = scores_df,
    scores = scores,
    n_total = n_total,
    results = results
  )

  # Generate visualizations if ggplot2 is available
  if(requireNamespace("ggplot2", quietly = TRUE)) {
    plots <- visualize_selection(report, scores)
    report$plots <- plots
  }

  # Print report
  print_selection_report(report)

  class(report) <- c("metabolite_selection_report", "list")

  return(report)
}

#' Find Natural Breaks in Score Distribution
#'
#' Identifies natural groupings in metabolite scores using Jenks natural breaks
#' or standard deviation-based thresholding.
#'
#' @param scores Numeric vector of metabolite scores (sorted descending)
#'
#' @return Integer, index of the first natural break
#' @keywords internal
find_natural_breaks <- function(scores) {
  if(requireNamespace("classInt", quietly = TRUE)) {
    # Use Jenks natural breaks algorithm
    breaks <- classInt::classIntervals(scores, n = 5, style = "jenks")
    # Find first break (boundary of highest quality group)
    first_break_idx <- which(scores <= breaks$brks[2])[1] - 1

    if(is.na(first_break_idx) || first_break_idx < 1) {
      first_break_idx <- min(20, length(scores))
    }
  } else {
    # Fallback: standard deviation method
    threshold <- mean(scores) + sd(scores)
    first_break_idx <- sum(scores >= threshold)

    if(first_break_idx < 5) {
      first_break_idx <- min(10, length(scores))
    }
  }

  return(first_break_idx)
}

#' Information Gain Analysis
#'
#' Calculates the point where adding more metabolites provides diminishing returns
#' based on information gain rate.
#'
#' @param scores Numeric vector of metabolite scores (sorted descending)
#'
#' @return Integer, optimal number of metabolites based on information gain
#' @keywords internal
information_gain_analysis <- function(scores) {
  if(length(scores) < 2) return(length(scores))

  # Calculate information gain for each additional metabolite
  cumsum_score <- cumsum(scores)
  info_gain <- numeric(length(scores) - 1)

  for(i in 2:length(scores)) {
    # Information gain = incremental score / cumulative score
    info_gain[i-1] <- scores[i] / cumsum_score[i-1]
  }

  # Find where information gain drops below 5%
  threshold_idx <- which(info_gain < 0.05)[1]

  if(is.na(threshold_idx)) {
    threshold_idx <- min(20, length(scores))
  }

  return(threshold_idx)
}

#' Sample Size Constraint
#'
#' Determines maximum number of features based on sample size to avoid overfitting.
#'
#' @param scores Numeric vector of metabolite scores
#' @param sample_size Integer, number of samples
#' @param results Results object (for inferring sample size if not provided)
#'
#' @return Integer, maximum recommended number of metabolites
#' @keywords internal
sample_size_constraint <- function(scores, sample_size, results = NULL) {
  # Infer sample size if not provided
  if(is.null(sample_size)) {
    if(!is.null(results$normalized_data$metabolome_norm)) {
      sample_size <- nrow(results$normalized_data$metabolome_norm)
    } else if(!is.null(results$normalized_data$metabolome)) {
      sample_size <- nrow(results$normalized_data$metabolome)
    } else {
      return(20)  # Default value
    }
  }

  # Apply rule of thumb for feature-to-sample ratio
  if(sample_size < 30) {
    max_features <- floor(sample_size / 5)
  } else if(sample_size < 100) {
    max_features <- floor(sample_size / 4)
  } else {
    max_features <- floor(sample_size / 3)
  }

  # Ensure minimum of 5 features
  max_features <- max(5, max_features)

  # Apply quality constraint
  score_threshold <- mean(scores) - 0.5 * sd(scores)
  quality_constraint <- sum(scores >= score_threshold)

  return(min(max_features, quality_constraint, length(scores)))
}

#' Independent Components Analysis
#'
#' Estimates number of independent metabolites based on correlation structure.
#'
#' @param results Results object from find_associated_metabolites
#'
#' @return Integer, estimated number of independent components
#' @keywords internal
independent_components <- function(results) {
  if(is.null(results$correlation_results)) {
    return(20)  # Default value
  }

  cor_data <- results$correlation_results

  # Count metabolites by correlation strength
  high_cor <- sum(abs(cor_data$pearson_cor) > 0.5 &
                    cor_data$pearson_fdr < 0.05, na.rm = TRUE)
  medium_cor <- sum(abs(cor_data$pearson_cor) > 0.3 &
                      cor_data$pearson_fdr < 0.1, na.rm = TRUE)

  # Recommend based on correlation patterns
  if(high_cor >= 10) {
    return(high_cor + floor(medium_cor * 0.3))
  } else if(high_cor >= 5) {
    return(high_cor + floor(medium_cor * 0.5))
  } else {
    return(min(30, high_cor + medium_cor))
  }
}

#' Network Modularity Analysis
#'
#' Determines optimal number based on network topology.
#'
#' @param results Results object from find_associated_metabolites
#'
#' @return Integer, recommended number based on network structure
#' @keywords internal
network_modularity <- function(results) {
  if(is.null(results$network_metrics)) {
    return(20)  # Default value
  }

  network_data <- results$network_metrics

  # Identify hub nodes (high degree)
  degree_threshold <- quantile(network_data$degree, 0.75)
  hub_metabolites <- sum(network_data$degree >= degree_threshold)

  # Identify bridge nodes (high betweenness)
  bridge_threshold <- quantile(network_data$betweenness, 0.70)
  bridge_metabolites <- sum(network_data$betweenness >= bridge_threshold)

  # Combined recommendation
  return(ceiling((hub_metabolites + bridge_metabolites) / 1.5))
}

#' Smoothness Analysis
#'
#' Finds the elbow point in score distribution using curvature analysis.
#'
#' @param scores Numeric vector of metabolite scores (sorted descending)
#'
#' @return Integer, index of the elbow point
#' @keywords internal
smoothness_analysis <- function(scores) {
  if(length(scores) <= 2) {
    return(length(scores))
  }

  # Calculate second derivative (curvature)
  second_diff <- diff(diff(scores))

  # Find point of maximum curvature
  max_curvature_idx <- which.min(second_diff) + 2

  # Apply smoothing if enough data points
  window_size <- min(5, floor(length(scores) / 10))
  if(window_size > 1) {
    smoothed <- stats::filter(scores, rep(1/window_size, window_size))
    smoothed[is.na(smoothed)] <- scores[is.na(smoothed)]

    # Find elbow in smoothed curve
    smooth_diff <- diff(diff(smoothed))
    smooth_idx <- which.min(smooth_diff[!is.na(smooth_diff)]) + 2

    # Average both estimates
    return(floor((max_curvature_idx + smooth_idx) / 2))
  }

  return(min(max_curvature_idx, 15))
}

#' Resource Adjustment
#'
#' Adjusts recommendations based on available resources and validation capability.
#'
#' @param base_n Integer, base recommendation
#' @param resource_level Character, resource level: "low", "medium", or "high"
#' @param validation Logical, whether validation cohort is available
#'
#' @return Integer, adjusted recommendation
#' @keywords internal
resource_adjustment <- function(base_n, resource_level, validation) {
  # Resource level multipliers
  resource_multiplier <- switch(resource_level,
                                "low" = 0.5,      # Limited resources
                                "medium" = 1.0,   # Standard resources
                                "high" = 1.5,     # Abundant resources
                                1.0               # Default
  )

  # Validation availability multiplier
  validation_multiplier <- ifelse(validation, 1.2, 0.9)

  # Calculate adjusted recommendation
  adjusted_n <- round(base_n * resource_multiplier * validation_multiplier)

  # Apply reasonable bounds
  adjusted_n <- max(5, min(50, adjusted_n))

  return(adjusted_n)
}

#' Compile Selection Report
#'
#' Creates comprehensive report structure from selection analysis.
#'
#' @param final_recommendation Integer, final recommended number
#' @param confidence Character, confidence level
#' @param recommendations List of method recommendations
#' @param weights List of method weights
#' @param base_recommendation Integer, base recommendation before adjustment
#' @param weighted_recommendation Numeric, weighted average recommendation
#' @param scores_df Data frame of all scores
#' @param scores Numeric vector of scores
#' @param n_total Integer, total number of metabolites
#' @param results Original results object
#'
#' @return List containing comprehensive selection report
#' @keywords internal
compile_selection_report <- function(final_recommendation, confidence,
                                     recommendations, weights,
                                     base_recommendation, weighted_recommendation,
                                     scores_df, scores, n_total, results) {

  report <- list(
    final_recommendation = final_recommendation,
    confidence = confidence,
    method_recommendations = recommendations,
    method_weights = weights,
    base_recommendation = base_recommendation,
    weighted_average = weighted_recommendation,

    # Tiered recommendations
    tiered_recommendations = list(
      core_set = min(5, final_recommendation),
      standard_set = final_recommendation,
      extended_set = round(min(final_recommendation * 1.5, 50)),
      maximum_set = round(min(final_recommendation * 2, 100))
    ),

    # Actual metabolite lists
    metabolites = list(
      core = head(scores_df, min(5, final_recommendation)),
      standard = head(scores_df, final_recommendation),
      extended = head(scores_df, round(min(final_recommendation * 1.5, 50)))
    ),

    # Summary statistics
    statistics = list(
      total_metabolites = n_total,
      above_mean_score = sum(scores > mean(scores)),
      significant_correlations = sum(results$correlation_results$pearson_fdr < 0.05,
                                     na.rm = TRUE),
      high_importance = sum(scores > quantile(scores, 0.75))
    )
  )

  return(report)
}

#' Visualize Selection Process
#'
#' Creates visualizations of the metabolite selection process.
#'
#' @param report Selection report object
#' @param scores Numeric vector of metabolite scores
#'
#' @return List containing ggplot objects
#' @keywords internal
visualize_selection <- function(report, scores) {

  # Create methods comparison data
  methods_df <- data.frame(
    Method = names(report$method_recommendations),
    Recommendation = unlist(report$method_recommendations),
    Weight = unlist(report$method_weights),
    stringsAsFactors = FALSE
  )

  # Methods comparison plot
  p1 <- ggplot2::ggplot(methods_df,
                        ggplot2::aes(x = reorder(Method, Recommendation),
                                     y = Recommendation)) +
    ggplot2::geom_col(ggplot2::aes(fill = Weight)) +
    ggplot2::geom_hline(yintercept = report$final_recommendation,
                        linetype = "dashed", color = "red", size = 1) +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient(low = "lightblue", high = "darkblue") +
    ggplot2::labs(
      title = "Metabolite Selection: Method Comparison",
      subtitle = paste("Final Recommendation:", report$final_recommendation,
                       "| Confidence:", report$confidence),
      x = "Method",
      y = "Recommended Number"
    ) +
    ggplot2::theme_minimal()

  # Score distribution plot
  scores_df <- data.frame(
    rank = 1:length(scores),
    score = scores,
    stringsAsFactors = FALSE
  )

  p2 <- ggplot2::ggplot(scores_df, ggplot2::aes(x = rank, y = score)) +
    ggplot2::geom_line(size = 1) +
    ggplot2::geom_point(
      data = scores_df[1:report$final_recommendation, ],
      color = "red", size = 2
    ) +
    ggplot2::geom_vline(xintercept = report$final_recommendation,
                        linetype = "dashed", color = "red") +
    ggplot2::geom_vline(xintercept = report$tiered_recommendations$extended_set,
                        linetype = "dotted", color = "orange") +
    ggplot2::annotate("text",
                      x = report$final_recommendation,
                      y = max(scores) * 0.9,
                      label = paste("Standard:", report$final_recommendation),
                      angle = 90, vjust = -0.5, color = "red") +
    ggplot2::annotate("text",
                      x = report$tiered_recommendations$extended_set,
                      y = max(scores) * 0.9,
                      label = paste("Extended:", report$tiered_recommendations$extended_set),
                      angle = 90, vjust = -0.5, color = "orange") +
    ggplot2::labs(
      title = "Score Distribution and Selection Cutoffs",
      x = "Metabolite Rank",
      y = "Integrated Score"
    ) +
    ggplot2::theme_minimal()

  # Combine plots if gridExtra is available
  if(requireNamespace("gridExtra", quietly = TRUE)) {
    gridExtra::grid.arrange(p1, p2, ncol = 1)
  }

  return(list(p1 = p1, p2 = p2))
}

#' Print Selection Report
#'
#' Prints formatted selection report to console.
#'
#' @param report Selection report object
#'
#' @keywords internal
print_selection_report <- function(report) {
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

  cat("\n📈 Statistical Summary:\n")
  cat(sprintf("   Total metabolites analyzed : %d\n",
              report$statistics$total_metabolites))
  cat(sprintf("   Above mean score          : %d\n",
              report$statistics$above_mean_score))
  cat(sprintf("   Significant correlations  : %d\n",
              report$statistics$significant_correlations))

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

#' Print Method for Selection Report
#'
#' @param x A metabolite_selection_report object
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method print metabolite_selection_report
print.metabolite_selection_report <- function(x, ...) {
  cat("Metabolite Selection Report\n")
  cat("---------------------------\n")
  cat(sprintf("Final Recommendation: %d metabolites\n", x$final_recommendation))
  cat(sprintf("Confidence Level: %s\n", x$confidence))
  cat(sprintf("Base Recommendation: %d\n", x$base_recommendation))
  cat("\nTiered Recommendations:\n")
  cat(sprintf("  Core Set: %d\n", x$tiered_recommendations$core_set))
  cat(sprintf("  Standard Set: %d\n", x$tiered_recommendations$standard_set))
  cat(sprintf("  Extended Set: %d\n", x$tiered_recommendations$extended_set))
  cat("\nUse plot() to visualize the selection process.\n")
  invisible(x)
}

#' Plot Method for Selection Report
#'
#' @param x A metabolite_selection_report object
#' @param ... Additional arguments (ignored)
#'
#' @export
#' @method plot metabolite_selection_report
plot.metabolite_selection_report <- function(x, ...) {
  if(!is.null(x$plots)) {
    if(requireNamespace("gridExtra", quietly = TRUE)) {
      gridExtra::grid.arrange(x$plots$p1, x$plots$p2, ncol = 1)
    } else {
      print(x$plots$p1)
      print(x$plots$p2)
    }
  } else {
    message("No plots available. Ensure ggplot2 is installed.")
  }
  invisible(x)
}
