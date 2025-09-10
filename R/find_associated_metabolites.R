#' MicrobiomeMetabolomeAnalysis: Cross-omics Association Analysis
#'
#' A comprehensive toolkit for identifying metabolites associated with specific microbes
#' through integrated multi-dimensional analysis including correlation, network topology,
#' differential expression, and machine learning approaches.
#'
#' _PACKAGE
#' @name MicrobiomeMetabolomeAnalysis
NULL

#' Find Associated Metabolites with Target Microbe
#'
#' Performs comprehensive analysis to identify metabolites associated with a specific microbe
#' using multiple analytical approaches including correlation, network analysis, differential
#' expression, and machine learning feature selection.
#'
#' @param phyloseq_obj A phyloseq object containing microbiome data
#' @param metabolome_mat A matrix or data.frame with metabolome data (samples as rows, metabolites as columns)
#' @param target_microbe_name Character string specifying the name of the target microbe
#' @param top_n Integer specifying the number of top metabolites to return (default: 20)
#' @param use_normalization Logical, whether to apply cross-omics normalization (default: TRUE)
#'
#' @return A list containing:
#' \describe{
#'   \item{top_metabolites}{Data frame of top-ranked metabolites with scores}
#'   \item{all_scores}{Complete scoring results for all metabolites}
#'   \item{correlation_results}{Detailed correlation analysis results}
#'   \item{network_metrics}{Network topology metrics for metabolites}
#'   \item{diff_results}{Differential expression analysis results}
#'   \item{ml_results}{Machine learning feature importance scores}
#'   \item{normalized_data}{Normalized data used in analysis}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load your data
#' # phyloseq_obj <- your_phyloseq_object
#' # metabolome_mat <- your_metabolome_matrix
#'
#' # Run analysis
#' results <- find_associated_metabolites(
#'   phyloseq_obj = phyloseq_obj,
#'   metabolome_mat = metabolome_mat,
#'   target_microbe_name = "Bacteroides",
#'   top_n = 30
#' )
#'
#' # View top metabolites
#' head(results$top_metabolites)
#' }
find_associated_metabolites <- function(phyloseq_obj,
                                        metabolome_mat,
                                        target_microbe_name,
                                        top_n = 20,
                                        use_normalization = TRUE) {

  cat("Step 1: Preparing data...\n")
  data_list <- prepare_data(phyloseq_obj, metabolome_mat, target_microbe_name)

  if(use_normalization) {
    cat("Step 2: Cross-omics normalization...\n")
    norm_data <- cross_omics_normalization(data_list$target_abundance,
                                           data_list$metabolome)
    analysis_microbe <- norm_data$microbe_norm
    analysis_metabolome <- norm_data$metabolome_norm
  } else {
    cat("Step 2: Skipping normalization...\n")
    analysis_microbe <- data_list$target_abundance
    analysis_metabolome <- data_list$metabolome
    norm_data <- data_list
  }

  cat("Step 3: Correlation analysis...\n")
  correlation_results <- comprehensive_correlation(analysis_microbe,
                                                   analysis_metabolome)

  cat("Step 4: Network analysis...\n")
  network_metrics <- metabolite_network_analysis(analysis_metabolome)

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

  class(results) <- "metabolite_association_results"

  cat(paste("\nAnalysis complete! Top", top_n, "metabolites identified.\n"))

  return(results)
}

#' MicrobiomeMetabolomeAnalysis: Cross-omics Association Analysis
#'
#' A comprehensive toolkit for identifying metabolites associated with specific microbes
#' through integrated multi-dimensional analysis including correlation, network topology,
#' differential expression, and machine learning approaches.
#'
#' _PACKAGE
#' @name MicrobiomeMetabolomeAnalysis
#' @import glmnet
#' @importFrom utils head
NULL

#' Find Associated Metabolites with Target Microbe
#'
#' Performs comprehensive analysis to identify metabolites associated with a specific microbe
#' using multiple analytical approaches including correlation, network analysis, differential
#' expression, and machine learning feature selection.
#'
#' @param phyloseq_obj A phyloseq object containing microbiome data
#' @param metabolome_mat A matrix or data.frame with metabolome data (samples as rows, metabolites as columns)
#' @param target_microbe_name Character string specifying the name of the target microbe
#' @param top_n Integer specifying the number of top metabolites to return (default: 20)
#' @param use_normalization Logical, whether to apply cross-omics normalization (default: TRUE)
#'
#' @return A list containing:
#' \describe{
#'   \item{top_metabolites}{Data frame of top-ranked metabolites with scores}
#'   \item{all_scores}{Complete scoring results for all metabolites}
#'   \item{correlation_results}{Detailed correlation analysis results}
#'   \item{network_metrics}{Network topology metrics for metabolites}
#'   \item{diff_results}{Differential expression analysis results}
#'   \item{ml_results}{Machine learning feature importance scores}
#'   \item{normalized_data}{Normalized data used in analysis}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load your data
#' # phyloseq_obj <- your_phyloseq_object
#' # metabolome_mat <- your_metabolome_matrix
#'
#' # Run analysis
#' results <- find_associated_metabolites.2(
#'   phyloseq_obj = phyloseq_obj,
#'   metabolome_mat = metabolome_mat,
#'   target_microbe_name = "Bacteroides",
#'   top_n = 30
#' )
#'
#' # View top metabolites
#' head(results$top_metabolites)
#' }

find_associated_metabolites.2 <- function(phyloseq_obj, metabolome_mat,
                                          target_microbe_name,
                                          top_n = 20) {

  cat("Step 1: Preparing data...\n")
  data_list <- prepare_data(phyloseq_obj, metabolome_mat, target_microbe_name)

  # cat("Step 2: Cross-omics normalization...\n")
  # norm_data <- cross_omics_normalization(data_list$target_abundance,
  #                                        data_list$metabolome)
  #
  cat("Step 3: Correlation analysis...\n")
  correlation_results <- comprehensive_correlation(data_list$target_abundance,
                                                   metabolome_mat)

  cat("Step 4: Network analysis...\n")



  network_metrics <- metabolite_network_analysis(data_list$metabolome)

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

  # 选择top代谢物
  top_metabolites <- head(final_scores, top_n)

  # 创建结果对象
  results <- list(
    top_metabolites = top_metabolites,
    all_scores = final_scores,
    correlation_results = correlation_results,
    network_metrics = network_metrics,
    diff_results = diff_results,
    ml_results = ml_results,
    normalized_data = data_list
  )

  cat(paste("\nAnalysis complete! Top", top_n, "metabolites identified.\n"))

  return(results)
}



#' Prepare Data for Analysis
#'
#' Extracts and aligns microbiome and metabolome data for analysis
#'
#' @param phyloseq_obj A phyloseq object
#' @param metabolome_mat Metabolome data matrix
#' @param target_microbe_name Name of target microbe
#'
#' @return A list containing aligned data
#' @export
prepare_data <- function(phyloseq_obj, metabolome_mat, target_microbe_name) {

  # Extract OTU table from phyloseq
  otu_table_data <- as.data.frame(otu_table(phyloseq_obj))

  # Ensure sample names are consistent
  common_samples <- intersect(rownames(metabolome_mat), colnames(otu_table_data))

  if(length(common_samples) == 0) {
    common_samples <- intersect(rownames(metabolome_mat), rownames(otu_table_data))
    otu_table_data <- t(otu_table_data)
  }

  if(length(common_samples) == 0) {
    stop("No common samples found between microbiome and metabolome data")
  }

  # Align samples
  otu_table_data <- otu_table_data[, common_samples]
  metabolome_mat <- metabolome_mat[common_samples, ]

  # Ensure metabolome_mat is numeric matrix
  if(!is.matrix(metabolome_mat)) {
    metabolome_mat <- as.matrix(metabolome_mat)
  }

  # Check and convert to numeric
  if(!is.numeric(metabolome_mat)) {
    metabolome_mat <- apply(metabolome_mat, 2, as.numeric)
  }

  # Extract target microbe abundance
  if(target_microbe_name %in% rownames(otu_table_data)) {
    target_microbe_abundance <- as.numeric(otu_table_data[target_microbe_name, ])
  } else {
    stop(paste("Target microbe", target_microbe_name, "not found in OTU table"))
  }

  return(list(
    otu_table = otu_table_data,
    metabolome = metabolome_mat,
    target_abundance = target_microbe_abundance,
    samples = common_samples
  ))
}

#' Cross-omics Data Normalization
#'
#' Normalizes microbiome and metabolome data for cross-omics comparison
#'
#' @param microbe_abundance Numeric vector of microbe abundance
#' @param metabolome_data Matrix of metabolome data
#'
#' @return A list containing normalized data
#' @export
cross_omics_normalization <- function(microbe_abundance, metabolome_data) {

  # Ensure inputs are numeric
  if(!is.numeric(microbe_abundance)) {
    microbe_abundance <- as.numeric(microbe_abundance)
  }

  if(!is.matrix(metabolome_data)) {
    metabolome_data <- as.matrix(metabolome_data)
  }

  # Remove metabolites with all NA values
  na_cols <- apply(metabolome_data, 2, function(x) all(is.na(x)))
  if(any(na_cols)) {
    warning(paste("Removing", sum(na_cols), "metabolites with all NA values"))
    metabolome_data <- metabolome_data[, !na_cols]
  }

  # CLR transform for microbiome data
  microbe_clr <- clr_transform(microbe_abundance + 1)

  # Log transform and scale metabolome data
  metabolome_log <- log2(metabolome_data + 1)
  metabolome_scaled <- scale(metabolome_log)

  # Handle any remaining NA/Inf values
  metabolome_scaled[!is.finite(metabolome_scaled)] <- 0

  # Quantile normalization for microbiome
  microbe_rank <- rank(microbe_clr) / length(microbe_clr)

  return(list(
    microbe_norm = microbe_clr,
    metabolome_norm = metabolome_scaled,
    microbe_quantile = microbe_rank
  ))
}

#' CLR Transformation
#'
#' Performs centered log-ratio transformation
#'
#' @param x Numeric vector
#'
#' @return Transformed vector
#' @keywords internal
clr_transform <- function(x) {
  x[x == 0] <- min(x[x > 0]) / 2  # Pseudocount for zeros
  log(x) - mean(log(x))
}

#' Comprehensive Correlation Analysis
#'
#' Calculates multiple correlation metrics between microbe and metabolites
#'
#' @param microbe_vec Numeric vector of microbe abundance
#' @param metabolome_mat Matrix of metabolome data
#'
#' @return Data frame with correlation results
#' @export
comprehensive_correlation <- function(microbe_vec, metabolome_mat) {

  n_metabolites <- ncol(metabolome_mat)
  correlation_results <- data.frame(
    metabolite = colnames(metabolome_mat),
    pearson_cor = numeric(n_metabolites),
    pearson_pval = numeric(n_metabolites),
    spearman_cor = numeric(n_metabolites),
    spearman_pval = numeric(n_metabolites),
    partial_cor = numeric(n_metabolites),
    stringsAsFactors = FALSE
  )

  for(i in 1:n_metabolites) {
    # Extract metabolite vector safely
    if(is.data.frame(metabolome_mat)) {
      metabolite_vec <- as.numeric(metabolome_mat[, i])
    } else if(is.matrix(metabolome_mat)) {
      metabolite_vec <- as.numeric(metabolome_mat[, i])
    } else {
      metabolite_vec <- as.numeric(unlist(metabolome_mat[, i]))
    }

    # Check for NA values
    if(all(is.na(metabolite_vec))) {
      warning(paste("Metabolite", colnames(metabolome_mat)[i], "contains all NA values"))
      correlation_results$pearson_cor[i] <- NA
      correlation_results$pearson_pval[i] <- NA
      correlation_results$spearman_cor[i] <- NA
      correlation_results$spearman_pval[i] <- NA
      next
    }

    # Pearson correlation
    pearson_test <- cor.test(microbe_vec, metabolite_vec, method = "pearson")
    correlation_results$pearson_cor[i] <- pearson_test$estimate
    correlation_results$pearson_pval[i] <- pearson_test$p.value

    # Spearman correlation
    spearman_test <- cor.test(microbe_vec, metabolite_vec, method = "spearman")
    correlation_results$spearman_cor[i] <- spearman_test$estimate
    correlation_results$spearman_pval[i] <- spearman_test$p.value

    # Partial correlation (if ppcor is available)
    if(n_metabolites > 2 && requireNamespace("ppcor", quietly = TRUE)) {
      other_metabolites <- metabolome_mat[, -i]
      if(ncol(as.matrix(other_metabolites)) > 1) {
        pc <- tryCatch({
          ppcor::pcor.test(microbe_vec, metabolite_vec,
                           other_metabolites[, 1:min(5, ncol(other_metabolites))])
        }, error = function(e) list(estimate = 0, p.value = 1))
        correlation_results$partial_cor[i] <- pc$estimate
      }
    }
  }

  # FDR correction
  correlation_results$pearson_fdr <- p.adjust(correlation_results$pearson_pval, method = "BH")
  correlation_results$spearman_fdr <- p.adjust(correlation_results$spearman_pval, method = "BH")

  return(correlation_results)
}

#' Metabolite Network Analysis
#'
#' Analyzes metabolite co-occurrence network
#'
#' @param metabolome_mat Matrix of metabolome data
#' @param correlation_threshold Threshold for including edges (default: 0.3)
#'
#' @return Data frame with network metrics
#' @export
metabolite_network_analysis <- function(metabolome_mat, correlation_threshold = 0.3) {

  # Calculate correlation matrix
  cor_matrix <- cor(metabolome_mat, use = "pairwise.complete.obs")

  # Handle NA/NaN values
  cor_matrix[!is.finite(cor_matrix)] <- 0

  # Apply threshold
  cor_matrix[abs(cor_matrix) < correlation_threshold] <- 0
  diag(cor_matrix) <- 0

  # Build graph
  g <- graph_from_adjacency_matrix(abs(cor_matrix), mode = "undirected", weighted = NULL)

  # Calculate network metrics
  network_metrics <- data.frame(
    metabolite = colnames(metabolome_mat),
    degree = degree(g),
    betweenness = betweenness(g, weights = NA),
    closeness = closeness(g, weights = NA),
    eigenvector = eigen_centrality(g, weights = NA)$vector,
    hub_score = hub_score(g, weights = NA)$vector,
    stringsAsFactors = FALSE
  )

  return(network_metrics)
}

#' Differential Analysis
#'
#' Performs differential expression analysis based on microbe abundance
#'
#' @param microbe_abundance Numeric vector of microbe abundance
#' @param metabolome_data Matrix of metabolome data
#' @param groups Optional grouping vector
#'
#' @return Data frame with differential analysis results
#' @export
differential_analysis <- function(microbe_abundance, metabolome_data, groups = NULL) {

  # Create groups if not provided
  if(is.null(groups)) {
    groups <- ifelse(microbe_abundance > median(microbe_abundance), "High", "Low")
  }

  # Use limma for differential analysis
  design <- model.matrix(~ groups)
  fit <- lmFit(t(metabolome_data), design)
  fit <- eBayes(fit)

  # Extract results
  diff_results <- topTable(fit, coef = 2, number = Inf)
  diff_results$metabolite <- rownames(diff_results)

  # Calculate fold change
  high_group <- metabolome_data[groups == "High", ]
  low_group <- metabolome_data[groups == "Low", ]

  # Handle single sample groups
  if(is.vector(high_group)) high_group <- matrix(high_group, nrow = 1)
  if(is.vector(low_group)) low_group <- matrix(low_group, nrow = 1)

  diff_results$mean_high <- colMeans(high_group)
  diff_results$mean_low <- colMeans(low_group)
  diff_results$fold_change <- diff_results$mean_high / (diff_results$mean_low + 1e-10)

  return(diff_results)
}

#' Machine Learning Feature Selection
#'
#' Uses multiple ML approaches to identify important metabolites
#'
#' @param microbe_abundance Numeric vector of microbe abundance
#' @param metabolome_data Matrix of metabolome data
#'
#' @return List of ML results
#' @export
ml_feature_selection <- function(microbe_abundance, metabolome_data) {

  results <- list()

  # Clean data for ML
  metabolome_data <- clean_data_for_ml(metabolome_data, microbe_abundance)

  # Ensure matrix format
  if(!is.matrix(metabolome_data)) {
    metabolome_data <- as.matrix(metabolome_data)
  }

  # 1. Random Forest
  tryCatch({
    rf_model <- randomForest(x = metabolome_data,
                             y = microbe_abundance,
                             ntree = 500,
                             importance = TRUE)

    importance_scores <- importance(rf_model)
    results$rf_importance <- data.frame(
      metabolite = rownames(importance_scores),
      importance = importance_scores[, "%IncMSE"],
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Random Forest failed: ", e$message)
    results$rf_importance <- NULL
  })

  # 2. LASSO regression
  tryCatch({
    lasso_cv <- cv.glmnet(as.matrix(metabolome_data),
                          microbe_abundance,
                          alpha = 1)

    lasso_coef <- coef(lasso_cv, s = "lambda.min")
    results$lasso_coef <- data.frame(
      metabolite = rownames(lasso_coef)[-1],
      coefficient = as.numeric(lasso_coef[-1, ]),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("LASSO failed: ", e$message)
    results$lasso_coef <- NULL
  })

  # 3. Elastic Net
  tryCatch({
    enet_cv <- cv.glmnet(as.matrix(metabolome_data),
                         microbe_abundance,
                         alpha = 0.5)

    enet_coef <- coef(enet_cv, s = "lambda.min")
    results$enet_coef <- data.frame(
      metabolite = rownames(enet_coef)[-1],
      coefficient = as.numeric(enet_coef[-1, ]),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Elastic Net failed: ", e$message)
    results$enet_coef <- NULL
  })

  # 4. PLS-DA (if mixOmics available)
  if(requireNamespace("mixOmics", quietly = TRUE)) {
    tryCatch({
      groups <- ifelse(microbe_abundance > median(microbe_abundance), "High", "Low")
      plsda_model <- mixOmics::plsda(metabolome_data, groups, ncomp = 2)

      vip_scores <- mixOmics::vip(plsda_model)
      results$plsda_vip <- data.frame(
        metabolite = rownames(vip_scores),
        vip_comp1 = vip_scores[, 1],
        vip_comp2 = vip_scores[, 2],
        stringsAsFactors = FALSE
      )
    }, error = function(e) {
      warning("PLS-DA failed: ", e$message)
      results$plsda_vip <- NULL
    })
  }

  return(results)
}

#' Clean Data for Machine Learning
#'
#' Removes NA values and zero variance features
#'
#' @param metabolome_data Matrix of metabolome data
#' @param microbe_abundance Numeric vector of microbe abundance
#'
#' @return Cleaned metabolome data matrix
#' @keywords internal
clean_data_for_ml <- function(metabolome_data, microbe_abundance) {

  # Remove zero variance columns
  zero_var <- apply(metabolome_data, 2, var, na.rm = TRUE) == 0
  if(any(zero_var)) {
    metabolome_data <- metabolome_data[, !zero_var]
  }

  # Handle NA values
  complete_cols <- complete.cases(t(metabolome_data))
  if(!all(complete_cols)) {
    warning(paste("Removing", sum(!complete_cols), "metabolites with missing values"))
    metabolome_data <- metabolome_data[, complete_cols]
  }

  # Remove features uncorrelated with target
  cors <- cor(metabolome_data, microbe_abundance, use = "pairwise.complete.obs")
  valid_features <- !is.na(cors) & is.finite(cors)
  metabolome_data <- metabolome_data[, valid_features]

  return(metabolome_data)
}

#' Integrated Scoring System
#'
#' Combines multiple metrics to rank metabolites
#'
#' @param correlation_results Correlation analysis results
#' @param network_metrics Network analysis results
#' @param diff_results Differential analysis results
#' @param ml_results Machine learning results
#' @param metabolome_data Original metabolome data
#' @param weights Named vector of weights for each score component
#'
#' @return Data frame with integrated scores
#' @export
integrated_scoring <- function(correlation_results, network_metrics,
                               diff_results, ml_results, metabolome_data,
                               weights = c(abundance = 0.15,
                                           correlation = 0.25,
                                           network = 0.15,
                                           differential = 0.20,
                                           ml = 0.25)) {

  # Initialize scores
  all_metabolites <- colnames(metabolome_data)
  scores <- data.frame(
    metabolite = all_metabolites,
    stringsAsFactors = FALSE
  )

  # 1. Abundance score
  abundance_score <- rank(colMeans(metabolome_data, na.rm = TRUE)) / length(all_metabolites)
  scores$abundance_score <- abundance_score[match(scores$metabolite, all_metabolites)]

  # 2. Correlation score
  scores <- merge(scores, correlation_results[, c("metabolite", "pearson_cor", "spearman_cor")],
                  by = "metabolite", all.x = TRUE)
  scores$correlation_score <- (abs(scores$pearson_cor) + abs(scores$spearman_cor)) / 2
  scores$correlation_score[is.na(scores$correlation_score)] <- 0

  # 3. Network score
  scores <- merge(scores, network_metrics[, c("metabolite", "degree", "betweenness")],
                  by = "metabolite", all.x = TRUE)
  scores$network_score <- (rank(scores$degree, na.last = "keep") +
                             rank(scores$betweenness, na.last = "keep")) / (2 * nrow(scores))
  scores$network_score[is.na(scores$network_score)] <- 0

  # 4. Differential score
  scores <- merge(scores, diff_results[, c("metabolite", "adj.P.Val", "fold_change")],
                  by = "metabolite", all.x = TRUE)
  scores$diff_score <- -log10(scores$adj.P.Val + 1e-10) * abs(log2(scores$fold_change + 1e-10))
  scores$diff_score[is.na(scores$diff_score) | !is.finite(scores$diff_score)] <- 0
  scores$diff_score <- rank(scores$diff_score) / nrow(scores)

  # 5. ML score
  scores$ml_score <- 0
  if(!is.null(ml_results$rf_importance)) {
    scores <- merge(scores, ml_results$rf_importance, by = "metabolite", all.x = TRUE)
    scores$ml_score <- rank(scores$importance, na.last = "keep") / nrow(scores)
    scores$ml_score[is.na(scores$ml_score)] <- 0
  }

  # Calculate integrated score
  scores$integrated_score <-
    weights["abundance"] * scores$abundance_score +
    weights["correlation"] * scores$correlation_score +
    weights["network"] * scores$network_score +
    weights["differential"] * scores$diff_score +
    weights["ml"] * scores$ml_score

  # Sort by integrated score
  scores <- scores[order(scores$integrated_score, decreasing = TRUE), ]
  rownames(scores) <- NULL

  return(scores)
}

#' Visualize Analysis Results
#'
#' Creates comprehensive visualizations of analysis results
#'
#' @param results Results object from find_associated_metabolites
#' @param top_n Number of top metabolites to visualize (default: 20)
#'
#' @export
#' @importFrom ggplot2 ggplot aes geom_point geom_bar geom_histogram geom_hline geom_vline
#' @importFrom ggplot2 scale_color_gradient2 theme_minimal labs coord_flip
visualize_results <- function(results, top_n = 20) {

  if(!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for visualization")
  }

  # Get top metabolites
  top_metabolites <- head(results$all_scores, top_n)

  # 1. Score heatmap
  if(requireNamespace("pheatmap", quietly = TRUE)) {
    score_matrix <- as.matrix(top_metabolites[, c("abundance_score", "correlation_score",
                                                  "network_score", "diff_score", "ml_score")])
    rownames(score_matrix) <- top_metabolites$metabolite

    pheatmap::pheatmap(t(score_matrix),
                       cluster_rows = FALSE,
                       cluster_cols = TRUE,
                       scale = "row",
                       main = "Metabolite Scores Heatmap",
                       color = colorRampPalette(c("navy", "white", "red"))(100))
  }

  # 2. Correlation scatter plot
  top_correlations <- results$correlation_results[
    results$correlation_results$metabolite %in% top_metabolites$metabolite, ]

  p1 <- ggplot(top_correlations, aes(x = pearson_cor, y = -log10(pearson_pval))) +
    geom_point(aes(size = abs(pearson_cor), color = pearson_cor)) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray") +
    theme_minimal() +
    labs(title = "Correlation Analysis Results",
         x = "Pearson Correlation",
         y = "-log10(P-value)")
  print(p1)

  # 3. Network degree distribution
  p2 <- ggplot(results$network_metrics, aes(x = degree)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    geom_vline(xintercept = median(results$network_metrics$degree),
               linetype = "dashed", color = "red") +
    theme_minimal() +
    labs(title = "Metabolite Network Degree Distribution",
         x = "Degree",
         y = "Count")
  print(p2)

  # 4. Integrated score bar plot
  p3 <- ggplot(top_metabolites[1:min(20, nrow(top_metabolites)), ],
               aes(x = reorder(metabolite, integrated_score), y = integrated_score)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top Metabolites by Integrated Score",
         x = "Metabolite",
         y = "Integrated Score")
  print(p3)

  invisible(list(p1 = p1, p2 = p2, p3 = p3))
}

#' Cross-validate Associations
#'
#' Performs k-fold cross-validation to assess association stability
#'
#' @param data_list Data list from prepare_data function
#' @param n_folds Number of folds for cross-validation (default: 5)
#'
#' @return Data frame with cross-validation results
#' @export
cross_validate_associations <- function(data_list, n_folds = 5) {

  n_samples <- length(data_list$samples)
  fold_size <- floor(n_samples / n_folds)

  cv_results <- list()

  for(i in 1:n_folds) {
    # 创建训练和测试集
    test_idx <- ((i-1)*fold_size + 1):(i*fold_size)
    train_idx <- setdiff(1:n_samples, test_idx)

    # 训练集数据
    train_microbe <- data_list$target_abundance[train_idx]
    train_metabolome <- data_list$metabolome[train_idx, ]

    # 测试集数据
    test_microbe <- data_list$target_abundance[test_idx]
    test_metabolome <- data_list$metabolome[test_idx, ]

    # 在训练集上找关联
    train_cor <- cor(train_microbe, train_metabolome)

    # 在测试集上验证
    test_cor <- cor(test_microbe, test_metabolome)

    cv_results[[i]] <- data.frame(
      metabolite = colnames(data_list$metabolome),
      train_cor = as.numeric(train_cor),
      test_cor = as.numeric(test_cor),
      fold = i
    )
  }

  # 汇总结果
  cv_summary <- do.call(rbind, cv_results) %>%
    group_by(metabolite) %>%
    summarise(
      mean_train_cor = mean(train_cor),
      sd_train_cor = sd(train_cor),
      mean_test_cor = mean(test_cor),
      sd_test_cor = sd(test_cor),
      stability = 1 - abs(mean_train_cor - mean_test_cor)
    ) %>%
    arrange(desc(stability))

  return(cv_summary)
}
