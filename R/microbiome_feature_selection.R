#' MicrobiomeFeatureSelection: Comprehensive Microbiome Feature Selection
#'
#' A comprehensive toolkit for identifying key microbial features through
#' integrated analysis including differential abundance, machine learning
#' importance, network topology, and abundance patterns.
#'
#' _PACKAGE
#' @name MicrobiomeFeatureSelection
#' @importFrom caret train trainControl varImp
#' @importFrom Hmisc rcorr
#' @importFrom utils write.csv head
NULL

#' Microbiome Feature Selection Analysis
#'
#' Performs comprehensive feature selection for microbiome data using multiple
#' approaches including differential abundance analysis, machine learning feature
#' importance, network analysis, and abundance patterns.
#'
#' @param physeq A phyloseq object containing the microbiome data
#' @param group_var Character string specifying the grouping variable in sample_data
#' @param prevalence_threshold Numeric, minimum prevalence threshold (default: 0.1)
#' @param detection_threshold Numeric, minimum detection threshold (default: 0.001)
#' @param diff_method Character, differential analysis method: "DESeq2" or "Wilcoxon" (default: "DESeq2")
#' @param ml_method Character, machine learning method: "rf", "gbm", or "boruta" (default: "rf")
#' @param cor_method Character, correlation method: "pearson" or "spearman" (default: "spearman")
#' @param cor_threshold Numeric, correlation threshold for network construction (default: 0.6)
#' @param p_threshold Numeric, p-value threshold for network edges (default: 0.05)
#' @param weights List, weights for each scoring component (default: balanced weights)
#' @param nfolds Integer, number of folds for cross-validation (default: 5)
#' @param ntree Integer, number of trees for random forest (default: 500)
#'
#' @return A list of class "microbiome_feature_selection" containing:
#' \describe{
#'   \item{final_ranking}{Data frame with summarized ranking results}
#'   \item{full_results}{Complete results with all metrics}
#'   \item{diff_results}{Differential abundance analysis results}
#'   \item{importance_results}{Machine learning importance scores}
#'   \item{network_results}{Network topology metrics}
#'   \item{abundance_results}{Abundance and prevalence statistics}
#'   \item{filtered_physeq}{Filtered phyloseq object used in analysis}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load phyloseq object
#' data(GlobalPatterns)
#'
#' # Run feature selection
#' results <- microbiome_feature_selection(
#'   physeq = GlobalPatterns,
#'   group_var = "SampleType",
#'   prevalence_threshold = 0.1,
#'   diff_method = "DESeq2",
#'   ml_method = "rf"
#' )
#'
#' # View top features
#' head(results$final_ranking)
#'
#' # Visualize results
#' plot(results)
#' }
microbiome_feature_selection <- function(physeq,
                                         group_var,
                                         prevalence_threshold = 0.1,
                                         detection_threshold = 0.001,
                                         diff_method = c("DESeq2", "Wilcoxon"),
                                         ml_method = c("rf", "gbm", "boruta"),
                                         cor_method = c("spearman", "pearson"),
                                         cor_threshold = 0.6,
                                         p_threshold = 0.05,
                                         weights = list(abundance = 0.3,
                                                        importance = 0.3,
                                                        differential = 0.2,
                                                        network = 0.2),
                                         nfolds = 5,
                                         ntree = 500) {

  # Validate inputs
  diff_method <- match.arg(diff_method)
  ml_method <- match.arg(ml_method)
  cor_method <- match.arg(cor_method)

  if(!inherits(physeq, "phyloseq")) {
    stop("Input must be a phyloseq object")
  }

  if(!group_var %in% sample_variables(physeq)) {
    stop(paste("Group variable", group_var, "not found in sample data"))
  }

  # Normalize weights
  weight_sum <- sum(unlist(weights))
  weights <- lapply(weights, function(x) x/weight_sum)

  cat("=====================================\n")
  cat("Microbiome Feature Selection Analysis\n")
  cat("=====================================\n\n")

  # Print basic information
  cat(sprintf("Samples: %d\n", nsamples(physeq)))
  cat(sprintf("OTUs/ASVs: %d\n", ntaxa(physeq)))
  cat(sprintf("Group variable: %s\n", group_var))
  cat(sprintf("Groups: %s\n\n", paste(unique(sample_data(physeq)[[group_var]]),
                                      collapse = ", ")))

  # Step 1: Data preparation
  cat("Step 1: Data preparation...\n")
  data_list <- prepare_phyloseq_data(physeq, group_var,
                                     prevalence_threshold,
                                     detection_threshold)

  # Step 2: Differential analysis
  cat("\nStep 2: Differential analysis...\n")
  diff_results <- perform_differential_analysis(data_list$physeq,
                                                group_var,
                                                diff_method)

  # Step 3: Machine learning importance
  cat("\nStep 3: Machine learning importance assessment...\n")
  importance_df <- assess_ml_importance(data_list$physeq_relative,
                                        group_var,
                                        ml_method,
                                        nfolds,
                                        ntree)

  # Step 4: Network analysis
  cat("\nStep 4: Network analysis...\n")
  network_df <- perform_network_analysis(data_list$physeq_relative,
                                         cor_method,
                                         cor_threshold,
                                         p_threshold)

  # Step 5: Calculate abundance
  cat("\nStep 5: Calculating abundance and prevalence...\n")
  abundance_df <- calculate_abundance_metrics(data_list$physeq)

  # Step 6: Integrated ranking
  cat("\nStep 6: Integrated weighted ranking...\n")
  ranking_results <- perform_weighted_ranking(diff_results,
                                              importance_df,
                                              network_df,
                                              abundance_df,
                                              weights)

  cat("\nAnalysis complete!\n")
  cat(sprintf("Significant OTUs (p<0.05): %d\n",
              sum(diff_results$pvalue < 0.05, na.rm = TRUE)))
  cat(sprintf("Highly connected OTUs (degree>5): %d\n",
              sum(network_df$Degree > 5)))

  # Create result object
  results <- list(
    final_ranking = ranking_results$summary_results,
    full_results = ranking_results$full_results,
    diff_results = diff_results,
    importance_results = importance_df,
    network_results = network_df,
    abundance_results = abundance_df,
    filtered_physeq = data_list$physeq,
    parameters = list(
      group_var = group_var,
      prevalence_threshold = prevalence_threshold,
      detection_threshold = detection_threshold,
      diff_method = diff_method,
      ml_method = ml_method,
      cor_method = cor_method,
      cor_threshold = cor_threshold,
      p_threshold = p_threshold,
      weights = weights
    )
  )

  class(results) <- c("microbiome_feature_selection", "list")

  return(results)
}

#' Prepare Data from Phyloseq Object
#'
#' Filters and prepares data from a phyloseq object for analysis.
#'
#' @param physeq A phyloseq object
#' @param group_var Character, grouping variable name
#' @param prevalence_threshold Numeric, minimum prevalence
#' @param detection_threshold Numeric, minimum detection
#'
#' @return List containing processed data
#' @keywords internal
prepare_phyloseq_data <- function(physeq, group_var,
                                  prevalence_threshold = 0.1,
                                  detection_threshold = 0.001) {

  cat("Preparing data from phyloseq object...\n")

  # Filter low abundance/prevalence OTUs
  if(prevalence_threshold > 0 || detection_threshold > 0) {
    if(requireNamespace("microbiome", quietly = TRUE)) {
      physeq_filtered <- microbiome::core(physeq,
                                          detection = detection_threshold,
                                          prevalence = prevalence_threshold)
    } else {
      # Manual filtering if microbiome package not available
      otu_tab <- as.matrix(otu_table(physeq))
      if(!taxa_are_rows(physeq)) {
        otu_tab <- t(otu_tab)
      }

      # Calculate prevalence
      prevalence <- rowSums(otu_tab > detection_threshold) / ncol(otu_tab)
      keep_taxa <- prevalence >= prevalence_threshold

      physeq_filtered <- prune_taxa(keep_taxa, physeq)
    }

    cat(sprintf("Retained %d OTUs/ASVs after filtering (original: %d)\n",
                ntaxa(physeq_filtered), ntaxa(physeq)))
  } else {
    physeq_filtered <- physeq
  }

  # Extract data
  otu_counts <- as.data.frame(otu_table(physeq_filtered))
  if(!taxa_are_rows(physeq_filtered)) {
    otu_counts <- t(otu_counts)
  }

  # Calculate relative abundance
  physeq_relative <- transform_sample_counts(physeq_filtered,
                                             function(x) x/sum(x))

  otu_relative <- as.data.frame(otu_table(physeq_relative))
  if(!taxa_are_rows(physeq_relative)) {
    otu_relative <- t(otu_relative)
  }

  # Extract metadata
  metadata <- as.data.frame(sample_data(physeq_filtered))

  return(list(
    physeq = physeq_filtered,
    physeq_relative = physeq_relative,
    otu_counts = t(otu_counts),
    otu_relative = t(otu_relative),
    metadata = metadata,
    group_var = group_var,
    taxonomy = tax_table(physeq_filtered)
  ))
}

#' Perform Differential Analysis
#'
#' Performs differential abundance analysis using specified method.
#'
#' @param physeq A phyloseq object
#' @param group_var Character, grouping variable
#' @param method Character, analysis method
#'
#' @return Data frame with differential analysis results
#' @keywords internal
perform_differential_analysis2 <- function(physeq, group_var,
                                           method = "DESeq2") {

  cat(sprintf("Using %s for differential analysis...\n", method))

  if(method == "DESeq2") {
    if(!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("DESeq2 package required for this method")
    }

    # Convert to DESeq2
    dds <- phyloseq_to_deseq2(physeq, as.formula(paste("~", group_var)))

    # Calculate size factors
    geoMeans <- apply(counts(dds), 1, function(row) {
      if (all(row == 0)) { 0 } else { exp(mean(log(row[row != 0]))) }
    })
    dds <- estimateSizeFactors(dds, geoMeans = geoMeans)

    # Run DESeq2
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    res <- results(dds)

    # Extract results
    diff_results <- data.frame(
      OTU = rownames(res),
      pvalue = res$pvalue,
      padj = res$padj,
      log2FoldChange = res$log2FoldChange,
      baseMean = res$baseMean,
      stringsAsFactors = FALSE
    )

  } else if(method == "Wilcoxon") {
    # Wilcoxon rank sum test
    metadata <- sample_data(physeq)
    groups <- unique(metadata[[group_var]])

    if(length(groups) != 2) {
      stop("Wilcoxon test requires exactly 2 groups")
    }

    # Get OTU table
    otu_tab <- as.data.frame(otu_table(physeq))
    if(!taxa_are_rows(physeq)) {
      otu_tab <- t(otu_tab)
    }

    # Group samples
    group1_samples <- sample_names(physeq)[metadata[[group_var]] == groups[1]]
    group2_samples <- sample_names(physeq)[metadata[[group_var]] == groups[2]]

    # Perform Wilcoxon test
    pvalues <- apply(otu_tab, 1, function(x) {
      tryCatch({
        wilcox.test(x[group1_samples], x[group2_samples])$p.value
      }, error = function(e) NA)
    })

    # Calculate fold change
    mean_group1 <- rowMeans(otu_tab[, group1_samples, drop = FALSE])
    mean_group2 <- rowMeans(otu_tab[, group2_samples, drop = FALSE])
    log2FC <- log2((mean_group2 + 1) / (mean_group1 + 1))

    diff_results <- data.frame(
      OTU = names(pvalues),
      pvalue = pvalues,
      padj = p.adjust(pvalues, method = "BH"),
      log2FoldChange = log2FC,
      baseMean = rowMeans(otu_tab),
      stringsAsFactors = FALSE
    )
  }

  # Remove NA values
  diff_results <- diff_results[!is.na(diff_results$pvalue), ]

  return(diff_results)
}



perform_differential_analysis <- function(physeq, group_var,
                                          method = "DESeq2",
                                          data_type = "auto") {

  cat(sprintf("Using %s for differential analysis...\n", method))

  # 自动检测数据类型
  if(data_type == "auto") {
    otu_tab_check <- as.matrix(otu_table(physeq))
    if(!taxa_are_rows(physeq)) {
      otu_tab_check <- t(otu_tab_check)
    }

    # 检查是否为相对丰度数据（总和接近1或100）
    sample_sums <- colSums(otu_tab_check)
    if(all(abs(sample_sums - 1) < 0.1) || all(abs(sample_sums - 100) < 10)) {
      data_type <- "relative"
      cat("Detected relative abundance data\n")
    } else if(all(otu_tab_check == floor(otu_tab_check))) {
      data_type <- "counts"
      cat("Detected count data\n")
    } else {
      data_type <- "relative"
      cat("Assuming relative abundance data\n")
    }
  }

  if(method == "DESeq2") {

    # 处理相对丰度数据
    if(data_type == "relative") {
      cat("Converting relative abundance to pseudo-counts for DESeq2...\n")

      # 获取OTU表
      otu_tab <- as.matrix(otu_table(physeq))
      if(!taxa_are_rows(physeq)) {
        otu_tab <- t(otu_tab)
      }

      # 转换相对丰度为伪计数
      # 方法1：乘以一个大的缩放因子（如10000或100000）
      scale_factor <- 10000

      # 如果是百分比（0-100），先转换为比例（0-1）
      if(max(otu_tab) > 1.5) {
        otu_tab <- otu_tab / 100
      }

      # 转换为伪计数并四舍五入
      pseudo_counts <- round(otu_tab * scale_factor)

      # 创建新的phyloseq对象
      if(taxa_are_rows(physeq)) {
        otu_table(physeq) <- otu_table(pseudo_counts, taxa_are_rows = TRUE)
      } else {
        otu_table(physeq) <- otu_table(t(pseudo_counts), taxa_are_rows = TRUE)
      }

      cat(sprintf("Scaled relative abundances by factor of %d\n", scale_factor))
    }

    if(!requireNamespace("DESeq2", quietly = TRUE)) {
      stop("DESeq2 package required for this method")
    }

    # Convert to DESeq2
    dds <- phyloseq_to_deseq2(physeq, as.formula(paste("~", group_var)))

    # Calculate size factors
    geoMeans <- apply(counts(dds), 1, function(row) {
      if (all(row == 0)) { 0 } else { exp(mean(log(row[row != 0]))) }
    })

    # 对于伪计数，使用更宽松的size factor估计
    if(data_type == "relative") {
      dds <- estimateSizeFactors(dds, type = "poscounts")
    } else {
      dds <- estimateSizeFactors(dds, geoMeans = geoMeans)
    }

    # Run DESeq2
    dds <- DESeq(dds, fitType = "local", quiet = TRUE)
    res <- results(dds)

    # Extract results
    diff_results <- data.frame(
      OTU = rownames(res),
      pvalue = res$pvalue,
      padj = res$padj,
      log2FoldChange = res$log2FoldChange,
      baseMean = res$baseMean,
      stringsAsFactors = FALSE
    )

  } else if(method == "Wilcoxon") {
    # Wilcoxon方法可以直接处理相对丰度数据

    metadata <- sample_data(physeq)
    groups <- unique(metadata[[group_var]])

    if(length(groups) != 2) {
      stop("Wilcoxon test requires exactly 2 groups")
    }

    # Get OTU table
    otu_tab <- as.data.frame(otu_table(physeq))
    if(!taxa_are_rows(physeq)) {
      otu_tab <- t(otu_tab)
    }

    # Group samples
    group1_samples <- sample_names(physeq)[metadata[[group_var]] == groups[1]]
    group2_samples <- sample_names(physeq)[metadata[[group_var]] == groups[2]]

    # Perform Wilcoxon test
    pvalues <- apply(otu_tab, 1, function(x) {
      tryCatch({
        wilcox.test(x[group1_samples], x[group2_samples])$p.value
      }, error = function(e) NA)
    })

    # Calculate fold change (适用于相对丰度)
    mean_group1 <- rowMeans(otu_tab[, group1_samples, drop = FALSE])
    mean_group2 <- rowMeans(otu_tab[, group2_samples, drop = FALSE])

    # 添加小的伪计数避免除零
    pseudocount <- min(otu_tab[otu_tab > 0]) / 10
    log2FC <- log2((mean_group2 + pseudocount) / (mean_group1 + pseudocount))

    diff_results <- data.frame(
      OTU = names(pvalues),
      pvalue = pvalues,
      padj = p.adjust(pvalues, method = "BH"),
      log2FoldChange = log2FC,
      baseMean = rowMeans(otu_tab),
      stringsAsFactors = FALSE
    )

  } else if(method == "ALDEx2") {
    # ALDEx2方法特别适合相对丰度数据
    if(!requireNamespace("ALDEx2", quietly = TRUE)) {
      cat("ALDEx2 not available, falling back to Wilcoxon test\n")
      method <- "Wilcoxon"
      return(perform_differential_analysis(physeq, group_var, "Wilcoxon", data_type))
    }

    otu_tab <- as.data.frame(otu_table(physeq))
    if(!taxa_are_rows(physeq)) {
      otu_tab <- t(otu_tab)
    }

    metadata <- sample_data(physeq)
    conditions <- metadata[[group_var]]

    # 如果是相对丰度，转换为伪计数
    if(data_type == "relative") {
      if(max(otu_tab) > 1.5) {
        otu_tab <- otu_tab / 100
      }
      otu_tab <- round(otu_tab * 10000)
    }

    # Run ALDEx2
    aldex_result <- ALDEx2::aldex(
      reads = otu_tab,
      conditions = conditions,
      test = "t",
      effect = TRUE,
      denom = "all"
    )

    diff_results <- data.frame(
      OTU = rownames(aldex_result),
      pvalue = aldex_result$we.ep,
      padj = aldex_result$we.eBH,
      log2FoldChange = aldex_result$effect,
      baseMean = rowMeans(otu_tab),
      stringsAsFactors = FALSE
    )

  } else if(method == "ANCOM-BC") {
    # ANCOM-BC适合处理组成数据
    if(!requireNamespace("ANCOMBC", quietly = TRUE)) {
      cat("ANCOMBC not available, falling back to Wilcoxon test\n")
      method <- "Wilcoxon"
      return(perform_differential_analysis(physeq, group_var, "Wilcoxon", data_type))
    }

    # Run ANCOM-BC
    ancom_result <- ANCOMBC::ancombc(
      phyloseq = physeq,
      formula = group_var,
      p_adj_method = "BH",
      zero_cut = 0.90,
      lib_cut = 1000,
      group = group_var,
      struc_zero = TRUE,
      neg_lb = TRUE,
      tol = 1e-5,
      max_iter = 100,
      conserve = TRUE,
      alpha = 0.05,
      global = FALSE
    )

    diff_results <- data.frame(
      OTU = rownames(ancom_result$res$beta),
      pvalue = ancom_result$res$p_val[,1],
      padj = ancom_result$res$q_val[,1],
      log2FoldChange = ancom_result$res$beta[,1] / log(2),
      baseMean = apply(otu_table(physeq), 1, mean),
      stringsAsFactors = FALSE
    )
  }

  # Remove NA values
  diff_results <- diff_results[!is.na(diff_results$pvalue), ]

  # 添加数据类型信息
  attr(diff_results, "data_type") <- data_type
  attr(diff_results, "method") <- method

  cat(sprintf("Differential analysis completed. Found %d features.\n", nrow(diff_results)))

  return(diff_results)
}

# 辅助函数：验证数据类型
check_data_type <- function(physeq) {
  otu_tab <- as.matrix(otu_table(physeq))
  if(!taxa_are_rows(physeq)) {
    otu_tab <- t(otu_tab)
  }

  # 检查样本总和
  sample_sums <- colSums(otu_tab)

  # 判断数据类型
  if(all(abs(sample_sums - 1) < 0.1)) {
    return("relative_proportion")  # 0-1比例
  } else if(all(abs(sample_sums - 100) < 10)) {
    return("relative_percentage")  # 0-100百分比
  } else if(all(otu_tab == floor(otu_tab))) {
    return("counts")  # 整数计数
  } else {
    return("unknown")  # 未知类型
  }
}
#' Assess Machine Learning Importance
#'
#' Evaluates feature importance using machine learning methods.
#'
#' @param physeq_relative A phyloseq object with relative abundances
#' @param group_var Character, grouping variable
#' @param method Character, ML method
#' @param nfolds Integer, number of CV folds
#' @param ntree Integer, number of trees for RF
#'
#' @return Data frame with importance scores
#' @keywords internal
assess_ml_importance <- function(physeq_relative, group_var,
                                 method = "rf", nfolds = 5, ntree = 500) {

  cat(sprintf("Using %s for importance assessment...\n", method))

  # Prepare data
  otu_data <- as.data.frame(t(otu_table(physeq_relative)))
  if(!taxa_are_rows(physeq_relative)) {
    otu_data <- as.data.frame(otu_table(physeq_relative))
  }

  metadata <- sample_data(physeq_relative)
  ml_data <- cbind(otu_data, Group = metadata[[group_var]])

  if(method == "rf") {
    if(!requireNamespace("caret", quietly = TRUE)) {
      stop("caret package required for machine learning analysis")
    }

    set.seed(123)

    # Cross-validation setup
    ctrl <- trainControl(
      method = "cv",
      number = nfolds,
      classProbs = TRUE,
      summaryFunction = defaultSummary
    )

    # Train model
    rf_model <- tryCatch({
      train(
        Group ~ .,
        data = ml_data,
        method = "rf",
        trControl = ctrl,
        ntree = ntree,
        importance = TRUE,
        metric = "Accuracy"
      )
    }, error = function(e) {
      warning("RF training failed, using fallback method")
      NULL
    })

    if(!is.null(rf_model)) {
      # Extract importance
      importance_scores <- varImp(rf_model, scale = TRUE)$importance

      importance_df <- data.frame(
        OTU = rownames(importance_scores),
        Importance = importance_scores[,1],
        stringsAsFactors = FALSE
      )

      cat(sprintf("Model accuracy: %.2f%%\n",
                  max(rf_model$results$Accuracy) * 100))
    } else {
      # Fallback: simple correlation-based importance
      importance_df <- data.frame(
        OTU = colnames(otu_data),
        Importance = rep(0.5, ncol(otu_data)),
        stringsAsFactors = FALSE
      )
    }

  } else if(method == "gbm") {
    if(!requireNamespace("gbm", quietly = TRUE)) {
      stop("gbm package required for this method")
    }

    ctrl <- trainControl(
      method = "cv",
      number = nfolds,
      classProbs = TRUE
    )

    gbm_model <- train(
      Group ~ .,
      data = ml_data,
      method = "gbm",
      trControl = ctrl,
      verbose = FALSE
    )

    importance_scores <- varImp(gbm_model, scale = TRUE)$importance

    importance_df <- data.frame(
      OTU = rownames(importance_scores),
      Importance = importance_scores[,1],
      stringsAsFactors = FALSE
    )

  } else if(method == "boruta") {
    if(!requireNamespace("Boruta", quietly = TRUE)) {
      stop("Boruta package required for this method")
    }

    set.seed(123)
    boruta_result <- Boruta::Boruta(Group ~ ., data = ml_data, doTrace = 0)

    importance_df <- data.frame(
      OTU = names(boruta_result$ImpHistory[,1]),
      Importance = apply(boruta_result$ImpHistory, 1, median, na.rm = TRUE),
      Decision = boruta_result$finalDecision,
      stringsAsFactors = FALSE
    )
  }

  # Normalize importance scores
  if(max(importance_df$Importance) > min(importance_df$Importance)) {
    importance_df$Importance_scaled <-
      (importance_df$Importance - min(importance_df$Importance)) /
      (max(importance_df$Importance) - min(importance_df$Importance))
  } else {
    importance_df$Importance_scaled <- 0.5
  }

  return(importance_df)
}

#' Perform Network Analysis
#'
#' Constructs and analyzes microbial co-occurrence network.
#'
#' @param physeq_relative A phyloseq object with relative abundances
#' @param cor_method Character, correlation method
#' @param cor_threshold Numeric, correlation threshold
#' @param p_threshold Numeric, p-value threshold
#'
#' @return Data frame with network metrics
#' @keywords internal
perform_network_analysis <- function(physeq_relative,
                                     cor_method = "spearman",
                                     cor_threshold = 0.6,
                                     p_threshold = 0.05) {

  cat("Constructing microbial correlation network...\n")

  # Get OTU table
  otu_data <- as.data.frame(otu_table(physeq_relative))
  if(taxa_are_rows(physeq_relative)) {
    otu_data <- t(otu_data)
  }

  # Calculate correlation matrix
  if(requireNamespace("Hmisc", quietly = TRUE)) {
    cor_matrix <- Hmisc::rcorr(as.matrix(otu_data), type = cor_method)
    cor_values <- cor_matrix$r
    cor_pvalues <- cor_matrix$P
  } else {
    # Fallback to base R
    cor_values <- cor(otu_data, method = cor_method)
    cor_pvalues <- matrix(0.05, nrow = ncol(otu_data), ncol = ncol(otu_data))
  }

  # Set diagonal
  diag(cor_values) <- 0
  diag(cor_pvalues) <- 1

  # Create adjacency matrix
  adj_matrix <- (abs(cor_values) >= cor_threshold) & (cor_pvalues < p_threshold)
  adj_matrix <- adj_matrix * 1

  # Calculate network metrics
  degree_df <- data.frame(
    OTU = colnames(otu_data),
    Degree = colSums(adj_matrix),
    Weighted_degree = colSums(abs(cor_values) * adj_matrix),
    stringsAsFactors = FALSE
  )

  # Build igraph object for additional metrics
  if(sum(adj_matrix) > 0 && requireNamespace("igraph", quietly = TRUE)) {
    g <- igraph::graph_from_adjacency_matrix(
      adj_matrix * abs(cor_values),
      mode = "undirected",
      weighted = TRUE,
      diag = FALSE
    )

    degree_df$Betweenness <- igraph::betweenness(g)
    degree_df$Closeness <- igraph::closeness(g)
    degree_df$Eigenvector <- igraph::eigen_centrality(g)$vector
    degree_df$Hub_score <- igraph::hub_score(g)$vector

    # Identify modules
    if(igraph::vcount(g) > 10) {
      communities <- igraph::cluster_fast_greedy(g)
      degree_df$Module <- igraph::membership(communities)
    } else {
      degree_df$Module <- 1
    }

  } else {
    degree_df$Betweenness <- 0
    degree_df$Closeness <- 0
    degree_df$Eigenvector <- 0
    degree_df$Hub_score <- 0
    degree_df$Module <- 1
  }

  cat(sprintf("Network contains %d edges (correlation threshold: %.2f, p-value threshold: %.3f)\n",
              sum(adj_matrix)/2, cor_threshold, p_threshold))

  return(degree_df)
}

#' Calculate Abundance Metrics
#'
#' Calculates abundance and prevalence statistics for OTUs.
#'
#' @param physeq A phyloseq object
#'
#' @return Data frame with abundance metrics
#' @keywords internal
calculate_abundance_metrics <- function(physeq) {

  # Get relative abundance
  physeq_rel <- transform_sample_counts(physeq, function(x) x/sum(x))

  # Get OTU table
  otu_tab <- as.data.frame(otu_table(physeq_rel))
  if(taxa_are_rows(physeq_rel)) {
    otu_tab <- t(otu_tab)
  }

  abundance_df <- data.frame(
    OTU = colnames(otu_tab),
    Mean_abundance = colMeans(otu_tab),
    Median_abundance = apply(otu_tab, 2, median),
    Max_abundance = apply(otu_tab, 2, max),
    SD_abundance = apply(otu_tab, 2, sd),
    Prevalence = colSums(otu_tab > 0) / nrow(otu_tab),
    stringsAsFactors = FALSE
  )

  # Add taxonomy information if available
  if(!is.null(tax_table(physeq, errorIfNULL = FALSE))) {
    tax_info <- as.data.frame(tax_table(physeq))
    last_rank <- tail(rank_names(physeq), 1)
    if(last_rank %in% colnames(tax_info)) {
      abundance_df$Taxonomy <- tax_info[abundance_df$OTU, last_rank]
    }
  }

  return(abundance_df)
}

#' Perform Weighted Ranking
#'
#' Integrates multiple metrics using weighted scoring.
#'
#' @param diff_results Differential analysis results
#' @param importance_df ML importance results
#' @param network_df Network analysis results
#' @param abundance_df Abundance metrics
#' @param weights List of weights
#'
#' @return List with ranking results
#' @keywords internal
perform_weighted_ranking <- function(diff_results, importance_df,
                                     network_df, abundance_df, weights) {

  cat("Performing integrated weighted ranking...\n")

  # Merge all results
  combined_df <- diff_results

  if(!is.null(importance_df)) {
    combined_df <- merge(combined_df, importance_df, by = "OTU", all = TRUE)
  }
  if(!is.null(network_df)) {
    combined_df <- merge(combined_df, network_df, by = "OTU", all = TRUE)
  }
  if(!is.null(abundance_df)) {
    combined_df <- merge(combined_df, abundance_df, by = "OTU", all = TRUE)
  }

  # Handle missing values
  numeric_cols <- sapply(combined_df, is.numeric)
  combined_df[, numeric_cols][is.na(combined_df[, numeric_cols])] <- 0

  # Normalization function
  normalize <- function(x) {
    if(length(unique(x)) == 1) return(rep(0.5, length(x)))
    if(max(x) == min(x)) return(rep(0.5, length(x)))
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }

  # Calculate dimension scores
  # Differential score
  combined_df$Diff_score_p <- normalize(1 - combined_df$pvalue)
  combined_df$Diff_score_fc <- normalize(abs(combined_df$log2FoldChange))
  combined_df$Diff_score <- (combined_df$Diff_score_p + combined_df$Diff_score_fc) / 2

  # Importance score
  if("Importance_scaled" %in% colnames(combined_df)) {
    combined_df$Importance_score <- combined_df$Importance_scaled
  } else {
    combined_df$Importance_score <- 0.5
  }

  # Network score
  if("Degree" %in% colnames(combined_df)) {
    combined_df$Network_score_degree <- normalize(combined_df$Degree)
    combined_df$Network_score_between <- normalize(combined_df$Betweenness)
    combined_df$Network_score <- (combined_df$Network_score_degree +
                                    combined_df$Network_score_between) / 2
  } else {
    combined_df$Network_score <- 0.5
  }

  # Abundance score
  if("Mean_abundance" %in% colnames(combined_df)) {
    combined_df$Abundance_score_mean <- normalize(combined_df$Mean_abundance)
    combined_df$Abundance_score_prev <- normalize(combined_df$Prevalence)
    combined_df$Abundance_score <- (combined_df$Abundance_score_mean * 0.7 +
                                      combined_df$Abundance_score_prev * 0.3)
  } else {
    combined_df$Abundance_score <- 0.5
  }

  # Calculate weighted score
  combined_df$Weighted_score <-
    weights$abundance * combined_df$Abundance_score +
    weights$importance * combined_df$Importance_score +
    weights$differential * combined_df$Diff_score +
    weights$network * combined_df$Network_score

  # Sort by weighted score
  combined_df <- combined_df[order(combined_df$Weighted_score, decreasing = TRUE), ]

  # Add rank
  combined_df$Rank <- 1:nrow(combined_df)

  # Select key columns for summary
  key_columns <- c("Rank", "OTU", "Weighted_score",
                   "Mean_abundance", "Prevalence",
                   "pvalue", "padj", "log2FoldChange",
                   "Importance", "Degree", "Betweenness")

  if("Taxonomy" %in% colnames(combined_df)) {
    key_columns <- c(key_columns, "Taxonomy")
  }

  # Ensure columns exist
  key_columns <- key_columns[key_columns %in% colnames(combined_df)]

  # Create summary table
  combined_df_clean <- combined_df[, key_columns]

  return(list(
    full_results = combined_df,
    summary_results = combined_df_clean
  ))
}

#' Visualize Microbiome Feature Selection Results
#'
#' Creates comprehensive visualizations of feature selection results.
#'
#' @param results A microbiome_feature_selection object
#' @param top_n Integer, number of top features to display (default: 20)
#'
#' @return A combined plot object
#' @export
visualize_results <- function(results, top_n = 20) {
  library(ggplot2)
  library(patchwork)
  library(reshape2)

  final_df <- results$full_results

  # 1. 综合得分条形图
  p1 <- ggplot(head(final_df, top_n),
               aes(x = reorder(OTU, Weighted_score), y = Weighted_score)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    coord_flip() +
    labs(title = paste("Top", top_n, "微生物特征综合得分"),
         x = "OTU/ASV", y = "加权得分") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 8))

  # 2. 各维度得分热图
  score_cols <- c("Abundance_score", "Importance_score",
                  "Diff_score", "Network_score")

  score_data <- final_df[1:min(top_n, nrow(final_df)), c("OTU", score_cols)]
  score_data_long <- melt(score_data, id.vars = "OTU",
                          variable.name = "Score_type",
                          value.name = "Score")

  # 重命名得分类型
  score_data_long$Score_type <- factor(score_data_long$Score_type,
                                       levels = score_cols,
                                       labels = c("丰度", "重要性", "差异", "网络"))

  p2 <- ggplot(score_data_long,
               aes(x = Score_type, y = factor(OTU, levels = rev(final_df$OTU[1:min(top_n, nrow(final_df))])),
                   fill = Score)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                         midpoint = 0.5, limits = c(0, 1)) +
    labs(title = "各维度得分热图",
         x = "得分类型", y = "OTU/ASV", fill = "得分") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 8))

  # 3. 多维度散点图
  p3 <- ggplot(final_df,
               aes(x = Importance_score, y = Diff_score,
                   size = Abundance_score, color = Network_score)) +
    geom_point(alpha = 0.6) +
    scale_color_gradient(low = "blue", high = "red", name = "网络得分") +
    scale_size_continuous(range = c(1, 8), name = "丰度得分") +
    labs(title = "多维度得分分布",
         x = "重要性得分", y = "差异得分") +
    theme_minimal() +
    theme(legend.position = "right")

  # 4. 添加前几个OTU的标签
  if(top_n <= 10) {
    top_otus <- head(final_df$OTU, 5)
    top_data <- final_df[final_df$OTU %in% top_otus, ]
    p3 <- p3 +
      geom_text(data = top_data,
                aes(label = OTU),
                hjust = -0.1, vjust = -0.1, size = 3)
  }

  # 组合图
  combined_plot <- (p1 | p2) / p3 +
    plot_annotation(title = "微生物特征综合分析结果",
                    theme = theme(plot.title = element_text(size = 16, hjust = 0.5)))

  return(combined_plot)
}
