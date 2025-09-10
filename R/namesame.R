
# library(EasyMultiOmics)
# ?count_to_mass
# ?find_associated_metabolites
# count_to_mass = find_associated_metabolites
# mass_to_count =find_associated_microbes
# count_selection = microbiome_feature_selection
# ms_feature_selection = mass_selection
# multi_omics_alignment #

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
#' @importFrom utils head
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
feature_selection <- function(results,
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


#' MicrobiomeFeatureSelection: Comprehensive Microbiome Feature Selection
#'
#' A comprehensive toolkit for identifying key microbial features through
#' integrated analysis including differential abundance, machine learning
#' importance, network topology, and abundance patterns.
#'
#' _PACKAGE
#' @name MicrobiomeFeatureSelection
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
#' # Example workflow:
#' ps <- ps.16s %>% subset_samples(Group %in% c("WT","KO"))
#' ps <- ps %>% tax_glom("Genus")
#'
#' res <- count_selection(ps,
#'                      group_var = "Group",
#'                      prevalence_threshold = 0.1,
#'                      detection_threshold = 0.001,
#'                      diff_method = "DESeq2",
#'                      ml_method = "rf",
#'                      cor_method = "spearman",
#'                      cor_threshold = 0.6,
#'                      p_threshold = 0.05,
#'                      weights = list(abundance = 0.3,
#'                                     importance = 0.3,
#'                                     differential = 0.2,
#'                                     network = 0.2))
#' top_features <- res$final_ranking %>% head(100)
#' }
#'
count_selection <- function(physeq,
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
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' \dontrun{
#' ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","OE"))
#' tab = ps2 %>% vegan_otu() %>% as.data.frame()
#' ps = ps.16s %>% subset_samples.wt("Group",c("WT","OE")) %>% tax_glom_wt(6)
#' results <- count_to_mass(
#'  phyloseq_obj = ps,
#'  metabolome_mat = tab,
#'  target_microbe_name = "Actinocorallia" ,
#'  top_n = 30
#')
#' }
count_to_mass <- function(phyloseq_obj,
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




#' MicrobeFinderFromMetabolite: find microbes associated with a target metabolite
#'
#' 给定一个代谢物（非靶向代谢组矩阵中的一列），在 phyloseq 微生物数据中寻找潜在互作的微生物特征。
#' 采用多维证据：相关（CLR/SparCC 可选）、网络中心性、差异丰度（代谢物高低分组）、机器学习
#' （用微生物预测该代谢物，得到微生物重要性），以及丰度/流行率优先。
#'
#' 依赖包可选：DESeq2（差异分析更稳），mixOmics（可做 PLS 回归 VIP），
#' igraph（网络）、glmnet / randomForest（机器学习）。没装会自动降级。
#'
#' @param phyloseq_obj phyloseq对象（微生物计数或相对丰度；建议原始计数）
#' @param metabolome_mat 矩阵/数据框（样本×代谢物）
#' @param target_metabolite_name 字符串：目标代谢物列名
#' @param top_n 返回Top微生物数
#' @param cor_method 相关方法："spearman"（默认，基于CLR）、"pearson"、"sparcc"（需 SpiecEasi）
#' @param prevalence_threshold 微生物流行率阈值（默认0.1）
#' @param detection_threshold 相对丰度检出阈值（默认1e-4；用于计算流行率）
#' @param use_clr 是否对微生物做CLR（默认TRUE；DESeq2差异仍用原始计数）
#' @param group_split 策略将代谢物分为High/Low：c("median","tertile","quantile_0.4_0.6")
#' @return list：top_microbes, all_scores, correlation_results, network_metrics,
#'               diff_results, ml_results, aligned_data（样本对齐/变换后的数据）
#' \dontrun{
#' ps2 = ps.ms %>% subset_samples.wt("Group",c("WT","OE"))
#' tab = ps2 %>% vegan_otu() %>%
#'  as.data.frame()
#'ps = ps.16s %>% subset_samples.wt("Group",c("WT","OE")) %>% tax_glom_wt(6)
#' res <- mass_to_count(
#'   phyloseq_obj = GlobalPatterns,
#'   metabolome_mat = metabolome_data,
#'   target_metabolite_name = "Metabolite_X",
#'   top_n = 15,
#'   cor_method = "spearman"
#' )
#'
#' head(res$top_features)
#' }
#'
#' @export

mass_to_count <- function(
    phyloseq_obj,
    metabolome_mat,
    target_metabolite_name,
    top_n = 30,
    cor_method = c("spearman","pearson","sparcc"),
    prevalence_threshold = 0.10,
    detection_threshold  = 1e-4,
    use_clr = TRUE,
    group_split = c("median","tertile","quantile_0.4_0.6")
){
  cor_method  <- match.arg(cor_method)
  group_split <- match.arg(group_split)

  # ---------- 1) 对齐样本 & 基础准备 ----------
  dat <- .prepare_data_micro_side(phyloseq_obj, metabolome_mat, target_metabolite_name,
                                  prevalence_threshold, detection_threshold, use_clr)
  X_clr  <- dat$micro_clr        # 样本×微生物（CLR后；做相关和ML）
  X_rel  <- dat$micro_rel        # 样本×微生物（相对丰度；做流行率/网络）
  X_cnt  <- dat$micro_counts     # 微生物原始计数（做DESeq2）
  y_vec  <- dat$metab_vec        # 目标代谢物（连续）
  samples <- dat$samples
  message(sprintf("Aligned samples: %d; microbes: %d; metabolite: %s",
                  length(samples), ncol(X_clr), target_metabolite_name))

  # ---------- 2) 相关分析 ----------
  cor_res <- .cor_micro_vs_metab(X = X_clr, y = y_vec, method = cor_method)

  # ---------- 3) 微生物共现网络（基于相对丰度或CLR） ----------
  net_res <- .network_microbes(X_rel)  # 用相对丰度更稳；也可改 X_clr

  # ---------- 4) 差异分析（以代谢物高/低组作为分组） ----------
  grp <- .split_by_metabolite(y_vec, strategy = group_split)
  diff_res <- .diff_micro_by_group(phyloseq_obj, samples, grp)  # 优先 DESeq2; fallback Wilcoxon

  # ---------- 5) 机器学习：用微生物（X_clr）预测代谢物（y） ----------
  ml_res <- .ml_importance_micro_to_metab(X = X_clr, y = y_vec)

  # ---------- 6) 丰度/流行率指标 ----------
  abn_res <- .abundance_prevalence(X_rel)

  # ---------- 7) 多指标集成打分（与您现有 integrated_scoring 风格一致） ----------
  scores <- .integrated_scoring_micro(
    cor_results = cor_res,
    network_metrics = net_res,
    diff_results = diff_res,
    ml_results = ml_res,
    abn_results = abn_res,
    X_rel = X_rel
  )

  top_tbl <- head(scores, top_n)

  out <- list(
    top_microbes        = top_tbl,
    all_scores          = scores,
    correlation_results = cor_res,
    network_metrics     = net_res,
    diff_results        = diff_res,
    ml_results          = ml_res,
    aligned_data        = list(X_clr = X_clr, X_rel = X_rel, X_counts = X_cnt,
                               y = y_vec, samples = samples)
  )
  class(out) <- "microbe_association_results"
  return(out)
}










#' MSFeatureSelection: Comprehensive Untargeted Metabolomics Feature Selection
#'
#' 针对非靶向代谢组（MS）数据的特征挑选流水线：缺失与批次/漂移校正、归一化、QC 稳定性过滤、
#' 差异分析（limma）、机器学习重要性（RF/GBM/Boruta）、相关网络、综合加权排名。
#'
#' _PACKAGE
#' @name MSFeatureSelection
NULL

# ========================= 1. 主函数 =========================

#' Untargeted Metabolomics Feature Selection (MS)
#'
#' @param ms_mat 矩阵/数据框：样本×代谢物（intensity），行名=样本名，列名=代谢物ID
#' @param sample_meta 数据框：包含分组/批次/QC 标识，行名=样本名
#' @param group_var 分组变量名（在 sample_meta 中）
#' @param batch_var 批次变量名（可选，用于ComBat）
#' @param qc_flag_var QC 标识变量名（TRUE/FALSE 或 0/1），用于漂移校正/稳定性过滤
#' @param prevalence_threshold 特征保留的样本占比阈值（默认0.7）
#' @param impute_method 缺失值填补：'halfmin'|'knn'（默认'halfmin'）
#' @param normalize_method 强度归一化：'TIC'|'median'|'quantile'（默认'median'）
#' @param combat 是否执行 ComBat 批次校正（默认 TRUE）
#' @param qc_drift_correction 是否执行 QC-based 漂移校正（默认 TRUE）
#' @param qc_rsd_cut 基于 QC 的 RSD 过滤阈值（默认0.3）
#' @param diff_method 差异分析方法：'limma'|'wilcoxon'（默认 'limma'）
#' @param ml_method ML 方法：'rf'|'gbm'|'boruta'（默认 'rf'）
#' @param cor_method 相关：'spearman'|'pearson'（默认 'spearman'）
#' @param cor_threshold 边阈值（默认 0.6）
#' @param p_threshold 相关 p 值阈值（默认 0.05）
#' @param weights 综合权重 list(abundance, importance, differential, network)
#' @param nfolds ML 交叉验证折数（默认 5）
#' @param ntree RF 森林数（默认 500）
#'
#' @return list: final_ranking / full_results / diff_results / importance_results /
#'   network_results / abundance_results / processed_ms（处理后矩阵）/ parameters
#' \dontrun{
#'tab = ps.ms %>%subset_samples.wt("Group",c("WT","KO")) %>% vegan_otu( ) %>%
#'  as.data.frame()
#'map = sample_data(ps.ms %>%subset_samples.wt("Group",c("WT","KO")) )
#'map$Group = as.factor(map$Group)
#'map$ID = NULL
#'res <- mass_selection(
#'  ms_mat         = tab,
#'  sample_meta    = map,
#'  group_var      = "Group",
#' batch_var      = NULL,
#'  qc_flag_var    = NULL,
#'  prevalence_threshold = 0.7,      # ≥70% 样本非缺失才保留
#'  impute_method  = "halfmin",      # 缺失填补：半最小值
#'  normalize_method = "median",     # 归一化：样本中位数
#'  combat         = F,           # 批次校正：ComBat
#' qc_drift_correction = TRUE,      # 用 QC 做 LOESS 漂移校正
#'  qc_rsd_cut     = 0.30,           # QC RSD ≤ 30% 保留
#'  diff_method    = "limma",        # 差异分析：limma（两组更稳）
#' ml_method      = "rf",           # 机器学习：随机森林
#' cor_method     = "spearman",     # 网络相关：斯皮尔曼
#'  cor_threshold  = 0.6,
#' p_threshold    = 0.05,
#' weights = list(
#'    abundance   = 0.3,   # 偏好高丰度与高流行率
#'   importance  = 0.3,   # ML 重要性
#'   differential= 0.2,   # 差异显著性与效应量
#'   network     = 0.2    # 网络中心性
#'  ),
#'  nfolds = 5,
#' ntree  = 300           # RF 森林数（示例里稍微降点速度更快）
#')
#' }
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
mass_selection <- function(ms_mat,
                           sample_meta,
                           group_var,
                           batch_var = NULL,
                           qc_flag_var = NULL,
                           prevalence_threshold = 0.7,
                           impute_method = c("halfmin","knn"),
                           normalize_method = c("median","TIC","quantile"),
                           combat = TRUE,
                           qc_drift_correction = TRUE,
                           qc_rsd_cut = 0.30,
                           diff_method = c("limma","wilcoxon"),
                           ml_method = c("rf","gbm","boruta"),
                           cor_method = c("spearman","pearson"),
                           cor_threshold = 0.6,
                           p_threshold = 0.05,
                           weights = list(abundance = 0.3,
                                          importance = 0.3,
                                          differential = 0.2,
                                          network = 0.2),
                           nfolds = 5,
                           ntree = 500) {

  impute_method   <- match.arg(impute_method)
  normalize_method<- match.arg(normalize_method)
  diff_method     <- match.arg(diff_method)
  ml_method       <- match.arg(ml_method)
  cor_method      <- match.arg(cor_method)

  stopifnot(is.matrix(ms_mat) || is.data.frame(ms_mat))
  ms_mat <- as.matrix(ms_mat)

  # ------- 样本对齐 -------
  common_samples <- intersect(rownames(ms_mat), rownames(sample_meta))
  if (length(common_samples) < 3) stop("样本行名需与元数据行名匹配，并且数量足够。")
  ms_mat <- ms_mat[common_samples, , drop = FALSE]
  sample_meta <- sample_meta[common_samples, , drop = FALSE]

  # 分组变量
  if (!group_var %in% colnames(sample_meta)) {
    stop(paste0("分组变量 ", group_var, " 不在 sample_meta 中"))
  }
  sample_meta[[group_var]] <- as.factor(sample_meta[[group_var]])

  # QC 标记/批次
  qc_vec <- NULL
  if (!is.null(qc_flag_var) && qc_flag_var %in% colnames(sample_meta)) {
    qc_vec <- as.logical(sample_meta[[qc_flag_var]])
    qc_vec[is.na(qc_vec)] <- FALSE
  }
  batch_vec <- NULL
  if (!is.null(batch_var) && batch_var %in% colnames(sample_meta)) {
    batch_vec <- as.factor(sample_meta[[batch_var]])
  }

  cat("==== MS feature selection start ====\n")
  cat(sprintf("Samples: %d, Features: %d\n", nrow(ms_mat), ncol(ms_mat)))
  cat(sprintf("Groups: %s\n", paste(levels(sample_meta[[group_var]]), collapse=", ")))
  if (!is.null(batch_vec)) cat(sprintf("Batches: %s\n", paste(levels(batch_vec), collapse=", ")))
  if (!is.null(qc_vec)) cat(sprintf("QC samples: %d\n", sum(qc_vec)))

  # ------- 1. 预处理（缺失过滤/补齐、归一化、批次/漂移） -------
  prep <- ms_preprocess(ms_mat, sample_meta,
                        prevalence_threshold = prevalence_threshold,
                        impute_method = impute_method,
                        normalize_method = normalize_method,
                        combat = combat,
                        batch_vec = batch_vec,
                        qc_flag = qc_vec,
                        qc_drift_correction = qc_drift_correction,
                        qc_rsd_cut = qc_rsd_cut)

  X <- prep$ms_proc           # 处理后的矩阵（样本×代谢物）
  sample_meta <- prep$sample_meta

  # ------- 2. 差异分析（limma / wilcoxon） -------
  cat("\n[Step 2] Differential analysis ...\n")
  diff_res <- ms_differential(X, sample_meta, group_var, method = diff_method)

  # ------- 3. 机器学习重要性 -------
  cat("\n[Step 3] Machine-learning importance ...\n")
  imp_res <- ms_ml_importance(X, sample_meta, group_var,
                              method = ml_method, nfolds = nfolds, ntree = ntree)

  # ------- 4. 相关网络 -------
  cat("\n[Step 4] Correlation network ...\n")
  net_res <- ms_network(X, cor_method = cor_method,
                        cor_threshold = cor_threshold, p_threshold = p_threshold)

  # ------- 5. 丰度/流行率（用于偏向高丰度） -------
  cat("\n[Step 5] Abundance metrics ...\n")
  abd_res <- ms_abundance_metrics(X)

  # ------- 6. 综合排名 -------
  cat("\n[Step 6] Integrated weighted ranking ...\n")
  rank_res <- ms_weighted_ranking(diff_res, imp_res, net_res, abd_res, weights)

  cat("\nDone.\n")
  cat(sprintf("Significant features (p<0.05): %d\n",
              sum(diff_res$pvalue < 0.05, na.rm = TRUE)))
  cat(sprintf("Network edges: %d\n", sum(net_res$Degree > 0)))

  res <- list(
    final_ranking = rank_res$summary_results,
    full_results  = rank_res$full_results,
    diff_results  = diff_res,
    importance_results = imp_res,
    network_results = net_res,
    abundance_results = abd_res,
    processed_ms  = X,
    parameters = list(
      group_var = group_var,
      prevalence_threshold = prevalence_threshold,
      impute_method = impute_method,
      normalize_method = normalize_method,
      combat = combat,
      qc_drift_correction = qc_drift_correction,
      qc_rsd_cut = qc_rsd_cut,
      diff_method = diff_method,
      ml_method = ml_method,
      cor_method = cor_method,
      cor_threshold = cor_threshold,
      p_threshold = p_threshold,
      weights = weights
    )
  )
  class(res) <- c("ms_feature_selection","list")
  return(res)
}




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
feature_selection <- function(results,
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

