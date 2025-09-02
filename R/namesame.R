
# library(EasyMultiOmics)
# ?count_to_mass
# ?find_associated_metabolites
# count_to_mass = find_associated_metabolites
# mass_to_count =find_associated_microbes
# count_selection = microbiome_feature_selection
# ms_feature_selection = mass_selection
# multi_omics_alignment #


#' MicrobiomeMetabolomeAnalysis: Cross-omics Association Analysis
#'
#' A comprehensive toolkit for identifying metabolites associated with specific microbes
#' through integrated multi-dimensional analysis including correlation, network topology,
#' differential expression, and machine learning approaches.
#'
#' _PACKAGE
#' @name MicrobiomeMetabolomeAnalysis
#' @import phyloseq
#' @import tidyverse
#' @import vegan
#' @import randomForest
#' @import glmnet
#' @import igraph
#' @import limma
#' @import mixOmics
#' @importFrom stats cor cor.test p.adjust median sd var lm predict prcomp dist
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






#' MicrobiomeMetabolomeAnalysis: Cross-omics Association Analysis
#'
#' A comprehensive toolkit for identifying metabolites associated with specific microbes
#' through integrated multi-dimensional analysis including correlation, network topology,
#' differential expression, and machine learning approaches.
#'
#' _PACKAGE
#' @name MicrobiomeMetabolomeAnalysis
#' @import phyloseq
#' @import tidyverse
#' @import vegan
#' @import randomForest
#' @import glmnet
#' @import igraph
#' @import limma
#' @import mixOmics
#' @importFrom stats cor cor.test p.adjust median sd var lm predict prcomp dist
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



#' MicrobiomeMetabolomeAnalysis: Cross-omics Association Analysis
#'
#' A comprehensive toolkit for identifying metabolites associated with specific microbes
#' through integrated multi-dimensional analysis including correlation, network topology,
#' differential expression, and machine learning approaches.
#'
#' _PACKAGE
#' @name MicrobiomeMetabolomeAnalysis
#' @import phyloseq
#' @import tidyverse
#' @import vegan
#' @import randomForest
#' @import glmnet
#' @import igraph
#' @import limma
#' @import mixOmics
#' @importFrom stats cor cor.test p.adjust median sd var lm predict prcomp dist
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



#' MSFeatureSelection: Comprehensive Untargeted Metabolomics Feature Selection
#'
#' 针对非靶向代谢组（MS）数据的特征挑选流水线：缺失与批次/漂移校正、归一化、QC 稳定性过滤、
#' 差异分析（limma）、机器学习重要性（RF/GBM/Boruta）、相关网络、综合加权排名。
#'
#' _PACKAGE
#' @name MSFeatureSelection
#' @importFrom stats median sd var prcomp quantile model.matrix loess predict wilcox.test
#' @importFrom utils write.csv
#' @importFrom limma lmFit eBayes topTable makeContrasts contrasts.fit
#' @importFrom sva ComBat
#' @importFrom Hmisc rcorr
#' @importFrom caret train trainControl varImp
#' @importFrom randomForest randomForest
#' @importFrom igraph graph_from_adjacency_matrix betweenness closeness eigen_centrality hub_score cluster_fast_greedy membership vcount
#' @import ggplot2
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

