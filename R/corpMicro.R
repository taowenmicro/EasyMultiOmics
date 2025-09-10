
#
# # 加载必要的包
# library(phyloseq)
# library(randomForest)
# library(caret)
# library(dplyr)
# library(ggplot2)
# library(pheatmap)
# library(corrplot)
# library(parallel)
# library(doParallel)

# 主要分析类
MicrobialInteractionAnalyzer <- function(ps_object) {

  # 验证phyloseq对象
  if (!inherits(ps_object, "phyloseq")) {
    stop("输入必须是phyloseq对象")
  }

  # 提取数据
  otu_table <- as.data.frame(otu_table(ps_object))
  sample_data <- as.data.frame(sample_data(ps_object))

  # 创建分析对象
  analyzer <- list(
    otu_matrix = otu_table,
    sample_groups = sample_data,
    otu_names = rownames(otu_table),
    interaction_matrix = NULL,
    interaction_names = NULL,
    model_results = NULL,
    feature_importance = NULL
  )

  class(analyzer) <- "MicrobialInteractionAnalyzer"
  return(analyzer)
}

# 创建OTU互作矩阵
create_interaction_matrix <- function(analyzer, group_column,
                                      interaction_type = "sum",
                                      min_abundance = 0.001,
                                      max_features = 1000) {

  cat("正在创建OTU互作矩阵...\n")

  # 转置矩阵 (样本为行，OTU为列)
  otu_data <- t(analyzer$otu_matrix)

  # 过滤低丰度OTU
  abundance_filter <- colMeans(otu_data) >= min_abundance
  otu_filtered <- otu_data[, abundance_filter]
  filtered_otu_names <- analyzer$otu_names[abundance_filter]

  cat(sprintf("过滤后保留 %d 个OTU (原始: %d)\n",
              ncol(otu_filtered), length(analyzer$otu_names)))

  # 限制特征数量以避免计算爆炸
  n_otus <- ncol(otu_filtered)
  n_combinations <- choose(n_otus, 2)

  if (n_combinations > max_features) {
    # 根据丰度和方差选择top OTUs
    abundance_rank <- rank(-colMeans(otu_filtered))
    variance_rank <- rank(-apply(otu_filtered, 2, var))
    combined_rank <- abundance_rank + variance_rank

    top_indices <- order(combined_rank)[1:floor(sqrt(2 * max_features))]
    otu_filtered <- otu_filtered[, top_indices]
    filtered_otu_names <- filtered_otu_names[top_indices]

    cat(sprintf("为控制计算复杂度，选择了前 %d 个OTU\n", ncol(otu_filtered)))
  }

  # 生成所有OTU两两组合
  n_otus_final <- ncol(otu_filtered)
  otu_combinations <- combn(n_otus_final, 2)
  n_interactions <- ncol(otu_combinations)

  cat(sprintf("生成 %d 个互作特征\n", n_interactions))

  # 创建互作矩阵
  interaction_matrix <- matrix(0, nrow = nrow(otu_filtered), ncol = n_interactions)
  interaction_names <- character(n_interactions)

  # 并行计算互作特征
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)

  interaction_data <- foreach(i = 1:n_interactions, .combine = cbind) %dopar% {
    otu1_idx <- otu_combinations[1, i]
    otu2_idx <- otu_combinations[2, i]

    otu1_abundance <- otu_filtered[, otu1_idx]
    otu2_abundance <- otu_filtered[, otu2_idx]

    if (interaction_type == "sum") {
      interaction_value <- otu1_abundance + otu2_abundance
    } else if (interaction_type == "product") {
      interaction_value <- otu1_abundance * otu2_abundance
    } else if (interaction_type == "ratio") {
      interaction_value <- log2((otu1_abundance + 1e-6) / (otu2_abundance + 1e-6))
    }

    interaction_value
  }

  stopCluster(cl)

  # 生成互作名称
  for (i in 1:n_interactions) {
    otu1_idx <- otu_combinations[1, i]
    otu2_idx <- otu_combinations[2, i]
    interaction_names[i] <- paste(filtered_otu_names[otu1_idx],
                                  filtered_otu_names[otu2_idx],
                                  sep = "_x_")
  }

  # 保存结果
  analyzer$interaction_matrix <- as.data.frame(interaction_data)
  colnames(analyzer$interaction_matrix) <- interaction_names
  analyzer$interaction_names <- interaction_names
  analyzer$filtered_otu_names <- filtered_otu_names
  analyzer$group_column <- group_column

  cat("互作矩阵创建完成！\n")
  return(analyzer)
}

# 机器学习分类和特征选择
perform_ml_analysis <- function(analyzer, test_size = 0.3,
                                feature_selection_method = "rf_importance",
                                top_n_features = 50) {

  cat("开始机器学习分析...\n")

  # 准备数据
  X <- analyzer$interaction_matrix
  y <- factor(analyzer$sample_groups[[analyzer$group_column]])

  # 检查分组
  if (length(levels(y)) != 2) {
    stop("目前只支持二分类问题")
  }

  cat(sprintf("样本数: %d, 特征数: %d, 分组: %s\n",
              nrow(X), ncol(X), paste(table(y), collapse = " vs ")))

  # 数据标准化
  X_scaled <- scale(X)

  # 特征选择
  if (feature_selection_method == "rf_importance") {
    # 使用随机森林重要性选择特征
    rf_temp <- randomForest(X_scaled, y, importance = TRUE, ntree = 100)
    importance_scores <- importance(rf_temp)[, "MeanDecreaseGini"]
    top_features_idx <- order(importance_scores, decreasing = TRUE)[1:min(top_n_features, length(importance_scores))]

  } else if (feature_selection_method == "univariate") {
    # 使用单变量统计检验
    p_values <- apply(X_scaled, 2, function(feature) {
      t.test(feature ~ y)$p.value
    })
    top_features_idx <- order(p_values)[1:min(top_n_features, length(p_values))]
  }

  X_selected <- X_scaled[, top_features_idx]
  selected_feature_names <- colnames(X)[top_features_idx]

  cat(sprintf("选择了前 %d 个重要特征\n", ncol(X_selected)))

  # 划分训练测试集
  set.seed(42)
  train_indices <- createDataPartition(y, p = 1 - test_size, list = FALSE)

  X_train <- X_selected[train_indices, ]
  X_test <- X_selected[-train_indices, ]
  y_train <- y[train_indices]
  y_test <- y[-train_indices]

  # 训练随机森林模型
  rf_model <- randomForest(X_train, y_train,
                           ntree = 1000,
                           importance = TRUE,
                           do.trace = 50)

  # 预测和评估
  y_pred <- predict(rf_model, X_test)
  accuracy <- mean(y_pred == y_test)

  # 交叉验证
  cv_scores <- cross_validate_model(X_selected, y, n_folds = 5)

  # 特征重要性
  feature_importance <- importance(rf_model)[, "MeanDecreaseGini"]
  names(feature_importance) <- selected_feature_names
  feature_importance <- sort(feature_importance, decreasing = TRUE)

  # 保存结果
  analyzer$model_results <- list(
    model = rf_model,
    accuracy = accuracy,
    cv_scores = cv_scores,
    predictions = y_pred,
    actual = y_test,
    selected_features = selected_feature_names
  )

  analyzer$feature_importance <- feature_importance

  cat(sprintf("模型准确率: %.3f\n", accuracy))
  cat(sprintf("交叉验证平均准确率: %.3f ± %.3f\n",
              mean(cv_scores), sd(cv_scores)))

  return(analyzer)
}

# 交叉验证函数
cross_validate_model <- function(X, y, n_folds = 5) {
  folds <- createFolds(y, k = n_folds)
  cv_scores <- numeric(n_folds)

  for (i in 1:n_folds) {
    train_idx <- unlist(folds[-i])
    test_idx <- folds[[i]]

    X_train_cv <- X[train_idx, ]
    X_test_cv <- X[test_idx, ]
    y_train_cv <- y[train_idx]
    y_test_cv <- y[test_idx]

    rf_cv <- randomForest(X_train_cv, y_train_cv, ntree = 200)
    y_pred_cv <- predict(rf_cv, X_test_cv)
    cv_scores[i] <- mean(y_pred_cv == y_test_cv)
  }

  return(cv_scores)
}

# 识别关键互作菌株对
identify_key_interactions <- function(analyzer, top_n = 20,
                                      significance_threshold = 0.01) {

  cat("识别关键互作菌株对...\n")

  if (is.null(analyzer$feature_importance)) {
    stop("请先运行机器学习分析")
  }

  # 获取top特征
  top_features <- head(analyzer$feature_importance, top_n)

  # 解析互作菌株对
  interaction_pairs <- data.frame(
    interaction = names(top_features),
    importance = as.numeric(top_features),
    stringsAsFactors = FALSE
  )

  # 分离菌株对名称
  interaction_pairs$otu1 <- sapply(strsplit(interaction_pairs$interaction, "_x_"), `[`, 1)
  interaction_pairs$otu2 <- sapply(strsplit(interaction_pairs$interaction, "_x_"), `[`, 2)

  # 计算统计显著性
  X <- analyzer$interaction_matrix[, names(top_features)]
  y <- factor(analyzer$sample_groups[[analyzer$group_column]])

  interaction_pairs$p_value <- apply(X, 2, function(feature) {
    t.test(feature ~ y)$p.value
  })

  # 计算效应大小 (Cohen's d)
  interaction_pairs$effect_size <- apply(X, 2, function(feature) {
    group1 <- feature[y == levels(y)[1]]
    group2 <- feature[y == levels(y)[2]]

    pooled_sd <- sqrt(((length(group1) - 1) * var(group1) +
                         (length(group2) - 1) * var(group2)) /
                        (length(group1) + length(group2) - 2))

    (mean(group1) - mean(group2)) / pooled_sd
  })

  # 过滤显著性特征
  significant_interactions <- interaction_pairs[interaction_pairs$p_value < significance_threshold, ]

  analyzer$top_interactions <- significant_interactions

  cat(sprintf("发现 %d 个显著的关键互作菌株对 (p < %.3f)\n",
              nrow(significant_interactions), significance_threshold))

  return(analyzer)
}

# 可视化结果

plot_feature_importance_enhanced <- function(analyzer) {

  if (is.null(analyzer$feature_importance)) return()

  top_20 <- head(analyzer$feature_importance, 20)

  # 解析互作对名称，简化显示
  interaction_df <- data.frame(
    interaction = names(top_20),
    importance = as.numeric(top_20),
    stringsAsFactors = FALSE
  )

  # 简化OTU名称用于显示
  interaction_df$display_name <- sapply(interaction_df$interaction, function(x) {
    parts <- strsplit(x, "_x_")[[1]]
    otu1 <- gsub("OTU_", "", parts[1])
    otu2 <- gsub("OTU_", "", parts[2])
    paste0("OTU", otu1, " × OTU", otu2)
  })

  # 按重要性排序
  interaction_df <- interaction_df[order(interaction_df$importance), ]
  interaction_df$display_name <- factor(interaction_df$interaction,
                                        levels = interaction_df$interaction)

  # 创建水平条形图
  p1 <- ggplot(interaction_df, aes(x = importance, y = display_name)) +
    geom_col(aes(fill = importance), alpha = 0.8, width = 0.7) +
    scale_fill_viridis_c(name = "重要性\n得分", option = "plasma") +
    labs(
      title = "🔬 Top 20 关键菌株互作特征",
      subtitle = "基于随机森林特征重要性排序",
      x = "特征重要性得分",
      y = "菌株互作对"
    )
  interactions <- analyzer$top_interactions
  interaction_stats <- data.frame(
    Category = c("总互作对", "显著互作", "强效应互作"),
    Count = c(
      nrow(interactions),
      sum(interactions$p_value < 0.01),
      sum(abs(interactions$effect_size) > 0.5)
    ),
    Color = c("#9b59b6", "#e74c3c", "#f39c12")
  )

  return(p1)
}


# res = plot_model_performance_enhanced (analyzer)
# res[[1]]
plot_model_performance_enhanced <- function(analyzer) {

  if (is.null(analyzer$model_results)) return()

  results <- analyzer$model_results

  # 2.1 美化的混淆矩阵
  cm <- table(Predicted = results$predictions, Actual = results$actual)
  cm_df <- as.data.frame(cm)

  # 计算百分比
  cm_df$Percentage <- round(cm_df$Freq / sum(cm_df$Freq) * 100, 1)

  p_cm <- ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = paste0(Freq, "\n(", Percentage, "%)")),
              size = 6, fontface = "bold", color = "white") +
    scale_fill_gradient(low = "#3498db", high = "#e74c3c",
                        name = "样本数") +
    labs(
      title = "🎯 分类结果混淆矩阵",
      subtitle = sprintf("总体准确率: %.1f%%", results$accuracy * 100),
      x = "实际分组",
      y = "预测分组"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray60"),
      axis.text = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    coord_equal()


  # 2.2 交叉验证性能图
  cv_data <- data.frame(
    Fold = 1:length(results$cv_scores),
    Accuracy = results$cv_scores,
    Mean_Accuracy = mean(results$cv_scores)
  )

  p_cv <- ggplot(cv_data, aes(x = Fold, y = Accuracy)) +
    geom_hline(yintercept = cv_data$Mean_Accuracy[1],
               linetype = "dashed", color = "#e74c3c", size = 1) +
    geom_line(color = "#3498db", size = 1.2) +
    geom_point(color = "#2980b9", size = 4, alpha = 0.8) +
    geom_point(color = "white", size = 2) +
    annotate("text", x = max(cv_data$Fold) * 0.7, y = cv_data$Mean_Accuracy[1] + 0.02,
             label = sprintf("平均: %.3f ± %.3f",
                             mean(results$cv_scores),
                             sd(results$cv_scores)),
             color = "#e74c3c", fontface = "bold") +
    scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
    labs(
      title = "📈 5折交叉验证性能",
      subtitle = "模型在不同数据子集上的稳定性评估",
      x = "交叉验证折数",
      y = "分类准确率"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray60"),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    )

  return(list(  p_cm,p_cv))
}


plot_interaction_heatmap <- function(analyzer) {

  interactions <- analyzer$top_interactions

  # 创建节点和边的数据框
  nodes <- unique(c(interactions$otu1, interactions$otu2))
  edges <- data.frame(
    from = interactions$otu1,
    to = interactions$otu2,
    weight = abs(interactions$effect_size),
    importance = interactions$importance,
    p_value = interactions$p_value
  )

  # 创建邻接矩阵
  adj_matrix <- matrix(0, nrow = length(nodes), ncol = length(nodes))
  rownames(adj_matrix) <- colnames(adj_matrix) <- nodes

  for (i in 1:nrow(edges)) {
    adj_matrix[edges$from[i], edges$to[i]] <- edges$weight[i]
    adj_matrix[edges$to[i], edges$from[i]] <- edges$weight[i]
  }

  # 将邻接矩阵转换为长格式以用于ggplot
  adj_long <- melt(adj_matrix, varnames = c("from", "to"), value.name = "weight")

  # 使用ggplot绘制热图
  p_network <- ggplot(adj_long, aes(x = to, y = from, fill = weight)) +
    geom_tile(color = "white", linewidth = 0.5) +
    scale_fill_gradient2(
      low = "#ECEFF4",
      mid = "#88CCEE",
      high = "#4C78A8",
      midpoint = median(adj_long$weight, na.rm = TRUE),
      name = "Effect"
    ) +
    labs(title = "Key Strain Interaction Network (Effect Size)") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 8),
      axis.title = element_text(size = 12),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black")
    )

  # 直接显示图表
  print(p_network)
}



plot_results <- function(analyzer,n = 20) {

  # 1. 特征重要性图
  if (!is.null(analyzer$feature_importance)) {
    top_20 <- head(analyzer$feature_importance, n)
    importance_data <- data.frame(
      Feature = names(top_20),
      Importance = as.numeric(top_20)
    )

    p_importance <- ggplot(importance_data, aes(x = reorder(Feature, Importance), y = Importance)) +
      geom_bar(stat = "identity", fill = "#4C78A8") +
      coord_flip() +
      labs(title = "Top 20 互作特征重要性", x = "", y = "重要性得分") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major.x = element_line(color = "grey80"),
        panel.grid.minor = element_blank()
      )

    print(p_importance)
  }

  # 2. 模型性能评估
  if (!is.null(analyzer$model_results)) {
    results <- analyzer$model_results

    # 混淆矩阵
    cm <- table(Predicted = results$predictions, Actual = results$actual)
    cm_data <- as.data.frame(cm)

    p_cm <- ggplot(cm_data, aes(x = Actual, y = Predicted, fill = Freq)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = Freq), color = "black", size = 5) +
      scale_fill_gradient(low = "#ECEFF4", high = "#4C78A8") +
      labs(title = "混淆矩阵") +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid = element_blank()
      )

    # 交叉验证结果
    cv_data <- data.frame(Fold = 1:length(results$cv_scores),
                          Accuracy = results$cv_scores)

    p_cv <- ggplot(cv_data, aes(x = Fold, y = Accuracy)) +
      geom_line(color = "#4C78A8", size = 1) +
      geom_point(color = "#4C78A8", size = 3) +
      geom_hline(yintercept = mean(results$cv_scores),
                 linetype = "dashed", color = "#F28E2B") +
      ylim(0, 1) +
      labs(
        title = "交叉验证准确率",
        subtitle = sprintf("平均准确率: %.3f ± %.3f",
                           mean(results$cv_scores),
                           sd(results$cv_scores))
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        panel.grid.major = element_line(color = "grey80"),
        panel.grid.minor = element_blank()
      )

    # 使用aplot合并两个图
    combined_plot <- aplot::plot_list(p_cm, p_cv, ncol = 2)
    print(combined_plot)
  }
  return(list(p_importance,combined_plot))

}


# 方法验证和比较
validate_method <- function(analyzer, known_interactions = NULL) {

  cat("\n=== 方法可行性评估 ===\n")

  # 1. 数据质量检查
  cat("1. 数据质量评估:\n")
  cat(sprintf("   - 样本数量: %d\n", nrow(analyzer$interaction_matrix)))
  cat(sprintf("   - 特征数量: %d\n", ncol(analyzer$interaction_matrix)))
  cat(sprintf("   - 分组平衡性: %s\n",
              paste(table(analyzer$sample_groups[[analyzer$group_column]]), collapse = " vs ")))

  # 2. 模型性能评估
  if (!is.null(analyzer$model_results)) {
    cat("\n2. 模型性能:\n")
    results <- analyzer$model_results
    cat(sprintf("   - 测试准确率: %.3f\n", results$accuracy))
    cat(sprintf("   - 交叉验证准确率: %.3f ± %.3f\n",
                mean(results$cv_scores), sd(results$cv_scores)))

    # 判断模型是否过拟合
    train_accuracy <- results$model$confusion[1,1] + results$model$confusion[2,2]
    train_accuracy <- train_accuracy / sum(results$model$confusion)

    if (train_accuracy - results$accuracy > 0.1) {
      cat("   - 警告: 可能存在过拟合\n")
    }
  }

  # 3. 特征稳定性检查
  if (!is.null(analyzer$feature_importance)) {
    cat("\n3. 特征选择稳定性:\n")
    top_features <- head(analyzer$feature_importance, 10)
    cat(sprintf("   - Top 10特征重要性范围: %.3f - %.3f\n",
                min(top_features), max(top_features)))

    # 检查特征分布
    feature_cv <- apply(analyzer$interaction_matrix[, names(top_features)], 2,
                        function(x) sd(x) / mean(x))
    cat(sprintf("   - 关键特征变异系数范围: %.3f - %.3f\n",
                min(feature_cv), max(feature_cv)))
  }

  # 4. 生物学意义评估
  cat("\n4. 生物学合理性:\n")
  if (!is.null(analyzer$top_interactions)) {
    sig_interactions <- analyzer$top_interactions
    cat(sprintf("   - 显著互作数量: %d\n", nrow(sig_interactions)))
    cat(sprintf("   - 平均效应大小: %.3f\n", mean(abs(sig_interactions$effect_size))))

    # 检查是否有已知互作
    if (!is.null(known_interactions)) {
      validate_known_interactions(analyzer, known_interactions)
    }
  }

  return(invisible(analyzer))
}

# 验证已知互作
validate_known_interactions <- function(analyzer, known_interactions) {

  predicted_pairs <- paste(analyzer$top_interactions$otu1,
                           analyzer$top_interactions$otu2, sep = "_")
  known_pairs <- paste(known_interactions$otu1, known_interactions$otu2, sep = "_")

  overlap <- intersect(predicted_pairs, known_pairs)

  cat(sprintf("   - 与已知互作重叠: %d/%d (%.1f%%)\n",
              length(overlap), length(known_pairs),
              100 * length(overlap) / length(known_pairs)))
}

# 与现有方法比较
compare_with_existing_methods <- function(analyzer) {

  cat("\n=== 与现有方法比较 ===\n")

  # 1. 与相关性分析比较
  cat("1. vs 传统相关性分析:\n")
  otu_data <- t(analyzer$otu_matrix)

  # 计算Spearman相关性
  cor_matrix <- cor(otu_data, method = "spearman")
  cor_matrix[upper.tri(cor_matrix, diag = TRUE)] <- NA
  significant_cors <- which(abs(cor_matrix) > 0.5 & !is.na(cor_matrix), arr.ind = TRUE)

  cat(sprintf("   - 传统方法发现显著相关: %d 对\n", nrow(significant_cors)))
  cat(sprintf("   - 本方法发现关键互作: %d 对\n",
              ifelse(is.null(analyzer$top_interactions), 0, nrow(analyzer$top_interactions))))

  # 2. 计算复杂度比较
  cat("\n2. 计算复杂度:\n")
  n_otus <- ncol(otu_data)
  cat(sprintf("   - 传统相关性分析: O(n²) = %d 计算\n", n_otus^2))
  cat(sprintf("   - 本方法特征数: %d\n",
              ifelse(is.null(analyzer$interaction_matrix), 0, ncol(analyzer$interaction_matrix))))

  # 3. 预测能力比较
  cat("\n3. 预测能力优势:\n")
  cat("   - 传统方法: 描述性，无预测能力\n")
  cat("   - 本方法: 具有分类预测能力，可用于新样本\n")

  return(invisible(analyzer))
}

# 主分析流程
run_interaction_analysis <- function(ps_object, group_column,
                                     interaction_type = "sum",
                                     min_abundance = 0.001,
                                     max_features = 1000,
                                     top_n_features = 50) {

  cat("开始菌株互作分析流程...\n")

  # 创建分析器
  analyzer <- MicrobialInteractionAnalyzer(ps_object)

  # 添加分组信息
  sample_info <- as.data.frame(sample_data(ps_object))
  analyzer$sample_groups <- sample_info

  # 创建互作矩阵
  analyzer <- create_interaction_matrix(analyzer, group_column,
                                        interaction_type, min_abundance, max_features)

  # 机器学习分析
  analyzer <- perform_ml_analysis(analyzer, top_n_features = top_n_features)

  # 识别关键互作
  analyzer <- identify_key_interactions(analyzer)

  # 方法验证
  analyzer <- validate_method(analyzer)

  # 与现有方法比较
  analyzer <- compare_with_existing_methods(analyzer)

  # 生成报告
  generate_analysis_report(analyzer)

  return(analyzer)
}

# 生成分析报告
generate_analysis_report <- function(analyzer) {

  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("菌株互作分析报告\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")


  if (!is.null(analyzer$model_results)) {
    cat(sprintf("\n【模型性能】准确率: %.3f\n", analyzer$model_results$accuracy))

    if (analyzer$model_results$accuracy > 0.8) {
      cat("→ 模型表现优秀，方法可行性高\n")
    } else if (analyzer$model_results$accuracy > 0.7) {
      cat("→ 模型表现良好，方法具有一定可行性\n")
    } else {
      cat("→ 模型表现一般，需要优化参数或增加样本\n")
    }
  }

  cat("\n【建议】\n")
  cat("1. 可以作为探索性分析工具使用\n")
  cat("2. 建议结合传统网络分析方法验证\n")
  cat("3. 关键互作对需要实验验证\n")
  cat("4. 可扩展到多分类和回归问题\n")

  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")
}


# ps_object = pst
# group_column = "Group"
# interaction_type = "sum"
# min_abundance = 0.0000001
# max_features = 20000  # 控制计算量
# top_n_features = 100000


# 主分析流程
run_interaction_analysis2 <- function(ps_object, group_column,
                                      interaction_type = "sum",
                                      min_abundance = 0.001,
                                      max_features = 1000,
                                      top_n_features = 50) {

  cat("开始菌株互作分析流程...\n")

  # 创建分析器
  analyzer <- MicrobialInteractionAnalyzer(ps_object)

  # 添加分组信息
  sample_info <- as.data.frame(sample_data(ps_object))
  analyzer$sample_groups <- sample_info

  # 创建互作矩阵
  analyzer <- create_interaction_matrix(analyzer, group_column,
                                        interaction_type, min_abundance, max_features)

  # 机器学习分析
  analyzer <- perform_ml_analysis(analyzer, top_n_features = top_n_features)

  # 识别关键互作
  analyzer <- identify_key_interactions(analyzer,top_n = top_n_features,significance_threshold = 0.01)

  # 方法验证
  analyzer <- validate_method(analyzer)

  # 与现有方法比较
  analyzer <- compare_with_existing_methods(analyzer)

  # 生成报告
  generate_analysis_report(analyzer)

  return(analyzer)
}

# 改进版互作函数 - 考虑组成性数据特点
create_advanced_interaction_matrix <- function(analyzer, group_column,
                                               interaction_types = c("sum", "ratio", "clr_product"),
                                               min_abundance = 0.001) {

  library(compositions)  # 用于CLR变换

  otu_data <- t(analyzer$otu_matrix)

  # CLR变换处理组成性数据
  otu_clr <- clr(otu_data + 1e-6)  # 添加伪计数避免0值

  # 过滤低丰度OTU
  abundance_filter <- colMeans(otu_data) >= min_abundance
  otu_filtered <- otu_data[, abundance_filter]
  otu_clr_filtered <- otu_clr[, abundance_filter]

  # 多种互作类型
  interaction_features <- list()

  for (type in interaction_types) {
    cat(sprintf("计算 %s 类型互作...\n", type))

    n_otus <- ncol(otu_filtered)
    combinations_matrix <- combn(n_otus, 2)

    type_features <- matrix(0, nrow = nrow(otu_filtered),
                            ncol = ncol(combinations_matrix))

    for (i in 1:ncol(combinations_matrix)) {
      idx1 <- combinations_matrix[1, i]
      idx2 <- combinations_matrix[2, i]

      if (type == "sum") {
        type_features[, i] <- otu_filtered[, idx1] + otu_filtered[, idx2]
      } else if (type == "ratio") {
        type_features[, i] <- log2((otu_filtered[, idx1] + 1e-6) /
                                     (otu_filtered[, idx2] + 1e-6))
      } else if (type == "clr_product") {
        type_features[, i] <- otu_clr_filtered[, idx1] * otu_clr_filtered[, idx2]
      }
    }

    interaction_features[[type]] <- type_features
  }

  # 合并所有互作类型
  all_interactions <- do.call(cbind, interaction_features)

  analyzer$interaction_matrix <- as.data.frame(all_interactions)
  analyzer$group_column <- group_column

  return(analyzer)
}

# 方法比较框架
compare_methods <- function(ps_object, group_column) {

  cat("\n=== 方法比较分析 ===\n")

  # 1. 传统差异分析 (DESeq2-like)
  cat("1. 运行传统差异分析...\n")
  deseq_results <- run_differential_analysis2(ps_object, group_column)

  # 2. 相关性网络分析
  cat("2. 运行相关性网络分析...\n")
  correlation_results <- run_correlation_analysis2(ps_object, group_column)

  # 3. 本方法
  cat("3. 运行互作特征学习...\n")
  interaction_results <- run_interaction_analysis2(ps_object, group_column)

  # 比较结果
  comparison <- list(
    differential = deseq_results,
    correlation = correlation_results,
    interaction = interaction_results
  )

  return(comparison)
}
# 定义函数：对 phyloseq 对象进行 CLR 变换
clr_transform_phyloseq <- function(ps_object) {
  # 验证输入是否为 phyloseq 对象
  if (!inherits(ps_object, "phyloseq")) {
    stop("输入必须是 phyloseq 对象")
  }

  # 提取 OTU 表（丰度矩阵）
  otu_data <- as(otu_table(ps_object), "matrix")

  # 检查 OTU 表方向（样本为行还是列）
  if (taxa_are_rows(ps_object)) {
    otu_data <- t(otu_data) # 转换为样本为行，OTU 为列
  }

  # 添加伪计数（pseudo-count）以避免零值问题
  pseudo_count <- 1e-6
  otu_data <- otu_data + pseudo_count

  # 应用 CLR 变换（使用 compositions 包）
  otu_clr <- compositions::clr(otu_data)

  # 将 CLR 变换后的数据转换回 OTU 表
  # 确保方向与原始 OTU 表一致
  if (taxa_are_rows(ps_object)) {
    otu_clr <- t(otu_clr) # 转回 OTU 为行，样本为列
  }
  otu_table_clr <- otu_table(otu_clr, taxa_are_rows = taxa_are_rows(ps_object))

  # 创建新的 phyloseq 对象，保留原始的样本元数据和分类信息
  ps_clr <- phyloseq(
    otu_table = otu_table_clr,
    sample_data = sample_data(ps_object), # 保留样本元数据
    tax_table = tax_table(ps_object)      # 保留分类信息（如物种注释）
  )

  # 如果原始对象有 phy_tree（系统发育树），也保留
  if (!is.null(phy_tree(ps_object, errorIfNULL = FALSE))) {
    ps_clr <- merge_phyloseq(ps_clr, phy_tree(ps_object))
  }

  cat("CLR 变换完成！新 phyloseq 对象已创建。\n")
  return(ps_clr)
}
