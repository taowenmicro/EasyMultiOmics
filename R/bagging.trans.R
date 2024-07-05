
bagging.trans <- function(ps=ps, top = 20, seed = 1010, k = 5) {
  set.seed(seed)

  # 数据准备
  ps.cs <- ps %>% filter_OTU_ps(top)
  map <- as.data.frame(phyloseq::sample_data(ps.cs))
  otutab <- as.data.frame(t(ggClusterNet::vegan_otu(ps.cs)))
  colnames(otutab) <- gsub("-", "_", colnames(otutab))
  test <- as.data.frame(t(otutab))
  test$OTUgroup <- factor(map$Group)

  # 初始化结果存储
  accuracy_values <- numeric(k)
  importance_list <- list()

  # 创建交叉验证的折
  folds <- createFolds(y = test$OTUgroup, k = k)

  # 进行交叉验证
  for (i in 1:k) {
    fold_test <- test[folds[[i]], ]
    fold_train <- test[-folds[[i]], ]

    # 训练 Bagging 模型
    a_bagging <- bagging(OTUgroup ~ ., data = fold_train, coob = TRUE)

    # 得到测试集的预测值
    pred <- predict(a_bagging, newdata = fold_test, type = 'class')

    # 计算准确率
    correct_predictions <- sum(pred == fold_test$OTUgroup)
    accuracy <- correct_predictions / nrow(fold_test)
    accuracy_values[i] <- accuracy

    # 提取特征重要性
    tree_importances <- lapply(a_bagging$mtrees, function(tree) {
      tree$btree$variable.importance
    })
    combined_importance <- do.call(cbind, tree_importances)
    avg_importance <- rowMeans(combined_importance, na.rm = TRUE)
    importance_list[[i]] <- avg_importance
  }
  #
  # tree$btree$variable.importance
  #
  # a_bagging$mtrees[[1]]$btree$variable.importance

  # 平均准确率
  mean_accuracy <- mean(accuracy_values)
  accuracy_result <- paste("Bagging Accuracy:", round(mean_accuracy, 3))
  print(accuracy_result)

  # 合并重要变量
  if (length(importance_list) > 0) {
    combined_importance <- do.call(cbind, importance_list)
    avg_importance <- rowMeans(combined_importance, na.rm = TRUE)
    importance_df <- data.frame(Feature = names(avg_importance), Importance = avg_importance)
    importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
  } else {
    importance_df <- data.frame(Feature = character(0), Importance = numeric(0))
  }
  importance_df
  # 返回结果
  list(Accuracy = accuracy_result, Importance = importance_df)
}

