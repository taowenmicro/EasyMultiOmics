#' @title Naive Bayes Classification of characteristic microorganisms
#' @description
#' naiveBayes, one of the machine learning methods, was used to screen for characteristic
#' microorganisms, and the model was evaluated using k-fold cross-validation.
#' @param ps A phyloseq format file used as an alternative for the input containing otu, tax, and map.
#' @param top The top microorganisms to consider.
#' @param seed The random seed for reproducibility.
#' @param k The number of folds for cross-validation.
#' @return A list object including the following components:
#' \item{Accuracy}{The average accuracy of the naiveBayes model.}
#' \item{Importance}{A data frame showing the feature importance ranked in descending order.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' library(ggClusterNet)
#' library(caret)
#' res = naiveBayes.micro(ps=ps, top = 20, seed = 1010, k = 5)
#' accuracy = res[[1]]
#' accuracy
#' importance = res[[2]]
#' importance
 naivebayes.micro <- function(ps= ps, top = 20, seed = 1010, k = 5) {
   set.seed(seed)
  # 数据准备
  ps.cs <- ps %>% filter_OTU_ps(top)
  map <- as.data.frame(phyloseq::sample_data(ps.cs))
  otutab <- as.data.frame(t(ggClusterNet::vegan_otu(ps.cs)))
  colnames(otutab) <- gsub("-", "_", colnames(otutab))
  test <- as.data.frame(t(otutab))
  test$OTUgroup <- factor(map$Group)
  # 确保因变量是具有两个水平的因子
  # 初始化结果存储
  accuracy_values <- numeric(k)
  feature_list <- list()
  # 创建交叉验证的折
  folds <- caret::createFolds(y = test$OTUgroup, k = k)
  # 进行交叉验证
  for (i in 1:k) {
    fold_test <- test[folds[[i]], ]
    fold_train <- test[-folds[[i]], ]

    # 训练 naiveBayes 模型
    a_nb <- e1071::naiveBayes(OTUgroup ~ ., data = fold_train)

    # 得到测试集的预测值
    pred <- predict(a_nb, newdata = fold_test)

    # 计算准确率
    correct_predictions <- sum(pred == fold_test$OTUgroup)
    accuracy <- correct_predictions / nrow(fold_test)
    accuracy_values[i] <- accuracy
   #  print(paste("Fold", i, "Accuracy:", round(accuracy, 3)))  # 输出每个折的准确率

    # 提取特征重要性
    feature_importance <- data.frame(Feature = names(a_nb$tables),
                                     Importance = sapply(a_nb$tables, function(table) mean(apply(table, 1, var))))
    feature_list[[i]] <- feature_importance
  }

  # 平均准确率
  mean_accuracy <- mean(accuracy_values)
  accuracy_result <- paste("naiveBayes Average Accuracy:", round(mean_accuracy, 3))
  print(accuracy_result)

  # 合并重要变量
  if (length(feature_list) > 0) {
    combined_importance <- Reduce(function(x, y) merge(x, y, by = "Feature", all = TRUE), feature_list)
    combined_importance[is.na(combined_importance)] <- 0
    combined_importance$AvgImportance <- rowMeans(combined_importance[,-1])
    importance_df <- combined_importance[, c("Feature", "AvgImportance")]
    importance_df <- importance_df[order(importance_df$AvgImportance, decreasing = TRUE), ]
  } else {
    importance_df <- data.frame(Feature = character(0), Importance = numeric(0))
  }
  importance_df
  # 返回结果
  list(Accuracy = accuracy_result, Importance = importance_df)
}

