#' @title Bagging screening of features in multi-omics data
#' @description
#' Bagging, one of the machine learning methods, was used to screen for characteristic
#' microorganisms, and the model was evaluated using k-fold cross-validation.
#' @param ps A phyloseq format file used as an alternative for the input containing otu, tax, and map.
#' @param top The top microorganisms to consider.
#' @param seed The random seed for reproducibility.
#' @param k The number of folds for cross-validation.
#' @return A list object including the following components:
#' \item{Accuracy}{The average accuracy of the bagging model.}
#' \item{Importance}{A data frame showing the feature importance ranked in descending order.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(ipred)
#' library(dplyr)
#' library(ggClusterNet)
#' library(caret)
#' res =bagging.omics(ps =  ps03, top = 100, seed = 1010, k = 5)
#' accuracy = res[[1]]
#' accuracy
#' importance = res[[2]]
#' importance

bagging.omics <- function(ps=ps, top = 20, seed = 1010, k = 5) {

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
  set.seed(seed)
  folds <- createFolds(y = test$OTUgroup, k = k)
  #i=1
  # 进行交叉验证
  for (i in 1:k) {
    fold_test <- test[folds[[i]], ]
    fold_train <- test[-folds[[i]], ]

    # 训练 Bagging 模型
    set.seed(seed)
    a_bagging <- bagging(OTUgroup ~ ., data = fold_train, coob = TRUE)

    # 得到测试集的预测值
    pred <- predict(a_bagging, newdata = fold_test, type = 'class')

    # 计算准确率
    correct_predictions <- sum(pred == fold_test$OTUgroup)
    accuracy <- correct_predictions / nrow(fold_test)
    accuracy_values[i] <- accuracy

    # 提取特征重要性
    importance <- caret::varImp(a_bagging)
    importance$ID <- rownames(importance)
    importance_list[[i]] <- importance
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
    combined_importance <- importance_list[[1]]
    for(i in 2:length(importance_list)) {
      combined_importance <- merge(combined_importance, importance_list[[i]], by = "ID", all = TRUE)
      # 注意：这里假设行名是默认的 "row.names"，如果实际情况不同，需要调整 by 参数
    }
    avg_importance <- rowMeans(combined_importance[,-which(colnames(combined_importance) == "ID")], na.rm = TRUE)
    importance_df <- data.frame(Feature=combined_importance$ID , Importance = avg_importance)
    importance_df <- importance_df[order(importance_df$Importance, decreasing = TRUE), ]
  } else {
    importance_df <- data.frame(Feature = character(0), Importance = numeric(0))
  }
  importance_df
  # 返回结果
  list(Accuracy = accuracy_result, Importance = importance_df)
}

