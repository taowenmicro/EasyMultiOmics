#' @title Neural network screening of characteristic genes
#' @description
#' A neural network, one of the machine learning methods, was used to screen for characteristic
#' genes, and the model was evaluated using k-fold cross-validation.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table,
#' transcriptome functional classification table, and sample metadata.
#' @param top The top genes to consider.
#' @param seed The random seed for reproducibility.
#' @param k The number of folds for cross-validation.
#' @return A list object including the following components:
#' \item{Accuracy}{The average accuracy of the neural network model.}
#' \item{Importance}{A data frame showing the feature importance ranked in descending order.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' library(dplyr)
#' library(ggClusterNet)
#' data(ps.trans)
#' ps =ps.trans %>% filter_OTU_ps(Top = 100)
#' res =nnet.trans(ps=ps, top = 100, seed = 1010, k = 5)
#' accuracy = res[[1]]
#' accuracy
#' importance = res[[2]]
#' importance

nnet.trans <- function(ps=ps, top = 20, seed = 1010, k = 5) {
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
  feature_importance_list <- list()
  # 创建交叉验证的折
  folds <- createFolds(y = test$OTUgroup, k = k)
  # 进行交叉验证
  for (i in 1:k) {
    fold_test <- test[folds[[i]], ]
    fold_train <- test[-folds[[i]], ]

    # 训练 nnet 模型
    a_nn <- nnet::nnet(OTUgroup ~ ., data = fold_train, size = 2, rang = 0.1, decay = 5e-4, MaxNWts = 84581, maxit = 200)

    # 得到测试集的预测值
    pred <- predict(a_nn, newdata = fold_test, type = 'class')

    # 计算准确率
    correct_predictions <- sum(pred == fold_test$OTUgroup)
    accuracy <- correct_predictions / nrow(fold_test)
    accuracy_values[i] <- accuracy
    #  print(paste("Fold", i, "Accuracy:", round(accuracy, 3)))  # 输出每个折的准确率

    # 提取特征重要性
    num_features <- ncol(fold_train) - 1  # 减去目标列
    input_to_hidden_weights <- a_nn$wts[1:(num_features * 2)]  # 输入到隐藏层的权重
    feature_importance <- rowSums(matrix(input_to_hidden_weights, nrow = num_features, byrow = TRUE))
    feature_importance_list[[i]] <- feature_importance
  }

  # 平均准确率
  mean_accuracy <- mean(accuracy_values)
  accuracy_result <- paste("nnet Average Accuracy:", round(mean_accuracy, 3))
  # print(accuracy_result)

  # 合并重要变量
  if (length(feature_importance_list) > 0) {
    combined_importance <- do.call(cbind, feature_importance_list)
    avg_importance <- rowMeans(combined_importance, na.rm = TRUE)
    importance_df <- data.frame(Feature = colnames(fold_train)[-ncol(fold_train)], Importance = avg_importance)
    importance_df <- importance_df[order(abs(importance_df$Importance), decreasing = TRUE), ]  # 使用绝对值排序
  } else {
    importance_df <- data.frame(Feature = character(0), Importance = numeric(0))
  }
  importance_df$Importance =  abs(importance_df$Importance)
  # 返回结果
  list(Accuracy = accuracy_result, Importance = importance_df)
}
