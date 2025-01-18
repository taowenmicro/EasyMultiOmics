#' @title k-Nearest Neighbors (k-NN) Classification for Microbial Community Data
#'
#' @description
#' The `knn.micro` function applies the k-Nearest Neighbors (k-NN) algorithm to classify microbial community data.
#' It uses a Centered Log-Ratio (CLR) transformation for compositional data and performs k-fold cross-validation to evaluate model accuracy.
#' Additionally, it computes feature importance by assessing the impact of each feature on classification accuracy.
#'
#' @param ps A `phyloseq` object containing microbial community data, including the OTU table and sample metadata.
#' @param seed An integer used to set the random seed for reproducibility. Default is `6358`.
#' @param k An integer specifying the number of folds for k-fold cross-validation. Default is `5`.
#' @param top An integer specifying the number of top OTUs to retain based on abundance. Default is `20`.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `Accuracy`: A string reporting the average k-NN classification accuracy across all folds.
#'   \item `Importance`: A data frame with feature importance metrics, showing the classification accuracy when using each feature individually.
#' }
#'
#' @examples
#' \dontrun{
#' res = knn.micro(ps =pst, top = 20,k =5 )
#' accuracy = res[[1]]
#' importance = res[[2]]
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#'
knn.micro <- function(ps =ps, seed = 6358, k = 5,top = 20) {
  set.seed(seed)

  # 数据准备
  ps.cs <- ps %>% filter_OTU_ps(top)
  tse <- ps.cs %>% makeTreeSummarizedExperimentFromPhyloseq()

  # 应用 CLR 变换
  tse <- mia::transformAssay(tse, assay.type = "counts", method = "clr",
                             MARGIN = "samples", pseudocount = 1)

  # 获取 assay 数据
  assay_data <- assay(tse, "clr")
  assay_data <- t(assay_data)

  # 转换为数据框
  df <- as.data.frame(assay_data)

  # 添加标签
  labels <- colData(tse)$Group
  labels <- as.factor(labels)
  df$OTUgroup <- labels

  # 初始化结果存储
  accuracy_values <- numeric(k)

  # 创建交叉验证的折
  folds <- createFolds(y = df$OTUgroup, k = k)

  # 进行交叉验证
  for (i in 1:k) {
    fold_test <- df[folds[[i]], ]
    fold_train <- df[-folds[[i]], ]

    # 训练 k-NN 模型
    a_kNN <- kknn(OTUgroup ~ .,
                  fold_train,
                  fold_test,
                  distance = 1,
                  kernel = "triangular")

    # 得到测试集的预测值
    fit <- fitted(a_kNN)

    # 计算准确率
    correct_predictions <- sum(fit == fold_test$OTUgroup)
    accuracy <- correct_predictions / nrow(fold_test)
    accuracy_values[i] <- accuracy
  }

  # 平均准确率
  mean_accuracy <- mean(accuracy_values)
  accuracy_result <- paste("k-NN Accuracy:", round(mean_accuracy, 3))
  print(accuracy_result)

  # 特征重要性评估：逐个特征训练模型并评估其影响
  feature_importance <- data.frame(Feature = names(df)[-ncol(df)], Importance = numeric(ncol(df) - 1))

  for (i in 1:(ncol(df) - 1)) {
    single_feature <- names(df)[i]
    accuracy_values <- numeric(k)

    for (j in 1:k) {
      fold_test <- df[folds[[j]], c(single_feature, "OTUgroup")]
      fold_train <- df[-folds[[j]], c(single_feature, "OTUgroup")]

      # 训练 k-NN 模型
      a_kNN <- kknn(OTUgroup ~ .,
                    fold_train,
                    fold_test,
                    distance = 1,
                    kernel = "triangular")

      # 得到测试集的预测值
      fit <- fitted(a_kNN)

      # 计算准确率
      correct_predictions <- sum(fit == fold_test$OTUgroup)
      accuracy <- correct_predictions / nrow(fold_test)
      accuracy_values[j] <- accuracy
    }

    # 平均准确率
    feature_importance$Importance[i] <- mean(accuracy_values)
  }

  # 返回结果
  list(Accuracy = accuracy_result, Importance = feature_importance)
}

