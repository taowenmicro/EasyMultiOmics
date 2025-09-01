#' @title lasso model screening of characteristic microorganisms
#' @description
#' lasso, one of the machine learning methods, was used to screen for characteristic
#' microorganisms, and the model was evaluated using k-fold cross-validation.
#' @param ps A phyloseq format file used as an alternative for the input containing otu, tax, and map.
#' @param top The top microorganisms to consider.
#' @param seed The random seed for reproducibility.
#' @param k The number of folds for cross-validation.
#' @return A list object including the following components:
#' \item{Accuracy}{The average accuracy of the lasso model.}
#' \item{Importance}{A data frame showing the feature importance ranked in descending order.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' library(ggClusterNet)
#' library(caret)
#' library(glmnet)
#' res =lasso.micro (ps =  ps, top = 100, seed = 1010, k = 5)
#' accuracy = res[[1]]
#' accuracy
#' importance = res[[2]]
#' importance
lasso.micro <- function(ps= ps, top = 20, seed = 1010, k = 5) {

  ps.cs <- ps %>% filter_OTU_ps(top)
  map <- as.data.frame(phyloseq::sample_data(ps.cs))
  otutab <- as.data.frame(t(ggClusterNet::vegan_otu(ps.cs)))
  colnames(otutab) <- gsub("-", "_", colnames(otutab))
  test <- as.data.frame(t(otutab))
  test$OTUgroup <- factor(map$Group)
  levels(test$OTUgroup) <- c(1:length(levels(test$OTUgroup))-1)
  # 初始化结果存储
  accuracy_values <- numeric(k)
  importance_list <- list()

  # 创建交叉验证的折
  set.seed(seed)
  folds <- createFolds(y = test$OTUgroup, k = k)
  # 进行交叉验证
  for (i in 1:k) {
    fold_test <- test[folds[[i]], ]
    fold_train <- test[-folds[[i]], ]
    # 创建模型矩阵
    x_train <- model.matrix(OTUgroup ~ . - 1, data = fold_train)
    y_train <- as.numeric(fold_train$OTUgroup)
    x_test <- model.matrix(OTUgroup ~ . - 1, data = fold_test)
    y_test <- as.numeric(fold_test$OTUgroup)
    # 训练 s模型
    if (length(levels(test$OTUgroup))>2) {
      lasso_model <- cv.glmnet(x_train, y_train, family = "multinomial", alpha = 1)
    }else{
      lasso_model <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1, type.measure = "class")
    }
    # 得到测试集的预测值
    pred <- predict(lasso_model, newx = x_test,s = "lambda.min",type="class")
    # 计算准确率
    correct_predictions <- sum(pred == y_test)
    accuracy <- correct_predictions / nrow(fold_test)
    accuracy_values[i] <- accuracy
    # 提取特征重要性
    if (length(levels(test$OTUgroup))>2) {
      coef_matrix <- coef(lasso_model,s = "lambda.min")
      importance <- Reduce("+", lapply(coef_matrix, function(c) abs(as.matrix(c)[-1, ])))
      names(importance) <- rownames(coef(lasso_model,s = "lambda.min")[[1]])[-1]}else{
        importance <- abs(as.vector(coef(lasso_model, s = "lambda.min"))[-1])
        names(importance) <- rownames(coef(lasso_model))[-1]}
    importance_list[[i]] <- importance
  }

  # 平均准确率
  mean_accuracy <- mean(accuracy_values)
  accuracy_result <- paste("LASSO Accuracy:", round(mean_accuracy, 3))
  print(accuracy_result)

  # 合并重要变量
  if (length(importance_list) > 0) {
    combined_importance <- do.call(cbind, importance_list)
    avg_importance <- rowMeans(combined_importance, na.rm = TRUE)
    importance_df <- data.frame(Feature = names(avg_importance), Importance = avg_importance)
    importance_df <- importance_df[order(abs(importance_df$Importance), decreasing = TRUE), ]
  } else {
    importance_df <- data.frame(Feature = character(0), Importance = numeric(0))
  }
  importance_df =  importance_df[order(importance_df$Importance,decreasing = TRUE),]
  # 返回结果
  list(Accuracy = accuracy_result, Importance = importance_df)
}


