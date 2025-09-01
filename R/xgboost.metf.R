#' @title XGBoost model screening of characteristic genes
#' @description
#' XGBoost, one of the machine learning methods, was used to screen for characteristic
#' genees, and the model was evaluated using k-fold cross-validation.
#' @param ps A phyloseq format file used as an alternative for the input containing metagenome functional composition table,
#' metagenome functional classification table, and sample metadata.
#' @param top The top genes to consider.
#' @param seed The random seed for reproducibility.
#' @param k The number of folds for cross-validation.
#' @param adjust A logical value for whether to adjust XGBoost model parameters,default F.
#' @return A list object including the following components:
#' \item{Accuracy}{The average accuracy of the XGBoost model.}
#' \item{Importance}{A data frame showing the feature importance ranked in descending order.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' library(ggClusterNet)
#' library(caret)
#' library(xgboost)
#' library(Ckmeans.1d.dp)
#' library(tidyverse)
#' library(mia)
#' ps =ps.kegg %>% filter_OTU_ps(Top = 1000)
#' res = xgboost.metf(ps =ps, top = 30)
#' accuracy = res[[1]]
#' accuracy
#' importance = res[[2]]
#' importance
xgboost.metf <- function(ps = ps, seed = 200, top = 20,k =5,adjust=F ) {
  set.seed(seed)

  # 数据准备
  ps.cs <- ps %>% filter_OTU_ps(top)
  tse <- ps.cs %>% makeTreeSummarizedExperimentFromPhyloseq()

  # 应用CLR变换
  tse <- mia::transformAssay(tse, assay.type = "counts", method = "clr",
                             MARGIN = "samples", pseudocount = 1)

  # 获取assay数据
  assay_data <- assay(tse, "clr")
  assay_data <- t(assay_data)

  # 转换为数据框
  df <- as.data.frame(assay_data)

  # 添加标签
  labels <- colData(tse)$Group
  labels <- as.factor(labels)
  df$Group <- labels

  # 指定训练控制
  train_control <- trainControl(method = "cv", number = k,
                                classProbs = TRUE,
                                savePredictions = "final",
                                allowParallel = TRUE)

  # 指定超参数调优网格
  if(adjust==TRUE){
    tune_grid <- expand.grid(nrounds = c(50, 100, 200),
                             max_depth = c(6, 8, 10),
                             colsample_bytree = c(0.6, 0.8, 1),
                             eta = c(0.1, 0.3),
                             gamma = 0,
                             min_child_weight = c(1,3, 4, 5),
                             subsample = c(0.6, 0.8))
  }else{
    tune_grid <- expand.grid(nrounds = 100,
                             max_depth = 6,
                             colsample_bytree = 1,
                             eta = 0.3,
                             gamma = 0,
                             min_child_weight = 1,
                             subsample = 1)
  }

  # 训练模型
  if (length(unique(labels))==2) {
    model <- train(x = df[, -ncol(df)],
                   y = df$Group,
                   method = "xgbTree",
                   objective = "binary:logistic",
                   trControl = train_control,
                   tuneGrid = tune_grid,
                   metric = "AUC",
                   verbosity = 0)
  }else{
    model <- train(x = df[, -ncol(df)],
                   y = df$Group,
                   method = "xgbTree",
                   objective = "multi:softprob",
                   num_class =length(unique(labels)),
                   trControl = train_control,
                   tuneGrid = tune_grid,
                   metric = "AUC",
                   verbosity = 0)}

  # 输出模型摘要
  print(summary(model))

  # 计算准确率
  pred <- predict(model, newdata = df[, -ncol(df)])
  conf_matrix <- confusionMatrix(pred, df$Group)
  accuracy <- conf_matrix$overall["Accuracy"]
  accuracy_result <- paste("xgboost accuracy:", round(accuracy, 3))
  print(accuracy_result)

  # 提取特征重要性
  importance <- varImp(model, scale = FALSE)
  #row.names(importance) =  gsub("OTU","",row.names(importance))
  # 返回结果
  list(Accuracy = accuracy_result, Importance = importance)
}
