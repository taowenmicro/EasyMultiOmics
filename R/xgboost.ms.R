xgboost.ms <- function(ps = ps, seed = 200, top = 20,k =5 ) {
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
  tune_grid <- expand.grid(nrounds = c(50, 100, 200),
                           max_depth = c(6, 8, 10),
                           colsample_bytree = c(0.6, 0.8, 1),
                           eta = c(0.1, 0.3),
                           gamma = 0,
                           min_child_weight = c(3, 4, 5),
                           subsample = c(0.6, 0.8))

  # 训练模型
  model <- train(x = df[, -ncol(df)],
                 y = df$Group,
                 method = "xgbTree",
                 objective = "binary:logistic",
                 trControl = train_control,
                 tuneGrid = tune_grid,
                 metric = "AUC",
                 verbosity = 0)

  # 输出模型摘要
 # print(summary(model))

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
