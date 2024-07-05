glm.omics <- function(ps = ps, k = 5) {
  # 数据准备
  map <- as.data.frame(phyloseq::sample_data(ps))
  otutab <- as.data.frame(t(ggClusterNet::vegan_otu(ps)))
  colnames(otutab) <- gsub("-", "_", colnames(otutab))
  test <- as.data.frame(t(otutab))
  test$group <- factor(map$Group)
  test$group

  colnames(test) <- paste("OTU", colnames(test), sep = "")

  colnames(test) <- gsub("-", "_", colnames(test))
  colnames(test) <- gsub("[/]", "_", colnames(test))
  colnames(test) <- gsub("[(]", "_", colnames(test))
  colnames(test) <- gsub("[)]", "_", colnames(test))
  colnames(test) <- gsub("[:]", "_", colnames(test))
  colnames(test) <- gsub("[[]", "", colnames(test))
  colnames(test) <- gsub("[]]", "_", colnames(test))
  colnames(test) <- gsub("[#]", "_", colnames(test))
  colnames(test) <- gsub("[+]", "_", colnames(test))
  colnames(test) <- gsub(" ", "_", colnames(test))

  test <- dplyr::select(test, OTUgroup, everything())
  train <- test

  folds <- createFolds(y = test$OTUgroup, k = k)

  fc <- as.numeric()
  mod_pre <- as.numeric()
  accuracy_values <- numeric(k)

  for (i in 1:k) {
    fold_test<-train[folds[[i]],]
    fold_train<-train[-folds[[i]],]
    model<-glm(OTUgroup~.,family='binomial',data=fold_train)
    model
    model_pre<-predict(model,type='response',newdata=fold_test)
    model_pre

    fc<-append(fc,fold_test$OTUgroup)
    mod_pre<-append(mod_pre,as.numeric(model_pre))
    # # 将概率转换为类标签并计算准确率
    # predicted_class <- ifelse(model_pre_fold > 0.5, "WT", "OE")
    # correct_predictions <- sum(predicted_class == fold_test$OTUgroup)
    # accuracy_values[i] <- correct_predictions / nrow(fold_test)
  }

  # 计算AUC
  pred <- prediction(mod_pre, fc)
  auc <- performance(pred, "auc")@y.values[[1]]
  auc_result <- paste("GLM AUC:", round(auc, 3))

  # 提取特征重要性
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  importance <- data.frame(
    Feature = rownames(coefficients),
    Coefficient = coefficients[, 1],
    Std_Error = coefficients[, 2],
    z_value = coefficients[, 3],
    p_value = coefficients[, 4]
  )

  # 按照绝对系数大小排序特征重要性
  importance <- importance[order(abs(importance$Coefficient), decreasing = TRUE), ]
  row.names(importance) =  gsub("OTU","",row.names(importance))
  importance$Feature = gsub("OTU","",importance$Feature)
  # 输出结果
  list(AUC = auc_result, Importance = importance)
}


