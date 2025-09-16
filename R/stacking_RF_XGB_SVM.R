# ====== 小工具：行归一化 + log1p ======

# library(xgboost)
# library(e1071)
# library(randomForest)
# library(pROC)
# # 假设 ps 是你的 phyloseq 对象，Group 为二分类标签
# fit <- stacking_RF_XGB_SVM(pst, outcome_col = "Group", k = 5)
#
# # 查看 CV 指标对比
# fit$summary_metrics
#
#
# # 预测（返回各基模型概率 + stacking 最终预测）
# pred <- fit$predict(pst)
# head(pred)


.preproc <- function(X) {
  X <- as.matrix(X)
  rs <- rowSums(X); rs[rs == 0] <- 1
  log1p(X / rs)
}



#' @export
stacking_RF_XGB_SVM <- function(ps, outcome_col = "Group", k = 5, seed = 123) {
  set.seed(seed)

  # 1) 提取数据
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- as.data.frame(otu, check.names = FALSE)

  meta <- as.data.frame(sample_data(ps))
  if (!outcome_col %in% names(meta)) stop("sample_data 没有 outcome 列：", outcome_col)
  y <- factor(meta[[outcome_col]])
  if (nlevels(y) != 2) stop("目前仅支持二分类")
  y_num <- as.numeric(y) - 1

  # 2) 简单预处理
  X <- .preproc(otu)
  colnames(X) <- colnames(otu)

  # 3) K 折交叉验证
  folds <- caret::createFolds(y, k = k, returnTrain = FALSE)

  # 保存 OOF 预测
  oof_rf   <- rep(NA, nrow(X))
  oof_xgb  <- rep(NA, nrow(X))
  oof_svm  <- rep(NA, nrow(X))

  # 逐折
  for (i in seq_along(folds)) {
    idx_te <- folds[[i]]
    idx_tr <- setdiff(seq_len(nrow(X)), idx_te)

    xtr <- X[idx_tr, , drop=FALSE]
    ytr <- y[idx_tr]
    ytr_num <- y_num[idx_tr]
    xte <- X[idx_te, , drop=FALSE]

    # --- RF ---
    rf_fit <- randomForest(x = xtr, y = ytr, ntree = 300)
    oof_rf[idx_te] <- predict(rf_fit, xte, type = "prob")[,2]

    # --- XGBoost ---
    dtrain <- xgb.DMatrix(data = xtr, label = ytr_num)
    dtest  <- xgb.DMatrix(data = xte)
    xgb_fit <- xgboost(dtrain, objective = "binary:logistic",
                       nrounds = 200, max_depth = 4, eta = 0.1,
                       subsample = 0.8, colsample_bytree = 0.8,
                       verbose = 0)
    oof_xgb[idx_te] <- predict(xgb_fit, dtest)

    # --- SVM (概率输出) ---
    svm_fit <- e1071::svm(xtr, ytr, probability = TRUE, kernel = "radial", cost = 1)
    svm_pred <- predict(svm_fit, xte, probability = TRUE)
    oof_svm[idx_te] <- attr(svm_pred, "probabilities")[,2]
  }

  # 4) 训练元学习器（逻辑回归）
  oof_df <- data.frame(rf = oof_rf, xgb = oof_xgb, svm = oof_svm)
  meta_fit <- glm(y_num ~ ., data = oof_df, family = binomial())

  # 5) 计算 OOF 指标
  auc_rf   <- auc(y, oof_rf)
  auc_xgb  <- auc(y, oof_xgb)
  auc_svm  <- auc(y, oof_svm)
  auc_meta <- auc(y, predict(meta_fit, type="response"))

  acc_rf   <- mean((oof_rf > 0.5) == (y_num == 1))
  acc_xgb  <- mean((oof_xgb > 0.5) == (y_num == 1))
  acc_svm  <- mean((oof_svm > 0.5) == (y_num == 1))
  acc_meta <- mean((predict(meta_fit, type="response") > 0.5) == (y_num == 1))

  summary_metrics <- data.frame(
    Model = c("RF","XGB","SVM","Stacking"),
    Accuracy = c(acc_rf, acc_xgb, acc_svm, acc_meta),
    AUC = c(auc_rf, auc_xgb, auc_svm, auc_meta)
  )

  # 6) 全量训练基模型（用于最终预测）
  rf_full   <- randomForest(x = X, y = y, ntree = 300)
  dtrain_all <- xgb.DMatrix(data = X, label = y_num)
  xgb_full <- xgboost(dtrain_all, objective = "binary:logistic",
                      nrounds = 200, max_depth = 4, eta = 0.1,
                      subsample = 0.8, colsample_bytree = 0.8,
                      verbose = 0)
  svm_full  <- svm(X, y, probability = TRUE, kernel = "radial", cost = 1)

  # 7) 返回结果对象
  list(
    summary_metrics = summary_metrics,
    meta_fit = meta_fit,
    rf_full = rf_full,
    xgb_full = xgb_full,
    svm_full = svm_full,
    levels = levels(y),
    predict = function(new_ps) {
      otu2 <- as(otu_table(new_ps), "matrix"); if (taxa_are_rows(new_ps)) otu2 <- t(otu2)
      X2 <- .preproc(otu2)
      # RF
      p_rf <- predict(rf_full, X2, type = "prob")[,2]
      # XGB
      dtest <- xgb.DMatrix(data = X2)
      p_xgb <- predict(xgb_full, dtest)
      # SVM
      svm_pred <- predict(svm_full, X2, probability = TRUE)
      p_svm <- attr(svm_pred, "probabilities")[,2]
      # Meta
      p_meta <- as.numeric(plogis(predict(meta_fit,
                                          newdata = data.frame(rf=p_rf, xgb=p_xgb, svm=p_svm))))
      cls <- ifelse(p_meta > 0.5, levels(y)[2], levels(y)[1])
      data.frame(SampleID = rownames(X2),
                 prob_RF = p_rf, prob_XGB = p_xgb, prob_SVM = p_svm,
                 prob_STACK = p_meta, PredClass = cls)
    }
  )
}
