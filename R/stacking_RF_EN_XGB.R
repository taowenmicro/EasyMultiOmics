
# # meta-learner
# fit <- stacking_RF_EN_XGB(pst, outcome_col="Group", k=5)
#
# # 查看 CV 性能比较
# fit$summary_metrics
#
# # 预测
# pred <- fit$predict(pst)
# head(pred)
#
# library(phyloseq)
# library(randomForest)
# library(xgboost)
# library(glmnet)
# library(caret)
# library(pROC)
# library(zCompositions)
# library(compositions)


# --- 预处理：零值替换 + CLR ---
.preproc_clr <- function(X) {
  X <- as.matrix(X)
  X[X < 0] <- 0
  rs <- rowSums(X)
  rs[rs == 0] <- 1e-6

  # 保存原始行名和列名
  rn <- rownames(X)
  cn <- colnames(X)

  # 零值替换（强制 matrix 输出）
  repl <- as.matrix(zCompositions::cmultRepl(X, method = "CZM", output = "p-counts"))

  # 确保维度一致（如果 cmultRepl 删除了列，就补回来）
  if (ncol(repl) != ncol(X)) {
    tmp <- matrix(0, nrow = nrow(X), ncol = ncol(X))
    colnames(tmp) <- cn
    rownames(tmp) <- rn
    common <- intersect(colnames(repl), cn)
    tmp[, common] <- repl[, common, drop = FALSE]
    repl <- tmp
  }

  # CLR 转换
  clrX <- compositions::clr(acomp(repl))
  clrX <- as.matrix(clrX)

  # 恢复行列名
  rownames(clrX) <- rn
  colnames(clrX) <- cn

  clrX
}


# --- 主函数 ---
#' @export
#'
stacking_RF_EN_XGB <- function(ps, outcome_col="Group", k=5, seed=123) {
  set.seed(seed)

  # 1. 提取数据
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- as.data.frame(otu, check.names=FALSE)
  meta <- as.data.frame(sample_data(ps))

  if (!outcome_col %in% names(meta)) stop("sample_data 中没有 outcome 列")
  y <- factor(meta[[outcome_col]])
  if (nlevels(y) != 2) stop("只支持二分类")
  y_num <- as.numeric(y) - 1

  # 2. 预处理：CLR
  X <- .preproc_clr(otu)
  # colnames(X) <- colnames(otu)

  # 3. K 折交叉验证
  folds <- caret::createFolds(y, k=k, returnTrain=FALSE)
  oof_rf <- oof_en <- oof_xgb <- rep(NA, nrow(X))

  for (i in seq_along(folds)) {
    idx_te <- folds[[i]]
    idx_tr <- setdiff(seq_len(nrow(X)), idx_te)
    xtr <- X[idx_tr,,drop=FALSE]; ytr <- y[idx_tr]; ytr_num <- y_num[idx_tr]
    xte <- X[idx_te,,drop=FALSE]

    # --- RF ---
    rf_fit <- randomForest(x=xtr, y=ytr, ntree=300)
    oof_rf[idx_te] <- predict(rf_fit, xte, type="prob")[,2]

    # --- Elastic Net ---
    en_fit <- cv.glmnet(as.matrix(xtr), ytr_num, family="binomial", alpha=0.5)
    oof_en[idx_te] <- as.numeric(predict(en_fit, newx=as.matrix(xte), s="lambda.min", type="response"))

    # --- XGB ---
    xgb_fit <- xgboost(data=as.matrix(xtr), label=ytr_num, objective="binary:logistic",
                       nrounds=200, max_depth=4, eta=0.1,
                       subsample=0.8, colsample_bytree=0.8, verbose=0)
    oof_xgb[idx_te] <- predict(xgb_fit, as.matrix(xte))
  }

  # 4. 元学习器（XGB）
  oof_df <- data.frame(rf=oof_rf, en=oof_en, xgb=oof_xgb)
  meta_fit <- xgboost(data=as.matrix(oof_df), label=y_num,
                      objective="binary:logistic",
                      nrounds=200, max_depth=3, eta=0.05, verbose=0)

  # 5. OOF 性能
  acc <- function(prob) mean((prob > 0.5) == (y_num==1))
  aucv <- function(prob) as.numeric(pROC::auc(y, prob))
  summary_metrics <- data.frame(
    Model=c("RF","ElasticNet","XGB","Stacking-XGB"),
    Accuracy=c(acc(oof_rf), acc(oof_en), acc(oof_xgb),
               acc(predict(meta_fit, as.matrix(oof_df)))),
    AUC=c(aucv(oof_rf), aucv(oof_en), aucv(oof_xgb),
          aucv(predict(meta_fit, as.matrix(oof_df))))
  )

  # 6. 全量训练基模型
  rf_full <- randomForest(x=X, y=y, ntree=300)
  en_full <- cv.glmnet(as.matrix(X), y_num, family="binomial", alpha=0.5)
  xgb_full <- xgboost(data=as.matrix(X), label=y_num,
                      objective="binary:logistic",
                      nrounds=200, max_depth=4, eta=0.1,
                      subsample=0.8, colsample_bytree=0.8, verbose=0)

  # 7. 预测函数（带特征对齐）
  predict_fun <- function(new_ps) {
    otu2 <- as(otu_table(new_ps), "matrix")
    if (taxa_are_rows(new_ps)) otu2 <- t(otu2)
    otu2 <- as.data.frame(otu2, check.names=FALSE)
    # 对齐特征
    X2 <- matrix(0, nrow=nrow(otu2), ncol=ncol(X), dimnames=list(rownames(otu2), colnames(X)))
    common <- intersect(colnames(otu2), colnames(X))
    X2[, common] <- as.matrix(otu2[, common, drop=FALSE])
    # CLR
    X2 <- .preproc_clr(X2)

    # 基模型预测
    p_rf <- predict(rf_full, X2, type="prob")[,2]
    p_en <- as.numeric(predict(en_full, newx=as.matrix(X2), s="lambda.min", type="response"))
    p_xgb <- predict(xgb_full, as.matrix(X2))

    # Meta 预测
    p_meta <- predict(meta_fit, as.matrix(data.frame(rf=p_rf, en=p_en, xgb=p_xgb)))
    cls <- ifelse(p_meta > 0.5, levels(y)[2], levels(y)[1])

    data.frame(SampleID=rownames(X2),
               prob_RF=p_rf, prob_EN=p_en, prob_XGB=p_xgb,
               prob_STACK=p_meta, PredClass=cls)
  }

  list(summary_metrics=summary_metrics,
       rf_full=rf_full, en_full=en_full, xgb_full=xgb_full,
       meta_fit=meta_fit, feature_names=colnames(X),
       predict=predict_fun)
}
