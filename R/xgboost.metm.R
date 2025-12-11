#' @title XGBoost model screening of characteristic microorganisms
#' @description
#' XGBoost, one of the machine learning methods, is used to screen for
#' characteristic microorganisms, and the model is evaluated using
#' k-fold cross-validation.
#'
#' @param ps A phyloseq object containing otu, tax, and sample data.
#' @param top The top microorganisms to keep (feature selection by abundance).
#' @param seed Random seed for reproducibility.
#' @param k Number of folds for cross-validation.
#' @param adjust Logical; whether to tune XGBoost hyperparameters (default FALSE).
#'
#' @return A list with components:
#' \item{Accuracy}{Overall accuracy of the final model.}
#' \item{Accuracy_text}{A character string summarizing the accuracy.}
#' \item{Importance}{A \code{varImp} object with feature importance.}
#' \item{Model}{The fitted \code{caret} model object.}
#' @export
xgboost.metm <- function(ps,
                         seed   = 200,
                         top    = 20,
                         k      = 5,
                         adjust = FALSE) {

  ## ---- 依赖检查 -------------------------------------------------------------
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be a phyloseq object.", call. = FALSE)
  }
  if (!requireNamespace("mia", quietly = TRUE)) {
    stop("Package 'mia' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Package 'SummarizedExperiment' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("xgboost", quietly = TRUE)) {
    stop("Package 'xgboost' is required but not installed.", call. = FALSE)
  }

  set.seed(seed)

  ## ---- 1. 筛 top OTU -------------------------------------------------------
  ps.cs <- ps %>% filter_OTU_ps(top)

  ## ---- 2. phyloseq -> TreeSummarizedExperiment -----------------------------
  # makeTreeSummarizedExperimentFromPhyloseq 已弃用，改用 convertFromPhyloseq
  tse <- mia::convertFromPhyloseq(ps.cs)

  ## ---- 3. CLR 变换 ---------------------------------------------------------
  tse <- mia::transformAssay(
    tse,
    assay.type  = "counts",
    method      = "clr",
    MARGIN      = "samples",
    pseudocount = 1
  )
  assay_data <- SummarizedExperiment::assay(tse, "clr")
  assay_data <- t(assay_data)

  ## ---- 4. 提取分组信息 -----------------------------------------------------
  cd <- SummarizedExperiment::colData(tse)

  if (!"Group" %in% colnames(cd)) {
    stop("`colData(tse)` 中未找到名为 'Group' 的分组信息，请检查样本分组列名。", call. = FALSE)
  }

  labels <- as.factor(cd$Group)

  df <- as.data.frame(assay_data)
  df$Group <- labels
  df$Group <- droplevels(df$Group)

  n_class <- length(levels(df$Group))

  ## ---- 5. caret 训练控制 ---------------------------------------------------
  # 注意：caret 里 AUC 对应的 metric 名叫 "ROC"，需要 twoClassSummary
  if (n_class == 2) {
    train_control <- caret::trainControl(
      method          = "cv",
      number          = k,
      classProbs      = TRUE,
      savePredictions = "final",
      summaryFunction = caret::twoClassSummary,
      allowParallel   = FALSE   # 如果你并行环境配置好，可以改 TRUE
    )
    metric <- "ROC"
  } else {
    train_control <- caret::trainControl(
      method          = "cv",
      number          = k,
      classProbs      = FALSE,
      savePredictions = "final",
      allowParallel   = FALSE
    )
    metric <- "Accuracy"
  }

  ## ---- 6. 超参数网格 -------------------------------------------------------
  if (isTRUE(adjust)) {
    tune_grid <- expand.grid(
      nrounds          = c(50, 100, 200),
      max_depth        = c(6, 8, 10),
      colsample_bytree = c(0.6, 0.8, 1),
      eta              = c(0.1, 0.3),
      gamma            = 0,
      min_child_weight = c(1, 3, 4, 5),
      subsample        = c(0.6, 0.8)
    )
  } else {
    tune_grid <- expand.grid(
      nrounds          = 100,
      max_depth        = 6,
      colsample_bytree = 1,
      eta              = 0.3,
      gamma            = 0,
      min_child_weight = 1,
      subsample        = 1
    )
  }

  ## ---- 7. 训练 XGBoost 模型（交给 caret 管） --------------------------------
  # 注意：caret::xgbTree 会根据因变量是否为因子自动选择合适 objective，
  # 一般不需要我们手工传 objective / num_class。
  model <- caret::train(
    x         = df[, -ncol(df)],
    y         = df$Group,
    method    = "xgbTree",
    trControl = train_control,
    tuneGrid  = tune_grid,
    metric    = metric,
    verbose   = 0
  )

  print(summary(model))

  ## ---- 8. 计算总体准确率 ----------------------------------------------------
  pred <- stats::predict(model, newdata = df[, -ncol(df)])
  cm   <- caret::confusionMatrix(pred, df$Group)
  acc  <- cm$overall["Accuracy"]
  acc_txt <- paste("xgboost accuracy:", round(acc, 3))
  print(acc_txt)

  ## ---- 9. 变量重要性 -------------------------------------------------------
  imp <- caret::varImp(model, scale = FALSE)

  ## ---- 10. 返回结果 ---------------------------------------------------------
  list(
    Accuracy      = unname(acc),
    Accuracy_text = acc_txt,
    Importance    = imp,
    Model         = model
  )
}
