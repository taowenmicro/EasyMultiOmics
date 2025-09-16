# BiocManager::install("neuralnet")
#
# mlp.metm (ps.cs, top = 20, seed = 1010, k = 5,
#           hidden = c(200, 100, 50)
#
#           )
# ps = pst
# top = 20
# seed = 1010
# k = 5
# hidden = c(20, 10)
# group = "Group"
#
#
# mlp_metm_neuralnet(pst, group = "Group",
#                    hidden = c(200, 100, 50),
#                                seed = 123,
#                                k = 5,
#                                threshold = 0.01)



#' @title MLP classification for microbiome data
#' @description
#' Train a multi-layer perceptron (MLP) using RSNNS on microbiome data (phyloseq object),
#' evaluate model accuracy by k-fold cross-validation, and compute feature importance
#' from input-to-hidden layer weights.
#' @param ps A phyloseq object containing otu, tax, and metadata.
#' @param group Group column name in sample_data (must be binary classification).
#' @param top Top N OTUs to keep (default 50).
#' @param hidden Number of hidden units in the MLP.
#' @param seed Random seed for reproducibility.
#' @param k Number of folds for cross-validation.
#' @return A list containing:
#' \item{Accuracy}{Average cross-validation accuracy.}
#' \item{Importance}{Data frame of feature importance.}
#' @export
#' @examples
#' res <- mlp_metm(ps, group = "Group", top = 50, hidden = 5, seed = 123, k = 5)
#' res$Accuracy
#' head(res$Importance)
mlp_metm_neuralnet <- function(ps,
                               group = "Group",
                               hidden = c(50, 20),
                               seed = 123,
                               k = 5,
                               threshold = 0.01) {
  set.seed(seed)

  # ---- data preparation ----
  # ps.cs <- ps %>% ggClusterNet::filter_OTU_ps(top)
  map <- as.data.frame(phyloseq::sample_data(ps))
  otu <- ggClusterNet::vegan_otu(ps)  %>% as.data.frame()
  y <- factor(map[[group]])

  if (length(levels(y)) != 2) {
    stop("Only binary classification is supported.")
  }

  dat <- data.frame(otu, y = as.numeric(y) - 1) # y must be numeric (0/1)
  folds <- caret::createFolds(y, k = k)

  accs <- numeric(k)
  importance_list <- list()

  # ---- cross validation ----
  for (i in seq_len(k)) {
    train <- dat[-folds[[i]], ]
    test  <- dat[folds[[i]], ]

    # build formula: y ~ x1 + x2 + ...
    fml <- as.formula(
      paste("y ~", paste(colnames(train)[-ncol(train)], collapse = " + "))
    )

    # train neuralnet
    model <- neuralnet::neuralnet(
      fml,
      data = train,
      hidden = hidden,
      linear.output = FALSE,
      threshold = threshold,
      stepmax = 5e5
    )

    # prediction
    preds <- predict(model, test[, -ncol(test)])
    pred_class <- ifelse(preds > 0.5, 1, 0)
    accs[i] <- mean(pred_class == test$y)

    # feature importance (absolute input-hidden weights)
    wts <- model$weights[[1]][[1]]  # 输入层 -> 第一隐藏层
    feat_imp <- apply(abs(wts[-1, , drop = FALSE]), 1, sum) # -1 去掉 bias
    names(feat_imp) <- colnames(train)[-ncol(train)]
    importance_list[[i]] <- feat_imp
  }

  # ---- average accuracy ----
  mean_acc <- mean(accs)

  # ---- aggregate feature importance ----
  if (length(importance_list) > 0) {
    imp_mat <- do.call(cbind, importance_list)
    avg_imp <- rowMeans(imp_mat, na.rm = TRUE)
    importance_df <- data.frame(
      Feature = names(avg_imp),
      Importance = avg_imp
    ) %>% dplyr::arrange(desc(Importance))
  } else {
    importance_df <- data.frame(Feature = character(0), Importance = numeric(0))
  }

  return(list(
    Accuracy = mean_acc,
    Importance = importance_df
  ))
}
