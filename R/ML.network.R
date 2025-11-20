#' @title Machine Learning-based Microbial Association Network
#' @description
#' Construct a machine learning-based association network from a `phyloseq` object.
#' Each taxon (OTU/ASV/species) is iteratively treated as the response variable,
#' with the remaining taxa as predictors. Various machine learning models are used
#' to estimate feature importance, followed by permutation testing to assess
#' significance and construct network edges.
#'
#' @param ps A `phyloseq` object containing OTU/ASV table.
#' @param method Character. Machine learning method to use. Options:
#'   \code{"rf"} (Random Forest, ranger),
#'   \code{"xgb"} (XGBoost),
#'   \code{"glmnet"} (ElasticNet regression),
#'   \code{"svm"} (Support Vector Machine, linear kernel),
#'   \code{"lgbm"} (LightGBM),
#'   \code{"mlp"} (Multilayer Perceptron, nnet),
#'   \code{"gam"} (Generalized Additive Model, mgcv).
#' @param ntree Integer. Number of trees or iterations for tree-based models.
#'   Default is \code{500}.
#' @param normalize Logical. Whether to normalize importance values. Default is \code{TRUE}.
#' @param scale.method Character. Normalization method for importance values.
#'   Options: \code{"l1"}, \code{"softmax"}, \code{"minmax"}, \code{"none"}.
#' @param n_perm Integer. Number of permutations for significance testing.
#'   Default is \code{100}. Set to \code{0} to skip permutation testing.
#' @param p.threshold Numeric. Significance threshold for permutation test.
#'   Default is \code{0.05}.
#' @param seed Integer. Random seed for reproducibility. Default is \code{123}.
#' @param n_cores Integer. Number of parallel cores. Default is \code{2}.
#'
#' @details
#' This function iterates over each taxon, treating it as the response variable
#' and fitting the selected machine learning model against all other taxa.
#' Variable importance is extracted and normalized (optional).
#' A permutation test is then conducted to obtain p-values and select significant
#' associations, forming the network edges.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{edges}{Data frame of significant associations (target, predictor, importance, pval).}
#'   \item{importance_all}{Data frame of all importance scores (target, predictor, importance).}
#'   \item{method}{The machine learning method used.}
#' }
#'
#' @author Wen Tao <\email{2018203048@@njau.edu.cn}>
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("GlobalPatterns")
#'
#' # Subset to 20 taxa for demonstration
#' ps_sub <- prune_taxa(taxa_names(GlobalPatterns)[1:20], GlobalPatterns)
#'
#' res <- ML.network(
#'   ps = ps_sub,
#'   method = "rf",
#'   ntree = 100,
#'   n_perm = 10,
#'   n_cores = 2
#' )
#'
#' head(res$edges)
#' head(res$importance_all)
#' }
#'
#' @export
ML.network <- function(ps,
                       method = c("rf","xgb","glmnet","svm","lgbm","mlp","gam"),
                       ntree = 500,
                       normalize = TRUE,
                       scale.method = c("l1","softmax","minmax","none"),
                       n_perm = 100,
                       p.threshold = 0.05,
                       seed = 123,
                       n_cores = 2) {
  set.seed(seed)
  scale.method <- match.arg(scale.method)
  method <- match.arg(method)

  # === 1) Extract OTU table ===
  otu <- as(otu_table(ps), "matrix")
  if (taxa_are_rows(ps)) otu <- t(otu)
  otu <- as.data.frame(otu, check.names = FALSE)

  taxa_names_vec <- colnames(otu)
  n_taxa <- ncol(otu)

  # --- helper: fit model + extract importance ---
  get_importance <- function(y, X, all_pred, method, ntree) {
    df <- data.frame(y = y, X, check.names = FALSE)

    if (method == "rf") {
      fit <- ranger::ranger(
        dependent.variable.name = "y",
        data = df,
        num.trees = ntree,
        importance = "permutation"
      )
      imp <- fit$variable.importance
    }
    if (method == "xgb") {
      dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
      fit <- xgboost::xgboost(
        data = dtrain,
        nrounds = 30,
        max_depth = 3,
        subsample = 0.8,
        colsample_bytree = 0.8,
        learning_rate = 0.2,
        objective = "reg:squarederror",
        verbose = 0,
        nthread = 1
      )
      imp_tab <- xgboost::xgb.importance(model = fit)
      imp <- setNames(imp_tab$Gain, imp_tab$Feature)
    }
    if (method == "glmnet") {
      fit <- glmnet::cv.glmnet(as.matrix(X), y, alpha = 0.5, nfolds = 5)
      coefs <- coef(fit, s = "lambda.min")
      imp <- setNames(abs(as.vector(coefs[-1])), rownames(coefs)[-1])
    }
    if (method == "svm") {
      fit <- e1071::svm(x = as.matrix(X), y = y, kernel = "linear", scale = TRUE)
      w <- t(fit$coefs) %*% fit$SV
      imp <- setNames(abs(as.vector(w)), colnames(X))
    }
    if (method == "lgbm") {
      dtrain <- lightgbm::lgb.Dataset(data = as.matrix(X), label = y)
      fit <- lightgbm::lgb.train(
        params = list(objective = "regression",
                      learning_rate = 0.1,
                      num_leaves = 31,
                      max_depth = -1,
                      min_data_in_leaf = 5,
                      verbose = -1),
        data = dtrain,
        nrounds = 200
      )
      imp_tab <- lightgbm::lgb.importance(fit)
      imp <- setNames(imp_tab$Gain, imp_tab$Feature)
    }
    if (method == "mlp") {
      fit <- nnet::nnet(x = scale(as.matrix(X)), y = y,
                        size = 5, linout = TRUE, trace = FALSE, maxit = 200)
      imp <- sapply(all_pred, function(v) {
        x_perm <- as.matrix(X)
        x_perm[, v] <- sample(x_perm[, v])
        y_pred_perm <- predict(fit, x_perm)
        cor(y, y_pred_perm, use = "complete.obs")
      })
      imp <- 1 - abs(imp)
    }
    if (method == "gam") {
      df <- data.frame(y = y, X, check.names = FALSE)
      form <- as.formula(paste("y ~", paste(sprintf("s(%s)", all_pred), collapse = " + ")))
      fit <- mgcv::gam(form, data = df)
      s_table <- summary(fit)$s.table
      imp <- setNames(s_table[, "F"], rownames(s_table))
    }

    imp <- imp[all_pred]; imp[is.na(imp)] <- 0
    imp
  }

  # === 2) Parallel loop with progress bar ===
  results <- pbapply::pblapply(seq_len(n_taxa), function(i) {
    y_name <- taxa_names_vec[i]
    y <- otu[[i]]
    X <- otu[, -i, drop = FALSE]
    all_pred <- colnames(X)

    # raw importance
    imp_raw <- get_importance(y, X, all_pred, method, ntree)
    imp_use <- if (isTRUE(normalize)) .normalize_imp(imp_raw, scale.method) else imp_raw

    imp_df <- data.frame(
      target = y_name,
      predictor = all_pred,
      importance = as.numeric(imp_use),
      stringsAsFactors = FALSE
    )

    # permutation test
    if (!is.null(n_perm) && n_perm > 0) {
      perm_mat <- replicate(n_perm, {
        y_perm <- sample(y)
        imp_p <- get_importance(y_perm, X, all_pred, method, ntree)
        if (isTRUE(normalize)) imp_p <- .normalize_imp(imp_p, scale.method)
        imp_p
      })
      if (is.null(dim(perm_mat))) perm_mat <- matrix(perm_mat, nrow = length(all_pred))
      rownames(perm_mat) <- all_pred

      pvals <- vapply(all_pred, function(v) {
        real <- imp_use[v]
        nulls <- perm_mat[v, ]
        if (all(is.na(nulls)) || length(unique(nulls[!is.na(nulls)])) <= 1) {
          rank(-c(real, nulls))[1] / (length(nulls) + 1)
        } else {
          mean(nulls >= real, na.rm = TRUE)
        }
      }, numeric(1))

      imp_df$pval <- pvals
      sig_df <- dplyr::filter(imp_df, pval < p.threshold)
    } else {
      sig_df <- imp_df
    }

    list(edges = sig_df, imp = imp_df)
  }, cl = n_cores)

  # === 3) Summarize results ===
  edges <- dplyr::bind_rows(lapply(results, `[[`, "edges"))
  importance_all <- dplyr::bind_rows(lapply(results, `[[`, "imp"))

  list(
    edges = edges,
    importance_all = importance_all,
    method = method
  )
}

#' @title Normalize importance scores
#' @description
#' Internal helper to normalize variable importance scores.
#'
#' @param imp Numeric vector of importance values.
#' @param method Character. Normalization method: \code{"l1"}, \code{"softmax"},
#'   \code{"minmax"}, or \code{"none"}.
#'
#' @return A numeric vector of normalized importance values.
#'
#' @keywords internal
.normalize_imp <- function(imp, method = c("l1","softmax","minmax","none")) {
  method <- match.arg(method)
  if (method == "none") return(imp)

  if (method == "l1") {
    imp_pos <- pmax(imp, 0)                     # 先把负数截为 0，避免负贡献
    s <- sum(imp_pos, na.rm = TRUE)
    if (s > 0) {
      return(imp_pos / s)                       # 和为 1
    } else {
      # 如果全是 0 或 NA，则退化为 softmax，避免全 0
      w <- exp(imp - max(imp, na.rm = TRUE))
      return(w / sum(w, na.rm = TRUE))
    }
  }

  if (method == "softmax") {
    w <- exp(imp - max(imp, na.rm = TRUE))
    return(w / sum(w, na.rm = TRUE))
  }

  if (method == "minmax") {
    mn <- min(imp, na.rm = TRUE); mx <- max(imp, na.rm = TRUE)
    if (mx > mn) {
      return((imp - mn) / (mx - mn))
    } else {
      return(rep(0, length(imp)))
    }
  }
}

