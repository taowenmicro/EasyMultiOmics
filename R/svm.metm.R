#' Support Vector Machine model for screening characteristic microorganisms
#'
#' This function fits a Support Vector Machine (SVM) classification model
#' on a phyloseq object using k-fold cross-validation. The OTU table is used
#' as predictors and the sample grouping variable \code{Group} is used as
#' the response. The function returns overall accuracy and feature importance.
#'
#' @param ps A \code{phyloseq} object containing an OTU table and sample data.
#'   The sample_data slot must contain a column named \code{Group} indicating
#'   the class labels.
#' @param k Integer; the number of folds for cross-validation. Default is 5.
#' @param seed Integer; random seed for reproducibility. Default is 123.
#' @param top Optional integer. If not \code{NULL}, the OTU table is first
#'   filtered to keep only the top abundant OTUs using \code{filter_OTU_ps(top)}.
#'
#' @return A list with the following elements:
#' \describe{
#'   \item{AUC}{A character string summarising the overall accuracy
#'   (kept for backward compatibility with previous naming).}
#'   \item{Importance}{A \code{varImp.train} object containing feature importance
#'   scores as computed by \code{caret::varImp}.}
#'   \item{Model}{The fitted \code{caret} \code{train} object using the
#'   \code{"svmRadial"} method.}
#' }
#'
#' @details
#' This function relies on the following packages:
#' \itemize{
#'   \item \pkg{phyloseq} for handling microbiome data structures.
#'   \item \pkg{ggClusterNet} for extracting OTU tables via \code{vegan_otu}.
#'   \item \pkg{caret} for cross-validation and model training.
#'   \item \pkg{kernlab} as the backend for the \code{"svmRadial"} method.
#' }
#' The function assumes that the grouping variable is stored in
#' \code{sample_data(ps)$Group}. Column names of the OTU predictors are
#' sanitized using \code{make.names} to avoid issues in model formulas.
#'
#' @examples
#' \dontrun{
#'   library(phyloseq)
#'   library(ggClusterNet)
#'   library(caret)
#'   library(kernlab)
#'
#'   res_svm <- svm_metm(ps = ps_example, k = 5, top = 100)
#'   res_svm$AUC
#'
#'   # Extract importance as a data frame
#'   imp_df <- res_svm$Importance$importance
#'   imp_df$Feature <- rownames(imp_df)
#'   rownames(imp_df) <- NULL
#'   head(imp_df)
#' }
#'
#' @export
svm_metm <- function(ps,
                     k    = 5,
                     seed = 123,
                     top  = NULL) {

  ## ---- basic checks ---------------------------------------------------------
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be a phyloseq object.", call. = FALSE)
  }

  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("ggClusterNet", quietly = TRUE)) {
    stop("Package 'ggClusterNet' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("caret", quietly = TRUE)) {
    stop("Package 'caret' is required but not installed.", call. = FALSE)
  }
  if (!requireNamespace("kernlab", quietly = TRUE)) {
    stop("Package 'kernlab' is required but not installed.", call. = FALSE)
  }

  set.seed(seed)

  ## ---- optional OTU filtering -----------------------------------------------
  if (!is.null(top)) {
    if (!is.numeric(top) || length(top) != 1L) {
      stop("`top` must be a single numeric value.", call. = FALSE)
    }
    # Assuming filter_OTU_ps(ps, top) is available in the same package
    ps <- filter_OTU_ps(ps, top = top)
  }

  ## ---- extract sample data and OTU table ------------------------------------
  map <- as.data.frame(phyloseq::sample_data(ps))

  if (!"Group" %in% colnames(map)) {
    stop("sample_data(ps) must contain a column named 'Group'.", call. = FALSE)
  }

  otu_mat <- ggClusterNet::vegan_otu(ps)   # OTU x sample
  otu_df  <- as.data.frame(t(otu_mat))     # sample x OTU

  # Build combined data frame: response + predictors
  df <- cbind(Group = map$Group, otu_df)
  df$Group <- as.factor(df$Group)

  # Sanitize predictor names
  colnames(df)[-1] <- make.names(colnames(df)[-1])

  ## ---- caret training control -----------------------------------------------
  n_class <- nlevels(df$Group)

  if (n_class == 2L) {
    # Binary classification: use ROC as metric
    train_control <- caret::trainControl(
      method          = "cv",
      number          = k,
      classProbs      = TRUE,
      summaryFunction = caret::twoClassSummary,
      savePredictions = "final",
      allowParallel   = FALSE
    )
    metric <- "ROC"
  } else {
    # Multiclass: use Accuracy
    train_control <- caret::trainControl(
      method          = "cv",
      number          = k,
      classProbs      = FALSE,
      savePredictions = "final",
      allowParallel   = FALSE
    )
    metric <- "Accuracy"
  }

  ## ---- train SVM model (svmRadial) ------------------------------------------
  model <- caret::train(
    Group ~ .,
    data      = df,
    method    = "svmRadial",
    trControl = train_control,
    metric    = metric
  )

  ## ---- overall accuracy -----------------------------------------------------
  pred <- stats::predict(model, newdata = df)
  cm   <- caret::confusionMatrix(pred, df$Group)
  acc  <- cm$overall["Accuracy"]
  acc_txt <- paste("SVM Average Accuracy:", round(acc, 3))

  ## ---- feature importance ----------------------------------------------------
  imp <- caret::varImp(model, scale = FALSE)

  ## ---- return ---------------------------------------------------------------
  list(
    AUC        = acc_txt,  # for backward compatibility with previous naming
    Importance = imp,
    Model      = model
  )
}
