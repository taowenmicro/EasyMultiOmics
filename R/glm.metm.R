#' @title Glm model screening of characteristic microorganisms
#' @description
#' Glm, one of the machine learning methods, was used to screen for characteristic
#' microorganisms, and the model was evaluated using k-fold cross-validation.
#' @param ps A phyloseq format file used as an alternative for the input containing otu, tax, and map.
#' @param k The number of folds for cross-validation.
#' @return A list object including the following components:
#' \item{AUC}{The average accuracy of the glm model.}
#' \item{Importance}{A data frame showing the feature importance ranked in descending order.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' library(ggClusterNet)
#' library(caret)
#' pst=subset_samples(ps.16s,Group %in% c("KO" ,"OE"))
#' res <- glm.metm(ps = pst%>% filter_OTU_ps(50), k = 5)
#' AUC = res[[1]]
#' AUC
#' importance = res[[2]]
#' importance

glm.metm =function (ps = ps, k = 5)
{
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
    fold_test <- train[folds[[i]], ]
    fold_train <- train[-folds[[i]], ]
    model <- glm(OTUgroup ~ ., family = "binomial", data = fold_train)
    model
    model_pre <- predict(model, type = "response", newdata = fold_test)
    model_pre
    fc <- append(fc, fold_test$OTUgroup)
    mod_pre <- append(mod_pre, as.numeric(model_pre))
  }
  pred <- prediction(mod_pre, fc)
  auc <- performance(pred, "auc")@y.values[[1]]
  auc_result <- paste("GLM AUC:", round(auc, 3))
  summary_model <- summary(model)
  coefficients <- summary_model$coefficients
  importance <- data.frame(Feature = rownames(coefficients),
                           Coefficient = coefficients[, 1], Std_Error = coefficients[,
                                                                                     2], z_value = coefficients[, 3], p_value = coefficients[,
                                                                                                                                             4])
  importance <- importance[order(abs(importance$Coefficient),
                                 decreasing = TRUE), ]
  row.names(importance) = gsub("OTU", "", row.names(importance))
  importance$Feature = gsub("OTU", "", importance$Feature)
  list(AUC = auc_result, Importance = importance)
}
