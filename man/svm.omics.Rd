% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/svm.omics.R
\name{svm.omics}
\alias{svm.omics}
\title{Svm model screening of features in multi-omics data}
\usage{
\method{svm}{omics}(ps = ps, k = 5)
}
\arguments{
\item{ps}{A phyloseq format file used as an alternative for the input containing otu, tax, and map.}

\item{k}{The number of folds for cross-validation.}
}
\value{
A list object including the following components:
\item{AUC}{The average accuracy of the svm model.}
\item{Importance}{A data frame showing the feature importance ranked in descending order.}
}
\description{
This function uses a Support Vector Machine (SVM) to classify microbial community samples based on OTU (Operational Taxonomic Unit) abundances. It performs k-fold cross-validation to evaluate the classification accuracy and ranks feature importance using Recursive Feature Elimination (RFE).
}
\examples{
library(dplyr)
library(ggClusterNet)
library(caret)
library(e1071)
res <- svm.micro(ps = ps.16s \%>\% filter_OTU_ps(20), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
