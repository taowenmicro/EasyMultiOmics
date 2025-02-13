% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/glm.metm.R
\name{glm.metm}
\alias{glm.metm}
\title{Glm model screening of characteristic microorganisms}
\usage{
glm.metm(ps = ps, k = 5)
}
\arguments{
\item{ps}{A phyloseq format file used as an alternative for the input containing otu, tax, and map.}

\item{k}{The number of folds for cross-validation.}
}
\value{
A list object including the following components:
\item{AUC}{The average accuracy of the glm model.}
\item{Importance}{A data frame showing the feature importance ranked in descending order.}
}
\description{
Glm, one of the machine learning methods, was used to screen for characteristic
microorganisms, and the model was evaluated using k-fold cross-validation.
}
\examples{
library(dplyr)
library(ggClusterNet)
library(caret)
pst=subset_samples(ps.16s,Group \%in\% c("KO" ,"OE"))
res <- glm.metm(ps = pst\%>\% filter_OTU_ps(50), k = 5)
AUC = res[[1]]
AUC
importance = res[[2]]
importance
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
