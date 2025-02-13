% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Roc.trans.R
\name{Roc.trans}
\alias{Roc.trans}
\title{Comparison of three machine methods (randomforest,SVM,GLM).}
\usage{
Roc.trans(
  otu = NULL,
  tax = NULL,
  map = NULL,
  tree = NULL,
  ps = NULL,
  group = "Group",
  repnum = 5
)
}
\arguments{
\item{otu}{OTU/ASV table;}

\item{tax}{taxonomy table}

\item{map}{Sample metadata;}

\item{ps}{phyloseq object of microbiome}

\item{repnum}{Modeling times}

\item{Group}{column name for groupID in map table.}
}
\value{
list contain ggplot object and table.
}
\description{
Comparison of three machine methods (randomforest,SVM,GLM).
}
\examples{
result = MicroRoc( ps = ps,group  = "Group")
#--提取roc曲线
result[[1]]
#提取AUC值
result[[2]]
}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Yong-Xin Liu \email{yxliu@genetics.ac.cn}
}
