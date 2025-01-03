% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomforest.micro.R
\name{MicroRoc}
\alias{MicroRoc}
\title{For cross-validation of microbiome data}
\usage{
MicroRoc(
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

\item{optimal}{important OTU number which selected}

\item{rfcv}{TURE or FELSE,whether need to do cross-validation}

\item{nrfcvnum}{Number of cross-validation}

\item{Group}{column name for groupID in map table.}
}
\value{
list contain ggplot object and table.

list contain ggplot object and table.
}
\description{
For cross-validation of microbiome data

Comparison of three machine methods (randomforest,SVM,GLM).
}
\examples{
# data form github
result = Micro.rfcv(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = ps_rela,group  = "Group",optimal = 20,nrfcvnum = 6)
prfcv = result[[1]]# plot rfcv
rfcvtable = result[[3]]# table rfcv
result = MicroRoc( ps = ps,group  = "Group")
#--提取roc曲线
result[[1]]
#提取AUC值
result[[2]]
}
\references{
Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}

Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Yong-Xin Liu \email{yxliu@genetics.ac.cn}
}
