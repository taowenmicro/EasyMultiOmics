% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluMicro.bar.trans.R
\name{cluMicro.bar.trans}
\alias{cluMicro.bar.trans}
\title{Plotting alpha diversity barplot for each group with anova statistics}
\usage{
cluMicro.bar.trans(
  dist = "bray",
  otu = NULL,
  tax = NULL,
  map = NULL,
  tree = NULL,
  j = "pathway",
  ps = ps,
  rep = 6,
  Top = 10,
  tran = TRUE,
  hcluter_method = "complete",
  Group = "Group",
  cuttree = 3
)
}
\arguments{
\item{otu}{OTU/ASV table;}

\item{tax}{taxonomy table}

\item{map}{Sample metadata;}

\item{rep}{Number of sample replicates measured in each group}

\item{Top}{Select the number of microorganisms with the highest relative abundance in all samples}

\item{hcluter_method}{hcluter method}

\item{Group}{column name for groupID in map table.}

\item{cuttree}{cut number}
}
\value{
list contain ggplot object and table.
}
\description{
Combine sample clustering and drawing stacked histograms of microbial categories
}
\details{
hcluter method is same an function hcluster
}
\examples{
# data form github
metadata = read.table("http://210.75.224.110/github/EasyAmplicon/data/metadata.tsv", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
otutab = read.table("http://210.75.224.110/github/EasyAmplicon/data/otutab.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
taxonomy = read.table("http://210.75.224.110/github/EasyAmplicon/data/taxonomy.txt", header=T, row.names=1, sep="\t", comment.char="", stringsAsFactors = F)
result <-  cluMicro.bar (dist = "bray",
                        otu = otutab,
                        tax = taxonomy,
                        map = metadata,
                        rep = 6 ,# 重复数量是6个
                        Top = 10, # 提取丰度前十的物种注释
                        tran = TRUE, # 转化为相对丰度值
                        hcluter_method = "complete",
                        cuttree = 3,
                        Group = "Group"
)
}
\references{
Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Yong-Xin Liu \email{yxliu@genetics.ac.cn}
}
