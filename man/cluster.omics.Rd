% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.omics.R
\name{cluster.omics}
\alias{cluster.omics}
\title{Cluster and Visualize Communities of multi-omics}
\usage{
\method{cluster}{omics}(
  ps = ps,
  hcluter_method = "complete",
  dist = "bray",
  cuttree = gnum,
  row_cluster = T,
  col_cluster = T
)
}
\arguments{
\item{ps}{A `phyloseq` object containing microbiome data.}

\item{hcluter_method}{A string specifying the hierarchical clustering method.
Default is `"complete"`. Options include `"single"`, `"average"`, and `"ward.D"`.}

\item{dist}{A string specifying the distance method for clustering.
Default is `"bray"`. Other options include methods available in `phyloseq::distance`.}

\item{cuttree}{An integer specifying the number of clusters for cutting the dendrogram.}

\item{row_cluster}{Logical. Whether to perform clustering on rows. Default is `TRUE`.}

\item{col_cluster}{Logical. Whether to perform clustering on columns. Default is `TRUE`.}
}
\value{
A list containing the following elements:
\item{p0}{A ggtree plot of the clustering dendrogram.}
\item{p1}{A heatmap showing relative abundances with optional row and column clustering.}
\item{p2}{A heatmap with sample IDs on the y-axis, showing relative abundances.}
\item{tem}{A data frame representing the distance matrix used for clustering.}
}
\description{
This function performs hierarchical clustering and visualizes microbial community composition
using phyloseq data. It provides various visualizations, including clustering dendrograms
and heatmaps, to explore relationships among samples.
}
\examples{
res = cluster.omics (ps= ps03,hcluter_method = "complete",dist = "bray",cuttree = 3,row_cluster = TRUE,col_cluster =  TRUE)
p4 = res[[1]]
p4
p4_1 = res[[2]]
p4_1
p4_2 = res[[3]]
p4_2
dat = res[[4]]# cluster distance
head(dat)
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
