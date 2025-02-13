% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggVen.Upset.metm.R
\name{ggVen.Upset.metm}
\alias{ggVen.Upset.metm}
\title{Plotting Ven and Upset plot for each group}
\usage{
ggVen.Upset.metm(
  otu = NULL,
  tax = NULL,
  map = NULL,
  tree = NULL,
  ps = NULL,
  group = "Group",
  N = 0.5
)
}
\arguments{
\item{group}{column name for groupID.}

\item{otutab}{rarefied OTU table, typical output of usearch -otutab_norm or -otutab_rare,}

\item{metadata}{matrix or dataframe, including sampleID and groupID;}

\item{rep}{repeat number of each group}
}
\value{
base plot object.
}
\description{
Input otutab and metadata
VennDiagram::Ven plot
UpSetR::Upset plot
}
\examples{
res = ggVen.Upset.metm(ps = ps.16s,group = "Group")
grid::grid.draw(res[[1]])
dat = res[[2]]
}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
