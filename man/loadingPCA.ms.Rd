% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loadingPCA.ms.R
\name{loadingPCA.ms}
\alias{loadingPCA.ms}
\title{PCA loading Matrix screening of characteristic metabolites}
\usage{
loadingPCA.ms(ps = ps, Top = 20)
}
\arguments{
\item{ps}{A phyloseq format file used as an alternative for the input containing metabolite composition table,
metabolite classification table, and sample metadata.}

\item{Top}{The top microorganisms to visualize.}
}
\value{
A list object containing the following components:
\item{p}{A PCA correlation plot of the selected number of characteristic metabolites and
 the correlation decreases from top to bottom.}
\item{index}{Data frame containing the PCA load matrix and relative abundance of all metabolites.}
}
\description{
This function conducts PCA analysis on metabolite data to
extract the loading Matrix, screen for characteristic metabolites and visualize them.
The importance of the variables is  sorted according to
the square value of the correlation between the variable and the PC1 axis.
}
\examples{
res = loadingPCA.ms(ps = ps.ms,Top = 20)
p = res[[1]]
p
dat = res[[2]]
dat
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
