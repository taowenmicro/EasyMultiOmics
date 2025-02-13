% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nullModel.R
\name{nullModel}
\alias{nullModel}
\title{Null Model Analysis for Microbial Communities}
\usage{
nullModel(
  otu = NULL,
  tax = NULL,
  map = NULL,
  tree = NULL,
  ps = NULL,
  group = "Group",
  dist.method = "bray",
  gamma.method = "total",
  transfer = "none",
  null.model = "ecosphere"
)
}
\arguments{
\item{otu}{A data frame or matrix containing OTU (Operational Taxonomic Unit) abundance data.
Rows represent taxa, and columns represent samples.}

\item{tax}{A data frame containing taxonomic annotation for each OTU.}

\item{map}{A data frame containing sample metadata, including group information.}

\item{tree}{A phylogenetic tree object (optional).}

\item{ps}{A `phyloseq` object containing OTU, taxonomic, and sample data.
If provided, this supersedes `otu`, `tax`, `map`, and `tree`.}

\item{group}{A character string specifying the grouping variable in the metadata.
Default is `"Group"`.}

\item{dist.method}{A character string specifying the distance method to calculate beta diversity.
Default is `"bray"`.}

\item{gamma.method}{A character string specifying the gamma diversity method.
Options are `"total"` or `"group"`. Default is `"total"`.}

\item{transfer}{Transformation method for the data (not currently implemented). Default is `"none"`.}

\item{null.model}{A character string specifying the null model to use.
Options are `"ecosphere"`, `"ecosim"`, or `"frequency"`. Default is `"ecosphere"`.}
}
\value{
A list containing:
\itemize{
  \item A data frame summarizing the observed and permuted similarities,
  including effect sizes and standardized effect sizes for each group.
  \item  A data frame containing the ratio of null model similarity to observed similarity for each group.
  \item `Results of an ANOVA test on the ratios between groups (if applicable).
}
}
\description{
This function performs null model analysis to compare the observed beta diversity
similarity to the expected similarity under various null model assumptions.
It provides effect sizes, standardized effect sizes (SES), and ratios for different
groups in the microbial community dataset.
}
\details{
The function performs the following steps:
\itemize{
  \item Splits the data by groups and calculates observed beta diversity similarity.
  \item Applies the specified null model (e.g., `ecosphere`, `ecosim`) to generate permuted similarity values.
  \item Compares observed and permuted similarity using t-tests, effect sizes, and standardized effect sizes.
  \item Calculates the ratio of permuted to observed similarity.
  \item If multiple groups are present, performs ANOVA on the similarity ratios between groups.
}
}
\examples{
\dontrun{
result <- nullModel(ps = psphy,group="Group",dist.method =  "bray",gamma.method = "total",transfer = "none",null.model = "ecosphere")
}

}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
