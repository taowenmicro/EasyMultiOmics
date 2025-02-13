% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/neutralModel.R
\name{neutralModel}
\alias{neutralModel}
\title{Neutral Model Fit for Microbial Communities}
\usage{
neutralModel(
  otu = NULL,
  tax = NULL,
  map = NULL,
  tree = NULL,
  ps = NULL,
  group = "Group",
  ncol = 3,
  nrow = 1
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

\item{ncol}{An integer specifying the number of columns for the arranged output plots. Default is 3.}

\item{nrow}{An integer specifying the number of rows for the arranged output plots. Default is 1.}
}
\value{
A list containing:
\itemize{
  \item  A combined plot showing neutral model fits for all groups.
  \item  A list of individual ggplot objects for each group.
  \item  A list of data frames containing observed occupancy and abundance for each group.
  \item  A list of data frames containing predicted frequencies and confidence intervals for each group.
}
}
\description{
This function implements a neutral community model to assess the relationship
between microbial taxa abundance in the metacommunity and their occurrence frequency
across samples. It fits the model for each group in the dataset and visualizes the results
using scatter plots with confidence intervals.
}
\examples{
\dontrun{
psphy = filter_taxa(ps16s, function(x) sum(x ) > 100, TRUE);psphy
map = sample_data(psphy)
n = map$Group \%>\% unique() \%>\%length()
result = neutralModel(ps = psphy,group  = "Group",ncol = n)

}

}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
