% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/buplot.micro.R
\name{buplot.micro}
\alias{buplot.micro}
\title{KEGG Enrichment Bubble Plot}
\usage{
buplot.micro(dt = kk@result, id = id)
}
\arguments{
\item{dt}{A data frame containing KEGG enrichment results, typically the output from `clusterProfiler::enrichKEGG`. It must include columns for `geneID`, `Description`, `pvalue`, and `p.adjust`.}

\item{id}{The column name in `dif` corresponding to the group or classification (e.g., "enriched", "depleted").}
}
\value{
A list containing two ggplot2 objects:
\itemize{
  \item Bubble plot with pathway descriptions labeled.
  \item Bubble plot without pathway descriptions labeled.
}
}
\description{
`buplot.micro` creates bubble plots for KEGG enrichment results. It visualizes the relationships between z-scores (effect size), p-values, and the pathway descriptions.
}
\details{
The function calculates a z-score for each pathway based on the proportion of "enriched" and "depleted" genes associated with the pathway, normalized by the pathway size. It then colors the bubbles by the combined z-score and p-value and optionally labels the pathways on the plot.

Steps:
\itemize{
  \item Extract `geneID` from the KEGG result and match it with the `dif` data frame.
  \item Calculate the z-score for each pathway as the difference in enriched and depleted counts normalized by the pathway size.
  \item Visualize the results using ggplot2.
}
}
\examples{
\dontrun{
result = buplot.micro(dt=dat1,id = id)
p1 = result[[1]]
p2 = result[[2]]
}

}
\author{
Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
