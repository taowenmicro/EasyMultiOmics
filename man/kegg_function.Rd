% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kegg_function.R
\name{kegg_function}
\alias{kegg_function}
\title{KEGG Pathway Analysis for Phyloseq Object}
\usage{
kegg_function(ps = ps)
}
\arguments{
\item{ps}{A `phyloseq` object containing OTU data to be transformed into KEGG pathways.}
}
\value{
A `phyloseq` object with transformed OTU table and KEGG pathway information.
}
\description{
This function downloads KEGG pathway information, transforms an OTU table in a
`phyloseq` object to KO identifiers, merges the KO identifiers with KEGG pathways,
and returns a `phyloseq` object with KEGG pathway data.
}
\details{
This function performs several steps:
\itemize{
  \item Downloads KEGG pathway data using the `clusterProfiler::download_KEGG` function.
  \item Transforms the OTU table in `ps` to KO identifiers using `vegan_otu()`.
  \item Merges KO identifiers with KEGG pathway names and IDs.
  \item Summarizes pathway data and creates a new `phyloseq` object containing summarized OTU data.
}
}
\examples{
# Assuming you have a phyloseq object `ps`
kegg_ps <- kegg_function(ps.kegg)
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
