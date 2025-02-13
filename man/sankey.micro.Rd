% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sankey.micro.R
\name{sankey.micro}
\alias{sankey.micro}
\title{Generate Taxonomic Sankey Diagram for Microbial Data}
\usage{
sankey.micro(ps = ps, rank = 6, Top = 50)
}
\arguments{
\item{ps}{A `phyloseq` object containing microbiome data.}

\item{rank}{A numeric value specifying the taxonomic rank for aggregation. Default is `6` (e.g., Genus level).}

\item{Top}{An integer specifying the number of top OTUs to include based on abundance. Default is `50`.}
}
\value{
A list containing:
\describe{
  \item{An interactive Sankey diagram generated with `networkD3::sankeyNetwork`.}
  \item{A data frame containing the Sankey diagram links and node data.}
}
}
\description{
This function creates Sankey diagrams to visualize microbial taxonomic transitions from higher to lower taxonomic levels (e.g., Kingdom to Genus).
It generates interactive Sankey diagrams for each group in the `phyloseq` object based on the specified taxonomic rank and top OTUs.
}
\details{
The function performs the following steps:
\itemize{
  \item Extracts OTU and taxonomic data for each group in the `phyloseq` object.
  \item Aggregates taxonomic information up to the specified rank and filters the top `Top` OTUs by abundance.
  \item Creates source-target relationships for taxonomic transitions between levels (e.g., Kingdom → Phylum, Phylum → Class).
  \item Calculates mean abundance values for each taxonomic level within groups.
  \item Generates an interactive Sankey diagram visualizing taxonomic transitions and their relative abundances.
}

The Sankey diagram is interactive and allows users to explore taxonomic transitions and their abundances between levels.
}
\examples{
\dontrun{
res = sankey.micro (ps = ps.16s,rank = 6,Top = 50)
p22 = res[[1]]
dat = res[[2]]
}

}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
