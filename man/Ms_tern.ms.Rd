% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Ms_tern.ms.R
\name{Ms_tern.ms}
\alias{Ms_tern.ms}
\title{Generate Ternary Plots for metabolites composition data}
\usage{
Ms_tern.ms(ps = ps, color = color)
}
\arguments{
\item{ps}{A `phyloseq` object containing metabolites composition data.}
}
\value{
A list containing the following elements:
\describe{
  \item{plot}{A list of ternary plots (`ggtern` objects) for different group combinations.}
  \item{dataplot}{A data frame with combined metabolites abundance and taxonomic information.}
  \item{groups}{A matrix showing all possible combinations of three groups used for ternary plots.}
}
}
\description{
This function generates ternary plots for visualizing metabolites composition data in three groups.
It uses relative abundance data from a `phyloseq` object and displays metabolites distribution across specified groupings.
}
\details{
This function performs the following steps:
\itemize{
  \item Standardizes the metabolites composition table to relative abundances using `phyloseq::transform_sample_counts`.
  \item Splits the OTU data into groups based on the metadata column `Group`.
  \item Computes the mean abundance for each metabolite within each group.
  \item Generates ternary plots for all possible combinations of three groups using the `ggtern` package.
}

The plots visualize the relative abundances of metabolites in three groups, with points sized by their mean abundance and colored by their taxonomic classification (mode).
}
\examples{
\dontrun{
ps1 = ps.ms \%>\% filter_OTU_ps(500)
res = Micro_tern.ms(ps1)
p15 = res[[1]]
dat =  res[[2]]
}

}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
