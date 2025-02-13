% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tax_glom_meta.R
\name{tax_glom_meta}
\alias{tax_glom_meta}
\title{Summarize Microbiome Data by Taxonomic Rank}
\usage{
tax_glom_meta(ps = ps, ranks = "Phylum")
}
\arguments{
\item{ps}{A `phyloseq` object containing OTU and taxonomic data.}

\item{ranks}{A character string or numeric value specifying the taxonomic rank
by which to group the OTU data. Default is "Phylum". If a numeric value is provided,
it will select the rank from the list of ranks in the `phyloseq` object.}
}
\value{
A `phyloseq` object containing the aggregated OTU data based on the specified taxonomic rank.
}
\description{
This function aggregates the OTU data in a `phyloseq` object by a specified taxonomic rank
(default is "Phylum") and returns a new `phyloseq` object with the summarized data.
}
\details{
The function handles missing taxonomic values (e.g., NA, empty strings, or "NA") and replaces them with "Unknown".
It also removes duplicate taxonomic names within the specified rank. The summarized OTU data is returned as a new
`phyloseq` object.
}
\examples{
# Example of usage with a phyloseq object `ps`
ps_summary <- tax_glom_meta(ps =ps.card, ranks ="Drug_Class")
head(ps_summary)

}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
