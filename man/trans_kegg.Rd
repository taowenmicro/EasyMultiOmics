% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_kegg.R
\name{trans_kegg}
\alias{trans_kegg}
\title{KEGG Pathway Transformation for Phyloseq Data}
\usage{
trans_kegg(ps.trans = ps.trans)
}
\arguments{
\item{ps.trans}{A Phyloseq object containing OTU table and associated taxonomic data.
The default is `ps.trans` which assumes that the Phyloseq object is provided externally.}
}
\value{
A transformed Phyloseq object with updated taxonomic information based on KEGG pathways.
The OTU table is aggregated by KEGG pathway IDs and the taxonomic table is updated to include pathway descriptions.
}
\description{
This function performs a transformation of the given Phyloseq object by mapping
the KEGG Orthology (KO) terms to KEGG pathways. It uses data from the KEGG
database and merges the KO terms with pathway descriptions. It processes the OTU
(Operational Taxonomic Unit) data and integrates it with KEGG pathway information,
returning a transformed Phyloseq object.
}
\examples{
\dontrun{
ps_transformed <- trans_kegg( ps.trans)
}
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
