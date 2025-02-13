% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/trans_rekegg.R
\name{trans_rekegg}
\alias{trans_rekegg}
\title{KEGG Reaction Transformation for Phyloseq Data}
\usage{
trans_rekegg(ps = ps.trans)
}
\arguments{
\item{ps}{A Phyloseq object containing OTU table and associated taxonomic data.
The default is `ps.trans`, which assumes that the Phyloseq object is provided externally.}
}
\value{
A transformed Phyloseq object with the OTU table aggregated by KEGG reaction
IDs and the taxonomic table updated with KEGG reaction descriptions.
}
\description{
This function transforms a Phyloseq object by mapping KEGG Orthology (KO) terms
to KEGG reactions. It uses data from the EasyMultiOmics database and merges KO
terms with reaction data. The function processes the OTU (Operational Taxonomic Unit)
table from a Phyloseq object, aggregates the data by KEGG reaction, and returns
a transformed Phyloseq object with KEGG reaction descriptions.
}
\examples{
\dontrun{
ps_transformed <- trans_rekegg(ps.trans)
}
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
