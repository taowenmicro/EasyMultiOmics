% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EdgerSuper2.trans.R
\name{EdgerSuper2.trans}
\alias{EdgerSuper2.trans}
\title{Differential Abundance Analysis Using EdgeR for Groups}
\usage{
EdgerSuper2.trans(
  otu = NULL,
  tax = NULL,
  map = NULL,
  tree = NULL,
  ps = NULL,
  group = "Group",
  pvalue = 0.05,
  lfc = 0,
  artGroup = NULL,
  method = "TMM",
  j = 2
)
}
\arguments{
\item{otu}{A data frame containing OTU counts. Optional if `ps` is provided.}

\item{tax}{A data frame containing taxonomic annotations. Optional if `ps` is provided.}

\item{map}{A data frame containing sample metadata. Optional if `ps` is provided.}

\item{tree}{A phylogenetic tree. Optional if `ps` is provided.}

\item{ps}{A `phyloseq` object containing microbiome data. If provided, overrides `otu`, `tax`, `map`, and `tree`.}

\item{group}{A string specifying the grouping variable in the sample metadata. Default is `"Group"`.}

\item{pvalue}{A numeric value specifying the significance threshold for adjusted p-values (FDR). Default is `0.05`.}

\item{lfc}{A numeric value specifying the log2 fold change cutoff for significance. Default is `0`.}

\item{artGroup}{A custom matrix specifying pairwise group comparisons. Optional.}

\item{method}{A string specifying the normalization method for EdgeR. Default is `"TMM"`.}

\item{j}{A string or integer specifying the taxonomic rank to perform the analysis on.
Can be a numeric rank (1-7), a taxonomic name (e.g., `"Phylum"`), or `"OTU"`. Default is `2`.}
}
\value{
A data frame containing the differential abundance results:
\describe{
  \item{logFC}{The log2 fold change of taxa between the compared groups.}
  \item{p}{The raw p-values for each taxon.}
  \item{padj}{The adjusted p-values (FDR).}
  \item{level}{The classification of taxa as `"enriched"`, `"depleted"`, or `"nosig"` based on the logFC and p-value thresholds.}
  \item{group}{The name of the group comparison (e.g., `"Group1-Group2"`).}
}
}
\description{
The `EdgerSuper2.micro` function performs differential abundance analysis between microbial groups using the `EdgeR` package.
It calculates log fold changes, p-values, and adjusted p-values (FDR) to identify significantly enriched or depleted taxa
between group comparisons. The function works on data stored in a `phyloseq` object or manually provided OTU, taxonomic, and metadata tables.
}
\details{
The function performs the following steps:
\itemize{
  \item Prepares OTU, taxonomic, and metadata tables from the `phyloseq` object or individual inputs.
  \item Aggregates data to the specified taxonomic rank (if applicable).
  \item Normalizes OTU counts using the specified EdgeR normalization method (`"TMM"` by default).
  \item Constructs a model matrix and fits a generalized linear model (GLM) to the data.
  \item Performs pairwise group comparisons using contrasts in the EdgeR framework.
  \item Calculates log fold changes, raw p-values, and FDR-adjusted p-values for each taxon.
  \item Labels taxa as `"enriched"`, `"depleted"`, or `"nosig"` based on thresholds.
  \item Outputs a combined data frame with differential abundance results and group comparisons.
}

The results can be visualized using additional plotting tools, such as volcano plots, which highlight the significant taxa.
}
\examples{
\dontrun{
res =  EdgerSuper2.trans (ps = ps.trans,group  = "Group",artGroup =NULL, j = "gene")
head(res)

}

}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
