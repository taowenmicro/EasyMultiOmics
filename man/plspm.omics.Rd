% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plspm.omics.R
\name{plspm.omics}
\alias{plspm.omics}
\title{Partial Least Squares Path Modeling (PLS-PM) for Omics Data}
\usage{
plspm.omics(
  ps.trans = ps.trans,
  ps.ms = ps.ms,
  ps.micro = ps.micro,
  ps.meta = ps.kegg,
  filter = 0.001
)
}
\value{
A list containing:
- `dat`: Summary results from the PLS-PM model.
- `cor`: Correlation matrix of the PLS-PM latent variables.
- `effect`: Estimated effects from the PLS-PM model.
}
\description{
This function performs Partial Least Squares Path Modeling (PLS-PM) on multiple omics datasets (e.g., transcriptomics, metabolomics, microbiome, and metagenome).
It first extracts the first principal component (PCA1) from each omics dataset, filters low-abundance taxa,
and then constructs a PLS-PM model to examine the relationships between different omics layers.
This function integrates four different omics datasets and analyzes their interdependencies using PLS-PM.
It follows these key steps:
1. Identifies common samples across all datasets.
2. Filters low-abundance taxa and extracts the first principal component (PCA1) for each dataset.
3. Constructs a PLS-PM model with a predefined inner and outer model.
4. Evaluates the model and returns key results, including correlation and effect size estimates.
}
\examples{
\dontrun{
library(vegan)
library(phyloseq)
library(dplyr)
library(ggClusterNet)
library(lavaan)
library(semPlot)
library(EasyMultiOmics)
data(ps.trans)
data(ps.ms)
data(ps.micro)
data(ps.kegg)
res <- plspm.omics(ps.trans=ps.trans,ps.ms=ps.ms,ps.micro=ps.micro,
ps.meta=ps.kegg,filter=0.05)
res
}
}
\author{
Tao Wen \email{2018203048@njau.edu.cn},
Peng-Hao Xie \email{2019103106@njau.edu.cn}.
}
