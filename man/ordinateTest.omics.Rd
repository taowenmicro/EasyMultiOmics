% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ordinateTest.omics.R
\name{ordinateTest.omics}
\alias{ordinateTest.omics}
\title{Beta diversity statistics by adonis/anosim/MRPP in all groups}
\usage{
ordinateTest.omics(
  ps = NULL,
  group = "Group",
  Micromet = "MRPP",
  dist = "bray"
)
}
\arguments{
\item{ps}{alternative input;}

\item{group}{group ID;}

\item{Micromet}{statistics default by adonis, alternative anosim or MRPP;}

\item{dist}{distance type, including "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski"  "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial"  "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co";}

\item{method}{DCA, CCA, RDA, NMDS, MDS, PCoA, PCA, LDA;}

\item{pvalue.cutoff}{Pvalue threshold, default in 0.05;}
}
\value{
stat table
}
\description{
Input phyloseq object, test method and distance type
}
\details{
By default, input phyloseq object include metadata and otutab
The available diversity indices include the following:
\itemize{
\item{most used indices: bray unifrac wunifrac}
\item{other used indices: dpcoa jsd manhattan euclidean canberra kulczynski jaccard gower altGower morisita horn mountford raup binomial chao cao w -1 c wb r I e t me j sor m -2 co}
}
}
\examples{
dat1 = ordinateTest.omics(ps = ps03, Micromet = "adonis", dist = "bray")
dat1
}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
