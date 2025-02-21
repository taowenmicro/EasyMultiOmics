% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pairtranstest.trans.R
\name{pairtranstest.trans}
\alias{pairtranstest.trans}
\title{Beta diversity statistics by adonis/anosim/MRPP in pair}
\usage{
pairtranstest.trans(ps = ps, Micromet = "adonis", dist = "bray")
}
\arguments{
\item{ps}{alternative input;}

\item{Micromet}{statistics default by adonis, alternative anosim or MRPP;}

\item{dist}{distance type, including "unifrac" "wunifrac" "dpcoa" "jsd" "manhattan" "euclidean"   "canberra" "bray" "kulczynski"  "jaccard" "gower" "altGower" "morisita" "horn" "mountford"  "raup" "binomial"  "chao"  "cao" "w"  "-1"  "c" "wb"  "r"   "I"  "e" "t" "me"   "j"  "sor"  "m"   "-2"  "co";}

\item{otu}{OTU/ASV table;}

\item{map}{Sample metadata;}

\item{tree}{tree/nwk file;}

\item{group}{group ID;}

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
data(ps.trans)
ps = ps.trans\%>\% filter_taxa(function(x) sum(x ) > 5 , TRUE)
otutab_rare=data.frame(phyloseq::otu_table(ps))
metadata=data.frame(phyloseq::sample_data(ps))
# Input otutab, metadata, and options group, method and distance
pairtranstest.trans(otu = otutab_rare, map = metadata, group = "Group", Micromet = "MRPP", dist = "bray")
# Input phyloseq format input, and options group, method and distance
pairtranstest.trans(ps = ps, group = "Group", Micromet = "MRPP", dist = "bray")
}
\author{
Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
}
