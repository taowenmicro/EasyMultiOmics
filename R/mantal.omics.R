#' @title Mantel Test for Multi-Omics Data
#'
#' @description
#' The `mantal.omics` function performs Mantel tests between two distance matrices derived from multi-omics data (e.g., microbial and metabolomics datasets). It also computes Procrustes analysis results to assess the similarity of ordination configurations between the datasets.
#'
#' @param ps01 A `phyloseq` object containing the first omics dataset (e.g., microbial data).
#' @param ps02 A `phyloseq` object containing the second omics dataset (e.g., metabolomics data).
#' @param method A character string specifying the correlation method for the Mantel test. Options are `"spearman"` (default) or `"pearson"`.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes Bray-Curtis distance matrices for the two input datasets (`ps01` and `ps02`) using the `vegan::vegdist` function.
#'   \item Performs Mantel tests to evaluate the correlation between the two distance matrices.
#'   \item Computes ordination results for both datasets using non-metric multidimensional scaling (NMDS) and performs Procrustes analysis to assess similarity.
#'   \item Returns a summary table with Mantel and Procrustes test results for all pairwise group comparisons.
#' }
#'
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' \dontrun{
#' ps_microbial <- ps01  # Replace with your microbial phyloseq object
#' ps_metabolomics <- ps02  # Replace with your metabolomics phyloseq object
#' # Perform Mantel and Procrustes tests
#' mantel_results <- mantal.omics(ps01 = ps_microbial, ps02 = ps_metabolomics, method = "spearman")
#' print(mantel_results)
#' }
#' @export
mantal.omics= function(ps01= ps01,ps02= ps02,method="spearman" ){
dist.01 = ps01 %>%
  scale_micro() %>%
  # tax_glom_wt(j) %>%
  # subset_taxa.wt("OTU",id.micro) %>%
  otu_table() %>% t() %>%
  vegan::vegdist(method="bray") %>%
  as.matrix()


dist.02 = ps02 %>%
  scale_micro() %>%
  # subset_taxa.wt("OTU",id.ms) %>%
  otu_table() %>%
  t() %>%
  vegan::vegdist(method="bray") %>%
  as.matrix()

# method =  "spearman"
# method =  "pearson"

dist.list = list(Bac =dist.01,
                 Fun = dist.02
)

sample_data(ps01)

gru = c("da1","da2")
id = combn(unique(gru),2)
names(dist.list) = gru


R_mantel = c()
p_mantel = c()
name = c()
R_pro <- c()
p_pro <- c()
plots = list()
i = 1
for (i in 1:dim(id)[2]) {


  dist1 = dist.list[[id[1,i]]]
  dist2 <- dist.list[[id[2,i]]]
  mt <- vegan::mantel(dist1,dist2,method = method)
  R_mantel[i] = mt$statistic
  p_mantel[i] = mt$signif

  name[i] = paste(id[1,i],"_VS_",id[2,i],sep = "")
  #--p

  mds.s <- vegan::monoMDS(dist1)
  mds.r <- vegan::monoMDS(dist2)
  pro.s.r <- vegan::protest(mds.s,mds.r)

  R_pro[i] <- pro.s.r$ss
  p_pro[i] <- pro.s.r$signif

}
dat = data.frame(name,R_mantel,p_mantel,R_pro,p_pro )
head(dat)
return(dat)
}
