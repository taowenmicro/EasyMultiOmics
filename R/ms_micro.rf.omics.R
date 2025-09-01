#' @title Random Forest Analysis for Microbiome Data Using NMDS
#' @description
#' The function first converts the OTU (Operational Taxonomic Unit) table from a phyloseq object into a
#' numerical matrix, computes a Bray-Curtis distance matrix, performs NMDS ordination, and then applies
#' Random Forest to determine the importance of environmental variables on the first NMDS axis.
#' A bar plot is generated to visualize the feature importance based on Random Forest results.
#'
#' @param ps A phyloseq object containing the microbiome data. This object should include OTU data and sample metadata.
#' @param env A data frame or vector containing the environmental variables to be used in the analysis.
#' @param seed An integer for setting the random seed to ensure reproducibility (default is 1).
#'
#' @return
#' A list containing two elements:
#' \item{df}{A data frame containing the importance of each feature (environmental variable) based on Random Forest analysis.}
#' \item{p}{A ggplot object representing a bar plot of feature importance.}
#' @examples
#' \dontrun{
#' ps.tem = ps.micro %>% filter_OTU_ps(500)
#' res = loadingPCA.ms(ps = ps.ms)
#' dat = res[[2]]
#' dat$id = row.names(dat)
#' id = dat %>% arrange(desc(PCone)) %>% head(15) %>% .$id
#' ftab = ps.ms %>% scale_micro() %>% subset_taxa(row.names(tax_table(ps.ms)) %in% c(id) ) %>% vegan_otu() %>%as.data.frame()
#' result <- ms_micro.rf.omics(ps = ps, env = env, seed = 42)
#' }
#' @export
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}

ms_micro.rf.omics <- function(
    ps = ps,
    env = env,
    seed = 1
){
  otu<- as.data.frame(t(ggClusterNet::vegan_otu(ps)))
  map <- design <- phyloseq::sample_data(ps)
  aa<-lapply(otu, function(x) as.numeric(as.character(x)))
  aa<-as.data.frame(aa)
  otu<-aa
  dist= vegan::vegdist(t(otu), method="bray")
  nmds= vegan::monoMDS(dist)
  index = as.data.frame(nmds$points)
  index<-cbind(design,index)
  NMDS1<- cbind(env,MDS1 = index$MDS1)
  set.seed(seed)
  bb= randomForest::randomForest(NMDS1$MDS1 ~., data=NMDS1,importance = TRUE, ntree = 500, nrep = 100, na.action=na.omit)
  bb
  df<-round(randomForest::importance(bb), 2)
  df<-as.data.frame(df)
  df<-cbind(df,rownames(df))
  p <- ggplot(df, aes(x =`%IncMSE` , y =reorder(`rownames(df)`,`%IncMSE`) )) +
    geom_bar(stat = "identity", width = 0.75,position = "dodge",colour="black",fill="#9ACD32",alpha=1) +
    labs(y="", x="", title = "",size=9)+
    theme_bw() +
    theme(axis.text=element_text(colour='black',size=9))
  return(list(df,p))

}
