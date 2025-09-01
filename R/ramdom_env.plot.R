#' @title Random Forest Analysis of Environmental Variables on NMDS Axis
#'
#' @description
#' The `ramdom_env.plot` function uses Random Forest analysis to identify the importance of environmental variables in explaining the variation along the first NMDS axis (MDS1) derived from Bray-Curtis distances of OTU data.
#'
#' @param ps A `phyloseq` object containing OTU/taxa abundance data and sample metadata.
#' @param env A data frame of environmental variables, where rows are samples and columns are variables.
#' @param seed An integer value to set the random seed for reproducibility. Default is `1`.
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `df`: A data frame summarizing the importance of each environmental variable, including `%IncMSE` (percent increase in mean squared error) and `IncNodePurity` (increase in node purity).
#'   \item `p`: A ggplot2 bar plot visualizing the importance of each environmental variable based on `%IncMSE`.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Converts OTU data to a data frame and calculates Bray-Curtis distances.
#'   \item Computes Non-Metric Multidimensional Scaling (NMDS) using Bray-Curtis distances.
#'   \item Combines the first NMDS axis (MDS1) with environmental variables to create a model dataset.
#'   \item Applies Random Forest regression to assess the importance of each environmental variable in explaining MDS1 variation.
#'   \item Outputs a bar plot of variable importance based on `%IncMSE`.
#' }
#'
#' @examples
#' \dontrun{
#' result <- ramdom_env.plot(ps = ps, env = env, seed = 42)
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @export
ramdom_env.plot <- function(
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
