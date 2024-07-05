


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
