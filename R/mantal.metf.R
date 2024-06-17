#
#
# #--门特尔检验-普氏分析#-----
#
# method = "spearman"
#
# head(sample_data(ps0))
# result <- mantal.micro(ps = ps0,method =  "spearman",group = "Group")
#
# combn(unique(map$Group),2)
# result[[2]]

#
#
# source("G:\\Shared_Folder\\Function_local\\R_function\\micro/matel_pro_plot.R")
#
# #--门特尔检验-普氏分析#-----
# library(vegan)
# size = combn(unique(map$Group),2) %>% dim()
# size
# result <- mantal.micro(ps = ps,method =  "spearman",group = "Group",
#                        ncol = size[2],
#                        nrow = 1
# )
# data <- result[[1]]
#
# p3_7 <- result[[2]] +  mytheme1
# p3_7
#
#
# FileName <- paste(betapath,"mantel_pro.csv", sep = "")
# write.csv(data,FileName)
# FileName1 <- paste(betapath,"/a2_","Mantel_Pro.pdf", sep = "")
# ggsave(FileName1 , p3_7, width = 8, height = 8)



mantal.meta <- function(ps = ps,
                         method =  "spearman",
                         group = "Group",
                         ncol = 5,
                         nrow = 2

){
  dist <-
    ggClusterNet::scale_micro(ps = ps,method = "rela") %>%
    ggClusterNet::vegan_otu()%>%
    vegan::vegdist(method="bray") %>%
    as.matrix()

  map = phyloseq::sample_data(ps)
  gru = map[,group][,1] %>% unlist() %>% as.vector()
  id = combn(unique(gru),2)

  R_mantel = c()
  p_mantel = c()
  name = c()
  R_pro <- c()
  p_pro <- c()
  plots = list()

  for (i in 1:dim(id)[2]) {

    id_dist <- row.names(map)[gru == id[1,i]]
    dist1 = dist[id_dist,id_dist]
    id_dist <- row.names(map)[gru == id[2,i]]
    id_dist = id_dist[1:nrow(dist1)]
    dist2 = dist[id_dist,id_dist]
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

    Y <- cbind(data.frame(pro.s.r$Yrot), data.frame(pro.s.r$X))
    X <- data.frame(pro.s.r$rotation)
    Y$ID <- rownames(Y)



    p1 <- ggplot(Y) +
      geom_segment(aes(x = X1, y = X2, xend = (X1 + MDS1)/2, yend = (X2 + MDS2)/2),
                   arrow = arrow(length = unit(0, 'cm')),
                   color = "#B2182B", size = 1) +
      geom_segment(aes(x = (X1 + MDS1)/2, y = (X2 + MDS2)/2, xend = MDS1, yend = MDS2),
                   arrow = arrow(length = unit(0, 'cm')),
                   color = "#56B4E9", size = 1) +
      geom_point(aes(X1, X2), fill = "#B2182B", size = 4, shape = 21) +
      geom_point(aes(MDS1, MDS2), fill = "#56B4E9", size = 4, shape = 21) +
      labs(title =  paste(id[1,i],"-",id[2,i]," ","Procrustes analysis:\n    M2 = ",
                          round(pro.s.r$ss,3),
                          ", p-value = ",
                          round(pro.s.r$signif,3),
                          "\nMantel test:\n    r = ",
                          round(R_mantel[i],3),
                          ", p-value =, ",
                          round(p_mantel[i],3),sep = "") )
    p1

    plots[[i]] = p1
  }

  dat = data.frame(name,R_mantel,p_mantel,R_pro,p_pro )
  pp  = ggpubr::ggarrange(plotlist = plots, common.legend = TRUE, legend="right",ncol = ncol,nrow = nrow)
  return(list(dat,pp))
}
