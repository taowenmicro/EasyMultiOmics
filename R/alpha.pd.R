# library(ape)
# library(picante)

alpha.pd = function(ps){
  com_2020 <- ps %>% vegan_otu() %>%
    as.data.frame()

  rooted <- phy_tree(ps)


  cover2020.pd<-pd(com_2020,rooted,include.root=F)

  map = sample_data(ps)
  head(map)

  data = cbind(map[,c("ID","Group")],pd = cover2020.pd[,1])
  head(data)
  colnames(data)[2] = "group"
  data$group = as.factor(data$group)
  return(data)
}
