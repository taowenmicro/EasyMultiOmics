
# result = sankey.m.Group(
#     ps = ps,
#     rank = 6,
#     Top = 50,
#
#
# )
#
# p = result[[1]]
# dat = result[[2]]
#
# FileName <-paste(snapath,"/sankey_Group.csv", sep = "")
# write.csv(dat,FileName,sep = "")
#
# saveNetwork(p,paste(snapath,"/sankey_Group.html", sep = ""))
# library(webshot)
# # webshot::install_phantomjs()
# # webshot(paste(snapath,"/sankey1.html", sep = "") ,paste(snapath,"/sankey1.png", sep = ""))
# webshot(paste(snapath,"/sankey_Group.html", sep = "") , paste(snapath,"/sankey_Group.pdf", sep = ""))
#




sankey.m.Group.micro  = function(
    ps = ps,
    rank = 6,# 参数目前不可修改
    Top = 50

  ){

  map = sample_data(ps) %>% as.tibble()
  map$ID = row.names(map)
  otu =  ps %>% vegan_otu() %>%
    as.data.frame()
  otu$ID = row.names(otu)

  tax = ps %>%
    scale_micro() %>%
    vegan_tax() %>%
    as.data.frame()
  tax$taxid = row.names(tax)
  head(tax)


  data = map %>% as.tibble() %>%
    inner_join(otu) %>%
    gather( taxa, count, starts_with("ASV_")) %>%
    inner_join(tax,by = c("taxa" = "taxid") )
  data


  tax = ps %>%
    ggClusterNet::tax_glom_wt(ranks = rank) %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    subset_taxa(!Genus  %in% c("Unassigned","Unknown")) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()
  head(tax)
  dim(tax)
  id2 = c("k","p","c","o","f","g")
  dat = NULL
  for (i in 1:5) {
    dat <- tax[,c(i,i+1)] %>% distinct(.keep_all = TRUE)
    colnames(dat) = c("source","target")
    dat$source = paste(id2[i],dat$source,sep = "_")
    dat$target = paste(id2[i+1],dat$target,sep = "_")
    if (i == 1) {
      dat2 = dat
    }

    dat2 = rbind(dat2,dat)

  }
  dim(dat2)
  # dat2 = dat2 %>% distinct(.keep_all = TRUE)



  otu = ps %>%
    ggClusterNet::tax_glom_wt(ranks = 6) %>%
    ggClusterNet::scale_micro() %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    subset_taxa(!Genus  %in% c("Unassigned","Unknown")) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>%
    as.data.frame()

  head(otu)
  otutax = cbind(otu,tax)

  id = rank.names(ps)[1:6]
  dat = NULL
  dat3 = NULL
  head(data)

  for (j in 1) {
    dat = data %>% group_by(Group,!!sym(rank_names(ps)[j])) %>%
      summarise_if(is.numeric,sum,na.rm = TRUE)
    #
    colnames(dat) = c("source","target","value")
    dat$target = paste(id2[j],dat$target,sep = "_")

    if (j == 1) {
      dat3 = dat
    }
    # dat3 = rbind(dat3,dat)
  }


  dat3 %>% tail()
  tem = data.frame(target = dat3$target,value = dat3$value)
  dat3$value = NULL
  head(dat2)
  dat2 = rbind(dat2,dat3)


  head( otutax)
  for (i in 1:6) {
    dat <- otutax %>%
      dplyr::group_by(!!sym(id[i])) %>%
      summarise_if(is.numeric,sum,na.rm = TRUE)

    dat = data.frame(Genus = dat[,1],value =rowSums(dat[,-1]) )
    colnames(dat) = c("target","value")
    dat$target = paste(id2[i],dat$target,sep = "_")
    if (i == 1) {
      dat3 = dat
    }

    dat3 = rbind(dat3,dat)
  }

  dat4 <- dat2 %>% left_join(dat3)
  sankey = dat4

  head( sankey )


  nodes <- data.frame(name = unique(c(as.character(sankey$source),as.character(sankey$target))),stringsAsFactors = FALSE)
  nodes$ID <- 0:(nrow(nodes)-1)
  sankey <- merge(sankey,nodes,by.x = "source",by.y = "name")
  sankey <- merge(sankey,nodes,by.x = "target",by.y = "name")
  colnames(sankey) <- c("X","Y","value","source","target")
  sankey <- subset(sankey,select = c("source","target","value"))
  nodes <- subset(nodes,select = c("name"))

  ColourScal='d3.scaleOrdinal() .range(["#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"])'

  sankey$energy_type <- sub(' .*', '', nodes[sankey$source + 1, 'name'])
  library(networkD3)
  p <- sankeyNetwork(Links = sankey, Nodes = nodes,
                     Source = "source",Target = "target",Value = "value",
                     NodeID = "name",
                     sinksRight=FALSE,
                     LinkGroup = 'energy_type',
                     colourScale= ColourScal,
                     # nodeWidth=40,
                     # fontSize=13
                     # nodePadding=20
  )

  return(list(p,sankey))

}
