
lefse.trans = function(
    ps = ps,
    group =  "Group",
    alpha = 0.05
){


  map= sample_data(ps)

  head(map)

  id.g = map$Group %>% unique() %>% as.character() %>% combn(2)

  aaa = id.g

  data_list =  list()

  for (i in 1:dim(aaa)[2]) {
    # i = 1
    Desep_group = aaa[,i]
    print( Desep_group)
    ps.cs = ps %>% subset_samples.wt("Group" ,id.g[,i])

    tablda = LDA_trans(ps = ps.cs ,
                       group = group,
                       Top = 0,
                       p.lvl = 0.05,
                       lda.lvl = 0,
                       seed = 11,
                       adjust.p = FALSE)

    dat = tablda[[2]]
    head(dat)
    dat$ID = row.names(dat)
    row.names(dat) = NULL


    dat2 <- dat

    tab.d6 = dat2 %>%
      dplyr::select(ID,Pvalues,LDAscore,class) %>%
      dplyr::filter(Pvalues < 0.05) %>%
      dplyr::rename(

        p = Pvalues
      )  %>%
      dplyr::mutate(group = "LEFse")

    head(tab.d6)

    res = tab.d6

    data_list[[i]]= res
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)

}











LDA_trans = function(ps = ps.cs,
                     group = "Group",
                     Top = 100,
                     p.lvl = 0.05,
                     lda.lvl = 2,
                     seed = 11,
                     adjust.p = F
){

  ps = ps %>%  filter_taxa(function(x) sum(x ) > 0 , TRUE)
  alltax = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_tax() %>%
    as.data.frame()
  alltax$OTU = row.names(alltax)

  # alltax$Kingdom = paste(alltax$Kingdom,sep = "_Rank_")
  # alltax$Phylum = paste(alltax$Kingdom,alltax$Phylum,sep = "_Rank_")
  # alltax$Class = paste(alltax$Phylum,alltax$Class,sep = "_Rank_")
  # alltax$Order = paste(alltax$Class,alltax$Order,sep = "_Rank_")
  # alltax$Family = paste(alltax$Order,alltax$Family,sep = "_Rank_")
  # alltax$Genus = paste(alltax$Family,alltax$Genus,sep = "_Rank_")
  # alltax$Species = paste(alltax$Genus,alltax$Species,sep = "_Rank_")
  #

  otu = ps %>%
    ggClusterNet::filter_OTU_ps(Top) %>%
    ggClusterNet::vegan_otu() %>%
    t() %>%
    as.data.frame()

  otu_tax = merge(otu,alltax,by = "row.names",all = F)
  head(otu_tax)

  rank <- otu_tax %>%
    dplyr::group_by(OTU) %>%
    dplyr::summarise_if(is.numeric, sum, na.rm = TRUE)

  colnames(rank)[1] = "id"

  #--LDA排序#--------
  data1 = as.data.frame(rank)
  row.names(data1) = data1$id
  data1$id = NULL

  #-构建phylose对象

  ps_G_graphlan = phyloseq::phyloseq(phyloseq::otu_table(as.matrix(data1),taxa_are_rows = TRUE),
                                     phyloseq::sample_data(ps)) %>%  filter_taxa(function(x) sum(x ) > 0, TRUE)
  ps_G_graphlan

  #----提取OTU表格






  otu = as.data.frame((ggClusterNet::vegan_otu(ps_G_graphlan)))
  otu[otu==0] <- 1
  otu = otu[ colMeans(otu) != 1]


  map = as.data.frame(phyloseq::sample_data(ps_G_graphlan))
  # otu = (otu_table)
  claslbl= map[,group] %>% as.vector() %>% .[[1]] %>% as.factor()
  # claslbl= map$Group %>% as.factor()
  set.seed(seed)
  #KW rank sum test

  rawpvalues <- apply(otu, 2, function(x) kruskal.test(x, claslbl)$p.value);
  #--得到计算后得到的p值
  ord.inx <- order(rawpvalues)
  rawpvalues <- rawpvalues[ord.inx]
  clapvalues <- p.adjust(rawpvalues, method ="fdr")

  # p.adjust
  wil_datadf <- as.data.frame(otu[,ord.inx])


  ldares <- MASS::lda(claslbl ~ .,data = wil_datadf)
  # ldares
  ldamean <- as.data.frame(t(ldares$means))
  ldamean
  class_no <<- length(unique(claslbl))
  ldamean$max <- apply(ldamean[,1:class_no],1,max);
  ldamean$min <- apply(ldamean[,1:class_no],1,min);
  #---计算LDA
  ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);
  head(ldamean)

  a = rep("A",length(ldamean$max))
  for (i in 1:length(ldamean$max)) {
    name =colnames(ldamean[,1:class_no])
    a[i] = name[ldamean[,1:class_no][i,] %in% ldamean$max[i]]
  }
  ldamean$class = a

  tem1 = row.names(ldamean)
  tem1 %>% as.character()
  ldamean$Pvalues <- signif(rawpvalues[match(row.names(ldamean),names(rawpvalues))],digits=5)
  ldamean$FDR <- signif(clapvalues,digits=5)
  resTable <- ldamean
  rawNms <- rownames(resTable);
  rownames(resTable) <- gsub("`", '', rawNms);


  if (adjust.p) {
    de.Num <- sum(clapvalues <= p.lvl & ldamean$LDAscore>=lda.lvl)

  } else {
    de.Num <- sum(rawpvalues <= p.lvl & ldamean$LDAscore>=lda.lvl)
  }

  if(de.Num == 0){
    current.msg <<- "No significant features were identified with given criteria.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
  }
  print(current.msg)
  # sort by p value
  ord.inx <- order(resTable$Pvalues, resTable$LDAscore)
  resTable <- resTable[ord.inx, ,drop=FALSE]
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))]

  # resTable %>% tail()
  ldamean$Pvalues[is.na(ldamean$Pvalues)] = 1
  if (adjust.p) {
    taxtree = resTable[clapvalues <=p.lvl & ldamean$LDAscore>=lda.lvl,]
  } else {
    # taxtree = resTable[ldamean$Pvalues <=p.lvl & ldamean$LDAscore>=lda.lvl,]
    taxtree = resTable[ldamean$Pvalues <=p.lvl,]
  }

  #-提取所需要的颜色
  colour = c('darkgreen','red',"blue","#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  selececol = colour[1:length(levels(as.factor(taxtree$class)))]
  names(selececol) = levels(as.factor(taxtree$class))
  A = rep("a",length(row.names(taxtree)))

  for (i in 1:length(row.names(taxtree))) {
    A[i] = selececol [taxtree$class[i]]
  }

  taxtree$color = A
  # taxtree <- taxtree[row.names(taxtree) != "k__Bacteria",]
  # node_ids <- p0$data
  # anno <- rep("white", nrow(p1$data))

  lefse_lists = data.frame(node=row.names(taxtree),
                           color=A,
                           Group = taxtree$class,
                           stringsAsFactors = FALSE )


  return(list(lefse_lists,taxtree))
}
