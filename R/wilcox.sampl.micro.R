
wilcox.sampl.micro = function(
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

    ASV_table = ps.cs %>%
    scale_micro(method = "sampling") %>%
    # subset_samples(Group %in% id.g[,i]) %>%
    vegan_otu() %>% t() %>%
    as.data.frame()

  groupings <- ps.cs %>% sample_data()
    # subset_samples(Group %in% id.g[,i]) %>%
  groupings$ID = row.names(groupings)

  map= sample_data(ps.cs)

  g = map[,group] %>%
    as.vector() %>% .[[1]] %>% as.factor()
  pvals <- apply(ASV_table, 1, function(x) wilcox.test(x ~ g, exact=FALSE)$p.value)
  dat <- pvals %>% as.data.frame()
  head(dat)
  colnames(dat) = "p"
  tab.d12 = dat %>%
    rownames_to_column(var = "id") %>%
    dplyr::select(id,p) %>%
    dplyr::filter( p < alpha) %>%
    dplyr::rename(
      OTU = id
      # p = p
    )  %>%
    dplyr::mutate(group = " wilcox.test.rare")

  # head(tab.d12)
  res = data.frame(micro = tab.d12$OTU,method = tab.d12$group,
                                   adjust.p = tab.d12$p)

  data_list[[i]]= res
  names( data_list)[i] = Desep_group %>% paste( collapse = "_")
  }
  return(data_list)

}

# write.table(dat, file = "wilcox.test.rare.txt", sep="\t", col.names = NA, quote=F)
