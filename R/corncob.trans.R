corncob.trans= function(
    ps = ps.cs,
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
    my_formula <- as.formula(paste("~",group,sep=" ", collapse = ""))
    # my_formula
    results <- corncob::differentialTest(formula= my_formula,
                                         phi.formula = my_formula,
                                         phi.formula_null = my_formula,
                                         formula_null = ~ 1,
                                         test="Wald", data= ps.cs ,
                                         boot=FALSE,
                                         fdr_cutoff = alpha)
    dat = results$p_fdr  %>% as.data.frame()
    head(dat)
    colnames(dat) = "p_fdr"
    dat$p = results$p

    tab.d3 = dat %>%
      rownames_to_column(var = "id") %>%
      dplyr::select(id,p_fdr) %>%
      dplyr::filter(p_fdr < alpha) %>%
      dplyr::rename(
        OTU = id,
        p = p_fdr
      )  %>%
      dplyr::mutate(group = "corncob")

    data_list[[i]]= tab.d3
    names( data_list) [i]= Desep_group %>% paste( collapse = "_")

  }
  return(data_list)

}
