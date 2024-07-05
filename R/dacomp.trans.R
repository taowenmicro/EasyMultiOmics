
dacomp.trans =  function(ps = ps0,
                         sd = 0.6,
                         q_BH =  0.1,
                         q_DSFDR =0.1){


  map1= sample_data(ps)
  head(map1)
  id.g = map1$Group %>% unique() %>% as.character() %>% combn(2)
  aaa = id.g
  data_list =  list()

  for (i in 1:dim(aaa)[2]) {
    # i = 1
    Desep_group = aaa[,i]
    print( Desep_group)
    ps.cs = ps %>% subset_samples.wt("Group" ,id.g[,i])
    otu = ps.cs  %>% vegan_otu()

    result.selected.references = dacomp.select_references(
      X = otu,
      median_SD_threshold = sd, #APPLICATION SPECIFIC
      verbose = F)

    print(result.selected.references)

    q_BH = q_DSFDR = 0.1
    map = ps.cs %>% sample_data()

    result.test = dacomp.test(X = otu, #counts data
                              y = map$Group, #phenotype in y argument
                              # obtained from dacomp.select_references(...):
                              ind_reference_taxa = result.selected.references,
                              test = DACOMP.TEST.NAME.WILCOXON, #constant, name of test
                              verbose = F,q = q_DSFDR) #multiplicity adjustment level

    rejected_BH = which(p.adjust(result.test$p.values.test,method = 'BH')<=q_BH)

    rejected_DSFDR = result.test$dsfdr_rejected
    ?p.adjust
    res = data.frame(id = colnames(otu), p = p.adjust(result.test$p.values.test,method = 'BH'))

    res$p[is.na(res$p)] = 1

    data_list[[i]]= res
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")}
  return(data_list)
}
