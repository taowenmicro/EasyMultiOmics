ancombc.micro = function(
    ps = ps,
    group = "Group",
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



    tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.cs)

    out = ancombc(data = tse,
                  assay_name = "counts",
                  tax_level = NULL,
                  phyloseq = NULL,
                  formula = group,# 微生物群落数据受到那些变量影响，用公式来写
                  p_adj_method = "holm",
                  prv_cut = 0.10,
                  lib_cut = 1000,
                  group = group,
                  struc_zero = TRUE,
                  neg_lb = TRUE,
                  tol = 1e-5,
                  max_iter = 100,
                  conserve = TRUE,
                  alpha = 0.05,
                  global = TRUE,
                  n_cl = 1,
                  verbose = TRUE)

    res = out$res
    res2 = out$res$diff_abn %>% as.data.frame()
    head(res2)
    colnames(res2)[3] = "lo"
    res2$`(Intercept)` = NULL
    res3 = res$q_val%>% as.data.frame()
    colnames(res3)[3] = "p"
    tem = res3 %>% dplyr::filter(p < alpha)

    tab.d1 = tem %>%
      as.data.frame() %>%
      dplyr::rename(
        OTU = taxon
      ) %>%
      dplyr::mutate(group = "ancombc")


    res =diff.tab = data.frame(micro = tab.d1$OTU,method = tab.d1$group,adjust.p = tab.d1$p %>% as.vector())
    data_list[[i]]= res
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)
}
