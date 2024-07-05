
ancombc2.trans = function(
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

    tse = mia::makeTreeSummarizedExperimentFromPhyloseq(ps.cs)
    out <- ancombc2(
      data = tse,
      assay_name = "counts",
      tax_level = NULL,
      fix_formula = group,
      p_adj_method = "fdr",
      lib_cut = 0,
      prv_cut = 0,
      group = group,
      struc_zero = TRUE,
      neg_lb = TRUE,
      alpha = 0.05,
      global = TRUE # multi group comparison will be deactivated automatically
    )

    str_detect("q_Group",colnames(out$res))
    # store the FDR adjusted results [test on v2.0.3]
    ancombc_result <- cbind.data.frame(taxid = out$res$taxon,
                                       p = out$res[,11],
                                       lfc = out$res[,3],
                                       w = out$res[,7]
    )
    head(ancombc_result)

    tab.d15 = ancombc_result %>%
      # rownames_to_column(var = "id") %>%
      dplyr::select(taxid,p) %>%
      dplyr::filter( p < alpha) %>%
      dplyr::rename(
        OTU = "taxid"
      )  %>%
      dplyr::mutate(group = "ancombc2")

    head(tab.d15)
    res =  data.frame(micro = tab.d15$OTU,method = tab.d15$group,
                      adjust.p = tab.d15$p)
    data_list[[i]]= res
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)
}



