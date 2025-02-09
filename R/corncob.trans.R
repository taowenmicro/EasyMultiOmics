#' @title  Perform differential analysis using Corncob
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting Corncob analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param alpha Significance level, default is 0.05.
#' @return A list containing the results of significant differences in gene abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat <- corncob.trans(ps = ps%>% filter_OTU_ps(50), group = "Group", alpha = 0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export
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
   # otu_table(ps.cs)
    otu_table(ps.cs)  =  round( otu_table(ps.cs))
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
