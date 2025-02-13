#' @title  Perform differential analysis using Wilcox-CLR
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting Wilcox-CLR analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param alpha Significance level, default is 0.05.
#' @return A list containing the results of significant differences in genes abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}.
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat = wilcox.clr.trans(ps = ps%>% filter_OTU_ps(50),group =  "Group",alpha = 0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export

wilcox.clr.trans = function(
    ps = ps,
    group =  "Group",
    alpha = 0.05 ){


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
      # subset_samples(Group %in% id.g[,i]) %>%
      vegan_otu() %>% t() %>%
      as.data.frame()
    groupings <-ps.cs %>%
      # subset_samples(Group %in% id.g[,i]) %>%
      sample_data()
    groupings$ID = row.names(groupings)


    #add pseudo count
    CLR_table <- data.frame(apply(ASV_table + 1, 2, function(x){log(x) - mean(log(x))}))
    ## get clr table
    #apply wilcox test to rarified table

    map =  sample_data(ps.cs)
    g = map[,group] %>%
      as.vector() %>% .[[1]] %>% as.factor()
    pvals <- apply(CLR_table, 1, function(x) wilcox.test(x ~ g, exact=FALSE)$p.value)

    dat <- pvals %>% as.data.frame()
    head(dat)
    colnames(dat) = "p"

    tab.d13 = dat %>%
      rownames_to_column(var = "id") %>%
      dplyr::select(id,p) %>%
      dplyr::filter( p < alpha) %>%
      dplyr::rename(
        OTU = id
        # p = p
      )  %>%
      dplyr::mutate(group = " wilcox.test.CLR")

    head(tab.d13)
    res = diff.tab = data.frame(micro = tab.d13$OTU,method = tab.d13$group,
                                adjust.p = tab.d13$p)

    data_list[[i]]= res

    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)

}
