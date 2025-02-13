#' @title  Perform differential analysis using Lefse
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting Lefse analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param alpha Significance level, default is 0.05.
#' @return A list containing the results of significant differences in genes abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat = lefse.trans(ps = ps%>% filter_OTU_ps(50),group =  "Group",alpha = 0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export
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

    tablda = LDA.trans(ps = ps.cs ,
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
