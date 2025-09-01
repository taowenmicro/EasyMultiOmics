#' @title  Perform differential analysis using LinDA
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting LinDA analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param alpha Significance level, default is 0.05.
#' @return A list containing the results of significant differences in genes abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}.
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat = linda.trans(ps = ps%>% filter_OTU_ps(50),group =  "Group",alpha = 0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export

linda.trans = function(
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
    meta <- as.data.frame(colData(tse)) %>% dplyr::select(group)
    my_formula <- as.formula(paste("~",group,sep=" ", collapse = ""))

    linda.res <- linda(
      as.data.frame(assay(tse)),
      meta,
      formula =paste0("~",group),
      alpha = 0.05,
      prev.filter = 0,
      mean.abund.filter = 0)
    linda_out <- linda.res$output[[1]]
    head(linda_out )
    tab.d3 = linda_out  %>% filter(padj < 0.05) %>% rownames_to_column("id") %>%
      dplyr::select(id,padj) %>%
      dplyr::rename(
        OTU = id,
        p = padj
      )  %>%
      dplyr::mutate(group = "linda")


    res =data.frame(micro = tab.d3$OTU,method = tab.d3$group,adjust.p = tab.d3$p)

    data_list[[i]]= res
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)
}
