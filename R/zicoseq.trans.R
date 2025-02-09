#' @title  Perform differential analysis using ZicoSeq
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting ZicoSeq analysis to detect significant differences in genes
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
#' dat = zicoseq.trans(ps = ps%>% filter_OTU_ps(50),group =  "Group",alpha = 0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export

zicoseq.trans = function(
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


    comm <- ps.cs %>%
      filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
      vegan_otu() %>% t()
    meta.dat <- ps.cs %>% sample_data()
    meta.dat
    meta.dat$id = row.names(meta.dat)
    meta.dat = meta.dat %>% as.tibble() %>% column_to_rownames("id")


    # meta.dat[, grp.name]


    ZicoSeq.obj <- ZicoSeq(meta.dat = meta.dat,
                           feature.dat = comm,
                           grp.name = 'Group',
                           # adj.name = 'Sex',
                           feature.dat.type = "count",
                           # Filter to remove rare taxa
                           # prev.filter = 0,
                           # mean.abund.filter = 0,
                           # max.abund.filter = 0,
                           # min.prop = 0,
                           # Winsorization to replace outliers
                           is.winsor = TRUE,
                           outlier.pct = 0.03,
                           winsor.end = 'top',
                           # Posterior sampling
                           is.post.sample = TRUE,
                           post.sample.no = 25,
                           # Use the square-root transformation
                           link.func = list(function (x) x^0.5),
                           stats.combine.func = max,
                           # Permutation-based multiple testing correction
                           perm.no = 99,  strata = NULL,
                           # Reference-based multiple stage normalization
                           ref.pct = 0.5, stage.no = 6, excl.pct = 0.2,
                           # Family-wise error rate control
                           is.fwer = TRUE, verbose = TRUE, return.feature.dat = TRUE)


    dat = data.frame(row.names = names(ZicoSeq.obj$p.raw),
                     id = names(ZicoSeq.obj$p.raw),
                     p.raw = ZicoSeq.obj$p.raw,
                     p.adj = ZicoSeq.obj$p.adj.fdr,
                     p.fwer = ZicoSeq.obj$p.adj.fwer
    )
    head(dat)

    tab.d16 = dat %>%
      # rownames_to_column(var = "id") %>%
      dplyr::select(id,p.adj) %>%
      dplyr::filter( p.adj < alpha) %>%
      dplyr::rename(
        p = "p.adj"
      )  %>%
      dplyr::mutate(group = "ZicoSeq")

    head(tab.d16)
    res =  data.frame(micro = tab.d16$id,method = tab.d16$group,
                      adjust.p = tab.d16$p)
    data_list[[i]]= res
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)
}



