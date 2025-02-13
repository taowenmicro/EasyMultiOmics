#' @title  Perform differential analysis using Limma-Voom-TMM/TMMwsp
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting Limma-Voom-TMM/TMMwsp analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param alpha Significance level, default is 0.05.
#' @param method Normalization method,"TMM" or "TMMwsp", default is "TMM".
#' @return A list containing the results of significant differences in genes abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat = limma.v.TMM.trans(ps = ps%>% filter_OTU_ps(50),group =  "Group",alpha = 0.05,method = "TMM")
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export
limma.v.TMM.trans = function(
    ps = ps,
    group =  "Group",
    alpha = 0.05,
    method = "TMM"# method = "TMMwsp"
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
    ASV_table = ps.cs %>%
      filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
      # subset_samples(Group %in% id.g[,i]) %>%
      vegan_otu() %>% t() %>%
      as.data.frame()
    groupings <-  ps.cs %>%
      # subset_samples(Group%in% id.g[,i]) %>%
      sample_data()
    groupings$ID = row.names(groupings)
    DGE_LIST <- edgeR::DGEList(ASV_table)
    ### do normalization
    ### Reference sample will be the sample with the highest read depth
    ### check if upper quartile method works for selecting reference
    Upper_Quartile_norm_test <- edgeR::calcNormFactors(DGE_LIST, method="upperquartile")
    summary_upper_quartile <- summary(Upper_Quartile_norm_test$samples$norm.factors)[3]
    # if(is.na(summary_upper_quartile) | is.infinite(summary_upper_quartile)){
    #   message("Upper Quartile reference selection failed will use find sample with largest sqrt(read_depth) to use as reference")
    #   Ref_col <- which.max(colSums(sqrt(ASV_table)))
    #   DGE_LIST_Norm <- edgeR::calcNormFactors(DGE_LIST, method = method, refColumn = Ref_col)
    #   fileConn<-file(args[[4]])
    #   writeLines(c("Used max square root read depth to determine reference sample"), fileConn)
    #   close(fileConn)
    #
    # }else{
    DGE_LIST_Norm <- edgeR::calcNormFactors(DGE_LIST, method=method)
    # }

    ## make matrix for testing
    # colnames(groupings) <- c("comparison")
    groupings = groupings %>% as.tibble() %>% as.data.frame()
    mm <- model.matrix(my_formula, groupings)

    voomvoom <- voom(DGE_LIST_Norm, mm, plot=FALSE)

    fit <- lmFit(voomvoom,mm)
    fit <- eBayes(fit)
    res <- topTable(fit, coef=2, n=nrow(DGE_LIST_Norm), sort.by="none")
    head(res)

    tab.d7 = res %>%
      rownames_to_column(var = "id") %>%
      dplyr::select(id,adj.P.Val) %>%
      dplyr::filter(adj.P.Val < 0.05) %>%
      dplyr::rename(
        OTU = id,
        p = adj.P.Val
      )  %>%
      dplyr::mutate(group = paste0("limma.voom.",method))

    data_list[[i]]= data.frame(micro = tab.d7$OTU,method = tab.d7$group,
                               adjust.p = tab.d7$p)
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)
}
