#' @title  Perform differential analysis using metagenomeSeq
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting metagenomeSeq analysis to detect significant differences in genes
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
#' dat = metaSeq.trans(ps = ps%>% filter_OTU_ps(50),group =  "Group",alpha = 0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export
metaSeq.trans = function(ps = ps,
                         group =  "Group",
                         alpha = 0.05){
  map= sample_data(ps)
  head(map)
  id.g = map$Group %>% unique() %>% as.character() %>% combn(2)
  aaa = id.g
  data_list2 =  list()

  for (i in 1:dim(aaa)[2]) {
    # i = 1
    Desep_group = aaa[,i]
    print( Desep_group)
    ps.cs = ps %>% subset_samples.wt("Group" ,id.g[,i])
    ASV_table = ps.cs %>%
      filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
      # subset_samples(Group %in% id.g[,i]) %>%
      vegan_otu() %>% t() %>%
      as.data.frame()
    groupings <- ps.cs %>%
      # subset_samples(Group%in% id.g[,i]) %>%
      sample_data()
    groupings$ID = row.names(groupings)

    data_list <- list()
    data_list[["counts"]] <- ASV_table
    data_list[["taxa"]] <- rownames(ASV_table)

    pheno <- AnnotatedDataFrame(groupings)
    pheno
    counts <- AnnotatedDataFrame(ASV_table)
    feature_data <- data.frame("ASV"=rownames(ASV_table),
                               "ASV2"=rownames(ASV_table))
    feature_data <- AnnotatedDataFrame(feature_data)
    rownames(feature_data) <- feature_data@data$ASV


    test_obj <- newMRexperiment(counts = data_list$counts,
                                phenoData = pheno,
                                featureData = feature_data)

    p <- cumNormStat(test_obj, pFlag = TRUE)

    test_obj_norm <- cumNorm(test_obj, p=p)
    fromula <- as.formula(paste(~1, group, sep=" + "))
    pd <- pData(test_obj_norm)
    mod <- model.matrix(fromula, data=pd)
    regres <- fitFeatureModel(test_obj_norm, mod)
    res_table <- MRfulltable(regres, number = length(rownames(ASV_table)))
    # head(res_table)
    tab.d10 = res_table %>%
      rownames_to_column(var = "id") %>%
      dplyr::select(id,adjPvalues) %>%
      dplyr::filter(adjPvalues < alpha) %>%
      dplyr::rename(
        OTU = id,
        p = adjPvalues
      )  %>%
      dplyr::mutate(group = "metagenomeSeq")

    # head(tab.d10)
    res = data.frame(micro = tab.d10$OTU,method = tab.d10$group,
                     adjust.p = tab.d10$p)

    data_list2[[i]]=  res

    names( data_list2)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list2)

}
