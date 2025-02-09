#' @title  Perform differential analysis using Maaslin2
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting Maaslin2 analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param alpha Significance level, default is 0.05.
#' @param rare Boolean indicating whether to apply rarefaction, default is TRUE.
#' @return A list containing the results of significant differences in genes abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}.
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat = maaslin2.trans(ps = ps%>% filter_OTU_ps(50),group =  "Group",alpha = 0.05,rare = FALSE)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export
maaslin2.trans = function(
    ps = ps,
    group =  "Group",
    alpha = 0.05,
    rare = TRUE
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


    if (rare == TRUE) {
      ASV_table = ps.cs %>%
        scale_micro(method = "sampling") %>%
        # subset_samples(Group %in% id.g[,i]) %>%
        vegan_otu() %>% t() %>%
        as.data.frame()

      tem = "Maaslin2.rare"
    } else if(rare == FALSE){
      ASV_table = ps.cs %>%
        filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
        # subset_samples(Group %in% id.g[,i]) %>%
        vegan_otu() %>% t() %>%
        as.data.frame()
      tem = "Maaslin2"
    }



    groupings <- ps.cs %>%
      # subset_samples(Group%in% id.g[,i]) %>%
      sample_data()
    groupings$ID = row.names(groupings)

    # 会生成一个文件夹
    # 一个稀释后的或未经稀释的特征表
    # 反正弦平方根变换(arcsine square-root transformation)。
    # 没有指定随机效应，并关闭了默认的标准化。该函数将线性模型拟合到指定样本分组上每个特征的转换丰度，
    # 使用Wald检验进行显著性检验，并输出BH FDR校正的p值
    ASV_table <- data.frame(t(ASV_table), check.rows = FALSE, check.names = FALSE,
                            stringsAsFactors = FALSE
    )

    row.names(groupings) = groupings$ID

    fit_data <- Maaslin2(
      ASV_table, groupings,"Maaslin2", transform = "AST",
      fixed_effects = group,
      standardize = FALSE, plot_heatmap = FALSE, plot_scatter = FALSE)

    dat = fit_data$results
    head(dat)
    tab.d9 = dat %>%
      # rownames_to_column(var = "id") %>%
      dplyr::select(feature,qval) %>%
      dplyr::filter(qval < alpha) %>%
      dplyr::rename(
        OTU = feature,
        p = qval
      )  %>%
      dplyr::mutate(group = tem)

    # head(tab.d9)
    res = diff.tab = data.frame(micro = tab.d9$OTU,method = tab.d9$group,
                                adjust.p = tab.d9$p)

    data_list[[i]]= res
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")

  }
  return(data_list)

}

