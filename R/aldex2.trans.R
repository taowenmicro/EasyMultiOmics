#' @title  Perform differential analysis using Aldex2
#' @description This function is used to analyze transcriptome functional composition
#'  data  by conducting Aldex2 analysis to detect significant differences in genes
#'  abundance between different groups.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param group Column name for groupID in map table(sample metadata).
#' @param test Statistical test method, default is "t" for Welch's t-test.
#' @param alpha Significance level, default is 0.05.
#' @return A list containing the results of significant differences in gene abundance between different groups.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' \dontrun{
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' dat <- aldex2.trans(ps = ps%>% filter_OTU_ps(50), group = "Group", test = "t", alpha = 0.05)
#' dat1 = dat$WT_OE
#' head(dat1)
#' dat2 = dat$WT_KO
#' head(dat2)
#' dat3 = dat$OE_KO
#' head(dat3)
#' }
#' @export
aldex2.trans = function(
    ps = ps0,
    group = "Group",
    test = "t",# 传承参数
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

    ASV_table = ps.cs%>%
      filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
      vegan_otu() %>% round() %>% t() %>%
      as.data.frame()
    groupings <- ps.cs %>%
      sample_data()
    groupings$ID = row.names(groupings)
    #mc.samples:一个整数。 估计基础分布时使用的蒙特卡洛样本数
    results <- aldex(reads=ASV_table,
                     conditions = groupings[,group] %>% as.vector() %>% .[[1]],
                     mc.samples = 128,
                     test=test,
                     effect=TRUE,
                     include.sample.summary = FALSE,
                     verbose=TRUE,
                     denom="all")

    # head(results)
    #使用Welchs t 检验矫正后的p值 作为显著差异的微生物
    tab.d1 = results %>%
      as.data.frame() %>%
      dplyr::filter(we.ep < alpha) %>% rownames_to_column(var = "id") %>%
      dplyr::select(id,we.ep) %>%
      dplyr::rename(
        OTU = id,
        p = we.ep
      ) %>%
      dplyr::mutate(group = "Aldex2")
    # 输出结果解读：we.ep：Welchs t 检验矫正后的p值，wi.ep 秩和检验p值，带有BH结尾的就是矫正后的P值
    # 估计效应量大小：rab.all 全部样本clr中值，带有处理的就是该组clr中值，带btw的就是处理见clr中位数，昂后win是clr最大差异中值
    # effect 效应量中值，overlap 效应量包含0的比例。
    # 显著差异的有哪些呢？可以使用Welchs t 检验矫正后的p值，wi.ep 秩和检验p值都小于0.05的
    # write.table(results, file="Aldex2.txt", quote=FALSE, sep='\t', col.names = NA)
    # res = list(tab.Aldex2 = results,
    #            diff.tab = data.frame(micro = tab.d1$OTU,method = tab.d1$group,adjust.p = tab.d1$p))
    str(results)
    data_list[[i]]= results
    names( data_list)[i] = Desep_group %>% paste( collapse = "_")
  }

  return(data_list)
}

