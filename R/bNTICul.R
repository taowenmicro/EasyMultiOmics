
#' @title Calculate βNTI for Phylogenetic Turnover
#'
#' @description
#' This function calculates the beta nearest taxon index (βNTI) for assessing phylogenetic turnover in microbial communities.
#' It compares observed βMNTD (beta mean nearest taxon distance) to null model expectations and standardizes the difference.
#'
#' @param otu A data frame or matrix containing OTU (Operational Taxonomic Unit) abundance data.
#' Rows represent taxa, and columns represent samples.
#' @param tax A data frame containing taxonomic annotation for each OTU.
#' @param map A data frame containing sample metadata, including group information.
#' @param tree A phylogenetic tree object.
#' @param ps A `phyloseq` object containing OTU, taxonomic, and sample data.
#' If provided, this supersedes `otu`, `tax`, `map`, and `tree`.
#' @param group A character string specifying the grouping variable in the metadata. Default is `"Group"`.
#' @param num An integer specifying the number of randomizations for the null model. Default is `99`.
#' @param thread An integer specifying the number of processor threads to use for parallel computation. Default is `1`.
#'
#' @return A list containing:
#' \itemize{
#'   \item `bNTI`: A data frame with calculated βNTI values, including observed βMNTD, mean and standard deviation of permuted βMNTD, and βNTI for each sample pair.
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Rarefies the input OTU table to standardize sample sequencing depth.
#'   \item Computes the observed βMNTD between sample pairs.
#'   \item Generates null communities by randomizing the phylogenetic tree and computes null model βMNTD values.
#'   \item Calculates βNTI as the standardized difference between observed and null model βMNTD.
#' }
#' This function supports parallel computation using multiple threads to reduce computation time for null model simulations.
#'
#' @examples
#' \dontrun{
#' result = bNTICul(ps = psphy,group  = "Group",num = 10,thread = 1)
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export
bNTICul = function(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = NULL,group  = "Group",num = 99,thread = 1){
  ps = inputMicro(otu,tax,map,tree,ps,group  = group)
  ps

  ps_sub <- ps

  # tree = phy_tree(ps)
  # tree
  #-------------调整map文件-----------------------------------------------------------------
  #添加一个ID列
  map = as.data.frame(sample_data(ps_sub))
  map$ID = row.names(map)
  sample_data(ps) = map


  #-----------准备OTU表格---------------------抽平-不设置抽平条数，默认按照最小序列数数目抽平
  set.seed(72)  # setting seed for reproducibility
  psrare = rarefy_even_depth(ps_sub)
  #检查序列数量
  sample_sums(psrare)
  # 标准化数据
  ps.norm = transform_sample_counts(psrare, function(x) x/sum(x))


  # ## -0--------------对样本序列数量可视化
  # read_count = data.frame("count" = colSums(otu_table(ps))) %>%
  #   rownames_to_column(var="ID") %>%
  #   inner_join(data.frame(sample_data(ps)), by="ID") %>%
  #   arrange(-count) %>%
  #   mutate(ID=factor(ID, levels=ID))
  #
  # # Now plot read count for each sample. The horizontal line represents a 2000 read threshold
  # ggplot(data=read_count, aes(x=ID, y=log10(count), fill=Group)) +
  #   geom_bar(stat="identity") +
  #   labs(x="Sample", y="Log10(Read count)") +
  #   geom_hline(yintercept=log10(10000)) +
  #   theme(text = element_text(size=16),
  #         axis.text.x = element_blank())
  # # Everything seems to be at or above 10000 total reads
  #
  # ps





  # 计算βMNTD对每个随机零模型群落
  bMNTD_null_func <- function(i, OTU.table, tree){
    tree$tip.label = sample(tree$tip.label)
    bMNTD_s = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
    A <- attr(bMNTD_s, "Size")
    B <- if (is.null(attr(bMNTD_s, "Labels"))) sequence(A) else attr(bMNTD_s, "Labels")
    if (isTRUE(attr(bMNTD_s, "Diag"))) attr(bMNTD_s, "Diag") <- FALSE
    if (isTRUE(attr(bMNTD_s, "Upper"))) attr(bMNTD_s, "Upper") <- FALSE
    bMNTD_s.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                            Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                            bMNTD = as.vector(bMNTD_s),
                            rep=i)
    return(bMNTD_s.df)
  }
  # 计算βNTI
  Phylo_turnover <- function(physeq, reps, nproc){
    # Extract OTU table
    OTU.table = t(otu_table(physeq))
    # Extract phylogenetic tree
    tree = phy_tree(physeq)

    # Get βMNTD between all communities
    bMNTD_o = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
    A <- attr(bMNTD_o, "Size")
    B <- if (is.null(attr(bMNTD_o, "Labels"))) sequence(A) else attr(bMNTD_o, "Labels")
    if (isTRUE(attr(bMNTD_o, "Diag"))) attr(bMNTD_o, "Diag") <- FALSE
    if (isTRUE(attr(bMNTD_o, "Upper"))) attr(bMNTD_o, "Upper") <- FALSE
    bMNTD_o.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                            Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                            bMNTD = as.vector(bMNTD_o))

    # Get βMNTD for randomized null communities
    rep.list = seq(1, reps)
    bMNTD_s.df.list = mclapply(rep.list, bMNTD_null_func, OTU.table=OTU.table, tree=tree, mc.cores=nproc)

    # Combine all data together and calculate βNTI for each sample pair
    bMNTD_s.df <- do.call("rbind", bMNTD_s.df.list)
    bMNTD_s.means.df = bMNTD_s.df %>%
      group_by(Sample_1, Sample_2) %>%
      dplyr::summarize(mean_bMNTD = mean(bMNTD),
                       sd_bMNTD = sd(bMNTD))

    bMNTD_o.df = inner_join(bMNTD_o.df, bMNTD_s.means.df, by=c("Sample_1", "Sample_2")) %>%
      mutate(bNTI = (bMNTD - mean_bMNTD)/sd_bMNTD)

    return(bMNTD_o.df)
  }







  #========这里一把单核就真实数据而言需要超过10个小时，跑999次，所以需要多核
  # 计算bnti，这里可以设置线程数量，是第三个参数，我们在linux下面可以设置，30个线程
  # 第二个参数设置迭代数量，这里文献一般999嘛。
  bNTI = Phylo_turnover(psrare, num, thread)

  return(list(bNTI))
}
