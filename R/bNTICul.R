
# result = bNTICul(ps = ps ,group  = "Group",num = 3,thread = 1)
#
# bNTI = result[[1]]
# head(bNTI)
#
# path = "./Result/bNTI"
# dir.create(path)
# filename = paste(path,"/bNTI.csv",sep = "")
# write.csv(bNTI, filename)


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
