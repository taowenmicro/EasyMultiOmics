#' @title Perform LDA analysis on metagenome functional composition data
#' @description
#' This function performs LDA analysis on metagenome functional composition data to
#' identify characteristic functions at different taxonomic levels.
#' @param ps A phyloseq format file used as an alternative for the input containing metagenome functional composition table, tax, and map.
#' @param group Column name for groupID in map table(sample metadata).
#' @param p.lvl p-value threshold for the characteristic functions.
#' @param lda.lvl  LDA score threshold for the characteristic functions.
#' @param seed The random seed for reproducibility.
#' @param adjust.p Logical indicating whether to adjust p-values.
#' @param Top The top functions to convert and visualize.
#' @returns A list object containing the following components:
#' \item{lefse_lists}{Data frame containing LDA analysis results.}
#' \item{taxtree}{Data frame containing taxonomic information and LDA results.}
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' ps =ps.cnps%>% filter_OTU_ps(Top = 1000)
#' tablda = LDA.metf(ps = ps,
#' Top = 100,
#' p.lvl = 0.05,
#' lda.lvl = 1,
#' seed = 11,
#' adjust.p = F)
#' dat1=tablda[[1]]
#' head(dat1)
#' dat2 = tablda[[2]]
#' head(dat2)
LDA.metf=function (ps = ps.cs, group = "Group", Top = 100, p.lvl = 0.05,
                   lda.lvl = 2, seed = 11, adjust.p = F)
{
  ps = ps %>% filter_taxa(function(x) sum(x) > 0, TRUE)
  alltax = ps %>% ggClusterNet::filter_OTU_ps(Top) %>% ggClusterNet::vegan_tax() %>%
    as.data.frame()
  alltax$OTU = row.names(alltax)
  otu = ps %>% ggClusterNet::filter_OTU_ps(Top) %>% ggClusterNet::vegan_otu() %>%
    t() %>% as.data.frame()
  otu_tax = merge(otu, alltax, by = "row.names", all = F)
  head(otu_tax)
  rank <- otu_tax %>% dplyr::group_by(OTU) %>% dplyr::summarise_if(is.numeric,
                                                                   sum, na.rm = TRUE)
  colnames(rank)[1] = "id"
  data1 = as.data.frame(rank)
  row.names(data1) = data1$id
  data1$id = NULL
  ps_G_graphlan = phyloseq::phyloseq(phyloseq::otu_table(as.matrix(data1),
                                                         taxa_are_rows = TRUE), phyloseq::sample_data(ps)) %>%
    filter_taxa(function(x) sum(x) > 0, TRUE)
  ps_G_graphlan
  otu = as.data.frame((ggClusterNet::vegan_otu(ps_G_graphlan)))
  otu[otu == 0] <- 1
  otu = otu[colMeans(otu) != 1]
  map = as.data.frame(phyloseq::sample_data(ps_G_graphlan))
  claslbl = map[, group] %>% as.vector() %>% .[[1]] %>% as.factor()
  set.seed(seed)
  rawpvalues <- apply(otu, 2, function(x) kruskal.test(x, claslbl)$p.value)
  ord.inx <- order(rawpvalues)
  rawpvalues <- rawpvalues[ord.inx]
  clapvalues <- p.adjust(rawpvalues, method = "fdr")
  wil_datadf <- as.data.frame(otu[, ord.inx])
  ldares <- MASS::lda(claslbl ~ ., data = wil_datadf)
  ldamean <- as.data.frame(t(ldares$means))
  ldamean
  class_no <<- length(unique(claslbl))
  ldamean$max <- apply(ldamean[, 1:class_no], 1, max)
  ldamean$min <- apply(ldamean[, 1:class_no], 1, min)
  ldamean$LDAscore <- signif(log10(1 + abs(ldamean$max - ldamean$min)/2),
                             digits = 3)
  head(ldamean)
  a = rep("A", length(ldamean$max))
  for (i in 1:length(ldamean$max)) {
    name = colnames(ldamean[, 1:class_no])
    a[i] = name[ldamean[, 1:class_no][i, ] %in% ldamean$max[i]]
  }
  ldamean$class = a
  tem1 = row.names(ldamean)
  tem1 %>% as.character()
  ldamean$Pvalues <- signif(rawpvalues[match(row.names(ldamean),
                                             names(rawpvalues))], digits = 5)
  ldamean$FDR <- signif(clapvalues, digits = 5)
  resTable <- ldamean
  rawNms <- rownames(resTable)
  rownames(resTable) <- gsub("`", "", rawNms)
  if (adjust.p) {
    de.Num <- sum(clapvalues <= p.lvl & ldamean$LDAscore >=
                    lda.lvl)
  }
  else {
    de.Num <- sum(rawpvalues <= p.lvl & ldamean$LDAscore >=
                    lda.lvl)
  }
  if (de.Num == 0) {
    current.msg <<- "No significant features were identified with given criteria."
  }
  else {
    current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
  }
  print(current.msg)
  ord.inx <- order(resTable$Pvalues, resTable$LDAscore)
  resTable <- resTable[ord.inx, , drop = FALSE]
  resTable <- resTable[, c(ncol(resTable), 1:(ncol(resTable) -
                                                1))]
  resTable <- resTable[, c(ncol(resTable), 1:(ncol(resTable) -
                                                1))]
  ldamean$Pvalues[is.na(ldamean$Pvalues)] = 1
  if (adjust.p) {
    taxtree = resTable[clapvalues <= p.lvl & ldamean$LDAscore >=
                         lda.lvl, ]
  }
  else {
    taxtree = resTable[ldamean$Pvalues <= p.lvl, ]
  }
  colour = c("darkgreen", "red", "blue", "#4DAF4A", "#984EA3",
             "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
  selececol = colour[1:length(levels(as.factor(taxtree$class)))]
  names(selececol) = levels(as.factor(taxtree$class))
  A = rep("a", length(row.names(taxtree)))
  for (i in 1:length(row.names(taxtree))) {
    A[i] = selececol[taxtree$class[i]]
  }
  taxtree$color = A
  lefse_lists = data.frame(node = row.names(taxtree), color = A,
                           Group = taxtree$class, stringsAsFactors = FALSE)
  return(list(lefse_lists, taxtree))
}
