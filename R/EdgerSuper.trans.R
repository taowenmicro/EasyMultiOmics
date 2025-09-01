#' @title Differential Abundance Analysis Using EdgeR for Metagenome Functional Data
#'
#' @description
#' This function performs differential abundance analysis for metagenome functional data using the `EdgeR` package.
#' It identifies differentially abundant taxa between groups in a `phyloseq` object, visualizes the results
#' with volcano plots, and outputs normalized abundance data with taxonomic annotations.
#'
#' @param otu A data frame containing metagenome functional counts. Optional if `ps` is provided.
#' @param tax A data frame containing taxonomic annotations. Optional if `ps` is provided.
#' @param map A data frame containing sample metadata. Optional if `ps` is provided.
#' @param tree A phylogenetic tree. Optional if `ps` is provided.
#' @param ps A `phyloseq` object containing metagenome functional data. If provided, overrides `otu`, `tax`, `map`, and `tree`.
#' @param group A string specifying the grouping variable in the sample metadata. Default is `"Group"`.
#' @param pvalue A numeric value specifying the significance threshold for adjusted p-values. Default is `0.05`.
#' @param lfc A numeric value specifying the log2 fold change cutoff. Default is `0`.
#' @param artGroup A matrix specifying custom group comparisons. Optional.
#' @param method A string specifying the normalization method for EdgeR. Default is `"TMM"`.
#' @param j A string or integer specifying the taxonomic rank to perform the analysis on. Metagenome functional data should select "meta".
#'
#' @return
#' A list containing:
#' \describe{
#'   \item{Plots}{A list of volcano plots for each group comparison, showing log fold changes and p-values.}
#'   \item{Results}{A data frame containing differential abundance results, including normalized abundances and taxonomic annotations.}
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Prepares metagenome functional data, taxonomic annotations, and metadata from a `phyloseq` object or individual input files.
#'   \item Aggregates data to the specified taxonomic rank (if applicable).
#'   \item Normalizes metagenome functional counts using the specified normalization method (`"TMM"` by default).
#'   \item Constructs a model matrix for group comparisons and estimates dispersions using EdgeR's generalized linear model (GLM) framework.
#'   \item Identifies differentially abundant taxa based on the specified `pvalue` and `lfc` thresholds.
#'   \item Generates volcano plots for each group comparison, highlighting enriched and depleted taxa.
#'   \item Appends taxonomic annotations and normalized abundances to the differential abundance results.
#' }
#'
#' The function supports both predefined and custom group comparisons (`artGroup`) and provides detailed visualizations for the results.
#'
#' @examples
#' \dontrun{
#' res = EdgerSuper.trans(ps = ps.kegg %>%
#'                        ggClusterNet::filter_OTU_ps(Top=500),
#'                        group = "Group",
#'                        artGroup = NULL,
#'                        j = "meta")
#' p25.1 = res[[1]][1]
#' p25.1
#' p25.2 = res[[1]][2]
#' p25.2
#' p25.3 = res[[1]][3]
#' p25.3
#' dat = res[[2]]
#' head(dat)
#' }
#'
#' @author Contact: Tao Wen \email{2018203048@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#'
#' @export

EdgerSuper.trans=function (otu = NULL, tax = NULL, map = NULL, tree = NULL, ps = NULL,
                          group = "Group", pvalue = 0.05, lfc = 0, artGroup = NULL,
                          method = "TMM", j = "meta")
{
  ps = ggClusterNet::inputMicro(otu, tax, map, tree, ps, group = group)
  ps
  if (j %in% c("OTU", "gene", "meta")) {
    ps = ps
  }
  else if (j %in% c(1:7)) {
    ps = ps %>% ggClusterNet::tax_glom_wt(ranks = j)
  }
  else if (j %in% c("Kingdom", "Phylum", "Class", "Order",
                    "Family", "Genus", "Species")) {
    ps = ps %>% ggClusterNet::tax_glom_wt(ranks = j)
  }
  else {
    ps = ps
    print("unknown j, checked please")
  }
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  colnames(sub_design) = "Group"
  Desep_group <- as.character(levels(as.factor(sub_design$Group)))
  Desep_group
  if (is.null(artGroup)) {
    aaa = combn(Desep_group, 2)
  }
  if (!is.null(artGroup)) {
    aaa = as.matrix(b)
  }
  otu_table = as.data.frame(ggClusterNet::vegan_otu(ps))
  count = as.matrix(otu_table)
  count <- t(count)
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  dim(sub_design)
  sub_design$SampleType = as.character(sub_design$Group)
  sub_design$SampleType <- as.factor(sub_design$Group)
  d = edgeR::DGEList(counts = count, group = sub_design$SampleType)
  d$samples
  d = edgeR::calcNormFactors(d, method = method)
  design.mat = model.matrix(~0 + d$samples$group)
  colnames(design.mat) = levels(sub_design$SampleType)
  d2 = edgeR::estimateGLMCommonDisp(d, design.mat)
  d2 = edgeR::estimateGLMTagwiseDisp(d2, design.mat)
  fit = edgeR::glmFit(d2, design.mat)
  plot_list <- list()
  for (i in 1:dim(aaa)[2]) {
    Desep_group = aaa[, i]
    print(Desep_group)
    group = paste(Desep_group[1], Desep_group[2], sep = "-")
    group
    BvsA <- limma::makeContrasts(contrasts = group, levels = design.mat)
    lrt = edgeR::glmLRT(fit, contrast = BvsA)
    de_lrt = edgeR::decideTestsDGE(lrt, adjust.method = "fdr",
                                   p.value = pvalue, lfc = lfc)
    summary(de_lrt)
    x = lrt$table
    x$sig = de_lrt
    head(x)
    row.names(count)[1:6]
    x <- cbind(x, padj = p.adjust(x$PValue, method = "fdr"))
    enriched = row.names(subset(x, sig == 1))
    depleted = row.names(subset(x, sig == -1))
    x$level = as.factor(ifelse(as.vector(x$sig) == 1, "enriched",
                               ifelse(as.vector(x$sig) == -1, "depleted", "nosig")))
    x = data.frame(row.names = row.names(x), logFC = x$logFC,
                   level = x$level, p = x$PValue)
    head(x)
    x1 = x %>% dplyr::filter(level %in% c("enriched", "depleted",
                                          "nosig"))
    head(x1)
    x1$Genus = row.names(x1)
    if (nrow(x1) <= 1) {
    }
    x2 <- x1 %>% dplyr::mutate(ord = logFC^2) %>% dplyr::filter(level !=
                                                                  "nosig") %>% dplyr::arrange(desc(ord)) %>% head(n = 5)
    x3 <- x1 %>% dplyr::mutate(ord = logFC^2) %>% dplyr::filter(level !=
                                                                  "nosig") %>% dplyr::arrange(desc(ord))
    head(x2)
    p <- ggplot(x1, aes(x = logFC, y = -log2(p), colour = level)) +
      geom_point() + geom_hline(yintercept = -log10(0.2),
                                linetype = 4, color = "black", size = 0.5) + geom_vline(xintercept = c(-1,
                                                                                                       1), linetype = 3, color = "black", size = 0.5) +
      ggrepel::geom_text_repel(data = x2, aes(x = logFC,
                                              y = -log2(p), label = Genus), size = 1) + scale_color_manual(values = c("blue2",
                                                                                                                      "red2", "gray30")) + ggtitle(group) + theme_bw()
    p
    plot_list[[i]] <- p
    colnames(x) = paste(group, colnames(x), sep = "")
    if (i == 1) {
      table = x
    }
    if (i != 1) {
      table = cbind(table, x)
    }
  }
  x = table
  count = as.matrix(count)
  norm = t(t(count)/colSums(count))
  dim(norm)
  norm1 = norm %>% t() %>% as.data.frame()
  library("tidyverse")
  head(norm1)
  iris.split <- split(norm1, as.factor(sub_design$SampleType))
  iris.apply <- lapply(iris.split, function(x) colMeans(x))
  iris.combine <- do.call(rbind, iris.apply)
  norm2 = t(iris.combine)
  str(norm2)
  norm2 = as.data.frame(norm2)
  head(norm2)
  x = cbind(x, norm2)
  head(x)
  return(list(plot_list, x))
}







