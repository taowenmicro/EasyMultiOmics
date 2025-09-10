#' @title Generate Manhattan Plots for Differential Abundance Analysis
#'
#' @description
#' The `edge_Manhattan.metm` function performs pairwise differential abundance analysis using the `edgeR` package.
#' It visualizes the results for each comparison in Manhattan plots, highlighting enriched, depleted, and nonsignificant taxa
#' with taxonomic annotations.
#'
#' @param ps A `phyloseq` object containing microbiome data (OTU table, taxonomic table, and sample metadata).
#' @param pvalue A numeric value specifying the significance threshold for adjusted p-values (FDR). Default is `0.05`.
#' @param lfc A numeric value specifying the log fold change threshold for significance. Default is `0`.
#'
#' @return
#' A list of Manhattan plots, one for each pairwise group comparison, showing:
#' \itemize{
#'   \item Log-transformed p-values (`-log(p)`) on the y-axis.
#'   \item OTUs (taxa) on the x-axis, arranged by Phylum.
#'   \item Different shapes and colors indicating significant enrichment, depletion, or nonsignificance.
#'   \item Phylum-specific colors and point sizes proportional to log-transformed counts per million (logCPM).
#' }
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Extracts OTU and sample data from the `phyloseq` object.
#'   \item Normalizes the OTU counts using the `edgeR` TMM method.
#'   \item Performs pairwise group comparisons using contrasts in a generalized linear model (GLM) framework.
#'   \item Identifies significantly enriched and depleted OTUs based on FDR-adjusted p-values and log fold changes.
#'   \item Annotates OTUs with taxonomic information.
#'   \item Generates Manhattan plots for each pairwise comparison.
#' }
#'
#' Top phyla are highlighted, while taxa from less abundant phyla are categorized as "Low Abundance."
#' Each Manhattan plot contains a horizontal line representing the minimum significant FDR value.
#'
#' @examples
#' \dontrun{
#' res = edge_Manhattan.metm(ps = ps.16s%>% ggClusterNet::filter_OTU_ps(500),pvalue = 0.05,lfc = 0)
#'  }
#'
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
#' @export


edge_Manhattan.metm=function (ps = ps, pvalue = 0.05, lfc = 0)
{
  count = ps %>% ggClusterNet::vegan_otu() %>% t()
  d = edgeR::DGEList(counts = count, group = phyloseq::sample_data(ps)$Group)
  d = edgeR::calcNormFactors(d)
  design.mat = model.matrix(~0 + d$samples$group)
  colnames(design.mat) = levels(as.factor(phyloseq::sample_data(ps)$Group))
  d2 = edgeR::estimateGLMCommonDisp(d, design.mat)
  d2 = edgeR::estimateGLMTagwiseDisp(d2, design.mat)
  fit = edgeR::glmFit(d2, design.mat)
  Desep_group <- as.character(levels(as.factor(phyloseq::sample_data(ps)$Group)))
  Desep_group
  aaa = combn(Desep_group, 2)
  plot_list <- list()
  # packageVersion("edgeR")
  for (i in 1:dim(aaa)[2]) {
    Desep_group = aaa[, i]
    print(Desep_group)
    group = paste(Desep_group[1], Desep_group[2], sep = "-")
    group
    BvsA <- limma::makeContrasts(contrasts = group, levels = design.mat)
    lrt = edgeR::glmLRT(fit, contrast = BvsA)
    de_lrt = limma::decideTests(lrt, adjust.method = "fdr",
                                   p.value = pvalue, lfc = lfc)
    summary(de_lrt)
    x = lrt$table
    x$sig = de_lrt
    x$sig = de_lrt
    head(x)
    row.names(count)[1:6]
    x <- cbind(x, padj = p.adjust(x$PValue, method = "fdr"))
    enriched = row.names(subset(x, sig == 1))
    depleted = row.names(subset(x, sig == -1))
    x$level = as.factor(ifelse(as.vector(x$sig) == 1, "enriched",
                               ifelse(as.vector(x$sig) == -1, "depleted", "nosig")))
    x$otu = rownames(x)
    x$neglogp = -log(x$PValue)
    tax = ps %>% ggClusterNet::vegan_tax() %>% as.data.frame()
    head(tax)
    x = cbind(x, tax)
    head(x)
    top_phylum = c("Bacteroidetes", "Firmicutes", "Planctomycetes",
                   "Proteobacteria", "Verrucomicrobia")
    x[!(x$Phylum %in% top_phylum), ]$Phylum = "Low Abundance"
    x$otu = factor(x$otu, levels = x$otu)
    x$level = factor(x$level, levels = c("enriched", "depleted",
                                         "nosig"))
    levels(x$Phylum) = c(top_phylum, "Low Abundance")
    x = x %>% arrange(Phylum)
    head(x)
    x$otu = factor(x$otu, levels = x$otu)
    FDR = min(x$neglogp[x$level == "depleted"])
    p = ggplot(x, aes(x = otu, y = neglogp, color = Phylum,
                      size = logCPM, shape = level)) + geom_point(alpha = 0.7) +
      geom_hline(yintercept = FDR, linetype = 2, color = "lightgrey") +
      scale_shape_manual(values = c(17, 25, 20)) + scale_size(breaks = c(5,
                                                                         10, 15)) + labs(x = "OTU", y = "-loge(P)") + theme(axis.ticks.x = element_blank(),
                                                                                                                            axis.text.x = element_blank(), legend.position = "top") +
      scale_color_manual(values = c(RColorBrewer::brewer.pal(9,
                                                             "Set1"))) + ggtitle(group)
    p
    plot_list[[i]] <- p
  }
  return(plot_list)
}
