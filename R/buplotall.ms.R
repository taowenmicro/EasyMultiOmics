#' @title Bubble Plot for Metabolic Pathway Enrichment Analysis
#'
#' @description
#' The `buplotall.ms` function performs pairwise comparisons of groups in a metabolomics dataset and generates bubble plots to visualize enriched metabolic pathways. It integrates statistical testing and pathway enrichment analysis, allowing users to explore significant pathways across groups.
#'
#' @param ps A `phyloseq` object containing metabolomics data, including the compound abundance table and sample metadata.
#' @param dif.method A character string specifying the statistical method for differential testing. Default is `"wilcox"` (Wilcoxon rank-sum test).
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item `plots`: A list of ggplot objects, where each plot corresponds to a pairwise comparison of groups and visualizes enriched metabolic pathways.
#'   \item `plotdata`: A list of data frames, where each data frame contains detailed results of the pathway enrichment analysis for a specific group comparison.
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Conducts pairwise group comparisons using the specified statistical method (`dif.method`).
#'   \item Identifies significantly altered compounds based on adjusted p-values.
#'   \item Maps significant compounds to metabolic pathways using a preloaded pathway database (`db.compound.pathway.bins`).
#'   \item Conducts pathway enrichment analysis using the `clusterProfiler::enricher` function.
#'   \item Generates bubble plots to visualize enriched pathways, where the x-axis represents the pathway coverage, the y-axis represents pathway enrichment significance (`-log10(p-value)`), and bubble size corresponds to the count of significant compounds in each pathway.
#' }
#'
#' @examples
#' \dontrun{
#' res= buplotall.ms(ps= ps.ms3,dif.method = "wilcox")
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export


buplotall.ms =function (ps, dif.method = "wilcox") {
  result1 = statSuper(ps = ps, group = "Group", artGroup = NULL,
                      method = "wilcox")
  ps %>% sample_data()
  rank_names(ps)
  group = sample_data(ps)$Group %>% unique()
  com = combn(group, 2)
  res = list()
  plot = list()
  # i=1
  for (i in 1:dim(com)[2]) {
    a = colnames(result1)[str_detect(colnames(result1), "[fdr]")]
    b = colnames(result1)[str_detect(colnames(result1), com[2,
                                                            i])]
    c = colnames(result1)[str_detect(colnames(result1), com[1,
                                                            i])]
    d = intersect(a, b)
    e = intersect(d, c)
    tem = e
    n.var = as.name(tem)
    tem2 = result1 %>% dplyr::filter(!!n.var < 0.1) %>% .$Row.names
    tem2 = tem2[!is.na(tem2)]
    id.tem = tem2 %>% strsplit("[,]") %>% sapply(`[`, 1)
    library(EasyMultiOmics.db)
    dat = db.compound.pathway.bins
    head(dat)
    colnames(dat)[1] = "compound"
    # 数据库格式调整
    head(dat)
    dat = dat%>%
      separate_rows(pathway, sep = "\\|")


    tax = ps.ms3 %>% phyloseq::tax_table() %>% as.data.frame()
    head(tax)
    allkeggid <- data.frame(ID = row.names(tax))
    allkeggid$ID = allkeggid$ID %>% strsplit("[,]") %>% sapply(`[`,
                                                               1)
    head(allkeggid)
    total = dat%>%
      filter(compound %in% c(paste0("cpd:",allkeggid$ID))) %>%  dplyr::select(2, 1)
    head(total)
    dim(total)
    x <- clusterProfiler::enricher(gene = c(paste0("cpd:",
                                                   id.tem)), TERM2GENE = total, minGSSize = 1, pvalueCutoff = 1,
                                   qvalueCutoff = 1)
    dat = as.data.frame(x@result)
    df <- dat %>% arrange(desc(Count)) %>% separate(GeneRatio,
           into = c("A", "B"), sep = "/") %>% mutate(A = as.numeric(A),
           B = as.numeric(B)) %>% mutate(count = A/B) %>% arrange(Count)
    head(df)
    row.names(df) = NULL
    df$Description <- factor(df$Description, levels = c(df$Description %>%
                                                          as.data.frame() %>% pull()))
    colnames(df)
    df$p.adjust
    df$ID2 = df$ID %>% strsplit("[|]") %>% sapply(`[`, 1)



    df$ID2 = gsub("path:", "", df$ID2)
    tab = db.compound.pathway.dis
    head(tab)
    df2 = df %>% left_join(tab, by = c(ID2 = "V1")) %>% distinct(ID2,
                                                                 .keep_all = TRUE) %>% arrange(desc(count))
    head(df2)
    res[[i]] = df2
    names(res)[i] = paste(com[2, i], com[1, i], sep = ".")
    tem = df2 %>% filter(p.adjust < 0.05)
    head(df2)
    if (dim(tem)[1] == 0) {
      laby = "Top 20 pathways"
      df3 = df2 %>% arrange(desc(p.adjust)) %>% head(20)
    }
    else {
      laby = "enriched pathways"
      df3 = df2 %>% filter(p.adjust < 0.05)
    }


    print("plot")
    head(df3)
    df3$color = df3$count + (-log10(df3$pvalue))
    p1 <- ggplot(df3, aes(x =count, y = -log10(pvalue)))  +
      geom_point(pch = 21,aes(fill=color,size = count ))  +
      ggrepel::geom_text_repel(aes(x =count, y = -log10(pvalue),label= V2)) +
      scale_fill_gradientn(colours =colorRampPalette(c("#F7F4F9","#FFFF33","#FF7F00"))(60)) +
      theme_bw()+
      labs(x = NULL)
    # p=df3 %>% ggplot(aes(x = count, y = reorder(V2, count))) +
    #   geom_bar(aes(x = count, y = reorder(V2, count), fill = pvalue),
    #            stat = "identity", width = 0.8) + guides(size = guide_legend(title = "Count")) +
    #   labs(x = NULL, y = laby) +
    #   theme(
    #     # axis.text.x = element_text(color = "black",
    #     #   angle = 0, hjust = 0.5, vjust = 0.5, margin = margin(b = 5)),
    #     axis.text.y = element_text(color = "black", angle = 0,
    #                                hjust = 1, vjust = 0.5), panel.background = element_rect(fill = NA,
    #                                                                                         color = NA), panel.grid.minor = element_line(size = 0.2,
    #                                                                                                                                      color = "#e5e5e5"), panel.grid.major = element_line(size = 0.2,
    #                                                                                                                                                                                          color = "#e5e5e5"), panel.border = element_rect(fill = NA,                                                                                                                                                            color = "black", size = 1, linetype = "solid"),
    #     legend.key = element_blank(), legend.title = element_text(color = "black",
    #                                                               size = 9), legend.text = element_text(color = "black",
    #                                                                                                     size = 8), legend.spacing.x = unit(0.1, "cm"),
    #     legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.5,
    #                                                                  "cm"), legend.background = element_blank(), legend.box = "horizontal",
    #     legend.box.background = element_rect(color = "black"),
    #     legend.position = c(1, 0), legend.justification = c(1,
    #                                                         0)) +
    #   scale_y_discrete(labels = function(y) str_wrap(y,
    #
    #                                                  width = 30))
    # p
    #
    plot[[i]] = p1

    names(plot)[i] = paste(com[2, i], com[1, i], "plot",
                           sep = ".")
  }
  return(list(plots = plot, plotdata = res))
}
