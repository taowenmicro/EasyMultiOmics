
#' This function performs an analysis comparing two groups based on microbiome (16S rRNA) and metabolomics data.
#' It identifies differentially abundant operational taxonomic units (OTUs) and metabolites, then compares
#' their association through statistical tests. The function visualizes the comparisons using a scatter plot.
#'
#' @param ps01 A phyloseq object containing microbiome data (default: `ps.16s`).
#' @param ps.ms A data frame or tibble containing metabolomics data (default: `ps.ms`).
#' @param group1 A string specifying the first group for comparison (default: `"micro"`).
#' @param group2 A string specifying the second group for comparison (default: `"ms"`).
#' @param r.threshold A numeric value specifying the correlation threshold for filtering the microbiome data (default: 0.6).
#' @param p.threshold A numeric value specifying the p-value threshold for significance in statistical tests (default: 0.1).
#' @param method A string specifying the statistical test method (e.g., `"spearman"`) (default: `"spearman"`).
#' @param top An integer specifying the number of top OTUs or metabolites to consider for analysis (default: 500).
#'
#' @return A list containing:
#' \item{dat}{A data frame with the merged OTU and metabolite data, with statistical results and category assignments.}
#' \item{plot}{A ggplot object visualizing the comparisons between the two groups.}
#' @examples
#' results <- quare.line.omics(ps.16s =ps.16s ,
#'                             ps.ms = ps.ms,
#'                             group1 = "microbe",
#'                             group2 = "metabolite",
#'                             r.threshold = 0.7,
#'                             p.threshold = 0.05,
#'                             method = "spearman",
#'                             top = 100)
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
quare.line.omics <- function(ps.16s=ps.16s,
                             ps.ms = ps.ms,
                             group1= "micro",
                             group2= "ms",
                             r.threshold = 0.6,
                             p.threshold = 0.1,
                             method = "spearman",
                             top = 500){

res = EdgerSuper.micro(ps = ps.16s %>% ggClusterNet::filter_OTU_ps(top),group  = "Group",
                       artGroup = NULL, j = "OTU")

da16 = res[[2]]
da16$id =  row.names(da16)
row.names(da16)= NULL
# 差异代谢物
dams = statSuper(ps = ps.ms,group  = "Group",artGroup = NULL,method = "wilcox")
head(dams)
head(da16)
# 提取分组
sub_design <- as.data.frame(sample_data(ps.16s))
Desep_group <- as.character(levels(as.factor(sub_design$Group)))
Desep_group
aaa = combn(Desep_group,2)

results <- list()

for (ii in 1:dim(aaa)[2]) {
  # ii = 1
  Desep_group = aaa[,ii]
  print( Desep_group)

  col_16 <- paste(paste(Desep_group[1], Desep_group[2], sep = "-"), c("logFC", "p"), sep = "")
  # id -NULL
  # 筛选ko富集的
  da16_1 =  da16 %>% dplyr::select(col_16,id) %>%
   #  dplyr::filter( .data[[col_16[1]]] >0 )  %>%
    mutate(group1=group1)
  head(da16_1)
  da16_1[is.na(da16_1)] = 0

  col_ms <- paste(Desep_group[1], Desep_group[2], c("log2_FC", "Pvalue"), sep = "_")
  head(dams)
  # 筛选ko富集的
  dams_1 =  dams%>% dplyr::select(col_ms,id) %>%
    #dplyr::filter( .data[[col_ms[1]]] >0 ) %>%
    mutate(group1=group2)
  dams_1[is.na(dams_1)] = 0
  dams_1[[col_ms[1]]] [dams_1[[col_ms[1]]] == "NA"] = 0
  dams_1[[col_ms[1]]] = as.numeric(dams_1[[col_ms[1]]])
 #  dams_1[[col_ms[1]]]=  - (dams_1[[col_ms[1]]])
  head(dams_1)

  dams_1[[col_ms[2]]] [dams_1[[col_ms[2]]] == "NA"] = 1
  dams_1[[col_ms[2]]] = as.numeric(dams_1[[col_ms[2]]])


  head(da16_1)
  da16_1$id2= 1:nrow(da16_1)
  dams_1$id2=1:nrow(dams_1)

  dat = inner_join(da16_1,dams_1, by ="id2")

  dat$group <- case_when(
    abs(dat[, col_16[1]]) >= 1 & abs(dat[, col_ms[1]]) >= 1 ~ "both_enrich",
    abs(dat[, col_16[1]]) < 1 & abs(dat[, col_ms[1]]) > 1 ~ paste(group2,"enrich",sep = "_"),
    abs(dat[, col_16[1]]) > 1 & abs(dat[, col_ms[1]]) < 1 ~ paste(group1,"enrich",sep = "_"),
    abs(dat[, col_16[1]]) < 1 & abs(dat[, col_ms[1]]) < 1 ~ "both_deplete"
  )

  # dat$group <- case_when(
  #   dat[, col_16[1]] >= 1 & dat[, col_ms[1]] >= 1 ~ "both_enrich",
  #   abs(dat[, col_16[1]]) < 1 & abs(dat[, col_ms[1]]) > 1 ~ paste(group2,"enrich",sep = "_"),
  #   abs(dat[, col_16[1]]) > 1 & abs(dat[, col_ms[1]]) < 1 ~ paste(group1,"enrich",sep = "_"),
  #   abs(dat[, col_16[1]]) < 1 & abs(dat[, col_ms[1]]) < 1 ~ "both_deplete"
  # )
   head(dat)

 p = ggplot(dat, aes( !!sym(col_16[1]), !!sym(col_ms[1]), color = group)) +
    geom_point(size = 1.2, alpha =0.5) +
   #guides(color = "none") +
   # scale_colour_manual(name = "", values = microbiome::alpha(mycolor, 0.7)) +
    geom_hline(yintercept = c(-1, 1), size = 0.7, color = "grey40", lty = "dashed") +
    geom_vline(xintercept = c(-1, 1), size = 0.7, color = "grey40", lty = "dashed") +
    scale_y_continuous(expand = expansion(add = c(0.5, 0.5)), limits = c(-5, 5)) +
   scale_x_continuous(expand = expansion(add = c(0.5, 0.5)), limits = c(-5, 5))+
    labs(x =paste(group1,"LogFC",sep="_") ,
         y = paste(group2,"LogFC",sep="_"),
         title = paste(Desep_group, collapse = "_"))+
    theme_bw()+
    theme(panel.grid=element_blank(),axis.text=element_blank())+
    theme(legend.position="bottom")

p

results[[paste(Desep_group, collapse = "_")]] <-dat
results[[paste(Desep_group, collapse = "_")]] <- p
}
return(results)
}
