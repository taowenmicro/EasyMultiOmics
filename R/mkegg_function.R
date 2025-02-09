#' @title Reaction-Based KEGG Pathway Analysis for Phyloseq Object
#' @description
#'
#' This function processes an OTU table from a `phyloseq` object, maps KO identifiers to KEGG reactions,
#' and creates a new `phyloseq` object that associates reactions with the OTU data.
#'
#' @param ps A `phyloseq` object containing OTU data to be transformed into KEGG reactions.
#' @param dat A data frame containing KO identifiers and associated reactions. If `NULL`, defaults to `EasyMultiOmics.db::db.ko.reaction`.
#' @param dat2 A data frame containing reactions and descriptions. If `NULL`, defaults to `EasyMultiOmics.db::db.reaction`.
#'
#' @return A `phyloseq` object with summarized OTU data and KEGG reaction information.
#' @details
#' This function performs several steps:
#' \itemize{
#'   \item Loads KO-reaction data from `EasyMultiOmics.db::db.ko.reaction` by default or the provided `dat`.
#'   \item Transforms the OTU table in `ps` to KO identifiers using `vegan_otu()`.
#'   \item Joins KO identifiers with KEGG reaction data.
#'   \item Summarizes the data by reaction and creates a new `phyloseq` object with summarized OTU data.
#'   \item Adds reaction descriptions from `EasyMultiOmics.db::db.reaction` or the provided `dat2`.
#' }
#' @examples
#' # Assuming you have a phyloseq object `ps`
#' reaction_ps <- reaction_function(ps)
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
mkegg_function= function(ps = ps) {

#匹配kegg基因名称
getOption("clusterProfiler.download.method")
R.utils::setOption( "clusterProfiler.download.method",'auto' )
#--基于Mkegg进行分析
Mkegg = clusterProfiler::download_KEGG('ko',keggType = "MKEGG")
PATH2ID <- Mkegg$KEGGPATHID2EXTID
PATH2NAME <- Mkegg$KEGGPATHID2NAME
head(PATH2NAME)
PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by='from')
colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRPTION')
head(PATH_ID_NAME )


ko = ps %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)


tem2 = ko %>% left_join(PATH_ID_NAME,by = c("KO" = "KO")) %>%
  tidyfst::filter_dt(KEGGID != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("KEGGID","DESCRPTION")) %>%
  as.data.frame()
head(tem2)
colnames(tem2)

otu = tem2[,sample_names(ps)]
tax = data.frame(MDESCRPTION = tem2$DESCRPTION,MDESCRPTION2 = tem2$DESCRPTION,row.names = tem2$KEGGID)
# tem2$DESCRPTION = gsub("[,]",".",tem2$DESCRPTION)
# tem2$DESCRPTION = gsub("[-]",".",tem2$DESCRPTION)
# tem2$DESCRPTION = gsub("[+]",".",tem2$DESCRPTION)
# tem2$DESCRPTION = gsub("[/]",".",tem2$DESCRPTION)
# tem2$DESCRPTION = gsub("[(]",".",tem2$DESCRPTION)
# tem2$DESCRPTION = gsub("[)]",".",tem2$DESCRPTION)
# tem2$DESCRPTION = gsub("[=>]",".",tem2$DESCRPTION)
# tem2$DESCRPTION = gsub("[ ]","_",tem2$DESCRPTION)
# tem2$DESCRPTION %>% unique()

row.names(otu) = tem2$KEGGID
head(otu)
ps1 = phyloseq(
  phyloseq::otu_table(as.matrix(otu),taxa_are_rows = TRUE),
 phyloseq:: tax_table(as.matrix(tax)),
 phyloseq:: sample_data(ps)
)
return(ps1)
}
