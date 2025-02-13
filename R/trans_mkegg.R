#' @title KEGG Pathway Transformation using Mkegg
#' @description
#' This function downloads and processes KEGG pathway data using the `Mkegg` option
#' for KEGG analysis from `clusterProfiler`. It maps KEGG Orthology (KO) terms to
#' KEGG pathways and integrates the data with OTU tables from a Phyloseq object.
#' The result is a transformed Phyloseq object with pathway-level information.
#'
#' @param ps.trans A Phyloseq object containing OTU table and taxonomic data.
#' The default is `ps.trans` which assumes that the Phyloseq object is provided externally.
#'
#' @return A transformed Phyloseq object where the OTU table is aggregated by KEGG pathway IDs
#' and taxonomic information is updated with KEGG pathway descriptions.
#'
#' @export
#' @examples
#' \dontrun{
#' ps_transformed <- trans_mkegg(ps.trans)
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
trans_mkegg = function(ps.trans =  ps.trans) {

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

ko =  ps.trans %>% vegan_otu() %>% t() %>%
  as.data.frame() %>% rownames_to_column("KO")
head(ko)


tem2 = ko %>% left_join(PATH_ID_NAME,by = c("KO" = "KO")) %>%
  tidyfst::filter_dt(KEGGID != "") %>%
  tidyfst::summarise_vars(is.numeric,sum,by =c("KEGGID","DESCRPTION")) %>%
  as.data.frame()
head(tem2)
colnames(tem2)

otu = tem2[,sample_names( ps.trans)]
tax = data.frame(MDESCRPTION = tem2$DESCRPTION,MDESCRPTION2 = tem2$DESCRPTION,
                 row.names = tem2$KEGGID)

row.names(otu) = tem2$KEGGID
dim(otu)

ps2 = phyloseq(
  phyloseq:: otu_table(as.matrix(otu),taxa_are_rows = TRUE),
  phyloseq::tax_table(as.matrix(tax)),
  phyloseq::sample_data( ps.trans)
)

ps2

return(ps2)
}
