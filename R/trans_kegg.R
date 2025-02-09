#' @title KEGG Pathway Transformation for Phyloseq Data
#' @description
#'
#' This function performs a transformation of the given Phyloseq object by mapping
#' the KEGG Orthology (KO) terms to KEGG pathways. It uses data from the KEGG
#' database and merges the KO terms with pathway descriptions. It processes the OTU
#' (Operational Taxonomic Unit) data and integrates it with KEGG pathway information,
#' returning a transformed Phyloseq object.
#'
#' @param ps.trans A Phyloseq object containing OTU table and associated taxonomic data.
#' The default is `ps.trans` which assumes that the Phyloseq object is provided externally.
#'
#' @return A transformed Phyloseq object with updated taxonomic information based on KEGG pathways.
#' The OTU table is aggregated by KEGG pathway IDs and the taxonomic table is updated to include pathway descriptions.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ps_transformed <- trans_kegg( ps.trans)
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}

trans_kegg = function(ps.trans =  ps.trans) {

  kegg <- clusterProfiler::download_KEGG('ko')
  PATH2ID <- kegg $KEGGPATHID2EXTID
  PATH2NAME <- kegg$KEGGPATHID2NAME
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

  colnames(tem2)

  otu = tem2[,sample_names( ps.trans)]
  row.names(otu) = tem2$KEGGID
  head(otu)
  head(tem2)
  tax = data.frame(ID = tem2$KEGGID,
                   row.names = tem2$KEGGID,
                   DESCRPTION = tem2$DESCRPTION,DESCRPTION2 = tem2$DESCRPTION)
  head(tax)

  dat = EasyMultiOmics.db::db.pathway.level
  # dat = read.delim("./pathway-levels.txt")
  head(dat)
  colnames(dat)[1] = "ID"
  dat$ID = gsub("map","ko",dat$ID)
  tax2 = tax %>% left_join(dat,by = "ID") %>% as.data.frame()
  head(tax2)
  row.names(tax2) = tax2$ID
  ps1 = phyloseq(
   phyloseq:: otu_table(as.matrix(otu),taxa_are_rows = TRUE),
   phyloseq::tax_table(as.matrix(tax2)),
   phyloseq::sample_data( ps.trans)
  )
  return(ps1)
}
