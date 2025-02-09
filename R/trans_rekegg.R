#' @title KEGG Reaction Transformation for Phyloseq Data
#' @description
#' This function transforms a Phyloseq object by mapping KEGG Orthology (KO) terms
#' to KEGG reactions. It uses data from the EasyMultiOmics database and merges KO
#' terms with reaction data. The function processes the OTU (Operational Taxonomic Unit)
#' table from a Phyloseq object, aggregates the data by KEGG reaction, and returns
#' a transformed Phyloseq object with KEGG reaction descriptions.
#'
#' @param ps A Phyloseq object containing OTU table and associated taxonomic data.
#' The default is `ps.trans`, which assumes that the Phyloseq object is provided externally.
#'
#' @return A transformed Phyloseq object with the OTU table aggregated by KEGG reaction
#' IDs and the taxonomic table updated with KEGG reaction descriptions.
#'
#' @export
#' @examples
#' \dontrun{
#' ps_transformed <- trans_rekegg(ps.trans)
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
trans_rekegg = function(ps =  ps.trans) {

  dat =EasyMultiOmics.db:: db.ko.reaction

  # colnames(dat) = c("ID","reaction")
  head(dat)

  ko = ps %>% vegan_otu() %>% t() %>%
    as.data.frame() %>% rownames_to_column("KO")
  head(ko)
  dat$ID = gsub("ko:","",dat$ID)
  dat$reaction = gsub("rn:","",dat$reaction)


  tem2 = ko %>% left_join(dat,by = c("KO" = "ID")) %>%
    tidyfst::filter_dt(reaction != "") %>%
    tidyfst::summarise_vars(is.numeric,sum,by =c("reaction")) %>%
    as.data.frame()

  colnames(tem2)

  otu = tem2[,sample_names(ps)]
  row.names(otu) = tem2$reaction
  head(otu)


  dat2 =  EasyMultiOmics.db:: db.reaction
  head(dat2)

  tax = data.frame(row.names = dat2$reaction,DESCRPTION = dat2$DESCRPTION,DESCRPTION2 = dat2$DESCRPTION)
  head(tax)

  ps3 = phyloseq(
    phyloseq:: otu_table(as.matrix(otu),taxa_are_rows = TRUE),
    phyloseq:: tax_table(as.matrix(tax)),
    phyloseq:: sample_data(ps)
  )

  ps3

  return(ps3)
}
