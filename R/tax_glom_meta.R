#' @title Summarize Microbiome Data by Taxonomic Rank
#' @description
#'
#' This function aggregates the OTU data in a `phyloseq` object by a specified taxonomic rank
#' (default is "Phylum") and returns a new `phyloseq` object with the summarized data.
#'
#' @param ps A `phyloseq` object containing OTU and taxonomic data.
#' @param ranks A character string or numeric value specifying the taxonomic rank
#'        by which to group the OTU data. Default is "Phylum". If a numeric value is provided,
#'        it will select the rank from the list of ranks in the `phyloseq` object.
#'
#' @return A `phyloseq` object containing the aggregated OTU data based on the specified taxonomic rank.
#'
#' @details
#' The function handles missing taxonomic values (e.g., NA, empty strings, or "NA") and replaces them with "Unknown".
#' It also removes duplicate taxonomic names within the specified rank. The summarized OTU data is returned as a new
#' `phyloseq` object.
#'
#' @examples
#' # Example of usage with a phyloseq object `ps`
#' ps_summary <- tax_glom_meta(ps =ps.card, ranks ="Drug_Class")
#' head(ps_summary)
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
tax_glom_meta <- function(ps = ps,ranks = "Phylum") {


  if (  is.numeric(ranks)) {
    ranks <- phyloseq::rank.names(ps)[ranks]
  }


  otu <- as.data.frame(t(vegan_otu(ps)))
  tax <- as.data.frame(vegan_tax(ps))

  # building group
  tax[[ranks]][is.na(tax[[ranks]])] = "Unknown"
  tax[[ranks]][tax[[ranks]] == ""] = "Unknown"
  tax[[ranks]][tax[[ranks]] == "NA"] = "Unknown"
  split <- split(otu,tax[[ranks]])
  #calculate sum by group
  apply <- lapply(split,function(x)colSums(x[]))
  # result chack
  otucon <- do.call(rbind,apply)

  taxcon <- tax[1:match(ranks,colnames(tax))]
  tax[[ranks]] = gsub("'","",tax[[ranks]])
  taxcon <- taxcon[!duplicated(tax[[ranks]]),]
  #-tax name with NA wound be repeated with unknown
  taxcon[[ranks]][is.na(taxcon[[ranks]])] = "unknown"
  row.names(taxcon) <- taxcon[[ranks]]

  row.names(otucon) %>% length()
  row.names(otucon) %>% unique() %>% length()

  # row.names(otucon) = gsub("-","_", row.names(otucon))
  row.names(otucon) = gsub("'","", row.names(otucon))
  row.names(taxcon) = gsub("'","", row.names(taxcon))
  # row.names(otucon[250:270,]) %>% unique()

  head(otucon)
  otucon = otucon %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>% dplyr::distinct(ID,.keep_all = TRUE) %>%
    column_to_rownames("ID")
  taxcon = taxcon %>%
    as.data.frame() %>%
    rownames_to_column("ID") %>% dplyr::distinct(ID,.keep_all = TRUE) %>%
    column_to_rownames("ID")


  pscon <- phyloseq::phyloseq(
    phyloseq::otu_table( as.matrix(otucon),taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(taxcon)),
    phyloseq::sample_data(ps)
  )

  return(pscon)
}
