#' KEGG Pathway Analysis for Phyloseq Object
#'
#' This function downloads KEGG pathway information, transforms an OTU table in a
#' `phyloseq` object to KO identifiers, merges the KO identifiers with KEGG pathways,
#' and returns a `phyloseq` object with KEGG pathway data.
#'
#' @param ps A `phyloseq` object containing OTU data to be transformed into KEGG pathways.
#'
#' @return A `phyloseq` object with transformed OTU table and KEGG pathway information.
#' @details
#' This function performs several steps:
#' \itemize{
#'   \item Downloads KEGG pathway data using the `clusterProfiler::download_KEGG` function.
#'   \item Transforms the OTU table in `ps` to KO identifiers using `vegan_otu()`.
#'   \item Merges KO identifiers with KEGG pathway names and IDs.
#'   \item Summarizes pathway data and creates a new `phyloseq` object containing summarized OTU data.
#' }
#' @examples
#' # Assuming you have a phyloseq object `ps`
#' kegg_ps <- kegg_function(ps.kegg)
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export

kegg_function <- function(ps = ps) {
  kegg <- clusterProfiler::download_KEGG('ko')
  PATH2ID <- kegg$KEGGPATHID2EXTID
  PATH2NAME <- kegg$KEGGPATHID2NAME

  # Merge PATH2ID and PATH2NAME
  PATH_ID_NAME <- merge(PATH2ID, PATH2NAME, by = 'from')
  colnames(PATH_ID_NAME) <- c('KEGGID', 'KO', 'DESCRIPTION')

  # Transform ps to ko
  ko <- ps %>%
    vegan_otu() %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("KO")

  # Join and summarize data
  tem2 <- ko %>%
    left_join(PATH_ID_NAME, by = c("KO" = "KO")) %>%
    tidyfst::filter_dt(KEGGID != "") %>%
    tidyfst::summarise_vars(is.numeric, sum, by = c("KEGGID", "DESCRIPTION")) %>%
    as.data.frame()

  # Prepare otu table
  otu <- tem2[, sample_names(ps)]
  row.names(otu) <- tem2$KEGGID

  # Prepare tax table
  tax <- data.frame(
    row.names = tem2$KEGGID,
    DESCRIPTION = tem2$DESCRIPTION,
    DESCRIPTION2 = tem2$DESCRIPTION
  )

  # Create phyloseq object
  ps1 <- phyloseq(
    phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(tax)),
    phyloseq::sample_data(ps)
  )

  return(ps1 )
}
