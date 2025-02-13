#' @title CNPS Gene Data for Phyloseq Object
#' @description
#' This function processes a `phyloseq` object by mapping the OTU data to CNPS gene information,
#' and returns a new `phyloseq` object with the CNPS gene data included.
#'
#' @param psko A `phyloseq` object (default is `ps.kegg`) containing OTU data to be mapped to CNPS gene information.
#'
#' @return A `phyloseq` object containing the OTU data, taxonomic data, and sample data, with CNPS gene information included.
#' @details
#' This function performs the following steps:
#' \itemize{
#'   \item Loads the CNPS gene data from `EasyMultiOmics.db::db.cnps`.
#'   \item Extracts OTU data from the `phyloseq` object and transforms it into a data frame.
#'   \item Joins the OTU data with the CNPS gene data using the KO identifier.
#'   \item Creates a unique ID for each sample by combining the CNPS gene group and KO number.
#'   \item Constructs a new `phyloseq` object, including the OTU table, taxonomy, and sample data, with CNPS gene information.
#' }
#' @examples
#' # Assuming you have a phyloseq object `ps.kegg`
#' ps_cnps<- cnps_ps(ps.kegg)
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}

cnps_ps = function(psko = ps.kegg ){
  dat = EasyMultiOmics.db::db.cnps
  head(dat)
  otu = psko %>%
    #scale_micro() %>%
    vegan_otu() %>% t() %>%
    as.data.frame()
  otu$ID = row.names(otu)
  head(otu)
  tab= otu %>% inner_join(dat,by = c("ID" = "K.number"))
  head(tab)

  # write.csv(C.tab,paste(otupath,"/C_cycle_imformation.csv",sep = ""),quote = FALSE)

  map=  psko %>% phyloseq::sample_data()
  head(map)
  tab = tab %>% dplyr::distinct(Group, ID, .keep_all = TRUE)
  otu = tab[,row.names(map)]
  head(otu)
  paste0(tab$Group,tab$ID) %>% unique()
  row.names(otu) = paste0(tab$Group,tab$ID)



  tax= tab[,setdiff(colnames(tab),row.names(map))]
  row.names(tax) = paste0(tab$Group,tab$ID)

  ps = phyloseq(
    phyloseq::otu_table(as.matrix(otu),taxa_are_rows = TRUE),
    phyloseq::tax_table(as.matrix(tax)),
    phyloseq::sample_data(psko)
  )

  return( ps)

}
