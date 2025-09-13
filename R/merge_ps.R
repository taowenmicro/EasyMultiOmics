#' This function merges two `phyloseq` objects (representing microbiome and metabolomics data) into a single `phyloseq` object.
#' It allows for scaling, filtering, and merging of OTU tables, taxonomic tables, and sample data. Additionally, the user can
#' specify whether to merge the taxonomic tables or just the sample group information.
#'
#' @param ps1 A `phyloseq` object containing microbiome data (e.g., 16S rRNA sequencing data).
#' @param ps2 A `phyloseq` object containing metabolomics data (e.g., ITS sequencing or functional data).
#' @param N1 An integer specifying the number of top OTUs to select from `ps1` (default: 100).
#' @param N2 An integer specifying the number of top OTUs to select from `ps2` (default: 100).
#' @param scale A logical value indicating whether to scale the OTU counts in the microbiome (`ps1`) and metabolomics (`ps2`) data (default: `TRUE`).
#' @param onlygroup A logical value indicating whether to merge only the group information in the taxonomic table, without merging the full taxonomic information (default: `FALSE`).
#' @param dat1.lab A character string specifying the label for the first dataset (default: `"bac"` for bacteria).
#' @param dat2.lab A character string specifying the label for the second dataset (default: `"fun"` for fungi or functional data).
#'
#' @return A `phyloseq` object containing the merged OTU table, taxonomic table, and sample data from `ps1` and `ps2`.
#'   If `onlygroup = TRUE`, the taxonomic table will only contain the group information, while the OTU table will be fully merged.
#'
#' @examples
#' merged_data <- merge_ps(ps1 = microbiome_data,
#'                         ps2 = metabolomics_data,
#'                         N1 = 100,
#'                         N2 = 100,
#'                         scale = TRUE,
#'                         onlygroup = FALSE,
#'                         dat1.lab = "microbe",
#'                         dat2.lab = "metabolite")
#' @export
#' @author
#' Tao Wen \email{taowen@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
merge_ps <- function(ps1 ,
                     ps2,
                     N1 = 100,
                     N2 = 100,
                     scale = TRUE,
                     onlygroup = FALSE,#不进行列合并，只用于区分不同域
                     dat1.lab = "bac",
                     dat2.lab = "fun") {

  if (scale == TRUE) {
    if (!is.null(ps16s)) {
      ps1  = phyloseq::transform_sample_counts(ps1, function(x) x / sum(x) )
    }
    if (!is.null(psITS)) {
      ps2  = phyloseq::transform_sample_counts(ps2, function(x) x / sum(x) )
    }
  }
  if (!is.null(ps1)) {
    # ps_16s = phyloseq::filter_taxa(ps16s, function(x) mean(x) > N16s, TRUE)#select OTUs according to  relative abundance
    ps_16s  =  filter_OTU_ps(ps = ps1,Top = N1)
    ###
    otu_table_16s = as.data.frame(t(vegan_otu(ps_16s)))
    row.names(otu_table_16s) = paste(dat1.lab,row.names(otu_table_16s),sep = "_")
    ## change the OTU name of bac and fungi OTU table
    tax_table_16s = as.data.frame(vegan_tax(ps_16s))
    #-- add a col marked the bac and fungi
    if ("filed" %in%colnames(tax_table_16s)) {

    } else{
      row.names(tax_table_16s) = paste(dat1.lab,row.names(tax_table_16s),sep = "_")
      tax_table_16s$filed = rep(dat1.lab,length(row.names(tax_table_16s)))
    }

  }
  if (!is.null(ps2)) {
    # ps_ITS = phyloseq::filter_taxa(psITS, function(x) mean(x) > NITS , TRUE)#select OTUs according to  relative abundance
    ps_ITS = filter_OTU_ps(ps = ps2,Top = N2)
    otu_table_ITS = as.data.frame(t(vegan_otu(ps_ITS)))
    row.names(otu_table_ITS) = paste(dat2.lab,row.names(otu_table_ITS ),sep = "_")
    tax_table_ITS = as.data.frame(vegan_tax(ps_ITS))
    row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
    tax_table_ITS$filed = rep(dat2.lab,length(row.names(tax_table_ITS)))

    if ("filed" %in%colnames(tax_table_ITS)) {
    } else{
      row.names(tax_table_ITS) = paste(dat2.lab,row.names(tax_table_ITS),sep = "_")
      tax_table_ITS$filed = rep(dat2.lab,length(row.names(tax_table_ITS)))
    }
  }


  if (!is.null(ps2) & !is.null(ps1) ) {
    ## merge OTU table of bac and fungi



    otu_table = rbind(otu_table_16s[,intersect(names(otu_table_ITS),names(otu_table_16s))],

                      otu_table_ITS[,intersect(names(otu_table_ITS),names(otu_table_16s))])

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s,tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table

    mapping = as.data.frame( phyloseq::sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))


  } else if(is.null(psITS) & !is.null(ps16s) ) {
    otu_table = rbind(otu_table_16s)

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_16s$filed,row.names = row.names(otu_table),id = row.names(otu_table)))
    }
    #on of map table as final map table
    mapping = as.data.frame(sample_data(ps_16s))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = TRUE),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))


  } else if (!is.null(ps2) & is.null(ps1)){
    otu_table = rbind(otu_table_ITS)

    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_ITS)
      dim(otu_table)
    } else if(onlygroup == TRUE){
      tax_table = data.frame(filed = c(tax_table_ITS$filed),row.names = row.names(otu_table),id = row.names(otu_table))
    }
    #on of map table as final map table
    mapping = as.data.frame( phyloseq::sample_data(psITS))
    head(mapping)
    # mapping$Group4 = "all_sample"
    # mapping$Group4 = as.factor(mapping$Group4)
    ##merge all abject of phyloseq
    pallps <-  phyloseq::phyloseq( phyloseq::otu_table(as.matrix(otu_table),taxa_are_rows = T),
                                   phyloseq::sample_data(mapping),
                                   phyloseq::tax_table(as.matrix(tax_table)))

  }

  tax = pallps %>% vegan_tax() %>%
    as.data.frame() %>% dplyr::select(filed,everything())
  phyloseq::tax_table(pallps) = as.matrix(tax)


  return(pallps)
}
