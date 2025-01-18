#' @title Merge 16S and ITS Microbiome Datasets
#'
#' @description
#' This function combines two phyloseq objects representing 16S rRNA and ITS datasets.
#' It normalizes the datasets (optional) and merges taxonomic and OTU tables. Users
#' can choose to retain only taxonomic groupings if needed.
#'
#' @param ps16s A phyloseq object containing 16S rRNA dataset. Default is `ps16`.
#' @param psITS A phyloseq object containing ITS dataset. Default is `psIT`.
#' @param N16s An integer specifying the number of top OTUs or taxa to retain from the 16S dataset. Default is 100.
#' @param NITS An integer specifying the number of top OTUs or taxa to retain from the ITS dataset. Default is 100.
#' @param scale A logical value indicating whether to normalize datasets to relative abundance. Default is `TRUE`.
#' @param onlygroup A logical value specifying whether to retain only taxonomic group information. Default is `FALSE`.
#' @param dat1.lab A character string used as a prefix for 16S dataset taxonomic IDs. Default is `"bac"`.
#' @param dat2.lab A character string used as a prefix for ITS dataset taxonomic IDs. Default is `"fun"`.
#'
#' @return
#' A merged phyloseq object containing OTU table, taxonomy table, and sample metadata.
#'
#' @examples
#' \dontrun{
#' merged_ps <- merge16S_ITS(ps16s = ps16, psITS = psIT, N16s = 50, NITS = 50, scale = TRUE)
#' }
#'
#'
#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#'
merge16S_ITS=function (ps16s = ps16, psITS = psIT, N16s = 100, NITS = 100,
          scale = TRUE, onlygroup = FALSE, dat1.lab = "bac", dat2.lab = "fun")
{
  if (scale == TRUE) {
    if (!is.null(ps16s)) {
      ps16s = phyloseq::transform_sample_counts(ps16s,
                                                function(x) x/sum(x))
    }
    if (!is.null(psITS)) {
      psITS = phyloseq::transform_sample_counts(psITS,
                                                function(x) x/sum(x))
    }
  }
  if (!is.null(ps16s)) {
    ps_16s = filter_OTU_ps(ps = ps16s, Top = N16s)
    otu_table_16s = as.data.frame(t(vegan_otu(ps_16s)))
    row.names(otu_table_16s) = paste(dat1.lab, row.names(otu_table_16s),
                                     sep = "_")
    tax_table_16s = as.data.frame(vegan_tax(ps_16s))
    row.names(tax_table_16s) = paste(dat1.lab, row.names(tax_table_16s),
                                     sep = "_")
    tax_table_16s$filed = rep(dat1.lab, length(row.names(tax_table_16s)))
  }
  if (!is.null(psITS)) {
    ps_ITS = filter_OTU_ps(ps = psITS, Top = NITS)
    otu_table_ITS = as.data.frame(t(vegan_otu(ps_ITS)))
    row.names(otu_table_ITS) = paste(dat2.lab, row.names(otu_table_ITS),
                                     sep = "_")
    tax_table_ITS = as.data.frame(vegan_tax(ps_ITS))
    row.names(tax_table_ITS) = paste(dat2.lab, row.names(tax_table_ITS),
                                     sep = "_")
    tax_table_ITS$filed = rep(dat2.lab, length(row.names(tax_table_ITS)))
  }
  if (!is.null(psITS) & !is.null(ps16s)) {
    otu_table = rbind(otu_table_16s, otu_table_ITS)
    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s, tax_table_ITS)
      dim(otu_table)
    }
    else if (onlygroup == TRUE) {
      tax_table = data.frame(filed = c(tax_table_16s$filed,
                                       tax_table_ITS$filed), row.names = row.names(otu_table),
                             id = row.names(otu_table))
    }
    mapping = as.data.frame(phyloseq::sample_data(ps_16s))
    head(mapping)
    pallps <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(otu_table),
                                                     taxa_are_rows = T), phyloseq::sample_data(mapping),
                                 phyloseq::tax_table(as.matrix(tax_table)))
  }
  else if (is.null(psITS) & !is.null(ps16s)) {
    otu_table = rbind(otu_table_16s)
    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_16s)
      dim(otu_table)
    }
    else if (onlygroup == TRUE) {
      tax_table = data.frame(filed = c(tax_table_16s$filed,
                                       row.names = row.names(otu_table), id = row.names(otu_table)))
    }
    mapping = as.data.frame(sample_data(ps_16s))
    head(mapping)
    pallps <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(otu_table),
                                                     taxa_are_rows = T), phyloseq::sample_data(mapping),
                                 phyloseq::tax_table(as.matrix(tax_table)))
  }
  else if (!is.null(psITS) & is.null(ps16s)) {
    otu_table = rbind(otu_table_ITS)
    if (onlygroup == FALSE) {
      tax_table = rbind(tax_table_ITS)
      dim(otu_table)
    }
    else if (onlygroup == TRUE) {
      tax_table = data.frame(filed = c(tax_table_ITS$filed),
                             row.names = row.names(otu_table), id = row.names(otu_table))
    }
    mapping = as.data.frame(phyloseq::sample_data(psITS))
    head(mapping)
    pallps <- phyloseq::phyloseq(phyloseq::otu_table(as.matrix(otu_table),
                                                     taxa_are_rows = T), phyloseq::sample_data(mapping),
                                 phyloseq::tax_table(as.matrix(tax_table)))
  }
  tax = pallps %>% vegan_tax() %>% as.data.frame() %>% dplyr::select(filed,
                                                                     everything())
  phyloseq::tax_table(pallps) = as.matrix(tax)
  return(pallps)
}
