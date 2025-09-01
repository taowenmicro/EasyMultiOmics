#' KEGG Reaction-Based Analysis for Phyloseq Object
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
#' This function performs the following steps:
#' \itemize{
#'   \item Loads KO-reaction data from `EasyMultiOmics.db::db.ko.reaction` by default or the provided `dat`.
#'   \item Transforms the OTU table in `ps` to KO identifiers using `vegan_otu()`.
#'   \item Joins KO identifiers with KEGG reaction data, mapping the KO to reactions.
#'   \item Summarizes the data by reaction and creates a new `phyloseq` object with summarized OTU data for each reaction.
#'   \item Adds reaction descriptions from `EasyMultiOmics.db::db.reaction` or the provided `dat2`.
#' }
#' @examples
#' # Assuming you have a phyloseq object `ps`
#' reaction_ps <- reaction_function(ps)

#' @export
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
reaction_function <- function(ps = ps,
                              dat =NULL,
                              dat2 = NULL){

dat = EasyMultiOmics.db::db.ko.reaction

colnames(dat) = c("ID","reaction")
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

dat2=  EasyMultiOmics.db::db.reaction
head(dat2)
colnames(dat2) = c("reaction","DESCRPTION")

tax = data.frame(row.names = dat2$reaction,DESCRPTION = dat2$DESCRPTION,DESCRPTION2 = dat2$DESCRPTION)
head(tax)

ps1 = phyloseq(
  phyloseq::otu_table(as.matrix(otu),taxa_are_rows = TRUE),
 phyloseq:: tax_table(as.matrix(tax)),
 phyloseq::sample_data(ps)
)

}
