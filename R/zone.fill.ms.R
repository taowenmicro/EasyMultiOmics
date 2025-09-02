#' @title Frequency filtering of metabolites in repeated samples
#' @description
#' This function mainly filters the frequency of metabolite occurrence in duplicate samples to improve data reliability.
#' By default, metabolites that appear in 2/3 of duplicate samples are retained,
#' otherwise their abundance in each repeated sample is filled to 0.
#' @param ps A phyloseq format file used as an alternative for the input containing metabolite composition table,
#' metabolite classification table, and sample metadata.
#' @param method The method to use for filling zones. Default is "repeat".
#' @param n Metabolite appearance frequency threshold. Default is 0.667(2/3).
#' @return The filled phyloseq format file after filtered based on the selected
#' metabolite appearance frequency threshold.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' ps.ms2 = zone.fill.ms(ps = ps.ms,method = "repeat")
#' @export
zone.fill.ms = function(ps ,method = "repeat",n = 0.667){
  if (method == "repeat") {
    map = ps %>% sample_data()
    ids = map$Group %>% unique()
    map$ID = row.names(map)
    for (i in 1:length(ids)) {
      id = map %>%
        as.tibble() %>%
        dplyr::filter(Group == ids[i]) %>%.$ID
      otu = ps %>% vegan_otu() %>% t() %>%
        as.data.frame() %>% #rownames_to_column("ID") %>%
        select(id) %>% as.matrix()
      head(otu)
      otu_01 <- otu
      otu_01 <- ifelse(otu_01>0,1,0)
      otu[rowSums(otu_01) < length(id)*n,] = 0

      if (FALSE) {
        otu1 = otu[rowSums(otu) > length(id)*n,]
      }
      if (i == 1 ) {
        otun = otu

      } else{
        otun = cbind(otun,otu)
      }

    }

  }
  otu_table(ps) = otu_table(as.matrix(otun),taxa_are_rows = TRUE)
  return(ps)
}
