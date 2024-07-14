
# ps.ms2 = zone.fill.ms (ps = ps.ms,method = "repeat")

zone.fill.ms = function(ps ,method = "repeat",n = 0.667){
  if (method == "repeat") {
    map = ps %>% sample_data()
    ids = map$Group %>% unique()
    for (i in 1:length(ids)) {

      map$ID = row.names(map)
      id = map %>%
        as.tibble() %>%
        filter(Group == ids[i]) %>%.$ID
      otu = ps %>% vegan_otu() %>% t() %>%
        as.data.frame() %>% #rownames_to_column("ID") %>%
        select(id) %>% as.matrix()
      head(otu)
      otu[rowSums(otu) < length(id)*n,] = 0

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
