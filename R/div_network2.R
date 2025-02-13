div_network2 =function (ps, group = "Group", flour = TRUE, N = 0.5)
{
  mapping = as.data.frame(phyloseq::sample_data(ps))
  mapping = mapping[, group]
  colnames(mapping[, group]) <- "Group"
  sample_data(ps) = mapping
  ps_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  ps_rela
  aa = vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count
  dim(count)
  sub_design <- as.data.frame(phyloseq::sample_data(ps))
  count[count > 0] <- 1
  count2 = as.data.frame(count)
  iris.split <- split(count2, as.factor(sub_design$Group))
  iris.apply <- lapply(iris.split, function(x) colSums(x[]))
  iris.combine <- do.call(rbind, iris.apply)
  ven2 = t(iris.combine)
  for (i in 1:length(unique(phyloseq::sample_data(ps)$Group))) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,
                                                                1]
    bb = as.data.frame(table(phyloseq::sample_data(ps)$Group))[i,
                                                               2]
    ven2[, aa] = ven2[, aa]/bb
  }
  ven2[ven2 < N] = 0
  ven2[ven2 >= N] = 1
  ven2 = as.data.frame(ven2)
  ven3 = as.list(ven2)
  ven2 = as.data.frame(ven2)
  if (flour == TRUE) {
    ven2 = ven2[rowSums(ven2) == dim(ven2)[2] | rowSums(ven2) ==
                  1, ]
  }
  tax_table = as.data.frame(vegan_tax(ps))
  otu_table = as.data.frame(t(vegan_otu(ps)))
  dim(otu_table)
  otu_net = merge(ven2, tax_table, by = "row.names", all = F)
  row.names(otu_net) = otu_net$Row.names
  otu_net$Row.names = NULL
  head(otu_net)
  OTU = as.matrix(otu_table)
  norm = t(t(OTU)/colSums(OTU))
  norm1 = norm %>% t() %>% as.data.frame()
  iris.split <- split(norm1, as.factor(mapping$Group))
  iris.apply <- lapply(iris.split, function(x) colSums(x))
  norm2 <- do.call(rbind, iris.apply) %>% t()
  head(norm2)
  colnames(norm2) = paste(colnames(norm2), "mean", sep = "")
  dim(otu_net)
  otu_net2 = merge(otu_net, norm2, by = "row.names", all = F)
  dim(otu_net2)
  head(otu_net2)
  colnames(otu_net2)[1] = "ID"
  sample_label = otu_net2[1:length(unique(mapping$Group)),
  ]
  sample_label$ID = unique(mapping$Group)
  point_label = rbind(otu_net2, sample_label)
  head(otu_net2)
  if (length(colnames(phyloseq::tax_table(ps))) == 6) {
    net_all = reshape2::melt(otu_net2, id = c("ID", "Kingdom",
                                              "Phylum", "Class", "Order", "Family", "Genus", paste(unique(mapping$Group),
                                                                                                   "mean", sep = "")))
  }
  if (length(colnames(phyloseq::tax_table(ps))) == 7) {
    net_all = reshape2::melt(otu_net2, id = c("ID", "Kingdom",
                                              "Phylum", "Class", "Order", "Family", "Genus", "Species",
                                              paste(unique(mapping$Group), "mean", sep = "")))
  }
  if (length(colnames(phyloseq::tax_table(ps))) != 7 | length(colnames(phyloseq::tax_table(ps))) !=
      6) {
    net_all = reshape2::melt(otu_net2, id = c( rank_names(ps),
                                              paste(unique(mapping$Group), "mean", sep = "")))
  }
  net_filter <- dplyr::filter(net_all, value != 0)
  colnames(net_all)
  net_fal = data.frame(source = net_filter$ID, target = net_filter$variable,
                       connect = rep("pp", nrow(net_filter)), value = net_filter$value,
                       Label = net_filter$ID)
  head(net_fal)
  return(list(net_fal, point_label, ven2))
}
