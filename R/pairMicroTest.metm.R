
#' @title Pairwise Microbial Community Statistical Testing
#' @description
#' The `pairMicroTest.metm` function performs pairwise statistical testing on microbial community data using selected metrics (Adonis, MRPP, or ANOSIM). It calculates statistical differences between groups in a `phyloseq` object based on a specified distance metric.
#' @param ps A `phyloseq` object containing microbial community data, including the OTU table and sample metadata.
#' @param Micromet A string specifying the statistical method to use. Options are `"adonis"` (default), `"MRPP"`, or `"anosim"`.
#' @param dist A string specifying the distance metric to use. Default is `"bray"` (Bray-Curtis distance).
#'
#' @return stat table
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' dat2 = pairMicroTest.metm(ps = ps.16s, Micromet = "MRPP", dist = "bray")
#' dat2
#' @export

pairMicroTest.metm=function (ps = ps, Micromet = "adonis", dist = "bray")
{
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  ps1_rela
  map = as.data.frame(phyloseq::sample_data(ps1_rela))
  otu = ps1_rela %>% ggClusterNet::vegan_otu() %>% as.data.frame()
  map$Group = as.factor(map$Group)
  aa = levels(map$Group)
  aaa = combn(aa, 2)
  aaa
  dim(aaa)[2]
  ID = rep("a", dim(aaa)[2])
  R = rep("a", dim(aaa)[2])
  P = rep("a", dim(aaa)[2])
  for (i in 1:dim(aaa)[2]) {
    print(i)
    Desep_group = aaa[, i]
    map = as.data.frame(phyloseq::sample_data(ps1_rela))
    map$ID = row.names(map)
    maps <- map[map$Group %in% Desep_group, ]
    row.names(maps) = maps$ID
    ps_sub = ps1_rela
    phyloseq::sample_data(ps_sub) = maps
    ps_sub = phyloseq::filter_taxa(ps_sub, function(x) sum(x) >
                                     0, TRUE)
    ps_sub
    map = as.data.frame(phyloseq::sample_data(ps_sub))
    gg = map$Group
    unif <- phyloseq::distance(ps_sub, method = dist)
    if (Micromet == "MRPP") {
      mrpp = vegan::mrpp(unif, map$Group)
      as1 = round(mrpp$delta, 3)
      R2 <- paste("MRPP.delta ", as1, sep = "")
      R2
      p_v = paste("p: ", round(mrpp$Pvalue, 3), sep = "")
      p_v
    }
    if (Micromet == "anosim") {
      dat.ano = vegan::anosim(unif, map$Group)
      a = round(dat.ano$statistic, 3)
      R2 <- paste("ANOSIM.r ", a, sep = "")
      R[i] = R2
      p_v = paste("p: ", round(dat.ano$signif, 3), sep = "")
      P[i] = p_v
    }
    if (Micromet == "adonis") {
      str(otu)
      str(dune.env)
      ado = vegan::adonis2(unif ~ map$Group, method = "bray",
                           by = NULL)
      R2 = paste("Adonis:R ", round(ado$R2[1], 3), sep = "")
      p_v = paste("p: ", ado$`Pr(>F)`[1], sep = "")
      R[i] = R2
      P[i] = p_v
    }
    ID[i] = paste(Desep_group[1], Desep_group[2], sep = "_VS_")
    P[i] = p_v
    R[i] = R2
  }
  result = data.frame(ID = ID, stat = R, p = P)
  return(result)
}
