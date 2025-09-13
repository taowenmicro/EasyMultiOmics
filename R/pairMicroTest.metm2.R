#' @title Pairwise Microbial Community Statistical Testing
#' @description
#' Perform pairwise statistical testing on microbial community data
#' (Adonis, MRPP, or ANOSIM) using a phyloseq object and a given distance metric.
#'
#' @param ps A `phyloseq` object containing OTU/ASV data and sample metadata.
#' @param Micromet Statistical method: `"adonis"` (default), `"MRPP"`, `"anosim"`.
#' @param dist Distance metric for community dissimilarity (default = `"bray"`).
#'
#' @return A data.frame with pairwise group comparisons, statistic, and p-value.
#' @author
#' Tao Wen \email{2018203048@@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' \dontrun{
#' pairMicroTest.metm2(ps, Micromet = "permdisp", dist = "bray")
#' pairMicroTest.metm2(ps, Micromet = "mantel", dist = "bray")
#' }
#' @export
pairMicroTest.metm2 = function (ps, Micromet = "adonis", dist = "bray")
{
  ps1_rela = phyloseq::transform_sample_counts(ps, function(x) x/sum(x))
  map = as.data.frame(phyloseq::sample_data(ps1_rela))
  otu = ps1_rela %>% ggClusterNet::vegan_otu() %>% as.data.frame()
  map$Group = as.factor(map$Group)
  aa = levels(map$Group)
  aaa = combn(aa, 2)

  ID = rep("a", ncol(aaa))
  R = rep("a", ncol(aaa))
  P = rep("a", ncol(aaa))

  for (i in 1:ncol(aaa)) {
    Desep_group = aaa[, i]
    map_sub = map[map$Group %in% Desep_group, ]
    ps_sub = phyloseq::prune_samples(rownames(map_sub), ps1_rela)
    ps_sub = phyloseq::filter_taxa(ps_sub, function(x) sum(x) > 0, TRUE)
    unif <- phyloseq::distance(ps_sub, method = dist)
    map_sub = as.data.frame(phyloseq::sample_data(ps_sub))

    if (Micromet == "MRPP") {
      mrpp = vegan::mrpp(unif, map_sub$Group)
      R2 <- paste("MRPP.delta ", round(mrpp$delta, 3))
      p_v <- paste("p: ", round(mrpp$Pvalue, 3))

    } else if (Micromet == "anosim") {
      dat.ano = vegan::anosim(unif, map_sub$Group)
      R2 <- paste("ANOSIM.r ", round(dat.ano$statistic, 3))
      p_v <- paste("p: ", round(dat.ano$signif, 3))

    } else if (Micromet == "adonis") {
      ado = vegan::adonis2(unif ~ map_sub$Group)
      R2 <- paste("Adonis:R ", round(ado$R2[1], 3))
      p_v <- paste("p: ", ado$`Pr(>F)`[1])

    } else if (Micromet == "permdisp") {
      mod = vegan::betadisper(unif, map_sub$Group)
      per = vegan::permutest(mod)
      R2 <- paste("PERMDISP.F ", round(per$tab[1,4], 3))
      p_v <- paste("p: ", per$tab[1,6])

    } else if (Micromet == "mantel") {
      # Mantel test 两组群体的子集距离矩阵相关性
      mat <- as.matrix(unif)
      group_dist = dist(model.matrix(~ map_sub$Group - 1))
      mt = vegan::mantel(unif, group_dist, permutations = 999)
      R2 <- paste("Mantel.r ", round(mt$statistic, 3))
      p_v <- paste("p: ", mt$signif)

    } else {
      stop("Unsupported Micromet method!")
    }

    ID[i] = paste(Desep_group[1], Desep_group[2], sep = "_VS_")
    R[i] = R2
    P[i] = p_v
  }

  result = data.frame(ID = ID, stat = R, p = P)
  return(result)
}
