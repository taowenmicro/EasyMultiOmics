MetalTast = function (env.dat, tabOTU, distance = "bray", method = "metal")
{
  for (j in 1:length(names(tabOTU))) {
    report = c()
    otu = tabOTU[[j]]
    name = names(tabOTU)[j]
    sle_env = vegan::decostand(env.dat, method = "standardize",
                               MARGIN = 2)
    env.std = vegan::decostand(sle_env, method = "standardize",
                               MARGIN = 2)
    BC.beta = vegan::vegdist(t(otu), method = "bray")
    JC.beta = vegan::vegdist(t(otu), method = "jaccard",
                             binary = T)
    if (method == "metal") {
      for (i in 1:ncol(env.std)) {
        envdis = vegan::vegdist(env.std[i], method = "euclidean",
                                na.rm = T)
        mantel.BC = vegan::mantel(envdis, BC.beta, na.rm = T)
        mantel.JC = vegan::mantel(envdis, JC.beta, na.rm = T)
        report = rbind(report, c(colnames(env.dat)[i],
                                 mantel.BC$statistic, mantel.BC$signif, mantel.JC$statistic,
                                 mantel.JC$signif))
      }
    }
    if (method == "Part.Mantel") {
      for (i in 1:ncol(env.std)) {
        envdis = dist(env.std[, i])
        envdis2 = dist(env.std[, -i])
        pmantel.BC = mantel.partial(BC.beta, envdis,
                                    envdis2, na.rm = T)
        pmantel.JC = mantel.partial(JC.beta, envdis,
                                    envdis2, na.rm = T)
        report = rbind(report, c(colnames(env.dat)[i],
                                 pmantel.BC$statistic, pmantel.BC$signif, pmantel.JC$statistic,
                                 pmantel.JC$signif))
      }
    }
    colnames(report) <- c("Envs", paste(rep(c("r", "p"),
                                            2), rep(c("BC", "JC"), each = 2), sep = "."))
    report = as.matrix(report)
    report = as.data.frame(report)
    report[, 2:dim(report)[2]] <- lapply(report[, 2:dim(report)[2]],
                                         as.character)
    report[, 2:dim(report)[2]] <- lapply(report[, 2:dim(report)[2]],
                                         as.numeric)
    if (distance == "bray") {
      report0 = report[1:3]
      colnames(report0)[-1] = paste(name, colnames(report0)[-1],
                                    sep = "")
    }
    else if (distance == "jcd") {
      report0 = report[c(1, c(4:5))]
      colnames(report0)[-1] = paste(name, colnames(report0)[-1],
                                    sep = "")
    }
    if (j == 1) {
      report1 = report0
    }
    if (j != 1) {
      report1 = report1 %>% inner_join(report0)
    }
  }
  return(report1)
}
