#' @title Perform Mantel and Partial Mantel Tests for Environmental and OTU Data
#'
#' @description
#' The `MetalTast` function calculates Mantel and Partial Mantel tests to assess the relationship between
#' environmental variables and microbial community composition. It supports multiple distance methods
#' (e.g., Bray-Curtis and Jaccard) and allows for partial Mantel tests to control for confounding variables.
#'
#' @param env.dat A data frame of environmental variables, where rows represent samples and columns represent variables.
#' @param tabOTU A list of OTU tables, where each element is a data frame/matrix of OTU abundances (rows = taxa, columns = samples).
#' @param distance A character string specifying the distance method to use for the microbial beta diversity matrix. Options are `"bray"` (Bray-Curtis) or `"jcd"` (Jaccard). Default is `"bray"`.
#' @param method A character string specifying the analysis type. Options are:
#'   \itemize{
#'     \item `"metal"`: Perform Mantel tests.
#'     \item `"Part.Mantel"`: Perform Partial Mantel tests to control for the influence of other variables.
#'   }
#' @return
#' A data frame summarizing the Mantel or Partial Mantel test results, with columns:
#' \itemize{
#'   \item `Envs`: Names of the environmental variables.
#'   \item `R.BC` and `P.BC`: Correlation coefficient (R) and significance (P) for Bray-Curtis distance.
#'   \item `R.JC` and `P.JC`: Correlation coefficient (R) and significance (P) for Jaccard distance (if applicable).
#' }
#'
#' @examples
#' \dontrun{
#' library(vegan)
#' library(dplyr)
#' result_partial <- MetalTast(env.dat = env.dat,
#'                              tabOTU = tabOTU,
#'                              distance = "bray",
#'                              method = "Part.Mantel")
#' print(result_partial)
#' }
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @export
MetalTast = function (env.dat, tabOTU, distance = "bray", method = "metal")
{
  # j=1
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
                             binary = TRUE)
    if (method == "metal") {
      i=1
      for (i in 1:ncol(env.std)) {
        envdis = vegan::vegdist(env.std[i], method = "euclidean",
                                na.rm = TRUE)
        mantel.BC = vegan::mantel(envdis, BC.beta, na.rm = TRUE)
        mantel.JC = vegan::mantel(envdis, JC.beta, na.rm = TRUE)
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
                                    envdis2, na.rm = TRUE)
        pmantel.JC = mantel.partial(JC.beta, envdis,
                                    envdis2, na.rm =TRUE)
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
      colnames(report0)[-1] = c("R","P")
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
      report1 = report1 %>% dplyr::inner_join(report0)

    }
  }
  return(report1)
}
