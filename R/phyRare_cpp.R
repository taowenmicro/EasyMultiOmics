#' @title Automatic Rarefaction for Phyloseq Objects
#' @description
#' Automatically performs rarefaction on a phyloseq object.
#' It tries to use the faster \code{phyRare_cpp()} if available;
#' otherwise, it falls back to the R version \code{phyRare()}.
#'
#' @param ps A \code{phyloseq} object containing OTU/ASV table and sample data.
#' @param N Integer. Rarefaction depth (number of sequences per sample). Default is 3000.
#' @param ... Additional parameters passed to \code{phyRare_cpp()} or \code{phyRare()}.
#' @return A rarefied \code{phyloseq} object.
#' @export
#' @examples
#' ps_rarefied <- phyRare_auto(ps, N = 3000)
phyRare_auto <- function(ps, N = 3000, ...) {
  if (exists("phyRare_cpp", mode = "function")) {
    message("Using phyRare_cpp() for rarefaction.")
    phyRare_cpp(ps, N = N, ...)
  } else if (exists("phyRare", mode = "function")) {
    message("phyRare_cpp() not found. Using phyRare() instead.")
    phyRare(ps, N = N, ...)
  } else {
    stop("Neither phyRare_cpp nor phyRare exists. Please define at least one.")
  }
}


#' @title R-based Rarefaction for Phyloseq Objects
#' @description
#' Performs rarefaction on a phyloseq object using R (vegan::rrarefy).
#' Supports parallel computation.
#'
#' @param ps A \code{phyloseq} object.
#' @param N Integer. Rarefaction depth. Default 3000.
#' @param cores Integer. Number of cores for parallel computation. Default 10.
#' @return A rarefied \code{phyloseq} object.
#' @export
#' @examples
#' ps_rarefied <- phyRare(ps, N = 3000, cores = 4)
phyRare <- function(ps, N = 3000, cores = 10) {

  vegan_otu <- function(physeq){
    OTU <- otu_table(physeq)
    if(taxa_are_rows(OTU)){
      OTU <- t(OTU)
    }
    as(OTU, "matrix")
  }

  otb <- vegan_otu(ps)

  if (cores > 1 && .Platform$OS.type == "windows") {
    cl <- makeCluster(cores)
    clusterEvalQ(cl, library(vegan))
    otb1 <- do.call(rbind, parallel::parLapply(cl, 1:nrow(otb), function(i){
      vegan::rrarefy(otb[i,,drop=FALSE], N)
    }))
    stopCluster(cl)
  } else if (cores > 1) {
    otb1 <- do.call(rbind, parallel::mclapply(1:nrow(otb), function(i){
      vegan::rrarefy(otb[i,,drop=FALSE], N)
    }, mc.cores = cores))
  } else {
    otb1 <- vegan::rrarefy(otb, N)
  }

  phyloseq(otu_table(otb1, taxa_are_rows = FALSE), sample_data(ps),tax_table(ps))
}


#' @title C++ Accelerated Rarefaction for Phyloseq Objects
#' @description
#' Performs rarefaction on a phyloseq object using a C++ backend via Rcpp.
#' This is faster than the R-based version and preserves OTU, sample, and taxonomic information.
#'
#' @param ps A \code{phyloseq} object.
#' @param N Integer. Rarefaction depth. Default 3000.
#' @param drop_insufficient Logical. Whether to drop samples with sequencing depth < N. Default TRUE.
#' @param rngseed Integer. Optional random seed for reproducibility.
#' @return A rarefied \code{phyloseq} object.
#' @export
#' @examples
#' ps_rarefied <- phyRare_cpp(ps, N = 3000, drop_insufficient = TRUE, rngseed = 123)
phyRare_cpp <- function(ps, N = 3000, drop_insufficient = TRUE, rngseed = NULL) {
  if (!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Please install the 'phyloseq' package.")
  }
  if (!is.null(rngseed)) set.seed(rngseed)

  # 提取 OTU 矩阵
  OTU <- phyloseq::otu_table(ps)
  taxa_rows <- phyloseq::taxa_are_rows(OTU)
  mat <- as.matrix(OTU)
  if (taxa_rows) mat <- t(mat)

  samp_names <- rownames(mat)
  taxa_names <- colnames(mat)

  if (drop_insufficient) {
    keep <- rowSums(mat, na.rm = TRUE) >= N
    if (!any(keep)) stop("No samples with depth >= N.")
    mat_sub <- mat[keep, , drop = FALSE]

    res <- rarefy_matrix(mat_sub, N)

    if (!is.null(samp_names)) rownames(res) <- samp_names[keep]
    if (!is.null(taxa_names)) colnames(res) <- taxa_names

    res_int <- matrix(as.integer(round(res)), nrow = nrow(res),
                      dimnames = dimnames(res))

    if (taxa_rows) {
      otu_tab <- phyloseq::otu_table(t(res_int), taxa_are_rows = TRUE)
    } else {
      otu_tab <- phyloseq::otu_table(res_int, taxa_are_rows = FALSE)
    }

    sd_sub <- phyloseq::sample_data(ps)[keep, , drop = FALSE]
    ps_out <- phyloseq::phyloseq(otu_tab, sd_sub)
  } else {
    res <- rarefy_matrix(mat, N)
    if (!is.null(samp_names)) rownames(res) <- samp_names
    if (!is.null(taxa_names)) colnames(res) <- taxa_names

    res_int <- matrix(as.integer(round(res)), nrow = nrow(res),
                      dimnames = dimnames(res))

    if (taxa_rows) {
      otu_tab <- phyloseq::otu_table(t(res_int), taxa_are_rows = TRUE)
    } else {
      otu_tab <- phyloseq::otu_table(res_int, taxa_are_rows = FALSE)
    }

    sd_all <- phyloseq::sample_data(ps)
    ps_out <- phyloseq::phyloseq(otu_tab, sd_all)

    warning("Some samples have depth < N and their OTU rows contain NA.")
  }
  if (!is.null(ps@tax_table)) {
    phyloseq::tax_table(ps_out) <- phyloseq::tax_table(ps)
  }
  return(ps_out)
}

# psRe = phyRare(ps = ps, N = 19141935,cores = 10)
# ps_rarefied <- phyloseq::rarefy_even_depth(ps, sample.size = 19141935, rngseed = 123)


# cppFunction('
# NumericMatrix rarefy_matrix(NumericMatrix mat, int N) {
#   int n_samples = mat.nrow();
#   int n_taxa = mat.ncol();
#   NumericMatrix result(n_samples, n_taxa);
#
#   for (int i = 0; i < n_samples; i++) {
#     double row_sum_d = 0;
#     for (int j = 0; j < n_taxa; j++) {
#       row_sum_d += mat(i, j);
#     }
#
#     int row_sum = (int) row_sum_d;
#
#     if (row_sum < N) {
#       // 深度不足，整行设为 NA
#       for (int j = 0; j < n_taxa; j++) {
#         result(i, j) = NA_REAL;
#       }
#     } else {
#       int remainingN = N;
#       int remainingTotal = row_sum;
#
#       for (int j = 0; j < n_taxa; j++) {
#         int count = (int) mat(i, j);
#         if (count == 0) {
#           result(i, j) = 0;
#         } else {
#           // Rf_rhyper(white, black, draws) 返回 double
#           int sampled = (int) Rf_rhyper(count, remainingTotal - count, remainingN);
#           result(i, j) = sampled;
#           remainingN -= sampled;
#           remainingTotal -= count;
#         }
#       }
#     }
#   }
#
#   return result;
# }
# ')



# # 假设 ps 是你的 phyloseq 对象
# ps_rarefied <- phyRare_auto(ps = ps.micro, N = 3000, drop_insufficient = TRUE, rngseed = 123)

