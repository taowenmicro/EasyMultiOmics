#' Fast ANCOM-II Analysis for Microbiome Data
#'
#' @description
#' This function performs a fast implementation of ANCOM-II (Analysis of Composition of Microbiomes)
#' for identifying differentially abundant taxa between two groups using pairwise log-ratio tests.
#'
#' @param ps A phyloseq object containing OTU table and sample metadata
#' @param group Character string specifying the grouping variable in sample metadata. Default is "Group"
#' @param alpha Numeric value for significance threshold. Default is 0.05
#' @param detect.cut Numeric value for detection cutoff. Taxa with W > detect.cut * (m-1) are considered significant. Default is 0.6
#' @param test Character string specifying the statistical test to use. Options are "t" (t-test) or "wilcox" (Wilcoxon test). Default is "t"
#' @param workers Integer specifying number of parallel workers. Default uses all available cores minus 1
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE
#'
#' @return A list containing:
#' \describe{
#'   \item{tab.ANCOM.fast}{Data frame with taxa_id, W statistic, mean p-value, and detection status}
#'   \item{diff.tab}{Data frame with significantly different taxa including method and adjusted p-values}
#' }
#'
#' @details
#' The function implements ANCOM-II using centered log-ratio (CLR) transformation followed by
#' pairwise log-ratio tests between all taxa pairs. The W statistic represents the number of
#' significant pairwise comparisons for each taxon. Taxa are considered differentially abundant
#' if W > detect.cut * (total_taxa - 1).
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- ancom.micro.fast(ps, group = "Treatment")
#'
#' # With custom parameters
#' result <- ancom.micro.fast(ps,
#'                           group = "Disease_Status",
#'                           alpha = 0.01,
#'                           detect.cut = 0.7,
#'                           test = "wilcox")
#'
#' # View results
#' head(result$diff.tab)
#' }
#'
#' @export
ancom.micro.fast <- function(ps,
                             group      = "Group",
                             alpha      = 0.05,
                             detect.cut = 0.6,
                             test       = c("t", "wilcox"),
                             workers    = max(1, parallel::detectCores()-1),
                             verbose    = TRUE) {

  test <- match.arg(test)

  # Data preparation
  mat <- as(phyloseq::otu_table(ps), "matrix")
  if (phyloseq::taxa_are_rows(ps)) {
    mat <- t(mat)
  }
  mat <- as.matrix(mat)

  meta <- as(phyloseq::sample_data(ps), "data.frame")
  common.samples <- intersect(rownames(meta), rownames(mat))
  if (length(common.samples) == 0L) stop("No overlapping samples between metadata and OTU table.")

  meta <- meta[common.samples, , drop = FALSE]
  mat  <- mat [common.samples, , drop = FALSE]

  g <- meta[[group]]
  if (is.null(g)) stop(sprintf("Grouping column '%s' not found in sample_data.", group))
  g <- droplevels(factor(g))

  if (length(levels(g)) != 2) {
    stop("This function only supports two-group comparisons")
  }

  keep <- !is.na(g)
  if (!all(keep)) {
    if (verbose) message("Removed ", sum(!keep), " samples with NA group.")
    g   <- droplevels(g[keep])
    mat <- mat[keep, , drop = FALSE]
  }

  # CLR transformation
  log_mat <- log(mat + 1)
  clr_mat <- log_mat - rowMeans(log_mat)

  taxa_names <- colnames(clr_mat)
  n_taxa <- ncol(clr_mat)
  if (n_taxa < 2L) stop("Need at least 2 taxa for pairwise log-ratio tests.")

  # Pre-compute group information
  g_levels <- levels(g)
  idx1 <- which(g == g_levels[1])
  idx2 <- which(g == g_levels[2])
  n1 <- length(idx1)
  n2 <- length(idx2)

  if (verbose) {
    message(sprintf("Computing W statistics for %d taxa...", n_taxa))
    start_time <- Sys.time()
  }

  # Core computation functions with p-value collection
  if (test == "t") {
    # Optimized t-test: compute all pairwise comparisons
    calc_W_ultra_fast <- function(taxa_indices) {
      n_indices <- length(taxa_indices)
      W_values <- integer(n_indices)
      mean_p_values <- numeric(n_indices)

      for (batch_idx in seq_along(taxa_indices)) {
        i <- taxa_indices[batch_idx]

        # Current taxa CLR values
        xi <- clr_mat[, i]

        # Calculate log-ratios with all other taxa
        other_indices <- setdiff(seq_len(n_taxa), i)
        n_others <- length(other_indices)

        if (n_others == 0) {
          W_values[batch_idx] <- 0L
          mean_p_values[batch_idx] <- 1.0
          next
        }

        # Pre-allocate result vector
        p_values <- numeric(n_others)

        # Batch computation: vectorized comparisons
        xi_g1 <- xi[idx1]
        xi_g2 <- xi[idx2]

        for (k in seq_len(n_others)) {
          j <- other_indices[k]

          # Calculate log-ratio
          xj <- clr_mat[, j]
          lr_g1 <- xi_g1 - xj[idx1]
          lr_g2 <- xi_g2 - xj[idx2]

          # Fast t-test calculation
          if (n1 < 2 || n2 < 2) {
            p_values[k] <- 1
            next
          }

          m1 <- mean(lr_g1)
          m2 <- mean(lr_g2)

          # Fast variance calculation
          v1 <- sum((lr_g1 - m1)^2) / (n1 - 1)
          v2 <- sum((lr_g2 - m2)^2) / (n2 - 1)

          if (v1 == 0 && v2 == 0) {
            p_values[k] <- ifelse(abs(m1 - m2) < .Machine$double.eps, 1, 0)
            next
          }

          se <- sqrt(v1/n1 + v2/n2)
          if (se == 0) {
            p_values[k] <- 1
            next
          }

          t_stat <- (m1 - m2) / se
          df <- (v1/n1 + v2/n2)^2 / ((v1/n1)^2/(n1-1) + (v2/n2)^2/(n2-1))

          p_values[k] <- 2 * pt(-abs(t_stat), df)
        }

        # Calculate W statistic and mean p-value
        W_values[batch_idx] <- sum(p_values < alpha, na.rm = TRUE)
        mean_p_values[batch_idx] <- mean(p_values, na.rm = TRUE)
      }

      return(list(W = W_values, mean_p = mean_p_values))
    }
  } else {
    # Wilcoxon version with p-value collection
    calc_W_wilcox <- function(taxa_indices) {
      n_indices <- length(taxa_indices)
      W_values <- integer(n_indices)
      mean_p_values <- numeric(n_indices)

      for (batch_idx in seq_along(taxa_indices)) {
        i <- taxa_indices[batch_idx]
        wi <- 0L
        xi <- clr_mat[, i]
        p_values <- numeric()

        for (j in seq_len(n_taxa)) {
          if (j == i) next
          lr <- xi - clr_mat[, j]

          p <- tryCatch({
            wilcox.test(lr[idx1], lr[idx2], exact = FALSE)$p.value
          }, error = function(e) 1)

          if (!is.na(p)) {
            p_values <- c(p_values, p)
            if (p < alpha) wi <- wi + 1L
          }
        }

        W_values[batch_idx] <- wi
        mean_p_values[batch_idx] <- if(length(p_values) > 0) mean(p_values) else 1.0
      }

      return(list(W = W_values, mean_p = mean_p_values))
    }
  }

  # Select computation function
  calc_W_func <- if (test == "t") calc_W_ultra_fast else calc_W_wilcox

  # Parallel strategy: use only for large datasets
  if (workers > 1 && n_taxa > 200) {
    # Reasonable chunking to avoid excessive communication overhead
    chunk_size <- max(20, ceiling(n_taxa / workers))
    taxa_chunks <- split(seq_len(n_taxa), ceiling(seq_len(n_taxa) / chunk_size))

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(workers)
      on.exit(parallel::stopCluster(cl), add = TRUE)

      # Export required variables and functions
      parallel::clusterExport(cl, c("clr_mat", "idx1", "idx2", "n1", "n2",
                                    "alpha", "n_taxa", "calc_W_func"),
                              envir = environment())

      results_chunks <- parallel::parLapply(cl, taxa_chunks, calc_W_func)
    } else {
      results_chunks <- parallel::mclapply(taxa_chunks, calc_W_func, mc.cores = workers)
    }

    # Merge results
    W.stat <- unlist(lapply(results_chunks, function(x) x$W), use.names = FALSE)
    mean_p.stat <- unlist(lapply(results_chunks, function(x) x$mean_p), use.names = FALSE)

  } else {
    # Single-thread processing
    results <- calc_W_func(seq_len(n_taxa))
    W.stat <- results$W
    mean_p.stat <- results$mean_p
  }

  if (verbose) {
    end_time <- Sys.time()
    message(sprintf("Computation completed in %.2f seconds",
                    as.numeric(end_time - start_time, units = "secs")))
  }

  # Significance determination
  cutoff <- detect.cut * (n_taxa - 1)
  detected <- W.stat > cutoff

  res <- data.frame(
    taxa_id  = taxa_names,
    W        = as.integer(W.stat),
    mean_p   = mean_p.stat,
    detected = detected,
    stringsAsFactors = FALSE
  )

  diff.tab <- if (any(detected)) {
    data.frame(
      micro    = taxa_names[detected],
      method   = "ANCOMII.fast",
      adjust.p = mean_p.stat[detected],
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(micro=character(0), method=character(0), adjust.p=numeric(0))
  }

  if (verbose) {
    message(sprintf(
      "ANCOMII.fast: taxa=%d, alpha=%.3f, detect.cut=%.2f -> detected=%d (cutoff=%.1f)",
      n_taxa, alpha, detect.cut, sum(detected), cutoff
    ))
  }

  list(tab.ANCOM.fast = res, diff.tab = diff.tab)
}
