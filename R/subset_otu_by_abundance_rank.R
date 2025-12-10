#' Subset phyloseq object by OTU abundance rank interval
#'
#' @param ps A phyloseq object.
#' @param from Integer, starting rank (1 = most abundant OTU).
#' @param to Integer, ending rank (inclusive).
#' @param remove_zero Logical, whether to remove OTUs with total abundance = 0
#'        before ranking. Default TRUE.
#'
#' @return A phyloseq object containing only OTUs whose total abundance ranks
#'         are between `from` and `to` (after sorting in descending order).
#'
#' @examples
#' # Get OTUs ranked 500 to 1000 by total abundance
#' ps_sub <- subset_otu_by_abundance_rank(ps.16s, from = 500, to = 1000)
subset_otu_by_abundance_rank <- function(ps,
                                         from = 1,
                                         to   = 500,
                                         remove_zero = TRUE) {
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object.")
  }
  if (!is.numeric(from) || !is.numeric(to)) {
    stop("from and to must be numeric.")
  }

  from <- as.integer(from)
  to   <- as.integer(to)
  if (from < 1 || to < 1) {
    stop("from and to must be >= 1.")
  }
  if (from > to) {
    stop("from must be <= to.")
  }

  # Extract OTU table and ensure taxa are rows
  otu <- as.matrix(phyloseq::otu_table(ps))
  if (!phyloseq::taxa_are_rows(ps)) {
    otu <- t(otu)
  }

  # Compute total abundance per OTU across all samples
  total_abund <- rowSums(otu, na.rm = TRUE)

  # Optionally remove zero-sum OTUs before ranking
  if (remove_zero) {
    non_zero_idx <- total_abund > 0
    total_abund  <- total_abund[non_zero_idx]
  }

  # If no OTUs remain, return an empty pruned phyloseq
  if (length(total_abund) == 0) {
    warning("No OTUs with non-zero abundance. Returning empty phyloseq.")
    return(phyloseq::prune_taxa(character(0), ps))
  }

  # Order OTUs by abundance (descending)
  ord <- order(total_abund, decreasing = TRUE)

  # Adjust bounds to available OTUs
  n_otus <- length(ord)
  if (from > n_otus) {
    stop("from (", from, ") is greater than the number of available OTUs (", n_otus, ").")
  }
  to <- min(to, n_otus)

  # Select OTU names in the desired rank interval
  ranked_otus <- names(total_abund)[ord]
  selected_otus <- ranked_otus[from:to]

  # Subset phyloseq object
  ps_sub <- phyloseq::prune_taxa(selected_otus, ps)

  return(ps_sub)
}
