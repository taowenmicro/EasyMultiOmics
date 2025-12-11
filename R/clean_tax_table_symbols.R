#' Clean taxonomy table characters in a phyloseq object
#'
#' This function cleans the `tax_table` slot of a `phyloseq` object by
#' replacing disallowed characters with a given replacement symbol.
#' By default, all characters that are not in `[A-Za-z0-9_.]` are replaced,
#' and the dash (`-`) is explicitly converted to the replacement symbol.
#' `NA` values are preserved.
#'
#' @param ps A \code{phyloseq} object.
#' @param replacement A single character string used to replace disallowed
#'   characters. Default is \code{"."}.
#' @param allowed_chars A single string defining the set of allowed characters
#'   inside a regular-expression character class (without the surrounding
#'   square brackets). The default \code{"A-Za-z0-9_."} allows letters,
#'   digits, underscore, and dot.
#' @param which_ranks Optional character vector of taxonomy ranks (column
#'   names of \code{tax_table(ps)}) to be cleaned. If \code{NULL} (default),
#'   all columns in the taxonomy table are cleaned.
#' @param verbose Logical; if \code{TRUE}, print a short summary of how many
#'   entries were modified. Default is \code{FALSE}.
#'
#' @return The same \code{phyloseq} object as input, with a cleaned
#'   \code{tax_table}.
#'
#' @examples
#' \dontrun{
#'   library(phyloseq)
#'   data("GlobalPatterns")
#'
#'   # Clean all taxonomy ranks
#'   gp_clean <- clean_tax_table_symbols(GlobalPatterns)
#'
#'   # Only clean Genus and Species columns
#'   gp_clean2 <- clean_tax_table_symbols(
#'     GlobalPatterns,
#'     which_ranks = c("Genus", "Species"),
#'     replacement = "_",
#'     verbose = TRUE
#'   )
#'
#'   head(as.data.frame(tax_table(gp_clean2)))
#' }
#'
#' @export
clean_tax_table_symbols <- function(ps,
                                    replacement   = ".",
                                    allowed_chars = "A-Za-z0-9_.",
                                    which_ranks   = NULL,
                                    verbose       = FALSE) {
  # ---- basic checks ---------------------------------------------------------
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be a phyloseq object.", call. = FALSE)
  }

  tax <- phyloseq::tax_table(ps)
  if (is.null(tax)) {
    stop("`ps` does not contain a tax_table slot.", call. = FALSE)
  }

  tax_mat <- as.matrix(tax)

  # determine which columns (ranks) to clean
  if (is.null(which_ranks)) {
    cols_to_clean <- colnames(tax_mat)
  } else {
    missing_ranks <- setdiff(which_ranks, colnames(tax_mat))
    if (length(missing_ranks) > 0L) {
      warning(
        "The following ranks were not found in tax_table(ps) and will ",
        "be ignored: ",
        paste(missing_ranks, collapse = ", ")
      )
    }
    cols_to_clean <- intersect(which_ranks, colnames(tax_mat))
    if (length(cols_to_clean) == 0L) {
      # nothing to do
      if (verbose) {
        message("No valid ranks to clean; returning input object unchanged.")
      }
      return(ps)
    }
  }

  # precompute pattern: replace any char NOT in allowed_chars
  pattern <- sprintf("[^%s]", allowed_chars)

  # helper to clean one character vector
  clean_vec <- function(x) {
    x_chr <- as.character(x)
    na_idx <- is.na(x_chr)

    # 1) replace "-" with replacement
    x_chr[!na_idx] <- gsub("-", replacement, x_chr[!na_idx], fixed = TRUE)

    # 2) replace any character outside allowed_chars with replacement
    x_chr[!na_idx] <- gsub(pattern, replacement, x_chr[!na_idx])

    x_chr
  }

  # keep a copy for change statistics (optional)
  original_mat <- tax_mat

  # clean selected columns
  for (col in cols_to_clean) {
    tax_mat[, col] <- clean_vec(tax_mat[, col])
  }

  # restore dimnames (as.matrix already keeps them, but be explicit)
  rownames(tax_mat) <- rownames(tax)
  colnames(tax_mat) <- colnames(tax)

  # write back to phyloseq object
  phyloseq::tax_table(ps) <- phyloseq::tax_table(tax_mat)

  # summary
  if (verbose) {
    changed <- original_mat != tax_mat
    changed[is.na(original_mat) & is.na(tax_mat)] <- FALSE
    n_changed <- sum(changed, na.rm = TRUE)

    message(
      "clean_tax_table_symbols: cleaned ",
      length(cols_to_clean), " rank(s); ",
      n_changed, " cell(s) were modified."
    )
  }

  ps
}
