# -------------------------------------------------------------------
# Internal helper: classify kingdom/domain strings into coarse types
# -------------------------------------------------------------------

#' Classify kingdom/domain into coarse organism types
#'
#' This is an internal helper used by \code{get_amplicon_type_ratio()}.
#'
#' @param x Character vector of kingdom/domain annotations.
#'
#' @return A character vector of the same length with values in
#'   \code{"bacteria"}, \code{"fungi"}, \code{"virus"}, \code{"protist"},
#'   or \code{NA}.
#' @keywords internal
.classify_kingdom_to_type <- function(x) {
  x <- gsub("^[a-z]__", "", x)
  x <- trimws(x)

  res <- rep(NA_character_, length(x))

  # bacteria / archaea
  res[grepl("bacteria|archaea", x, ignore.case = TRUE)] <- "bacteria"

  # fungi
  res[grepl("fungi|ascomycota|basidiomycota", x, ignore.case = TRUE)] <- "fungi"

  # virus
  res[grepl("virus|viridae|phage", x, ignore.case = TRUE)] <- "virus"

  # protist / protozoa
  res[grepl("protist|protozoa|ciliophora|apicomplexa|amoebozoa",
            x, ignore.case = TRUE)] <- "protist"

  res
}

# -------------------------------------------------------------------
# Estimate abundance ratio of coarse organism types (Bac/Fun/Prot/Vir)
# -------------------------------------------------------------------

#' Estimate abundance ratio of coarse amplicon types
#'
#' This helper computes the total read abundance of different organism
#' types (bacteria, fungi, protist, virus) based on the taxonomy table
#' in a phyloseq object, and returns their relative contributions as
#' integer percentages (rounded).
#'
#' If taxonomy is missing or no types can be classified, it returns a
#' single-element vector \code{c(amplicon = 100)} as a generic label.
#'
#' @param ps A phyloseq object.
#'
#' @return An integer named vector giving percentage abundance of each
#'   type; names are a subset of \code{c("bacteria", "fungi", "protist", "virus")}
#'   or \code{"amplicon"}.
#' @keywords internal
get_amplicon_type_ratio <- function(ps) {
  otu <- phyloseq::otu_table(ps)
  otu_mat <- as(otu, "matrix")

  if (!phyloseq::taxa_are_rows(otu)) {
    otu_mat <- t(otu_mat)
  }

  tt <- phyloseq::tax_table(ps)
  if (is.null(tt)) {
    return(structure(100L, names = "amplicon"))
  }
  tt <- as.matrix(tt)

  col_lower <- tolower(colnames(tt))
  kingdom_cols <- which(col_lower %in% c("kingdom", "k", "domain", "superkingdom"))

  if (length(kingdom_cols) == 0) {
    kingdom_cols <- seq_len(ncol(tt))
  }

  kingdom_raw <- apply(tt[, kingdom_cols, drop = FALSE], 1L, function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0) return(NA_character_)
    v[1]
  })

  taxon_type <- .classify_kingdom_to_type(kingdom_raw)

  type_order <- c("bacteria", "fungi", "protist", "virus")
  abund_vec  <- numeric(length(type_order))
  names(abund_vec) <- type_order

  for (tp in type_order) {
    sel <- which(taxon_type == tp)
    if (length(sel) > 0) {
      abund_vec[tp] <- sum(otu_mat[sel, , drop = FALSE])
    }
  }

  total_abund <- sum(abund_vec)
  if (total_abund <= 0) {
    return(structure(100L, names = "amplicon"))
  }

  pct <- round(abund_vec / total_abund * 100)
  pct <- pct[pct > 0]

  pct
}

# -------------------------------------------------------------------
# Build human-readable label for amplicon type ratios
# -------------------------------------------------------------------

#' Build amplicon type label from percentage vector
#'
#' Given a named percentage vector for organism types (typically the
#' output of \code{get_amplicon_type_ratio()}), construct a compact label.
#'
#' Examples:
#' \itemize{
#'   \item \code{c(bacteria = 100)} -> \code{"bacteria"}
#'   \item \code{c(bacteria = 60, fungi = 40)} -> \code{"multi_Bac60-Fun40"}
#' }
#'
#' @param type_pct Integer named vector of percentages.
#'
#' @return A single character scalar label.
#' @keywords internal
build_amplicon_type_label_with_ratio <- function(type_pct) {
  if (length(type_pct) == 1) {
    nm <- names(type_pct)
    if (identical(nm, "amplicon")) {
      return("amplicon_unclassified")
    }
    return(nm)
  }

  short_map <- c(
    bacteria = "Bac",
    fungi    = "Fun",
    protist  = "Pro",
    virus    = "Vir"
  )

  type_order <- c("bacteria", "fungi", "protist", "virus")
  keep <- intersect(type_order, names(type_pct))
  type_pct <- type_pct[keep]

  if (length(type_pct) == 0) {
    return("amplicon_unclassified")
  }

  parts <- paste0(short_map[names(type_pct)], type_pct)

  paste0("multi_", paste(parts, collapse = "-"))
}

# -------------------------------------------------------------------
# Create result directory for amplicon data
# -------------------------------------------------------------------

#' Create result directory for amplicon data
#'
#' This function creates a result directory for amplicon-based analyses.
#' The directory name encodes a prefix (e.g. \code{"amplicon"}), a
#' type label derived from organism-type ratios (e.g.
#' \code{"multi_Bac60-Fun40"}), and the current date (and optionally
#' time).
#'
#' @param ps A phyloseq object.
#' @param base_dir Base directory for all results (default: \code{"./result"}).
#' @param prefix Prefix for the folder name (default: \code{"amplicon"}).
#' @param include_time Logical; if \code{TRUE}, append time (HHMMSS) to
#'   the folder name.
#'
#' @return A character scalar giving the path of the created directory.
#' @export
create_amplicon_result_dir <- function(ps,
                                       base_dir     = "./result",
                                       prefix       = "amplicon",
                                       include_time = FALSE) {
  type_pct <- get_amplicon_type_ratio(ps)
  type_label <- build_amplicon_type_label_with_ratio(type_pct)

  date_part <- format(Sys.Date(), "%Y%m%d")

  if (isTRUE(include_time)) {
    time_part <- format(Sys.time(), "%H%M%S")
    datetime_label <- paste(date_part, time_part, sep = "_")
  } else {
    datetime_label <- date_part
  }

  folder_name <- paste(prefix, type_label, datetime_label, sep = "_")
  out_path    <- file.path(base_dir, folder_name)

  if (exists("dir_create")) {
    dir_create(out_path)
  } else if (requireNamespace("fs", quietly = TRUE)) {
    fs::dir_create(out_path, recurse = TRUE)
  } else {
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }

  out_path
}

# -------------------------------------------------------------------
# Infer sequencing type (amplicon vs metagenome)
# -------------------------------------------------------------------

#' Infer sequencing type (amplicon vs. metagenome) from a phyloseq object
#'
#' This function uses simple heuristics to guess whether a phyloseq
#' object likely comes from amplicon sequencing (e.g. 16S/ITS) or
#' metagenomic sequencing.
#'
#' Heuristics include:
#' \itemize{
#'   \item Library size distribution (\code{sample_sums})
#'   \item Feature ID patterns (ASV-like long A/C/G/T sequences)
#'   \item Number and composition of kingdoms in \code{tax_table}
#'   \item Presence of functional annotation columns (KO, KEGG, eggNOG, etc.)
#' }
#'
#' @param ps A phyloseq object.
#' @param verbose Logical; if \code{TRUE}, print a short summary.
#'
#' @return An object of class \code{"psSeqType"} with elements:
#'   \describe{
#'     \item{label}{One of \code{"amplicon"}, \code{"metagenome"},
#'       or \code{"ambiguous"}.}
#'     \item{amp_score}{Score supporting "amplicon".}
#'     \item{meta_score}{Score supporting "metagenome".}
#'     \item{metrics}{List of intermediate metrics used for scoring.}
#'     \item{call}{The original function call.}
#'   }
#'
#' @examples
#' \dontrun{
#'   res <- infer_sequencing_type(ps)
#'   res$label
#' }
#' @export
infer_sequencing_type <- function(ps, verbose = TRUE) {
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object.")
  }

  otutab    <- phyloseq::otu_table(ps)
  taxa_ids  <- phyloseq::taxa_names(ps)
  depths    <- phyloseq::sample_sums(ps)

  tax <- phyloseq::tax_table(ps, errorIfNULL = FALSE)
  if (!is.null(tax)) {
    tax <- as.matrix(tax)
  }

  n_samples  <- length(depths)
  n_features <- length(taxa_ids)

  depth_median <- stats::median(depths)
  depth_log10  <- if (depth_median > 0) unname(log10(depth_median)) else NA_real_

  is_seq_like <- function(x) {
    grepl("^[ACGT]+$", x)
  }
  id_len        <- nchar(taxa_ids)
  seq_like      <- is_seq_like(taxa_ids) & id_len >= 100
  prop_seq_like <- if (length(seq_like) > 0) mean(seq_like) else NA_real_

  n_kingdoms    <- NA_integer_
  kingdom_tab   <- NULL
  major_kingdoms <- NULL

  if (!is.null(tax)) {
    king_col <- which(colnames(tax) %in%
                        c("Kingdom", "kingdom", "Domain", "domain", "Superkingdom"))
    if (length(king_col) > 0) {
      king <- tax[, king_col[1]]
      king[king == "" | is.na(king)] <- NA
      king_tab       <- sort(table(king), decreasing = TRUE)
      kingdom_tab    <- king_tab
      n_kingdoms     <- length(king_tab)
      if (n_kingdoms > 0) {
        major_kingdoms <- names(king_tab)[seq_len(min(5L, n_kingdoms))]
      }
    }
  }

  has_fun_cols <- FALSE
  fun_cols     <- character(0)

  if (!is.null(tax)) {
    fun_candidates <- c(
      "KO", "KEGG", "KEGG_ID", "Kegg", "KO_ID",
      "eggNOG", "COG", "EC", "PFAM", "pfam", "gene", "Gene"
    )
    fun_cols     <- intersect(colnames(tax), fun_candidates)
    has_fun_cols <- length(fun_cols) > 0
  }

  amp_score  <- 0L
  meta_score <- 0L

  if (!is.na(depth_log10)) {
    if (depth_log10 <= 5.2) amp_score  <- amp_score  + 1L
    if (depth_log10 >= 5.8) meta_score <- meta_score + 1L
  }

  if (!is.na(prop_seq_like)) {
    if (prop_seq_like >= 0.3) amp_score  <- amp_score  + 2L
    if (prop_seq_like <= 0.05) meta_score <- meta_score + 1L
  }

  if (!is.na(n_kingdoms)) {
    if (n_kingdoms >= 3L) {
      meta_score <- meta_score + 2L
    } else if (n_kingdoms <= 2L) {
      amp_score <- amp_score + 1L
    }

    if (!is.null(kingdom_tab) && length(kingdom_tab) > 0) {
      main_king <- names(kingdom_tab)[1]
      main_prop <- kingdom_tab[1] / sum(kingdom_tab)

      if (!is.na(main_prop) &&
          main_prop >= 0.95 &&
          main_king %in% c("Bacteria", "Fungi")) {
        amp_score <- amp_score + 1L
      }

      if (any(names(kingdom_tab) %in% c("Viruses", "Virus", "Eukaryota", "Eukarya"))) {
        meta_score <- meta_score + 1L
      }
      if (any(names(kingdom_tab) %in% c("Archaea", "Archaeota"))) {
        meta_score <- meta_score + 1L
      }
    }
  }

  if (has_fun_cols) {
    meta_score <- meta_score + 2L
  }

  label <- "ambiguous"
  if (amp_score - meta_score >= 2L) {
    label <- "amplicon"
  } else if (meta_score - amp_score >= 2L) {
    label <- "metagenome"
  }

  metrics <- list(
    n_samples        = n_samples,
    n_features       = n_features,
    depth_median     = depth_median,
    depth_log10      = depth_log10,
    prop_seq_like_id = prop_seq_like,
    n_kingdoms       = n_kingdoms,
    kingdom_table    = kingdom_tab,
    major_kingdoms   = major_kingdoms,
    has_fun_cols     = has_fun_cols,
    fun_cols         = fun_cols
  )

  res <- list(
    label      = label,
    amp_score  = amp_score,
    meta_score = meta_score,
    metrics    = metrics,
    call       = match.call()
  )
  class(res) <- "psSeqType"

  if (verbose) {
    print(res)
  }

  invisible(res)
}

# -------------------------------------------------------------------
# Print method for psSeqType
# -------------------------------------------------------------------

#' Print method for \code{psSeqType} objects
#'
#' @param x A \code{psSeqType} object.
#' @param ... Additional arguments (currently ignored).
#'
#' @export
#' @method print psSeqType
print.psSeqType <- function(x, ...) {
  cat("Inferred sequencing type:", x$label, "\n")
  cat("  amp_score  :", x$amp_score, "\n")
  cat("  meta_score :", x$meta_score, "\n\n")

  m <- x$metrics

  cat("Basic metrics:\n")
  cat("  samples       :", m$n_samples, "\n")
  cat("  features      :", m$n_features, "\n")
  cat("  median depth  :", format(m$depth_median, digits = 3), "\n")
  cat("  log10(depth)  :", format(m$depth_log10, digits = 3), "\n")
  cat("  seq-like IDs  :", sprintf("%.1f%%", 100 * m$prop_seq_like_id), "\n")

  cat("\nTaxonomy metrics:\n")
  cat("  n_kingdoms    :", m$n_kingdoms, "\n")
  if (!is.null(m$major_kingdoms)) {
    cat("  major kingdoms:", paste(m$major_kingdoms, collapse = ", "), "\n")
  }

  cat("\nFunctional annotation:\n")
  cat("  has_fun_cols  :", m$has_fun_cols, "\n")
  if (m$has_fun_cols) {
    cat("  fun_cols      :", paste(m$fun_cols, collapse = ", "), "\n")
  }
  invisible(x)
}

# -------------------------------------------------------------------
# Auto-create result directory based on sequencing type
# -------------------------------------------------------------------

#' Automatically create result directory based on sequencing type
#'
#' If the phyloseq object is inferred as metagenome, a folder like
#' \code{"metagenome_YYYYMMDD"} (or with time) will be created.
#' If inferred as amplicon (or ambiguous), it will call
#' \code{create_amplicon_result_dir()} to create a folder such as
#' \code{"amplicon_multi_Bac60-Fun40_YYYYMMDD"}.
#'
#' @param ps A phyloseq object.
#' @param base_dir Base directory for all results (default: \code{"./result"}).
#' @param include_time Logical; if \code{TRUE}, append time (HHMMSS) to
#'   folder name.
#' @param metagenome_prefix Prefix for metagenome folders
#'   (default: \code{"metagenome"}).
#' @param amplicon_prefix Prefix for amplicon folders
#'   (default: \code{"amplicon"}).
#' @param verbose Logical; if \code{TRUE}, print inferred sequencing type.
#'
#' @return A character scalar giving the path of the created directory.
#' @export
create_omics_result_dir_auto <- function(ps,
                                         base_dir          = "./result",
                                         include_time      = FALSE,
                                         metagenome_prefix = "metagenome",
                                         amplicon_prefix   = "amplicon",
                                         verbose           = TRUE) {
  seq_type <- infer_sequencing_type(ps, verbose = FALSE)
  if (verbose) {
    message("Inferred sequencing type: ", seq_type$label,
            " (amp_score = ", seq_type$amp_score,
            ", meta_score = ", seq_type$meta_score, ")")
  }

  if (identical(seq_type$label, "metagenome")) {
    date_part <- format(Sys.Date(), "%Y%m%d")
    if (isTRUE(include_time)) {
      time_part <- format(Sys.time(), "%H%M%S")
      folder_name <- paste(metagenome_prefix, date_part, time_part, sep = "_")
    } else {
      folder_name <- paste(metagenome_prefix, date_part, sep = "_")
    }

    out_path <- file.path(base_dir, folder_name)

    if (exists("dir_create")) {
      dir_create(out_path)
    } else if (requireNamespace("fs", quietly = TRUE)) {
      fs::dir_create(out_path, recurse = TRUE)
    } else {
      dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
    }

    return(out_path)
  }

  out_path <- create_amplicon_result_dir(
    ps           = ps,
    base_dir     = base_dir,
    prefix       = amplicon_prefix,
    include_time = include_time
  )

  out_path
}
