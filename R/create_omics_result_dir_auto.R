# utils_metabolome.R
# Utility functions for metabolomics-style phyloseq objects
# Internal helpers start with '.' and are not exported by default.
# Public functions (with roxygen2 docs) can be exported in NAMESPACE.

# -------------------------------------------------------------------
# Internal: detect whether a taxonomy table looks like metabolomics
# -------------------------------------------------------------------

.is_metabolome_tax <- function(tt) {
  #' Internal helper: detect metabolomics-style taxonomy table
  #'
  #' @param tt A taxonomy matrix (typically from phyloseq::tax_table(ps)).
  #' @return Logical scalar, TRUE if this looks like metabolomics annotations.

  if (is.null(tt) || ncol(tt) == 0) return(FALSE)
  tt <- as.matrix(tt)
  col_lower <- tolower(colnames(tt))

  # Typical metabolite annotation columns
  has_hmdb        <- any(grepl("hmdb", col_lower))
  has_cas         <- any(col_lower %in% c("cas", "cas.id", "cas_id"))
  has_kegg_comp   <- any(grepl("kegg", col_lower) & grepl("compound|cmpd", col_lower))
  has_formula     <- any(col_lower %in% c("formula", "mf", "molecular_formula"))
  has_super_class <- any(col_lower %in% c("super.class", "superclass"))
  has_type_col    <- any(col_lower == "type")  # POS / NEG

  metabolite_score <- sum(c(
    has_hmdb,
    has_cas,
    has_kegg_comp,
    has_formula,
    has_super_class,
    has_type_col
  ))

  # Typical biological taxonomy columns (amplicon-like)
  has_kingdom <- any(col_lower %in% c("kingdom", "k", "domain", "superkingdom"))
  has_phylum  <- any(col_lower == "phylum")

  # Simple rule:
  # - at least 2 metabolite-style columns
  # - and no kingdom/phylum
  if (metabolite_score >= 2 && !has_kingdom && !has_phylum) {
    return(TRUE)
  }
  FALSE
}

# -------------------------------------------------------------------
# Public: guess phyloseq data type (amplicon vs metabolome)
# -------------------------------------------------------------------

#' Guess the data type of a phyloseq object (amplicon vs metabolome)
#'
#' This function inspects the taxonomy table of a phyloseq object and
#' tries to guess whether it represents amplicon-style data (with
#' biological taxonomy such as Kingdom / Phylum) or metabolomics-style
#' data (with HMDB / CAS / KEGG compound IDs, etc.).
#'
#' @param ps A \code{phyloseq} object.
#'
#' @return A character scalar, one of:
#'   \itemize{
#'     \item \code{"metabolome"} – looks like metabolomics data;
#'     \item \code{"amplicon"}   – looks like amplicon microbiome data;
#'     \item \code{"unknown"}    – cannot be reliably classified.
#'   }
#'
#' @examples
#' \dontrun{
#'   guess_phyloseq_data_type(ps)
#' }
#'
#' @export
guess_phyloseq_data_type <- function(ps) {
  if (!inherits(ps, "phyloseq")) {
    stop("`ps` must be a phyloseq object.")
  }

  tt <- phyloseq::tax_table(ps)
  if (is.null(tt)) return("unknown")
  tt <- as.matrix(tt)
  col_lower <- tolower(colnames(tt))

  # 1) Metabolomics-style?
  if (.is_metabolome_tax(tt)) {
    return("metabolome")
  }

  # 2) Amplicon-style? (biological taxonomy)
  has_kingdom <- any(col_lower %in% c("kingdom", "k", "domain", "superkingdom"))
  has_phylum  <- any(col_lower == "phylum")
  has_genus   <- any(col_lower == "genus")

  if (has_kingdom || has_phylum || has_genus) {
    return("amplicon")
  }

  # 3) Unknown
  "unknown"
}

# -------------------------------------------------------------------
# Internal: build metabolome type label from tax_table
# -------------------------------------------------------------------
.build_metabolome_type_label <- function(ps) {
  #' Internal helper: build a simple metabolome type label
  #'
  #' Only keeps polarity info (POS / NEG), e.g. "POS", "NEG", "POS-NEG".
  #' If no "type" column is found, returns NA and the caller may decide
  #' to omit it in folder name.
  #'
  #' @param ps A phyloseq object.
  #' @return A character scalar like "POS", "NEG", "POS-NEG", or NA.

  tt <- phyloseq::tax_table(ps)
  tt <- as.matrix(tt)
  col_lower <- tolower(colnames(tt))

  type_col <- which(col_lower == "type")
  if (length(type_col) != 1) {
    return(NA_character_)
  }

  types <- tt[, type_col]
  types <- toupper(trimws(types))
  types <- unique(na.omit(types))

  if (length(types) == 0) {
    return(NA_character_)
  }

  # 只保留极性：POS / NEG / POS-NEG
  paste(types, collapse = "-")
}

# -------------------------------------------------------------------
# Public: create metabolome result directory
# -------------------------------------------------------------------

#' Create result directory for metabolomics-style phyloseq objects
#'
#' This function checks whether a given \code{phyloseq} object looks
#' like metabolomics data (via \code{guess_phyloseq_data_type()}) and
#' if so, constructs a result directory with a structured name and
#' creates it on disk.
#'
#' Directory name pattern (without time):
#'   \verb{<prefix>_<type_label>_<YYYYMMDD>}
#'
#' Example:
#'   \verb{metabolome_meta_POS_20251129}
#'
#' @param ps A \code{phyloseq} object that should represent metabolomics data.
#' @param base_dir Base directory in which to create the result folder.
#'   Default is \code{"./result"}.
#' @param prefix Prefix used in the folder name. Default is \code{"metabolome"}.
#' @param include_time Logical; if \code{TRUE}, append \code{HHMMSS} to
#'   the date part to ensure uniqueness.
#'
#' @return A character scalar: the path of the created directory.
#'
#' @examples
#' \dontrun{
#'   out_dir <- create_metabolome_result_dir(ps.meta, base_dir = "./result")
#' }
#'
#' @export
create_metabolome_result_dir <- function(ps,
                                         base_dir     = "./result",
                                         prefix       = "metabolome",
                                         include_time = FALSE) {
  data_type <- guess_phyloseq_data_type(ps)
  if (!identical(data_type, "metabolome")) {
    stop("Current phyloseq object does not look like metabolomics data ",
         "(guess_phyloseq_data_type(ps) != 'metabolome').")
  }

  # 1. Build type label (POS/NEG + Super.Class)
  type_label <- .build_metabolome_type_label(ps)

  # 2. Date part, e.g. 20251129
  date_part <- format(Sys.Date(), "%Y%m%d")

  # 3. Optional time part, e.g. 103522
  if (isTRUE(include_time)) {
    time_part      <- format(Sys.time(), "%H%M%S")
    datetime_label <- paste(date_part, time_part, sep = "_")
  } else {
    datetime_label <- date_part
  }

  # 4. Final folder name:
  #    without time: metabolome_meta_POS_20251129
  #    with time:    metabolome_meta_POS_20251129_103522
  folder_name <- paste(prefix, type_label, datetime_label, sep = "_")
  out_path    <- file.path(base_dir, folder_name)

  # 5. Create directory, prefer user-defined dir_create() if available
  if (exists("dir_create")) {
    dir_create(out_path)
  } else if (requireNamespace("fs", quietly = TRUE)) {
    fs::dir_create(out_path, recurse = TRUE)
  } else {
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  }

  out_path
}


