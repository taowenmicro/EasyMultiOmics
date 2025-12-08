#' Infer a matrix name from a call expression
#'
#' This helper tries to infer a "pretty" name for a matrix argument based on how
#' it was passed to a function. For example:
#' \itemize{
#'   \item \code{mat[[1]]} uses \code{names(mat)[1]} if available.
#'   \item \code{mat[["KO"]]} becomes \code{"KO"}.
#'   \item \code{A} becomes \code{"A"}.
#' }
#'
#' This function is primarily intended for internal use in other functions
#' (e.g. \code{compare_corr_two()}, \code{align_corr_pair()}, \code{get_network_nodes()}),
#' so that output objects can carry meaningful labels.
#'
#' @param expr_chr Character string of the expression (typically from
#'   \code{deparse(substitute(x))}).
#' @param env Environment in which to look up objects (usually
#'   \code{parent.frame()}).
#'
#' @return A single character string representing the inferred name.
#'
#' @keywords internal
.get_matrix_name <- function(expr_chr, env = parent.frame()) {
  expr_chr <- expr_chr[1]

  ## case 1: mat[["KO"]]
  m_named <- regexec('^([A-Za-z0-9_.]+)\\[\\["([^"]+)"\\]\\]$', expr_chr)
  r_named <- regmatches(expr_chr, m_named)[[1]]
  if (length(r_named) == 3) {
    inner_nm <- r_named[3]
    return(inner_nm)
  }

  ## case 2: mat[[1]]
  m_idx <- regexec('^([A-Za-z0-9_.]+)\\[\\[(\\d+)\\]\\]$', expr_chr)
  r_idx <- regmatches(expr_chr, m_idx)[[1]]
  if (length(r_idx) == 3) {
    base_obj <- r_idx[2]
    idx      <- as.integer(r_idx[3])
    if (exists(base_obj, envir = env)) {
      obj <- get(base_obj, envir = env)
      if (!is.null(names(obj))) {
        nm <- names(obj)
        if (idx <= length(nm) && !is.na(nm[idx]) && nzchar(nm[idx])) {
          return(nm[idx])
        }
      }
    }
    return(expr_chr)
  }

  ## default: use raw expression string
  expr_chr
}


#' Extract nodes that participate in a correlation network
#'
#' Given a square correlation (or adjacency) matrix, this function returns all
#' node IDs that participate in at least one non-zero, non-\code{NA} edge.
#' The matrix is assumed to represent a simple undirected network with the same
#' node set in rows and columns.
#'
#' If \code{as_list = FALSE} (default), a character vector of node names is
#' returned, with an attribute \code{"mat_name"} storing an inferred matrix name
#' based on the call (e.g. \code{"WT"}, \code{"OE"}).
#'
#' If \code{as_list = TRUE}, the result is wrapped as a one-element list, whose
#' name is the inferred matrix name.
#'
#' @param mat A square numeric matrix with row and column names representing
#'   node IDs. Row names and column names must be identical and in the same
#'   order.
#' @param as_list Logical; if \code{TRUE}, return a named list
#'   \code{list(<matrix_name> = node_ids)}; if \code{FALSE}, return a character
#'   vector with attribute \code{"mat_name"}.
#'
#' @return Either a character vector of node IDs (with attribute
#'   \code{"mat_name"}) or a one-element list if \code{as_list = TRUE}.
#'
#' @examples
#' mat <- matrix(c(1, 0.5, 0, 0.5, 1, 0.2, 0, 0.2, 1),
#'               nrow = 3, byrow = TRUE,
#'               dimnames = list(c("A","B","C"), c("A","B","C")))
#' get_network_nodes(mat)
#'
#' @export
get_network_nodes <- function(mat, as_list = FALSE) {
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix must have rownames and colnames.")
  }
  if (!identical(rownames(mat), colnames(mat))) {
    stop("Rownames and colnames must be identical.")
  }

  ## infer matrix name (e.g. "WT", "OE", "A", etc.)
  expr_mat <- deparse(substitute(mat))
  mname    <- .get_matrix_name(expr_mat, env = parent.frame())

  nodes_all <- rownames(mat)

  ## non-NA and non-zero entries are treated as edges
  edge_mask <- !is.na(mat) & (mat != 0)
  ## remove diagonal
  diag(edge_mask) <- FALSE

  ## whether each row/column participates in at least one edge
  involved_row <- apply(edge_mask, 1, any)
  involved_col <- apply(edge_mask, 2, any)
  involved <- involved_row | involved_col

  nodes <- nodes_all[involved]

  if (as_list) {
    out <- list(nodes)
    names(out) <- mname
    return(out)
  } else {
    attr(nodes, "mat_name") <- mname
    return(nodes)
  }
}


#' Align two correlation matrices to a common variable set
#'
#' This function takes two square matrices (e.g. correlation matrices) and
#' aligns them to the same set of variables by either intersection or union of
#' their row/column names. The resulting matrices share identical row and
#' column names and ordering, making them directly comparable.
#'
#' Optionally, a name-standardization function \code{std_names()} can be
#' applied to the row/column names of both matrices (e.g. to strip or normalize
#' IDs). This function is assumed to be defined elsewhere in the package or
#' analysis code.
#'
#' The output is a list whose first two elements are the aligned matrices,
#' named according to the inferred names of the input arguments (e.g.
#' \code{"WT"} and \code{"OE"} if you called
#' \code{align_corr_pair(mat[["WT"]], mat[["OE"]])}). The third element
#' \code{$vars} contains the variable names used.
#'
#' @param A A square numeric matrix with row and column names.
#' @param B A square numeric matrix with row and column names.
#' @param mode Character, either \code{"intersect"} (default) for the
#'   intersection of variable names, or \code{"union"} for the union.
#' @param fill Value used to fill missing entries when \code{mode = "union"}.
#' @param name_std Logical; if \code{TRUE}, apply \code{std_names()} to
#'   row/column names of both matrices before alignment.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{[[<name_of_A>]]} Aligned version of matrix \code{A}.
#'   \item \code{[[<name_of_B>]]} Aligned version of matrix \code{B}.
#'   \item \code{$vars} Character vector of variable names used.
#' }
#'
#' @examples
#' # Suppose A and B are correlation matrices over overlapping sets of taxa/metabolites:
#' # res <- align_corr_pair(A, B, mode = "intersect")
#' # names(res) might be c("WT", "OE", "vars") if A and B were passed as mat[["WT"]], mat[["OE"]].
#'
#' @export
align_corr_pair <- function(A, B,
                            mode = c("intersect","union"),
                            fill = NA,
                            name_std = TRUE) {
  mode <- match.arg(mode)
  stopifnot(is.matrix(A), is.matrix(B))

  ## infer names for A and B
  exprA <- deparse(substitute(A))
  exprB <- deparse(substitute(B))
  aname <- .get_matrix_name(exprA, env = parent.frame())
  bname <- .get_matrix_name(exprB, env = parent.frame())

  rnA <- rownames(A); rnB <- rownames(B)
  cnA <- colnames(A); cnB <- colnames(B)

  if (name_std) {
    rnA <- std_names(rnA); cnA <- std_names(cnA)
    rnB <- std_names(rnB); cnB <- std_names(cnB)
    rownames(A) <- rnA; colnames(A) <- cnA
    rownames(B) <- rnB; colnames(B) <- cnB
  }

  ## basic checks
  if (!identical(sort(rnA), sort(cnA)))
    warning("Matrix A: row/column name sets differ (may not be strictly symmetric).")
  if (!identical(sort(rnB), sort(cnB)))
    warning("Matrix B: row/column name sets differ (may not be strictly symmetric).")
  if (any(duplicated(rnA)) || any(duplicated(cnA)))
    stop("Matrix A has duplicated row/column names.")
  if (any(duplicated(rnB)) || any(duplicated(cnB)))
    stop("Matrix B has duplicated row/column names.")

  S <- switch(
    mode,
    intersect = intersect(rnA, rnB),
    union     = union(rnA, rnB)
  )
  if (length(S) < 2) stop("Number of shared variables < 2; cannot compare.")

  ## expand/crop a matrix to S, preserving symmetry
  expand_to <- function(M, S) {
    nm <- rownames(M)
    M2 <- matrix(fill, nrow = length(S), ncol = length(S),
                 dimnames = list(S, S))
    keep <- intersect(S, nm)
    if (length(keep))
      M2[keep, keep] <- M[keep, keep, drop = FALSE]
    ## enforce symmetry using upper triangle
    M2[lower.tri(M2)] <- t(M2)[lower.tri(M2)]
    diag(M2) <- 1
    M2
  }

  A2 <- expand_to(A, S)
  B2 <- expand_to(B, S)

  out <- list()
  out[[aname]] <- A2
  out[[bname]] <- B2
  out$vars     <- S
  out
}


#' Compare two correlation matrices and their derived networks
#'
#' This function compares two correlation matrices (\code{A} and \code{B})
#' that share the same set of variables (node IDs) and differ only in the
#' correlation values. It provides:
#' \itemize{
#'   \item Per-edge differences in correlation (\code{r_B - r_A}).
#'   \item Flags for equality within a tolerance, sign flips, increases and decreases.
#'   \item Network-level comparison based on a correlation threshold
#'         (\code{edge_thr}) and edge rule (\code{"abs"}, \code{"pos"},
#'         or \code{"neg"}).
#'   \item Counts of shared and condition-specific edges.
#'   \item Optional comparison of "strong" correlations using \code{thr}.
#' }
#'
#' The input matrices must:
#' \itemize{
#'   \item Be square numeric matrices.
#'   \item Have identical row and column names.
#'   \item Represent the same variables in the same order.
#' }
#'
#' The output includes a summary list with counts, a per-edge data frame
#' (\code{$edges}) containing the correlations in both matrices and various
#' flags, as well as convenience subsets (e.g. \code{$changed},
#' \code{$flipped}, \code{$shared_edges}, \code{$<Aname>_only_edges},
#' \code{$<Bname>_only_edges}). The correlation columns in \code{edges}
#' are named \code{r_<Aname>} and \code{r_<Bname>} based on how the matrices
#' were passed (e.g. \code{mat[["WT"]]}, \code{mat[["OE"]]}).
#'
#' @param A A square numeric correlation matrix with row and column names.
#' @param B A square numeric correlation matrix with row and column names.
#' @param tol Numeric scalar; tolerance within which correlations are treated
#'   as "equal".
#' @param delta Numeric scalar; \code{|Δr|} threshold to classify edges as
#'   increased or decreased (\code{r_B - r_A > delta} or
#'   \code{r_B - r_A < -delta}).
#' @param thr Optional numeric; absolute correlation threshold for classifying
#'   "strong" correlations when constructing the \code{$strong} component of
#'   the output.
#' @param edge_thr Numeric; threshold used to define network edges when
#'   comparing networks (e.g. \code{0.3}, \code{0.5}, \code{0.7}).
#' @param edge_rule Character; one of \code{"abs"}, \code{"pos"} or
#'   \code{"neg"}:
#'   \itemize{
#'     \item \code{"abs"}: treat edges as present when \code{|r| >= edge_thr}.
#'     \item \code{"pos"}: treat edges as present when \code{r >= edge_thr}.
#'     \item \code{"neg"}: treat edges as present when \code{r <= -edge_thr}.
#'   }
#' @param same_sign_only Logical; if \code{TRUE}, shared edges must have the
#'   same sign in \code{A} and \code{B} to be counted as shared.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{$summary}: Named numeric counts (total edges, changed,
#'         sign flips, increased, decreased, NA patterns, and counts of
#'         edges in each network and shared/unique edges under the chosen
#'         \code{edge_thr}).
#'   \item \code{$edges}: Data frame with one row per unique edge (upper
#'         triangle), including correlations in both matrices, difference,
#'         flags, and network membership.
#'   \item \code{$changed}: Subset of \code{edges} where correlations differ
#'         by more than \code{tol} and are not NA in only one matrix.
#'   \item \code{$flipped}: Subset of \code{edges} where correlation sign
#'         differs between \code{A} and \code{B}.
#'   \item \code{$increased}: Subset of \code{edges} with
#'         \code{r_B - r_A > delta}.
#'   \item \code{$decreased}: Subset of \code{edges} with
#'         \code{r_B - r_A < -delta}.
#'   \item \code{$shared_edges}: Subset of \code{edges} present in both
#'         networks under the specified \code{edge_thr} and \code{edge_rule}.
#'   \item \code{$<Aname>_only_edges}: Subset of \code{edges} present only
#'         in network \code{A}.
#'   \item \code{$<Bname>_only_edges}: Subset of \code{edges} present only
#'         in network \code{B}.
#'   \item \code{$strong}: Optional list with strong-correlation comparisons
#'         if \code{thr} is provided.
#'   \item \code{$top_changes}: Edges with the largest absolute changes
#'         in correlation.
#'   \item \code{$params}: List of parameters used in the comparison
#'         (names of matrices, thresholds, etc.).
#' }
#'
#' @examples
#' # Assume cor_WT and cor_OE are aligned correlation matrices
#' # with identical row/column names:
#' # res <- compare_corr_two(cor_WT, cor_OE, edge_thr = 0.6, edge_rule = "abs")
#' # head(res$shared_edges)
#'
#' @export
compare_corr_two <- function(
    A, B,
    tol = 1e-6,
    delta = 0.1,
    thr = NULL,
    edge_thr = 0.3,
    edge_rule = c("abs", "pos", "neg"),
    same_sign_only = FALSE
) {
  edge_rule <- match.arg(edge_rule)
  stopifnot(is.matrix(A), is.matrix(B))
  stopifnot(nrow(A) == ncol(A), nrow(B) == ncol(B))
  stopifnot(identical(rownames(A), colnames(A)))
  stopifnot(identical(rownames(B), colnames(B)))
  stopifnot(identical(rownames(A), rownames(B)))

  ## infer names of A and B
  exprA <- deparse(substitute(A))
  exprB <- deparse(substitute(B))
  aname <- .get_matrix_name(exprA, env = parent.frame())
  bname <- .get_matrix_name(exprB, env = parent.frame())

  ## correlation column names
  rA_col <- paste0("r_", aname)
  rB_col <- paste0("r_", bname)

  rn  <- rownames(A)
  idx <- which(upper.tri(A), arr.ind = TRUE)
  v1 <- rn[idx[, 1]]
  v2 <- rn[idx[, 2]]
  rA <- A[idx]
  rB <- B[idx]

  ## NA handling
  both_na <- is.na(rA) & is.na(rB)
  one_na  <- xor(is.na(rA), is.na(rB))

  ## numeric differences
  equal      <- (!one_na) & (both_na | (abs(rA - rB) <= tol))
  diff       <- rB - rA
  abs_change <- abs(diff)
  sA <- sign(rA)
  sB <- sign(rB)
  sign_flip <- (!one_na & !both_na) & (sA * sB == -1)

  increased <- (!one_na & !both_na) & (diff >  delta)
  decreased <- (!one_na & !both_na) & (diff < -delta)

  ## edge mask according to edge_rule
  edge_mask <- function(r, rule, thr) {
    if (rule == "abs") return(!is.na(r) & abs(r) >= thr)
    if (rule == "pos") return(!is.na(r) & r >=  thr)
    if (rule == "neg") return(!is.na(r) & r <= -thr)
  }
  A_edge <- edge_mask(rA, edge_rule, edge_thr)
  B_edge <- edge_mask(rB, edge_rule, edge_thr)

  ## shared vs unique edges
  shared_raw <- A_edge & B_edge
  if (same_sign_only) {
    shared <- shared_raw & (sA == sB) & !is.na(sA) & !is.na(sB)
  } else {
    shared <- shared_raw
  }
  A_only <- A_edge & !B_edge
  B_only <- B_edge & !A_edge

  ## per-edge dataframe
  edges <- data.frame(
    var1 = v1,
    var2 = v2,
    r_A  = rA,
    r_B  = rB,
    diff = diff,
    abs_change = abs_change,
    equal      = equal,
    increased  = increased,
    decreased  = decreased,
    sign_flip  = sign_flip,
    one_na     = one_na,
    both_na    = both_na,
    A_edge     = A_edge,
    B_edge     = B_edge,
    shared     = shared,
    A_only     = A_only,
    B_only     = B_only,
    stringsAsFactors = FALSE
  )
  ## rename r_A / r_B to r_<name>
  names(edges)[names(edges) == "r_A"] <- rA_col
  names(edges)[names(edges) == "r_B"] <- rB_col

  ## strong correlation comparison (optional)
  strong <- NULL
  if (!is.null(thr)) {
    A_strong <- (!is.na(rA)) & (abs(rA) >= thr)
    B_strong <- (!is.na(rB)) & (abs(rB) >= thr)

    strong <- list()
    strong[[paste0(aname, "_only")]] <- subset(
      edges,  A_strong & !B_strong,
      select = c("var1", "var2", rA_col, rB_col, "diff")
    )
    strong[[paste0(bname, "_only")]] <- subset(
      edges, !A_strong &  B_strong,
      select = c("var1", "var2", rA_col, rB_col, "diff")
    )
    strong[["both"]] <- subset(
      edges,  A_strong &  B_strong,
      select = c("var1", "var2", rA_col, rB_col, "diff")
    )
  }

  ## summary counts
  summary_counts <- list(
    total_edges   = nrow(edges),
    equal         = sum(equal,    na.rm = TRUE),
    changed       = sum(!equal & !one_na, na.rm = TRUE),
    sign_flip     = sum(sign_flip, na.rm = TRUE),
    increased     = sum(increased, na.rm = TRUE),
    decreased     = sum(decreased, na.rm = TRUE),
    one_na        = sum(one_na,    na.rm = TRUE),
    both_na       = sum(both_na,   na.rm = TRUE)
  )

  ## legacy names (compatible with older code)
  summary_counts$A_edge       <- sum(A_edge,   na.rm = TRUE)
  summary_counts$B_edge       <- sum(B_edge,   na.rm = TRUE)
  summary_counts$shared_edges <- sum(shared,   na.rm = TRUE)
  summary_counts$A_only_edges <- sum(A_only,   na.rm = TRUE)
  summary_counts$B_only_edges <- sum(B_only,   na.rm = TRUE)

  ## extended names with condition labels
  summary_counts[[paste0(aname, "_edges")]]       <- summary_counts$A_edge
  summary_counts[[paste0(bname, "_edges")]]       <- summary_counts$B_edge
  summary_counts[[paste0("shared_", aname, "_", bname, "_edges")]] <- summary_counts$shared_edges
  summary_counts[[paste0(aname, "_only_edges")]]  <- summary_counts$A_only_edges
  summary_counts[[paste0(bname, "_only_edges")]]  <- summary_counts$B_only_edges

  ## largest changes
  top_changes <- edges[order(-edges$abs_change), ][1:min(50, nrow(edges)), ]

  ## output
  out <- list(
    summary   = summary_counts,
    edges     = edges,
    changed   = subset(edges, !equal & !one_na),
    flipped   = subset(edges, sign_flip),
    increased = subset(edges, increased),
    decreased = subset(edges, decreased),
    strong    = strong,
    top_changes = top_changes,
    shared_edges = subset(edges, shared,
                          select = c("var1", "var2", rA_col, rB_col, "diff")),
    params = list(
      A_name = aname,
      B_name = bname,
      edge_thr = edge_thr,
      edge_rule = edge_rule,
      same_sign_only = same_sign_only
    )
  )

  ## attach condition-specific "only" edge tables
  out[[paste0(aname, "_only_edges")]] <- subset(
    edges, A_only,
    select = c("var1", "var2", rA_col, rB_col, "diff")
  )
  out[[paste0(bname, "_only_edges")]] <- subset(
    edges, B_only,
    select = c("var1", "var2", rA_col, rB_col, "diff")
  )

  out
}
