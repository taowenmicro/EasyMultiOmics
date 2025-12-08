#' Summarise taxonomy of nodes involved in a set of edges
#'
#' This helper summarises the taxonomic composition (phylum/genus) of nodes
#' involved in a given edge set (e.g. \code{res_net$shared_edges} or a
#' condition-specific edge set from \code{\link{compare_corr_two}}).
#'
#' @param edges A data frame of edges, typically one of the edge tables from
#'   \code{\link{compare_corr_two}} (e.g. \code{res_net$shared_edges} or
#'   \code{res_net[["WT_only_edges"]]}).
#' @param tax_table A taxonomy table. Can be a data frame derived from
#'   \code{phyloseq::tax_table(ps) \%>\% as.data.frame()} or any table with
#'   per-node taxonomy.
#' @param var1_col Name of the column in \code{edges} giving the first node
#'   ID of each edge (default \code{"var1"}).
#' @param var2_col Name of the column in \code{edges} giving the second node
#'   ID of each edge (default \code{"var2"}).
#' @param id_col Name of the column in \code{tax_table} that stores node IDs
#'   matching \code{var1_col} / \code{var2_col}. If \code{NULL}, row names of
#'   \code{tax_table} are used as node IDs.
#' @param phylum_col Name of the phylum column in \code{tax_table}.
#' @param genus_col Name of the genus column in \code{tax_table}.
#' @param set_name Character label for this edge set (e.g. \code{"shared"},
#'   \code{"WT"}, \code{"OE"}). Used in the output summaries.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{set_name}: The label used for this set.
#'   \item \code{nodes}: Character vector of unique node IDs appearing in
#'         the edge set.
#'   \item \code{node_tax}: Subset of \code{tax_table} corresponding to
#'         these nodes.
#'   \item \code{phylum_summary}: Data frame with counts and relative
#'         frequencies per phylum (\code{set}, \code{Phylum}, \code{n_nodes},
#'         \code{rel}).
#'   \item \code{genus_summary}: Data frame with counts and relative
#'         frequencies per genus (\code{set}, \code{Genus}, \code{n_nodes},
#'         \code{rel}).
#' }
#'
#' @examples
#' # edges <- res_net$shared_edges
#' # tax   <- phyloseq::tax_table(ps.micro) %>% as.data.frame()
#' # out   <- summarise_edge_nodes_tax(edges, tax, id_col = "geneid",
#' #                                   phylum_col = "Phylum", genus_col = "Genus",
#' #                                   set_name = "shared")
#'
#' @export
summarise_edge_nodes_tax <- function(
    edges,
    tax_table,
    var1_col   = "var1",
    var2_col   = "var2",
    id_col     = "geneid",
    phylum_col = "Phylum",
    genus_col  = "Genus",
    set_name   = "edge_set"
) {
  if (is.null(edges) || nrow(edges) == 0) {
    warning("`edges` is empty. Returning empty summaries.")
    empty_phylum <- data.frame(
      set    = character(0),
      Phylum = character(0),
      n_nodes = integer(0),
      rel    = numeric(0),
      stringsAsFactors = FALSE
    )
    empty_genus <- data.frame(
      set    = character(0),
      Genus  = character(0),
      n_nodes = integer(0),
      rel    = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(
      set_name       = set_name,
      nodes          = character(0),
      node_tax       = tax_table[0, , drop = FALSE],
      phylum_summary = empty_phylum,
      genus_summary  = empty_genus
    ))
  }

  # All unique node IDs in this edge set
  nodes <- unique(c(edges[[var1_col]], edges[[var2_col]]))

  tax_df <- as.data.frame(tax_table)

  # Prepare ID column
  if (!is.null(id_col)) {
    if (!id_col %in% colnames(tax_df)) {
      stop("`id_col` not found in `tax_table`.")
    }
    tax_df$NodeID <- tax_df[[id_col]]
  } else {
    if (is.null(rownames(tax_df))) {
      stop("tax_table has no rownames and `id_col` is NULL.")
    }
    tax_df$NodeID <- rownames(tax_df)
  }

  # Only nodes that appear in taxonomy table are considered (microbes)
  node_tax <- tax_df %>%
    dplyr::filter(NodeID %in% nodes)

  # Phylum-level summary
  phylum_summary <- node_tax %>%
    dplyr::filter(!is.na(.data[[phylum_col]]), .data[[phylum_col]] != "") %>%
    dplyr::count(.data[[phylum_col]], name = "n_nodes") %>%
    dplyr::mutate(
      set = set_name,
      rel = n_nodes / sum(n_nodes)
    ) %>%
    dplyr::rename(Phylum = .data[[phylum_col]]) %>%
    dplyr::select(set, Phylum, n_nodes, rel) %>%
    dplyr::arrange(dplyr::desc(n_nodes))

  # Genus-level summary
  genus_summary <- node_tax %>%
    dplyr::filter(!is.na(.data[[genus_col]]), .data[[genus_col]] != "") %>%
    dplyr::count(.data[[genus_col]], name = "n_nodes") %>%
    dplyr::mutate(
      set = set_name,
      rel = n_nodes / sum(n_nodes)
    ) %>%
    dplyr::rename(Genus = .data[[genus_col]]) %>%
    dplyr::select(set, Genus, n_nodes, rel) %>%
    dplyr::arrange(dplyr::desc(n_nodes))

  list(
    set_name       = set_name,
    nodes          = nodes,
    node_tax       = node_tax,
    phylum_summary = phylum_summary,
    genus_summary  = genus_summary
  )
}


#' Compare phylum/genus composition of shared / condition-specific edge sets
#'
#' This function compares the taxonomic composition (phylum/genus) of nodes
#' involved in three edge sets derived from \code{\link{compare_corr_two}}:
#' \itemize{
#'   \item \code{shared_edges}: edges present in both networks.
#'   \item edges unique to the first matrix (A-only).
#'   \item edges unique to the second matrix (B-only).
#' }
#'
#' It supports both the new naming scheme of \code{compare_corr_two()}:
#' \itemize{
#'   \item \code{res_net$shared_edges}
#'   \item \code{res_net[[paste0(A_name, "_only_edges")]]}
#'   \item \code{res_net[[paste0(B_name, "_only_edges")]]}
#' }
#' and the legacy generic names:
#' \itemize{
#'   \item \code{res_net$A_only_edges}
#'   \item \code{res_net$B_only_edges}.
#' }
#'
#' The output includes separate summaries for each edge set and two wide tables
#' comparing phylum- and genus-level node counts across the sets.
#'
#' @param res_net A list returned by \code{\link{compare_corr_two}}.
#' @param tax_table Taxonomy table.
#' @param var1_col Name of the column in the edge tables for the first node
#'   (default \code{"var1"}).
#' @param var2_col Name of the column in the edge tables for the second node
#'   (default \code{"var2"}).
#' @param id_col Name of the column in \code{tax_table} matching node IDs
#'   (e.g. ASV or gene IDs).
#' @param phylum_col Name of the phylum column in \code{tax_table}.
#' @param genus_col Name of the genus column in \code{tax_table}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{shared}, \code{A_only}, \code{B_only}: Each is the output
#'         of \code{\link{summarise_edge_nodes_tax}} for the corresponding
#'         edge set.
#'   \item \code{phylum_compare}: Wide table of phylum-level counts and
#'         relative frequencies across the sets (\code{n_nodes_*}, \code{rel_*}).
#'   \item \code{genus_compare}: Wide table of genus-level counts and
#'         relative frequencies across the sets.
#' }
#'
#' @examples
#' # res_net <- compare_corr_two(mat[["WT"]], mat[["OE"]], edge_thr = 0.6)
#' # tax     <- phyloseq::tax_table(ps.micro) %>% as.data.frame()
#' # cmp_tax <- compare_shared_A_B_tax(res_net, tax,
#' #                                   id_col = "geneid",
#' #                                   phylum_col = "Phylum",
#' #                                   genus_col  = "Genus")
#'
#' @export
compare_shared_A_B_tax <- function(
    res_net,
    tax_table,
    var1_col   = "var1",
    var2_col   = "var2",
    id_col     = "geneid",
    phylum_col = "Phylum",
    genus_col  = "Genus"
) {
  if (!("shared_edges" %in% names(res_net))) {
    stop("`res_net` must contain $shared_edges.")
  }

  # Infer A/B names from compare_corr_two params if available
  A_name <- if (!is.null(res_net$params$A_name)) res_net$params$A_name else "A"
  B_name <- if (!is.null(res_net$params$B_name)) res_net$params$B_name else "B"

  # Locate A-only and B-only edge sets:
  #   1) Try legacy generic names: A_only_edges / B_only_edges
  #   2) Otherwise, use dynamic names: <A_name>_only_edges / <B_name>_only_edges
  if ("A_only_edges" %in% names(res_net)) {
    A_edges <- res_net$A_only_edges
  } else {
    dynA <- paste0(A_name, "_only_edges")
    if (!dynA %in% names(res_net)) {
      stop("Cannot find A-only edges: neither `A_only_edges` nor `", dynA, "` exist in res_net.")
    }
    A_edges <- res_net[[dynA]]
  }

  if ("B_only_edges" %in% names(res_net)) {
    B_edges <- res_net$B_only_edges
  } else {
    dynB <- paste0(B_name, "_only_edges")
    if (!dynB %in% names(res_net)) {
      stop("Cannot find B-only edges: neither `B_only_edges` nor `", dynB, "` exist in res_net.")
    }
    B_edges <- res_net[[dynB]]
  }

  # Summaries for each edge set (label sets as "shared", A_name, B_name)
  res_shared <- summarise_edge_nodes_tax(
    edges      = res_net$shared_edges,
    tax_table  = tax_table,
    var1_col   = var1_col,
    var2_col   = var2_col,
    id_col     = id_col,
    phylum_col = phylum_col,
    genus_col  = genus_col,
    set_name   = "shared"
  )

  res_Aonly <- summarise_edge_nodes_tax(
    edges      = A_edges,
    tax_table  = tax_table,
    var1_col   = var1_col,
    var2_col   = var2_col,
    id_col     = id_col,
    phylum_col = phylum_col,
    genus_col  = genus_col,
    set_name   = A_name
  )

  res_Bonly <- summarise_edge_nodes_tax(
    edges      = B_edges,
    tax_table  = tax_table,
    var1_col   = var1_col,
    var2_col   = var2_col,
    id_col     = id_col,
    phylum_col = phylum_col,
    genus_col  = genus_col,
    set_name   = B_name
  )

  # Phylum-level comparison: long -> wide
  phylum_long <- dplyr::bind_rows(
    res_shared$phylum_summary,
    res_Aonly$phylum_summary,
    res_Bonly$phylum_summary
  )

  phylum_compare <- phylum_long %>%
    tidyr::pivot_wider(
      names_from  = set,
      values_from = c(n_nodes, rel),
      values_fill = 0
    ) %>%
    dplyr::mutate(
      total_n = rowSums(dplyr::across(dplyr::starts_with("n_nodes_")), na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(total_n)) %>%
    dplyr::select(-total_n)

  # Genus-level comparison: long -> wide
  genus_long <- dplyr::bind_rows(
    res_shared$genus_summary,
    res_Aonly$genus_summary,
    res_Bonly$genus_summary
  )

  genus_compare <- genus_long %>%
    tidyr::pivot_wider(
      names_from  = set,
      values_from = c(n_nodes, rel),
      values_fill = 0
    ) %>%
    dplyr::mutate(
      total_n = rowSums(dplyr::across(dplyr::starts_with("n_nodes_")), na.rm = TRUE)
    ) %>%
    dplyr::arrange(dplyr::desc(total_n)) %>%
    dplyr::select(-total_n)

  list(
    shared         = res_shared,
    A_only         = res_Aonly,
    B_only         = res_Bonly,
    phylum_compare = phylum_compare,
    genus_compare  = genus_compare
  )
}


#' Summarise within- and between-phylum structure of microbe–microbe edges
#'
#' This helper analyses the within- and between-phylum structure of edges,
#' typically microbe–microbe edges from a network comparison (e.g. one of the
#' sets from \code{\link{compare_corr_two}}).
#'
#' Only edges where both endpoints can be matched to a phylum in
#' \code{tax_table} are considered as "microbe–microbe" edges.
#'
#' @param edges Data frame of edges (e.g. \code{res_net$shared_edges} or a
#'   condition-specific edge set).
#' @param tax_table Taxonomy table.
#' @param var1_col Name of the column in \code{edges} storing the first node
#'   ID (default \code{"var1"}).
#' @param var2_col Name of the column in \code{edges} storing the second node
#'   ID (default \code{"var2"}).
#' @param id_col Name of the ID column in \code{tax_table} matching node IDs.
#'   If \code{NULL}, row names of \code{tax_table} are used as IDs.
#' @param phylum_col Name of the phylum column in \code{tax_table}.
#' @param set_name Label for this edge set (e.g. \code{"shared"},
#'   \code{"WT"}, \code{"OE"}).
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{set_name}: Edge-set label.
#'   \item \code{edge_n_all}: Total number of edges in \code{edges}.
#'   \item \code{edge_n_mm}: Number of microbe–microbe edges (both endpoints
#'         have a phylum annotation).
#'   \item \code{within_between_summary}: Table with counts and relative
#'         proportions of within-phylum vs between-phylum edges.
#'   \item \code{within_phylum}: For within-phylum edges, counts per phylum.
#'   \item \code{between_pairs}: For between-phylum edges, counts per
#'         phylum–phylum pair (order-invariant).
#' }
#'
#' @export
summarise_edge_phylum_structure <- function(
    edges,
    tax_table,
    var1_col   = "var1",
    var2_col   = "var2",
    id_col     = "geneid",
    phylum_col = "Phylum",
    set_name   = "edge_set"
) {
  edge_n_all <- nrow(edges)

  if (edge_n_all == 0) {
    warning("`edges` is empty. Returning empty summaries.")
    empty_within_between <- data.frame(
      set     = set_name,
      type    = character(0),
      n_edges = integer(0),
      rel     = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(
      set_name               = set_name,
      edge_n_all             = 0L,
      edge_n_mm              = 0L,
      within_between_summary = empty_within_between,
      within_phylum          = data.frame(),
      between_pairs          = data.frame()
    ))
  }

  tax_df <- as.data.frame(tax_table)

  # Prepare NodeID
  if (!is.null(id_col)) {
    if (!id_col %in% colnames(tax_df)) {
      stop("`id_col` not found in `tax_table`.")
    }
    tax_df$NodeID <- tax_df[[id_col]]
  } else {
    if (is.null(rownames(tax_df))) {
      stop("tax_table has no rownames and `id_col` is NULL.")
    }
    tax_df$NodeID <- rownames(tax_df)
  }

  # Retain only NodeID and phylum
  tax_df <- tax_df[, c("NodeID", phylum_col), drop = FALSE]
  colnames(tax_df)[colnames(tax_df) == phylum_col] <- "Phylum"

  # Attach phylum to each endpoint
  edge_tax <- edges %>%
    dplyr::transmute(
      node1 = .data[[var1_col]],
      node2 = .data[[var2_col]]
    ) %>%
    dplyr::left_join(tax_df, by = c("node1" = "NodeID")) %>%
    dplyr::rename(Phylum1 = Phylum) %>%
    dplyr::left_join(tax_df, by = c("node2" = "NodeID")) %>%
    dplyr::rename(Phylum2 = Phylum)

  # Only microbe–microbe edges where both endpoints have phylum
  edge_mm <- edge_tax %>%
    dplyr::filter(!is.na(Phylum1), !is.na(Phylum2))

  edge_n_mm <- nrow(edge_mm)

  if (edge_n_mm == 0) {
    warning("No microbe–microbe edges with phylum annotation in this set.")
    empty_within_between <- data.frame(
      set     = set_name,
      type    = character(0),
      n_edges = integer(0),
      rel     = numeric(0),
      stringsAsFactors = FALSE
    )
    return(list(
      set_name               = set_name,
      edge_n_all             = edge_n_all,
      edge_n_mm              = 0L,
      within_between_summary = empty_within_between,
      within_phylum          = data.frame(),
      between_pairs          = data.frame()
    ))
  }

  # Mark within- vs between-phylum edges
  edge_mm <- edge_mm %>%
    dplyr::mutate(
      type = ifelse(Phylum1 == Phylum2, "within", "between")
    )

  # 1) Overall within- vs between-phylum proportions
  within_between_summary <- edge_mm %>%
    dplyr::count(type, name = "n_edges") %>%
    dplyr::mutate(
      set = set_name,
      rel = n_edges / sum(n_edges)
    ) %>%
    dplyr::select(set, type, n_edges, rel)

  # 2) Within-phylum edges by phylum
  within_phylum <- edge_mm %>%
    dplyr::filter(type == "within") %>%
    dplyr::count(Phylum1, name = "n_edges") %>%
    dplyr::mutate(
      set        = set_name,
      rel_within = n_edges / sum(n_edges)
    ) %>%
    dplyr::rename(Phylum = Phylum1) %>%
    dplyr::select(set, Phylum, n_edges, rel_within) %>%
    dplyr::arrange(dplyr::desc(n_edges))

  # 3) Between-phylum edges by unordered phylum–phylum pair
  between_pairs <- edge_mm %>%
    dplyr::filter(type == "between") %>%
    dplyr::mutate(
      Phy_a = pmin(Phylum1, Phylum2),
      Phy_b = pmax(Phylum1, Phylum2)
    ) %>%
    dplyr::count(Phy_a, Phy_b, name = "n_edges") %>%
    dplyr::mutate(
      set         = set_name,
      rel_between = n_edges / sum(n_edges)
    ) %>%
    dplyr::select(set, Phylum1 = Phy_a, Phylum2 = Phy_b, n_edges, rel_between) %>%
    dplyr::arrange(dplyr::desc(n_edges))

  list(
    set_name               = set_name,
    edge_n_all             = edge_n_all,
    edge_n_mm              = edge_n_mm,
    within_between_summary = within_between_summary,
    within_phylum          = within_phylum,
    between_pairs          = between_pairs
  )
}


#' Compare within- and between-phylum structure for shared / A-only / B-only edges
#'
#' This function compares the within- and between-phylum edge structure for
#' three edge sets derived from \code{\link{compare_corr_two}}:
#' \itemize{
#'   \item shared edges (\code{res_net$shared_edges}),
#'   \item edges unique to matrix A,
#'   \item edges unique to matrix B.
#' }
#'
#' It supports both the dynamic naming scheme (e.g.
#' \code{res_net[["WT_only_edges"]]}, \code{res_net[["OE_only_edges"]]}) and
#' the legacy \code{A_only_edges} / \code{B_only_edges} names if present.
#'
#' @param res_net The result list from \code{\link{compare_corr_two}}.
#' @param tax_table Taxonomy table.
#' @param var1_col Name of the first-node column in edge tables (default
#'   \code{"var1"}).
#' @param var2_col Name of the second-node column in edge tables (default
#'   \code{"var2"}).
#' @param id_col Name of ID column in \code{tax_table} matching node IDs.
#' @param phylum_col Name of phylum column in \code{tax_table}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{shared}, \code{A_only}, \code{B_only}: Each is the result of
#'         \code{\link{summarise_edge_phylum_structure}} for that edge set.
#'   \item \code{within_between_compare}: Wide table comparing within- vs
#'         between-phylum edge counts/proportions across the three sets.
#' }
#'
#' @export
compare_phylum_structure_three <- function(
    res_net,
    tax_table,
    var1_col   = "var1",
    var2_col   = "var2",
    id_col     = "geneid",
    phylum_col = "Phylum"
) {
  if (!("shared_edges" %in% names(res_net))) {
    stop("`res_net` must contain $shared_edges.")
  }

  # Infer A/B names
  A_name <- if (!is.null(res_net$params$A_name)) res_net$params$A_name else "A"
  B_name <- if (!is.null(res_net$params$B_name)) res_net$params$B_name else "B"

  # Locate A-only and B-only edge sets with backward compatibility
  if ("A_only_edges" %in% names(res_net)) {
    A_edges <- res_net$A_only_edges
  } else {
    dynA <- paste0(A_name, "_only_edges")
    if (!dynA %in% names(res_net)) {
      stop("Cannot find A-only edges: neither `A_only_edges` nor `", dynA, "` exist in res_net.")
    }
    A_edges <- res_net[[dynA]]
  }

  if ("B_only_edges" %in% names(res_net)) {
    B_edges <- res_net$B_only_edges
  } else {
    dynB <- paste0(B_name, "_only_edges")
    if (!dynB %in% names(res_net)) {
      stop("Cannot find B-only edges: neither `B_only_edges` nor `", dynB, "` exist in res_net.")
    }
    B_edges <- res_net[[dynB]]
  }

  shared_res <- summarise_edge_phylum_structure(
    edges      = res_net$shared_edges,
    tax_table  = tax_table,
    var1_col   = var1_col,
    var2_col   = var2_col,
    id_col     = id_col,
    phylum_col = phylum_col,
    set_name   = "shared"
  )

  A_res <- summarise_edge_phylum_structure(
    edges      = A_edges,
    tax_table  = tax_table,
    var1_col   = var1_col,
    var2_col   = var2_col,
    id_col     = id_col,
    phylum_col = phylum_col,
    set_name   = A_name
  )

  B_res <- summarise_edge_phylum_structure(
    edges      = B_edges,
    tax_table  = tax_table,
    var1_col   = var1_col,
    var2_col   = var2_col,
    id_col     = id_col,
    phylum_col = phylum_col,
    set_name   = B_name
  )

  within_between_compare <- dplyr::bind_rows(
    shared_res$within_between_summary,
    A_res$within_between_summary,
    B_res$within_between_summary
  ) %>%
    tidyr::pivot_wider(
      names_from  = type,
      values_from = c(n_edges, rel),
      values_fill = 0
    )

  list(
    shared                 = shared_res,
    A_only                 = A_res,
    B_only                 = B_res,
    within_between_compare = within_between_compare
  )
}


#' Get shared nodes across multiple correlation networks
#' 专门用在 align_corr_pair 的结果上
#'
#' @param align_res 输出来自 align_corr_pair()，或者是一个仅包含相关性矩阵的 list
#' @param tax_table 可选，微生物的 taxonomy 表
#' @param id_col tax_table 中与网络节点名对应的列名；若为 NULL，则用 rownames(tax_table)
#' @param phylum_col 门所在列名
#' @param genus_col 属所在列名
#'
#' @return list:
#'   - shared_edges: data.frame，共享边
#'   - shared_nodes: 字符向量，共享网络中的所有节点
#'   - node_tax: 这些节点中能在 tax_table 里匹配到的微生物注释
#'   - phylum_summary: 按门统计的节点数
#'   - genus_summary: 按属统计的节点数
get_shared_nodes_from_align <- function(
    align_res,
    tax_table  = NULL,
    id_col     = NULL,
    phylum_col = "Phylum",
    genus_col  = "Genus"
) {
  # 1. 从 align_corr_pair 输出中抽出真正的矩阵 ------------------------
  if (is.list(align_res) && "mats" %in% names(align_res)) {
    ab <- align_res$mats              # 典型情况：align_res$mats 里是 G/B 等矩阵
  } else {
    ab <- align_res                   # 否则假定 align_res 本身就是 list(G=, B=, ...)
  }

  if (!is.list(ab) || length(ab) < 2) {
    stop("`align_res` 里面找不到至少两个相关矩阵，请检查 align_corr_pair 的输出结构。")
  }

  # 2. 只保留矩阵元素 ---------------------------------------------------
  is_mat <- vapply(ab, is.matrix, logical(1))
  ab <- ab[is_mat]

  if (length(ab) < 2) {
    stop("`align_res` 中矩阵元素少于两个，无法计算共享网络。")
  }

  net_names <- names(ab)
  if (is.null(net_names)) {
    net_names <- paste0("net", seq_along(ab))
  }

  # 3. 检查维度和节点名是否一致 ----------------------------------------
  ref_mat <- ab[[1]]
  ref_nodes_row <- rownames(ref_mat)
  ref_nodes_col <- colnames(ref_mat)

  if (is.null(ref_nodes_row) || is.null(ref_nodes_col)) {
    stop("相关性矩阵必须有 rownames 和 colnames。")
  }
  if (!identical(ref_nodes_row, ref_nodes_col)) {
    stop("第一个矩阵的行名和列名不一致，请确认是方阵，且 dimnames 对齐。")
  }

  for (i in seq_along(ab)[-1]) {
    mat <- ab[[i]]
    if (!all(dim(mat) == dim(ref_mat))) {
      stop("All matrices in `ab` must have the same dimensions.\n",
           "⚠️ 你现在传进去的似乎是原始 `mats`，请改为传 align_corr_pair 的结果：\n",
           "   res <- align_corr_pair(mats$G, mats$B, ...)\n",
           "   get_shared_nodes_from_align(res, tax_table = ...)")
    }
    if (!identical(rownames(mat), ref_nodes_row) ||
        !identical(colnames(mat), ref_nodes_col)) {
      stop("所有矩阵的 rownames / colnames 必须完全一致（同一批节点、同一顺序）。")
    }
  }

  n     <- nrow(ref_mat)
  nodes <- ref_nodes_row

  # 4. 计算“共享边”：所有网络中都非 NA 且 != 0 --------------------------
  edge_mask <- matrix(TRUE, nrow = n, ncol = n,
                      dimnames = list(nodes, nodes))

  for (i in seq_along(ab)) {
    mat <- ab[[i]]
    m_edge <- !is.na(mat) & (mat != 0)
    edge_mask <- edge_mask & m_edge
  }

  # 去掉对角线，只保留上三角，避免重复
  edge_mask[lower.tri(edge_mask, diag = TRUE)] <- FALSE

  idx <- which(edge_mask, arr.ind = TRUE)

  if (nrow(idx) == 0) {
    shared_edges <- data.frame(
      from = character(0),
      to   = character(0),
      stringsAsFactors = FALSE
    )
    return(list(
      shared_edges   = shared_edges,
      shared_nodes   = character(0),
      node_tax       = NULL,
      phylum_summary = NULL,
      genus_summary  = NULL
    ))
  }

  shared_edges <- data.frame(
    from = nodes[idx[, 1]],
    to   = nodes[idx[, 2]],
    stringsAsFactors = FALSE
  )

  shared_nodes <- sort(unique(c(shared_edges$from, shared_edges$to)))

  # 5. 不给 taxonomy 就到此为止 -----------------------------------------
  if (is.null(tax_table)) {
    return(list(
      shared_edges   = shared_edges,
      shared_nodes   = shared_nodes,
      node_tax       = NULL,
      phylum_summary = NULL,
      genus_summary  = NULL
    ))
  }

  # 6. 对接 taxonomy，做门/属统计 ---------------------------------------
  tax_df <- as.data.frame(tax_table)

  if (!is.null(id_col)) {
    if (!id_col %in% colnames(tax_df)) {
      stop("`id_col` not found in `tax_table`.")
    }
    tax_df$NodeID <- tax_df[[id_col]]
  } else {
    if (is.null(rownames(tax_df))) {
      stop("tax_table 没有 rownames，请指定 `id_col`。")
    }
    tax_df$NodeID <- rownames(tax_df)
  }

  node_tax <- tax_df %>%
    dplyr::filter(NodeID %in% shared_nodes)

  phylum_summary <- node_tax %>%
    dplyr::filter(!is.na(.data[[phylum_col]]), .data[[phylum_col]] != "") %>%
    dplyr::count(.data[[phylum_col]], name = "n_nodes") %>%
    dplyr::arrange(dplyr::desc(n_nodes))

  genus_summary <- node_tax %>%
    dplyr::filter(!is.na(.data[[genus_col]]), .data[[genus_col]] != "") %>%
    dplyr::count(.data[[genus_col]], name = "n_nodes") %>%
    dplyr::arrange(dplyr::desc(n_nodes))

  list(
    shared_edges   = shared_edges,
    shared_nodes   = shared_nodes,
    node_tax       = node_tax,
    phylum_summary = phylum_summary,
    genus_summary  = genus_summary
  )
}

