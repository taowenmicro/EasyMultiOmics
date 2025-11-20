#' @title Guess microbe type (fungi / bacteria / protist) from a phyloseq object
#'
#' @description
#' 根据 phyloseq 对象中的 taxonomy 表（Kingdom / Domain / Phylum 等），
#' 粗略判断这是一个真菌扩增子数据集、细菌/古菌数据，还是以原生动物为主的数据。
#'
#' @param ps A \code{phyloseq} object.
#' @param kingdom_cols Possible column names for kingdom/domain,
#'   e.g. \code{c("Kingdom", "kingdom", "Domain", "domain")}.
#' @param min_major_frac 最主要类群所占比例阈值（默认 0.5），
#'   低于该值则返回 "mixed" 或 "unknown".
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{type}: "fungi", "bacteria", "protist", "mixed", or "unknown".
#'   \item \code{kingdom_counts}: table of kingdom/domain counts (if available).
#'   \item \code{phylum_counts}: table of phylum counts (if available).
#'   \item \code{reason}: a character string explaining the decision.
#' }
#'
#' @examples
#' \dontrun{
#' res <- guess_ps_microbe_type(ps.16s)
#' res$type
#' res$reason
#' }
#'
#' @export
guess_ps_microbe_type <- function(
    ps,
    kingdom_cols   = c("Kingdom", "kingdom", "Domain", "domain"),
    min_major_frac = 0.5
) {
  ## ------- 这里改掉 is.phyloseq -------
  if (!inherits(ps, "phyloseq")) {
    stop("ps must be a phyloseq object.")
  }
  
  tax <- phyloseq::tax_table(ps) %>%
    as.data.frame(stringsAsFactors = FALSE)
  
  # ---- 1. 找 Kingdom / Domain 列 ----
  k_col <- intersect(kingdom_cols, colnames(tax))
  kingdom_counts <- NULL
  phylum_counts  <- NULL
  
  if (length(k_col) > 0) {
    k_col <- k_col[1]
    k_vec <- tax[[k_col]]
    k_vec <- as.character(k_vec)
    k_vec[is.na(k_vec)] <- "Unknown"
    kingdom_counts <- sort(table(k_vec), decreasing = TRUE)
    
    total    <- sum(kingdom_counts)
    top_k    <- names(kingdom_counts)[1]
    top_frac <- kingdom_counts[1] / total
    top_k_upper <- toupper(top_k)
    
    # 先按 Kingdom 判断
    if (grepl("FUNGI", top_k_upper)) {
      if (top_frac >= min_major_frac) {
        return(list(
          type           = "fungi",
          kingdom_counts = kingdom_counts,
          phylum_counts  = phylum_counts,
          reason = paste0(
            "Major kingdom is '", top_k,
            "' with proportion ", round(top_frac, 3),
            " ≥ ", min_major_frac, ", classified as fungal amplicon."
          )
        ))
      } else {
        type <- "mixed"
        reason <- paste0(
          "Kingdom '", top_k, "' is dominant but proportion (",
          round(top_frac, 3), ") < ", min_major_frac,
          ", classified as mixed."
        )
      }
    } else if (grepl("BACTERIA|BACTERIOD|EUBACTERIA", top_k_upper)) {
      if (top_frac >= min_major_frac) {
        return(list(
          type           = "bacteria",
          kingdom_counts = kingdom_counts,
          phylum_counts  = phylum_counts,
          reason = paste0(
            "Major kingdom is '", top_k,
            "' with proportion ", round(top_frac, 3),
            " ≥ ", min_major_frac, ", classified as bacterial amplicon."
          )
        ))
      } else {
        type <- "mixed"
        reason <- paste0(
          "Kingdom '", top_k, "' is dominant but proportion (",
          round(top_frac, 3), ") < ", min_major_frac,
          ", classified as mixed."
        )
      }
    } else if (grepl("ARCHAEA", top_k_upper)) {
      if (top_frac >= min_major_frac) {
        return(list(
          type           = "bacteria",  # 古菌也归在这类里
          kingdom_counts = kingdom_counts,
          phylum_counts  = phylum_counts,
          reason = paste0(
            "Major kingdom is '", top_k,
            "' (Archaea) with proportion ", round(top_frac, 3),
            " ≥ ", min_major_frac, ", classified as archaeal/bacterial amplicon."
          )
        ))
      } else {
        type <- "mixed"
        reason <- paste0(
          "Kingdom '", top_k, "' (Archaea) proportion (",
          round(top_frac, 3), ") < ", min_major_frac,
          ", classified as mixed."
        )
      }
    } else if (grepl("EUKARYA|EUKARYOTA", top_k_upper)) {
      # 真核再看门水平
      type   <- NA_character_
      reason <- NA_character_
    } else {
      type   <- NA_character_
      reason <- NA_character_
    }
  } else {
    type   <- NA_character_
    reason <- "No Kingdom/Domain column found, try using Phylum-level patterns."
  }
  
  # ---- 2. 如果 Kingdom 判不清，再看 Phylum ----
  phy_cols <- intersect(c("Phylum", "phylum"), colnames(tax))
  if (length(phy_cols) > 0) {
    phy_col <- phy_cols[1]
    p_vec   <- tax[[phy_col]] %>% as.character()
    p_vec[is.na(p_vec)] <- "Unknown"
    phylum_counts <- sort(table(p_vec), decreasing = TRUE)
    
    total_p    <- sum(phylum_counts)
    top_p      <- names(phylum_counts)[1]
    top_p_frac <- phylum_counts[1] / total_p
    top_p_upper <- toupper(top_p)
    
    fungi_phyla <- c(
      "ASCOMYCOTA", "BASIDIOMYCOTA", "MUCOROMYCOTA",
      "MORTIERELLAMYCOTA", "CHYTRIDIOMYCOTA", "GLOMEROMYCOTA"
    )
    
    protist_phyla <- c(
      "CILIOPHORA", "APICOMPLEXA", "CERCOZOA", "CHLOROPHYTA",
      "DINOPHYTA", "DINOPHYCEAE", "EUGLENOZOA", "HAPTOPHYTA",
      "OCHROPHYTA", "CRYPTOPHYTA", "RADIOLARIA", "AMOEBOZOA"
    )
    
    if (top_p_upper %in% fungi_phyla && top_p_frac >= min_major_frac) {
      type2 <- "fungi"
      reason2 <- paste0(
        "Major phylum is '", top_p,
        "' with proportion ", round(top_p_frac, 3),
        " ≥ ", min_major_frac, ", classified as fungal amplicon."
      )
    } else if (top_p_upper %in% protist_phyla && top_p_frac >= min_major_frac) {
      type2 <- "protist"
      reason2 <- paste0(
        "Major phylum is '", top_p,
        "' with proportion ", round(top_p_frac, 3),
        " ≥ ", min_major_frac, ", classified as protist-dominated amplicon."
      )
    } else {
      if (is.na(type)) {
        type2 <- "unknown"
        reason2 <- paste0(
          "No clear kingdom-based type; top phylum '", top_p,
          "' proportion ", round(top_p_frac, 3),
          " < ", min_major_frac, " or not in predefined fungi/protist sets."
        )
      } else {
        type2   <- type
        reason2 <- reason
      }
    }
  } else {
    if (is.na(type)) {
      type2   <- "unknown"
      reason2 <- "No Kingdom/Domain or Phylum information that matches known patterns."
    } else {
      type2   <- type
      reason2 <- reason
    }
  }
  
  return(list(
    type           = type2,
    kingdom_counts = kingdom_counts,
    phylum_counts  = phylum_counts,
    reason         = reason2
  ))
}


#' @title Subset phyloseq object by dominant kingdom
#'
#' @description
#' 使用 guess_ps_microbe_type 的结果，自动选择数量最多的 Kingdom，
#' 然后用 subset_taxa.wt() 按该 Kingdom 过滤 OTU/ASV。
#' - 如果是细菌数据，会选出 Bacteria 对应的所有条目；
#' - 真菌数据类似；
#' - 如果是混合数据，则选计数最多的那个 Kingdom；
#' - 如果没有 Kingdom 信息，则直接返回原始 ps。
#'
#' @param ps A phyloseq object.
#' @param kingdom_col tax_table 中 Kingdom 对应的列名，默认 "Kingdom"。
#' @param min_major_frac 传给 guess_ps_microbe_type，用于类型判断（这里只用于理由说明）。
#'
#' @return 一个新的 phyloseq 对象（按主 Kingdom 过滤后的子集）。
#'
#' @export
subset_main_kingdom_ps <- function(
    ps,
    kingdom_col   = "Kingdom",
    min_major_frac = 0.5
) {
  # 调用前面写的猜类型函数
  info <- guess_ps_microbe_type(
    ps,
    kingdom_cols   = c(kingdom_col, "Kingdom", "kingdom", "Domain", "domain"),
    min_major_frac = min_major_frac
  )
  
  kc <- info$kingdom_counts
  
  # 如果 kingdom 信息都没有，就不筛选，直接返回
  if (is.null(kc) || length(kc) == 0) {
    warning("No kingdom/domain information found; returning original phyloseq object.")
    return(ps)
  }
  
  # 取数量最多的 Kingdom 名称（原始字符串，不做大小写改动）
  main_kingdom <- names(kc)[1]
  
  message(
    "guess_ps_microbe_type: type = ", info$type,
    " ; using Kingdom = '", main_kingdom,
    "' (n = ", kc[1], ") as dominant group."
  )
  
  # 用你自己的 subset_taxa.wt 做筛选
  ps_sub <- ps %>% subset_taxa.wt(kingdom_col, main_kingdom)
  
  return(ps_sub)
}
