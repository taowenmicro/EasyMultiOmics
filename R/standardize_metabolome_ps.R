#' 基于语义理解的智能代谢组学数据标准化
#' 不依赖穷举列表，而是通过语义分析和规则推理识别列名
#'

# ========== 1. 文本预处理 ==========

#' 标准化字符串（保留语义信息）
normalize_text <- function(x) {
  x %>%
    tolower() %>%
    # 统一分隔符为空格
    gsub("[._-]", " ", .) %>%
    # 去除括号内容但保留括号前的词
    gsub("\\s*\\([^)]*\\)", "", .) %>%
    # 去除多余空格
    gsub("\\s+", " ", .) %>%
    trimws()
}

#' 提取关键词（词干）
extract_keywords <- function(text) {
  words <- unlist(strsplit(normalize_text(text), " "))
  # 过滤停用词和超短词
  words <- words[nchar(words) >= 2]
  words <- words[!words %in% c("the", "and", "or", "of", "in", "for", "to", "a", "an")]
  return(words)
}


#' Get semantic rules for column name matching
#'
#' This function returns a list of semantic rules used for matching
#' column names in metabolomics data. The rules include mandatory,
#' optional, and excluded keywords for different types of columns such
#' as "Metabolite", "Formula", and "KEGG_ID".
#'
#' @return A list of semantic rules, each containing:
#'   - `must_contain`: A list of word groups that must appear in the column name.
#'   - `should_contain`: A list of words that can improve the match.
#'   - `must_not_contain`: A list of words that should not appear.
#'   - `weight`: A numeric weight for the importance of the rule.
#' @examples
#' get_semantic_rules()
#' @export
get_semantic_rules <- function() {
  list(
    # Define rules for different types of columns (e.g., Metabolite, Formula, KEGG_ID, etc.)
    "metab_id" = list(
      must_contain = list(c("id")),  # At least one word group
      should_contain = c("metabolite", "metab", "compound", "feature"),
      must_not_contain = c("name", "description"),
      weight = 10
    ),
    "Metabolite" = list(
      must_contain = list(c("name"), c("metabolite")),  # At least one group with "name" or "metabolite"
      should_contain = c("compound", "annotation"),
      must_not_contain = c("id", "number", "score", "error"),
      weight = 10
    ),
    "Formula" = list(
      must_contain = list(c("formula")),
      should_contain = c("chemical", "molecular"),
      must_not_contain = c(),
      weight = 15
    ),
    "m.z" = list(
      must_contain = list(c("m.z", "m/z", "mass")),
      should_contain = c("precursor", "observed", "measured"),
      must_not_contain = c("error", "difference", "delta", "accuracy"),
      weight = 10
    ),
    "Retention.time" = list(
      must_contain = list(c("Retention.time", "retention", "time")),
      should_contain = c("chromatographic"),
      must_not_contain = c("error", "difference", "delta", "start", "end", "window"),
      weight = 10
    ),
    "Mode" = list(
      must_contain = list(c("mode", "polarity")),
      should_contain = c("ion", "ionization", "detection"),
      must_not_contain = c("score", "error"),
      weight = 8
    ),
    "KEGG_ID" = list(
      must_contain = list(c("kegg")),
      should_contain = c("id", "compound"),
      must_not_contain = c("pathway", "map", "module"),
      weight = 12
    ),
    "HMDB_ID" = list(
      must_contain = list(c("hmdb")),
      should_contain = c("id", "number"),
      must_not_contain = c(),
      weight = 12
    ),
    "CAS" = list(
      must_contain = list(c("cas")),
      should_contain = c("number", "registry", "rn"),
      must_not_contain = c("score", "match", "similarity", "fragmentation", "theoretical"),
      weight = 12
    ),
    "PUBCHEM_ID" = list(
      must_contain = list(c("pubchem", "cid")),
      should_contain = c("id", "compound"),
      must_not_contain = c(),
      weight = 12
    ),
    "CHEBI_ID" = list(
      must_contain = list(c("chebi")),
      should_contain = c("id", "number"),
      must_not_contain = c(),
      weight = 12
    ),
    "INCHIKEY" = list(
      must_contain = list(c("inchi")),
      should_contain = c("key", "standard"),
      must_not_contain = c(),
      weight = 12
    ),
    "SMILES" = list(
      must_contain = list(c("smiles")),
      should_contain = c("canonical", "isomeric"),
      must_not_contain = c(),
      weight = 12
    ),
    "Super_Class" = list(
      must_contain = list(c("super", "main", "primary", "kingdom")),
      should_contain = c("class"),
      must_not_contain = c("sub", "secondary", "minor"),
      weight = 8
    ),
    "Class" = list(
      must_contain = list(c("class")),
      should_contain = c("taxonomy", "classification"),
      must_not_contain = c("super", "sub", "main", "primary", "secondary"),
      weight = 7
    ),
    "Sub_Class" = list(
      must_contain = list(c("sub", "minor", "secondary")),
      should_contain = c("class"),
      must_not_contain = c("super", "main", "primary"),
      weight = 8
    ),
    "Pathway" = list(
      must_contain = list(c("pathway")),
      should_contain = c("metabolic", "biological", "kegg"),
      must_not_contain = c(),
      weight = 10
    ),
    "Frag_Score" = list(
      must_contain = list(c("fragmentation", "fragment", "frag", "ms2", "msms", "spectrum")),
      should_contain = c("score", "quality"),
      must_not_contain = c("theoretical", "expected", "reference", "library"),
      weight = 10
    ),
    "Theoretical_Frag_Score" = list(
      must_contain = list(c("theoretical", "expected", "reference")),
      should_contain = c("fragmentation", "fragment", "frag", "score"),
      must_not_contain = c(),
      weight = 12
    ),
    "Score" = list(
      must_contain = list(c("score", "confidence", "match", "similarity")),
      should_contain = c("identification", "annotation"),
      must_not_contain = c("fragmentation", "fragment", "frag", "ms2", "msms",
                           "spectrum", "mass", "error", "ppm", "theoretical"),
      weight = 6
    ),
    "Mass_Error_PPM" = list(
      must_contain = list(c("ppm", "mass"), c("error", "accuracy")),
      should_contain = c("delta", "difference"),
      must_not_contain = c("score", "fragmentation"),
      weight = 10
    ),
    "Level" = list(
      must_contain = list(c("level")),
      should_contain = c("confidence", "identification", "annotation", "id"),
      must_not_contain = c(),
      weight = 8
    ),
    "Adducts" = list(
      must_contain = list(c("adduct")),
      should_contain = c("ion", "type", "form"),
      must_not_contain = c(),
      weight = 10
    ),
    "Intensity" = list(
      must_contain = list(c("intensity", "abundance", "signal")),
      should_contain = c("peak", "ion"),
      must_not_contain = c("relative", "normalized", "ratio"),
      weight = 8
    ),
    "Concentration" = list(
      must_contain = list(c("concentration", "conc", "amount", "quantity")),
      should_contain = c("measured", "calculated"),
      must_not_contain = c(),
      weight = 8
    ),
    "Description" = list(
      must_contain = list(c("description", "comment", "note", "remark")),
      should_contain = c(),
      must_not_contain = c(),
      weight = 5
    ),
    "Synonym" = list(
      must_contain = list(c("synonym", "alias", "alternative")),
      should_contain = c("name"),
      must_not_contain = c(),
      weight = 8
    )
  )
}




#' 计算列名与规则的匹配分数
calculate_match_score <- function(col_name, rule) {
  keywords <- extract_keywords(col_name)

  if (length(keywords) == 0) return(0)

  score <- 0

  # 1. 检查必须包含的词（任一组满足即可）
  must_match <- FALSE
  for (must_group in rule$must_contain) {
    if (any(must_group %in% keywords)) {
      must_match <- TRUE
      score <- score + rule$weight * 2  # 必须词权重最高
      break
    }
  }

  if (!must_match) return(0)  # 不满足必须条件，直接返回0

  # 2. 检查排除词（如果包含任何排除词，直接返回0）
  if (any(rule$must_not_contain %in% keywords)) {
    return(0)
  }

  # 3. 检查建议包含的词（加分项）
  should_matches <- sum(rule$should_contain %in% keywords)
  score <- score + should_matches * rule$weight * 0.5

  # 4. 惩罚过长的列名（倾向于简洁的匹配）
  length_penalty <- max(0, (length(keywords) - 3) * 0.5)
  score <- score - length_penalty

  return(max(score, 0))
}

#' 为所有列名找到最佳匹配
smart_column_mapping <- function(col_names) {
  rules <- get_semantic_rules()

  # 计算每个列名与每个规则的匹配分数
  score_matrix <- matrix(0, nrow = length(col_names), ncol = length(rules),
                         dimnames = list(col_names, names(rules)))

  for (i in seq_along(col_names)) {
    for (j in seq_along(rules)) {
      score_matrix[i, j] <- calculate_match_score(col_names[i], rules[[j]])
    }
  }

  # 为每个规则找到最佳匹配的列（贪心算法）
  mapping <- list()
  used_cols <- character(0)

  # 按规则权重排序，优先匹配重要的规则
  rule_order <- names(rules)[order(sapply(rules, function(r) r$weight), decreasing = TRUE)]

  for (rule_name in rule_order) {
    available_cols <- setdiff(col_names, used_cols)
    if (length(available_cols) == 0) break

    scores <- score_matrix[available_cols, rule_name, drop = FALSE]
    max_score <- max(scores)

    if (max_score > 0) {
      best_col <- rownames(scores)[which.max(scores)]
      mapping[[best_col]] <- rule_name
      used_cols <- c(used_cols, best_col)
    }
  }

  return(mapping)
}





standardize_metabolome_tax <- function(tax, verbose = FALSE) {

  # 0. 保存原始行名
  original_rownames <- rownames(tax)
  has_meaningful_rownames <- !all(original_rownames == as.character(1:nrow(tax)))

  # 1. 创建新数据框
  tax_std <- tax

  # 2. 如果有有意义的行名但没有name列，先添加
  if (has_meaningful_rownames) {
    name_cols <- grep("name|metabolite|compound", colnames(tax_std),
                      ignore.case = TRUE, value = TRUE)
    name_cols <- name_cols[!grepl("id|file|path", name_cols, ignore.case = TRUE)]

    if (length(name_cols) == 0) {
      tax_std <- data.frame(
        Metabolite_Name_from_rownames = original_rownames,
        tax_std,
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
    }
  }

  # 3. 智能匹配列名
  mapping <- smart_column_mapping(colnames(tax_std))

  if (verbose) {
    cat("\n========== 智能列名识别结果 ==========\n\n")
    if (length(mapping) > 0) {
      for (i in seq_along(mapping)) {
        cat(sprintf("%-40s  -->  %s\n", names(mapping)[i], mapping[[i]]))
      }
    } else {
      cat("未识别到任何标准列名\n")
    }
  }

  # 4. 应用重命名
  for (old_name in names(mapping)) {
    new_name <- mapping[[old_name]]
    colnames(tax_std)[colnames(tax_std) == old_name] <- new_name
  }

  # 5. 处理从行名来的Metabolite_Name
  if ("Metabolite_Name_from_rownames" %in% colnames(tax_std)) {
    if (!"Metabolite" %in% colnames(tax_std) || all(is.na(tax_std$Metabolite_Name))) {
      tax_std$Metabolite_Name <- tax_std$Metabolite_Name_from_rownames
    }
    tax_std$Metabolite_Name_from_rownames <- NULL
  }

  # 6. 保留原始行名
  rownames(tax_std) <- original_rownames

  # 7. 确保关键列存在
  essential_cols <- c("metab_id", "Metabolite", "Formula", "KEGG_ID")
  for (col in essential_cols) {
    if (!col %in% colnames(tax_std)) {
      tax_std[[col]] <- NA
    }
  }

  # 8. 智能列重排序
  priority_cols <- c(
    "metab_id", "Metabolite", "Formula", "m.z", "Retention.time", "Mode",
    "KEGG_ID", "HMDB_ID", "PUBCHEM_ID", "CHEBI_ID", "CAS",
    "Super_Class", "Class", "Sub_Class", "Pathway",
    "Level", "Score", "Frag_Score", "Theoretical_Frag_Score",
    "Adducts", "Mass_Error_PPM", "INCHIKEY", "SMILES"
  )

  existing_priority <- intersect(priority_cols, colnames(tax_std))
  other_cols <- setdiff(colnames(tax_std), existing_priority)

  tax_std <- tax_std[, c(existing_priority, other_cols)]

  # 9. 添加元数据
  attr(tax_std, "standardization_info") <- list(
    original_ncol = ncol(tax),
    standardized_ncol = ncol(tax_std),
    matched_columns = mapping,
    timestamp = Sys.time(),
    method = "semantic_rules"
  )

  return(tax_std)
}


#' 直接标准化phyloseq对象
standardize_metabolome_ps <- function(ps, verbose = TRUE) {

  tax_df <- as.data.frame(phyloseq::tax_table(ps))

  if (verbose) {
    cat("原始taxonomy维度:", nrow(tax_df), "行 x", ncol(tax_df), "列\n")
    cat("原始列名:", paste(colnames(tax_df), collapse = ", "), "\n\n")
  }

  tax_std <- standardize_metabolome_tax(tax_df, verbose = verbose)

  if (verbose) {
    cat("\n标准化后维度:", nrow(tax_std), "行 x", ncol(tax_std), "列\n")
    cat("标准化列名:", paste(colnames(tax_std), collapse = ", "), "\n")
  }

  phyloseq::tax_table(ps) <- phyloseq::tax_table(as.matrix(tax_std))

  return(ps)
}

