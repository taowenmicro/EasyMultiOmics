#' Reorder Significance Letters Based on Group Means
#'
#' This function reorders group significance letters (e.g., from multiple comparisons)
#' according to descending order of group means, while preserving the original
#' significance grouping structure (i.e., groups that share a letter remain grouped).
#'
#' @param data A data frame containing group labels, means, and significance letters.
#' @param group_col The column name (as string) for the group identifiers (e.g., strain names).
#' @param letter_col The column name (as string) containing the significance letters (e.g., "a", "b", "bc", etc.).
#' @param mean_col The column name (as string) containing the group means.
#'
#' @return A data frame with the same structure as the input, but with reordered significance letters
#'         (i.e., high-mean groups receive earlier letters like "a").
#'
#' @details
#' This function is useful after statistical tests such as Tukey's HSD or Dunn's test,
#' where shared letters indicate non-significant differences between groups.
#' The function does not alter which groups are considered significantly different —
#' it only reassigns letter labels to reflect the relative magnitude of group means.
#'
#' @importFrom dplyr mutate select
#' @importFrom stringr str_split
#' @importFrom purrr map_dbl map_chr
#'
#' @examples
#' library(dplyr)
#' library(stringr)
#' library(purrr)
#'
#' data <- data.frame(
#'   groups = c("c", "bc", "b", "a"),
#'   group = c("Rhi", "CK", "CF", "WT"),
#'   mean = c(243000, 58500, 46500, 43000),
#'   stringsAsFactors = FALSE
#' )
#'
#' ord_sig(data)
#'
#' @export

ord_sig <- function(data, group_col = "group", letter_col = "groups", mean_col = "mean") {

  df <- data

  # 每个组的字母
  group_letters <- df %>%
    mutate(letters = str_split(.[[letter_col]], "")) %>%
    select(group = !!sym(group_col), mean = !!sym(mean_col), letters)

  # 创建字母到组的映射（反向）
  letter_to_groups <- list()
  for (i in seq_len(nrow(group_letters))) {
    gr <- group_letters$group[i]
    ls <- group_letters$letters[[i]]
    for (ltr in ls) {
      letter_to_groups[[ltr]] <- unique(c(letter_to_groups[[ltr]], gr))
    }
  }

  # 计算每个字母组的平均值（即该字母所覆盖的所有 group 的平均值）
  letter_means <- map_dbl(letter_to_groups, function(grps) {
    mean(df[[mean_col]][df[[group_col]] %in% grps])
  })

  # 按照 mean 从大到小排序字母
  sorted_letters <- names(sort(letter_means, decreasing = TRUE))

  # 分配新的字母（a, b, c...）
  new_letter_map <- setNames(letters[seq_along(sorted_letters)], sorted_letters)

  # 替换原始 groups 中的字母
  df[[letter_col]] <- map_chr(df[[letter_col]], function(x) {
    x_split <- str_split(x, "")[[1]]
    x_new <- sort(unique(new_letter_map[x_split]))
    paste(x_new, collapse = "")
  })

  return(df)
}

