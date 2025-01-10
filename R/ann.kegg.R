#' @title Retrieve KEGG metabolite information based on input metabolite IDs
#' @description
#' This function matches the input metabolite ID to allMetabolites in the KEGG database by formatting it,
#' and finally outputs its storage ID in the KEGG database.
#' @param id The metabolite ID for which KEGG information is to be retrieved.
#' @return A data frame with KEGG IDs and matched metabolite names for input metabolite IDs.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' tax= ps.ms %>% vegan_tax() %>%as.data.frame()
#' id = tax$Metabolite
#' tax2 = ann.kegg(id)
#' head(tax2)
#' @export
ann.kegg = function(id){
  mk = db.ms.kegg

  mk$allMetabolites = str_to_lower(mk$allMetabolites)

  id2 = gsub("[ ][0-9]","",id)
  id2 = str_to_lower(id2)
  A = c()
  B = c()
  for (i in 1:length(id2)) {
    a = id2[i]
    tem = mk %>% dplyr::filter(str_detect(allMetabolites , stringr::fixed(a)))
    tem
    if (dim(tem)[1] ==0) {
      A[i] = ""
      A[i] = ""
    } else{
      # A[i] =  str_c(tem$keggID,collapse = ";")
      # B[i] = str_c(tem$allMetabolites,collapse = "#")
      A[i] =  tem$keggID[1]
      B[i] =tem$allMetabolites[1]

    }
  }


  tax = data.frame(ID = id,keggID = A,matchname = B)
  return(tax)
}

#' @title Another matching mode for retrieving KEGG metabolite information based on input metabolite IDs
#' @description
#' This function matches the input metabolite ID to common names in the KEGG database by formatting it,
#' and finally outputs its storage ID in the KEGG database.
#' @param id The metabolite ID for which KEGG information is to be retrieved.
#' @return A data frame with KEGG IDs and matched metabolite names for input metabolite IDs.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' tax= ps.ms %>% vegan_tax() %>%as.data.frame()
#' id = tax$Metabolite
#' tax2 = ann.kegg2(id)
#' head(tax2)
#' @export
ann.kegg2 = function(id,repath){
  mk =db.ms.kegg
  head(mk)
  mk$allMetabolites = str_to_lower(mk$allMetabolites)

  id2 = gsub("[ ][0-9]","",id)
  id2 = str_to_lower(id2)
  A = c()
  B = c()
  #

  compoundM = mk[,1:2]
  colnames(compoundM) = c("KEGG compound", "common names")
  rownames(compoundM) = NULL


  for (i in 1:length(id)) {
    match= id2[i]
    if (length(match) >= 1) {
      target_matrix = compoundM
      target_column = compoundM[, 2]
      matchM = match_KEGG(match, target_column, target_matrix)

      if (!is.null(matchM[[1]])) {
        A[i] = matchM[[1]][1,1]
        B[i] = matchM[[1]][1,2]
      } else{
        A[i] = ""
        B[i] = ""
      }

    } else{

    }
  }

  tax = data.frame(ID = id,keggID = A,matchname = B)
  return(tax)
}

#' @title  Match KEGG metabolites based on input metabolite(s) ID, KEGG database
#' @description
#' The match_KEGG function performs custom matching operations for retrieving
#' KEGG metabolite information based on the input match criteria, target column, and target matrix.
#' @param match The metabolite(s) to be matched.
#' @param target_column The target column for matching in the metabolite database.
#' @param target_matrix The matrix containing  the metabolite database.
#' @return A list of matched KEGG metabolites based on the input metabolite(s) ID, KEGG database.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' library(dplyr)
#' tax= ps.ms %>% vegan_tax() %>%as.data.frame()
#' id = tax$Metabolite
#' tax2 = ann.kegg2(id)
#' head(tax2)
#' @export
match_KEGG = function(match, target_column, target_matrix) {
  match = tolower(match)
  match = unique(match)
  matchM = vector (mode = "list", length = length(match))
  names (matchM) = match

  for (key_word in match) {
    key_wordS = unlist(strsplit(key_word, " "))

    if (length(key_wordS) > 1) {
      main_index = grep(key_wordS[1], tolower(target_column))
      if (length(main_index) == 0) {
        to_print = paste("Could not match", key_word, sep = " ")
        message(to_print)
      } else {
        secondary_index = c()
        index_to_intersect = main_index
        for (p in 2:length(key_wordS)) {
          indexS = grep(key_wordS[p], tolower(target_column))
          if (length(intersect(main_index, indexS)) == 0) {
            indexS = NULL
          }
          if (length(indexS) > 0) {
            secondary_index = c(secondary_index, indexS)
            secondary_index = unique(secondary_index)
            index_to_intersect = intersect(index_to_intersect, indexS)
          }
        }
        if (length(secondary_index) > 0) {
          if (length(index_to_intersect) > 0) {
            index = index_to_intersect
          } else {
            intersection = intersect(main_index, secondary_index)
            if (length(intersection) > 0) {
              index = intersection
            } else {
              index = main_index
            }
          }
        } else {
          index = main_index
        }
        key_wordLines = target_matrix[index, ]
        matchM[[key_word]] = key_wordLines
      }
    } else {
      key_word = key_wordS
      index = grep(key_word, tolower(target_column))
      if (length(index) == 0) {
        to_print = paste("Could not match", key_word,
                         sep = " ")
        warning(to_print)
      } else {
        key_wordLines = target_matrix[index, ]
        matchM[[key_word]] = key_wordLines
      }
    }
  }
  return(matchM)
}


