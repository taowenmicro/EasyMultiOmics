#' @title  Retrieve metabolite information from the HMDB database
#' @description
#' The ann.HMDB function retrieves relevant information about a metabolite based
#'  on its ID from the HMDB (Human Metabolome Database) stored in the db.metabolites database.
#' @param id The metabolite ID for which information is to be retrieved.
#' @return A data frame containing the metabolite information including HMDB ID, KEGG ID,
#' Super_class, Class, Sub_class and related details from the database.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' library(dplyr)
#' tax= ps.ms %>% vegan_tax() %>%as.data.frame()
#' id = tax$Metabolite
#' tax1 = ann.HMDB (id = id)
#' head(tax1)
#' @export

ann.HMDB = function(id = id){
  db = db.metabolites
  id2 = gsub("[ ][0-9]","",id)
  id2 = str_to_lower(id2)
  db$Name2 = str_to_lower(db$Name)
  tem2 = data.frame(id = id,id2 = id2,id.org = id)
  tem3 = tem2 %>% dplyr::left_join(db,by = c(id2 = "Name2"))
  return(tem3)
}
