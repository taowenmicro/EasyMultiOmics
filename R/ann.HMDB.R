


# library(EasyMultiOmics.db)
ann.HMDB = function(id = id){
  db = db.metabolites
  id2 = gsub("[ ][0-9]","",id)
  id2 = str_to_lower(id2)
  db$Name2 = str_to_lower(db$Name)
  # tem = db %>% filter(Name2 %in% id2)
  # head(tem)
  tem2 = data.frame(id = id,id2 = id2,id.org = id)
  tem3 = tem2 %>% left_join(db,by = c( "id2" = "Name2" ))
  return(tem3)
}
