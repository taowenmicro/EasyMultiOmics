

#  (6 >gnum & gnum > 2)
distance.micro = function(
    ps = ps,
    group = "Group"
  ){

  map = as.data.frame(sample_data(ps))
  gro = map[,group] %>% unique()
  colnames(gro) = "group"
  conbgroup = combn(gro$group,2)
  # 计算包括终点均值的所有样品bray距离
  bray_curtis = vegan::vegdist(vegan_otu(ps), method = "bray")
  bray_curtis = as.matrix(bray_curtis)

  for (i in 1:dim(conbgroup)[2]) {
    a = conbgroup[,i]
    map = as.data.frame(sample_data(ps))
    head(map)

    chose1 = map[as.matrix(map[,group]) %>% as.vector() == a[1],] %>% row.names()
    chose2 = map[as.matrix(map[,group]) %>% as.vector() == a[2],] %>% row.names()

    dat = data.frame(group = paste(a[1],a[2],sep = "_VS_"), Distance =bray_curtis[chose1,chose2] %>% as.dist() %>% as.vector() )
    head(dat)

    if (i == 1) {
      table = dat
    }

    if (i != 1) {
      table = rbind(table,dat)
    }
  }

  head(table)
  table$id = 1:dim(table)[1]
  data <- table %>% dplyr::select(id,everything())


  result = EasyStat::MuiKwWlx(data = data,num = c(3))
  result1 = EasyStat::FacetMuiPlotresultBox(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1 )
  p1_1 = result1[[1]]
  p1_1
  res = EasyStat::FacetMuiPlotresultBar(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1)
  p1_2 = res[[1]]
  p1_2
  res = EasyStat::FacetMuiPlotReBoxBar(data = data,num = c(3),result = result,sig_show ="abc",ncol = 1)
  p1_3 = res[[1]]
  p1_3

  return(list(p1_1,p1_2,p1_3,distance.data = data))
}






