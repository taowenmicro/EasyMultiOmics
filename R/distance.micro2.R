distance.micro2 = function(ps = NULL, group = "Group") {
  map = as.data.frame(sample_data(ps))

  # 分组向量
  gro = data.frame(group = unique(map[[group]]))
  conbgroup = combn(gro$group, 2)

  # Bray-Curtis 距离矩阵
  bray_curtis = vegan::vegdist(vegan_otu(ps), method = "bray") %>% as.matrix()

  table = NULL
  for (i in 1:ncol(conbgroup)) {
    a = conbgroup[, i]

    chose1 = row.names(map[map[[group]] == a[1], ])
    chose2 = row.names(map[map[[group]] == a[2], ])

    if (length(chose1) == 0 | length(chose2) == 0) next  # 跳过空组

    dat = data.frame(
      group = paste(a[1], a[2], sep = "_VS_"),
      Distance = bray_curtis[chose1, chose2] %>% as.dist() %>% as.vector()
    )

    table = rbind(table, dat)
  }

  if (is.null(table)) stop("No valid group comparisons found.")

  # 整理结果
  table$id = 1:nrow(table)
  data <- table %>% dplyr::select(id, everything())

  result = MuiKwWlx2(data = data, num = c(3))
  result1 = FacetMuiPlotresultBox(data = data, num = c(3), result = result, sig_show ="abc", ncol = 1)
  p1_1 = result1[[1]] +theme_nature()
  res = FacetMuiPlotresultBar(data = data, num = c(3), result = result, sig_show ="abc", ncol = 1)
  p1_2 = res[[1]]+theme_nature()
  res = FacetMuiPlotReBoxBar(data = data, num = c(3), result = result, sig_show ="abc", ncol = 1)
  p1_3 = res[[1]]+theme_nature()

  return(list(p1_1, p1_2, p1_3, distance.data = data))
}
