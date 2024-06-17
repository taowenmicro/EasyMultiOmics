

# res = Ven.Upset(ps =  ps,
#     group = "Group",
#     N = 0.5,
#     size = 3)
#
# res[[1]]
# res[[2]]


Ven.Upset.metf = function(
    otu = NULL,
    tax = NULL,
    map = NULL,
    tree = NULL,
    ps = NULL,
    group = "Group",
    N = 0.5,
    size = 3
){
  library(ggVennDiagram)
  ps =  ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  aa =  ggClusterNet::vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count

  sub_design <- as.data.frame(phyloseq::sample_data(ps))


  #
  # pick_val_num <- rep/2
  count[count > 0] <- 1
  count2 = as.data.frame(count )
  aa = sub_design[,"Group"]
  colnames(aa) = "Vengroup"

  #-div group
  iris.split <- split(count2,as.factor(aa$Vengroup))
  #数据分组计算平均值
  iris.apply <- lapply(iris.split,function(x)colSums(x[]))
  # 组合结果
  iris.combine <- do.call(rbind,iris.apply)
  ven2 = t(iris.combine)
  for (i in 1:length(unique(phyloseq::sample_data(ps)[,"Group"]))) {
    aa <- as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,1]
    bb =  as.data.frame(table(phyloseq::sample_data(ps)[,"Group"]))[i,2]
    ven2[,aa] = ven2[,aa]/bb
  }
  # N = 0.5
  ven2[ven2 < N]  = 0
  ven2[ven2 >=N]  = 1
  ven2 = as.data.frame(ven2)

  ven3 = list()
  for (i in 1:ncol(ven2)) {
    ven3[[i]] = row.names(ven2)[ven2[,i] == 1]
  }
  names(ven3) = colnames(ven2)
  # ven3 = as.list(ven2)
  p = ggVennDiagram(
    ven3, label_alpha = 0,
    category.names = names(ven3)
  ) +
    ggplot2::scale_fill_gradientn(colours = c("#E41A1C","#377EB8" , "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999"))



  # install.packages("ggupset")
  # library(tidyverse)
  # library(ggupset)
  A = list()
  for (i in 1:nrow(ven2)) {
    tem = ven2[i,] %>% t() %>%.[,1]
    A[[i]] = tem[tem == 1] %>% names()
  }

  dat = tibble(ID = row.names(ven2),A)
  head(dat)

  main_plot <-ggplot(data = dat,aes(x=A)) +
    geom_bar() +
    ggupset::scale_x_upset(n_intersections = 20) +
    theme_classic()

  side_plot <- dat %>%
    select(A) %>%
    unnest(cols = A) %>%
    dplyr::count(A) %>%
    mutate(A = fct_reorder(as.factor(A), n)) %>%
    ggplot(aes(y = n, x = A)) +
    geom_col() +
    coord_flip() +
    scale_y_reverse() +
    xlab("") + ylab("") +
    theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(),
          panel.background=element_blank(),
          panel.grid=element_blank()
    )



  p2 = cowplot::plot_grid(
    cowplot::plot_grid(NULL, side_plot + theme(plot.margin = unit(c(1, -5, -5, 1), "pt")), ncol = 1, rel_heights = c(size, 1)),
    main_plot, nrow = 1, rel_widths = c(1, 3)
  )
  return(list(p,p2,ven2,main_plot,side_plot))

}

