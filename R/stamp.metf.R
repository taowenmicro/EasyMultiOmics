
stamp.metf <- function(ps_sub = ps,Top = 20,method = "rela",
                       test.method = "t.test"){

  #--门水平合并
  data = ps_sub %>%
    ggClusterNet::scale_micro(method = method) %>%
    ggClusterNet::filter_OTU_ps(Top = 200) %>%
    ggClusterNet::vegan_otu() %>%
    as.data.frame()
  tem = colnames(data)
  data$ID = row.names(data)
  map= phyloseq::sample_data(ps)
  map$ID = row.names(map)
  data <- data %>%
    dplyr::inner_join(as.tibble(map),by = "ID")
  data$Group = as.factor(data$Group)

  if (test.method == "t.test") {
    diff <- data[,tem] %>%
      # dplyr::select_if(is.numeric) %>%
      purrr::map_df(~ broom::tidy(t.test(. ~ Group,data = data)), .id = 'var')




  } else if(test.method == "wilcox.test"){
    diff <- data[,tem] %>%
      # dplyr::select_if(is.numeric) %>%
      purrr::map_df(~ broom::tidy(wilcox.test(. ~ Group,data = data)), .id = 'var')

  }

  diff$p.value[is.nan(diff$p.value)] = 1
  diff$p.value <- p.adjust(diff$p.value,"bonferroni")
  tem = diff$p.value [diff$p.value < 0.05] %>% length()
  if (tem > 30) {
    diff <- diff %>%
      dplyr::filter(p.value < 0.05) %>%
      head(30)
  } else {
    diff <- diff %>%
      # filter(p.value < 0.05) %>%
      head(30)
  }

  # diff <- diff %>% filter(p.value < 0.05)

  # diff1$p.value <- p.adjust(diff1$p.value,"bonferroni")
  # diff1 <- diff1 %>% filter(p.value < 0.05)

  abun.bar <- data[,c(diff$var,"Group")] %>%
    tidyr::gather(variable,value,-Group) %>%
    dplyr::group_by(variable,Group) %>%
    dplyr::summarise(Mean = mean(value))


  diff.mean <- diff[,c("var","estimate","conf.low","conf.high","p.value")]
  diff.mean$Group <- c(ifelse(diff.mean$estimate >0,levels(data$Group)[1],
                              levels(data$Group)[2]))
  diff.mean <- diff.mean[order(diff.mean$estimate,decreasing = TRUE),]


  cbbPalette <- c("#E69F00", "#56B4E9")
  abun.bar$variable <- factor(abun.bar$variable,levels = rev(diff.mean$var))


  p1 <- ggplot(abun.bar,aes(variable,Mean,fill = Group)) +
    scale_x_discrete(limits = levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Mean proportion (%)") +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          legend.title=element_blank(),
          # legend.text=element_text(size=12,face = "bold",colour = "black",
          #                          margin = margin(r = 20)),
          legend.position = c(-1,-0.1),
          legend.direction = "horizontal",
          legend.key.width = unit(0.8,"cm"),
          legend.key.height = unit(0.5,"cm"))

  p1

  for (i in 1:(nrow(diff.mean) - 1))
    p1 <- p1 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

  p1
  p1 <- p1 +
    geom_bar(stat = "identity",position = "dodge",width = 0.7,colour = "black") +
    scale_fill_manual(values=cbbPalette) + theme(legend.position = "bottom")
  p1

  diff.mean$var <- factor(diff.mean$var,levels = levels(abun.bar$variable))
  diff.mean$p.value <- signif(diff.mean$p.value,3)
  diff.mean$p.value <- as.character(diff.mean$p.value)





  p2 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
    theme(panel.background = element_rect(fill = 'transparent'),
          panel.grid = element_blank(),
          axis.ticks.length = unit(0.4,"lines"),
          axis.ticks = element_line(color='black'),
          axis.line = element_line(colour = "black"),
          axis.title.x=element_text(colour='black', size=12,face = "bold"),
          axis.text=element_text(colour='black',size=10,face = "bold"),
          axis.text.y = element_blank(),
          legend.position = "none",
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          plot.title = element_text(size = 15,face = "bold",colour = "black",hjust = 0.5)) +
    scale_x_discrete(limits = levels(diff.mean$var)) +
    coord_flip() +
    xlab("") +
    ylab("Difference in mean proportions (%)") +
    labs(title="95% confidence intervals")

  for (i in 1:(nrow(diff.mean) - 1))
    p2 <- p2 + annotate('rect', xmin = i+0.5, xmax = i+1.5, ymin = -Inf, ymax = Inf,
                        fill = ifelse(i %% 2 == 0, 'white', 'gray95'))

  p2 <- p2 +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high),
                  position = position_dodge(0.8), width = 0.5, size = 0.5) +
    geom_point(shape = 21,size = 3) +
    scale_fill_manual(values=cbbPalette) +
    geom_hline(aes(yintercept = 0), linetype = 'dashed', color = 'black')

  p3 <- ggplot(diff.mean,aes(var,estimate,fill = Group)) +
    geom_text(aes(y = 0,x = var),label = diff.mean$p.value,
              hjust = 0,fontface = "bold",inherit.aes = FALSE,size = 3) +
    geom_text(aes(x = nrow(diff.mean)/2 +0.5,y = 0.85),label = "P-value (corrected)",
              srt = 90,fontface = "bold",size = 5) +
    coord_flip() +
    ylim(c(0,1)) +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank())


  library(patchwork)
  p <- p1 + p2 + p3 + plot_layout(widths = c(4,6,2))
  p
  return(list(p,abun.bar,diff.mean))
}
