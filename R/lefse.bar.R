#' @title Visualization of characteristic microorganisms at all taxonomic levels
#' @description
#' This function generates a bar plot for Lefse analysis results, displaying the LDA scores
#' of all classification levels. It arranges the bars based on the LDA scores and colors them
#' according to the Group.
#' @param taxtree LDA.micro funtion output or a data frame containing taxonomic information and LDA results.
#' @return A ggplot object representing the bar plot of the Lefse analysis.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#' @examples
#' tablda = LDA.micro(ps = ps.16s,
#' Top = 100,
#' p.lvl = 0.05,
#' lda.lvl = 1,
#' seed = 11,
#' adjust.p = F)
#' p <- lefse_bar(taxtree = tablda[[2]])
#' p
lefse_bar = function(taxtree = tablda[[2]]){
  # taxtree = tablda[[2]]
  taxtree$ID = row.names(taxtree)
  head(taxtree)
  taxtree$ID = gsub("_Rank_",";",taxtree$ID)
  taxtree <- taxtree %>%
    arrange(class,LDAscore)
  taxtree$ID = factor(taxtree$ID,levels=taxtree$ID)
  taxtree$class = factor(taxtree$class,levels = unique(taxtree$class))

  pbar <- ggplot(taxtree) + geom_bar(aes(y =ID, x = LDAscore,fill = class),stat = "identity") +
    scale_fill_manual(values = unique(taxtree$color))  +
    scale_x_continuous(limits = c(0,max(taxtree$LDAscore)*1.2))
  return(pbar)
}
#' @title Visualization of characteristic microorganisms only at OTU/ASV level
#' @description
#' This function generates a customized bar plot for Lefse analysis results, displaying the LDA scores
#' of only asv/otu classification levels. It arranges the bars based on the LDA scores
#' and colors them according to the Group.
#' @param ps phyloseq object of microbiome,including the Classified information of Characteristic microorganisms.
#' @param taxtree LDA.micro funtion output or a data frame containing taxonomic information and LDA results.
#' @return A ggplot object representing the bar plot of the Lefse analysis,only OTU/ASV level.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#' @examples
#' tablda = LDA.micro(ps = ps.16s,
#' Top = 100,
#' p.lvl = 0.05,
#' lda.lvl = 1,
#' seed = 11,
#' adjust.p = F)
#' p <- lefse_bar2(ps = ps.16s,taxtree = tablda[[2]])
#' p
lefse_bar2 = function(ps = pst,taxtree = tablda[[2]]){
  taxtree = tablda[[2]]
  taxtree$ID = row.names(taxtree)
  head(taxtree)
  taxtree$ID = gsub("_Rank_",";",taxtree$ID)

  taxtree$rankid = sapply(strsplit(taxtree$ID, "__"), `[`, 1)
  tem = taxtree %>% dplyr::filter(rankid %in% c("st"))
  tem$ASV = sapply(strsplit(tem$ID, "st__"), `[`, 2)
  row.names(tem) = tem$Genus
  tem$ID = tem$ASV
  tem <- tem %>%
    arrange(class,LDAscore)
  tem$ID = factor(tem$ID,levels=tem$ID)
  tem$class = factor(tem$class,levels = unique(tem$class))
  head(tem)

  tax.t = ps %>% subset_taxa.wt("OTU",as.character(tem$ID)) %>% vegan_tax() %>% as.data.frame() %>%
    rownames_to_column("ID")
  head(tax.t)
  tem = tem %>% left_join(tax.t,by = c("ASV"= "ID"))
  tem$ID2 = paste0(tem$Genus,"(",tem$ID,")")
  if (unique(tem$class) %>% length() == 2) {
    # a = unique(tem$class)
    # tem2 = tem %>% dplyr::filter(class == a[2]) %>% mutate(LDAscore2 = LDAscore * -1)
    # tem$LDAscore[tem$class == a[2]] = tem2$LDAscore2
    # pbar <- ggplot(tem) + geom_bar(aes(y =ID, x = LDAscore,fill = class),stat = "identity") +
    #   scale_fill_manual(values = unique(taxtree$color)) + theme_classic()
    tem$just = 1
    a = unique(tem$class) %>% as.character()
    tem2 = tem %>% dplyr::filter(class == a[2]) %>% mutate(LDAscore2 = LDAscore * -1,
                                                    just2 = just *0
    )
    tem$LDAscore[tem$class == a[2]] = tem2$LDAscore2
    tem$just[tem$class == a[2]] = tem2$just2
    head(tem)

    pbar <- ggplot(tem) + geom_bar(aes(y =ID, x = LDAscore,fill = class),stat = "identity") +
      scale_fill_manual(values = unique(taxtree$color))  +
      geom_text(aes(y =ID, x =0,label = ID2),hjust = tem$just) +
      theme_classic()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()
      )

  } else{
    pbar <- ggplot(tem) + geom_bar(aes(y =ID2, x = LDAscore,fill = class),stat = "identity") +
      scale_fill_manual(values = unique(taxtree$color)) +
      scale_x_continuous(limits = c(0,max(taxtree$LDAscore)*1.2)) +
      theme_classic()+
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title.y = element_blank(),
            axis.line.y = element_blank()
      )
  }
  return(pbar)
}
