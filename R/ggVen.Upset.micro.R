# 绘制韦恩图和upset图表
#
# This is the first function named 'VenUpset'
# which draw Ven plot and Upset plot with otutab and metadata, and reture  a base plot object
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'
#   Generate doc:              'Ctrl + Shift + Alt + R'

#' @title Plotting Ven and Upset plot for each group
#' @description Input otutab and metadata
#' VennDiagram::Ven plot
#' UpSetR::Upset plot
#' @param otutab rarefied OTU table, typical output of usearch -otutab_norm or -otutab_rare,
#' @param metadata matrix or dataframe, including sampleID and groupID;
#' @param group column name for groupID.
#' @param rep repeat number of each group
#' @return base plot object.
#' @author Contact: Yong-Xin Liu \email{metagenome@@126.com}
#' @references
#'
#' Zhang, J., Zhang, N., Liu, Y.X., Zhang, X., Hu, B., Qin, Y., Xu, H., Wang, H., Guo, X., Qian, J., et al. (2018).
#' Root microbiota shift in rice correlates with resident time in the field and developmental stage.
#' Sci China Life Sci 61, DOI: \url{https://doi.org/10.1007/s11427-018-9284-4}
#'
#' @seealso Ven-Upset
#' @examples
#' # Set four parameters: otutab, metadata, group and rep
#' VenUpset(otu = otutab,map = metadata,group = "genotype",rep = 6)
#' @export

#代码测试


# #清空内存
# rm(list=ls())
# load("../data/otutab.rda")
# load("../data/metadata.rda")
#
# otu = otutab
# map = metadata
# result = VenUpset(otu = otutab,map = metadata,group = "genotype",rep = 6)
# grid.draw(T)
#
# ## upset 我不会直接输出一个图形对象，我现在先在这里把数据提出来在出图
# upset(result[[2]], sets = colnames(result[[2]]),
#       number.angles = 30, point.size = 2, line.size = 1,
#       mainbar.y.label = "OTU", sets.x.label = "OTU Per Treatment",
#       text.scale = c(2, 2, 2,2, 2, 2),mb.ratio = c(0.7, 0.3),order.by = "freq",keep.order = TRUE,
#       queries = list(list(query = intersects, params =
#                             list(colnames(result[[2]])), color = "red", active = T),
#                      list(query = intersects, params =
#                             list(colnames(result[[2]])), color = "red", active = T),
#                      list(query = intersects, params =
#                             list(colnames(result[[2]])), color = "red", active = T)))

# otu = "./otu.txt"
# tax = "./tax.txt"
# map = "./map.txt"
# group = "Group"
# rep = 9
# path  =flowpath



# require(ggplotify)
# g = as.ggplot(p)


# otutab
ggVen.Upset.micro = function(otu = NULL,
                    tax = NULL,
                    map = NULL,
                    tree = NULL,
                    ps = NULL,
                    group = "Group",
                    # path = path,
                    N = 0.5
                    ){

  # library(UpSetR)
  # library (VennDiagram)

  ps =  ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)

  aa =  ggClusterNet::vegan_otu(ps)
  otu_table = as.data.frame(t(aa))
  count = aa
  countA = count

  sub_design <- as.data.frame(phyloseq::sample_data(ps))

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


  ven2[ven2 < N]  = 0
  ven2[ven2 >=N]  = 1
  ven2 = as.data.frame(ven2)



  #
  ven3 = as.list(ven2)

  # ven_pick = get.venn.partitions(ven3)

  for (i in 1:ncol(ven2)) {


    ven3[[i]] <-  row.names(ven2[ven2[i] == 1,])

  }


  if (length(names(ven3)) == 2) {
    # filename3 = paste(path,"ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    # pdf(file=filename3,width = 12, height = 8)
    T<- VennDiagram::venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue"), #填充颜色
                    col=c('red',"blue"), #圈线颜色
                    cat.col=c('red',"blue"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    grid.draw(T)
    dev.off();


  } else if (length(names(ven3)) == 3) {
    # filename3 = paste(path,"ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    # pdf(file=filename3,width = 18, height = 15)
    T<- VennDiagram::venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue","yellow"), #填充颜色
                    col=c('red',"blue","yellow"), #圈线颜色
                    cat.col=c('red',"blue","yellow"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid::grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    grid::grid.draw(T)
    dev.off()
    grid::grid.draw(T)
  } else if (length(names(ven3)) == 4) {
    # filename3 = paste(path,"ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    # pdf(file=filename3,width = 18, height = 15)
    T<-VennDiagram::venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue","yellow","#7ad2f6"), #填充颜色
                    col=c('red',"blue","yellow","#7ad2f6"), #圈线颜色
                    cat.col=c('red',"blue","yellow","#7ad2f6"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid::grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    grid::grid.draw(T)
    dev.off()
    grid::grid.draw(T)
  }else if (length(names(ven3)) == 5) {
    # filename3 = paste(path,"ven_",paste(names(ven3),sep = "",collapse="-"),".pdf",sep = "",collapse="_")
    # pdf(file=filename3,width = 12, height = 12)
    T<- VennDiagram::venn.diagram(ven3,
                    filename=NULL,
                    lwd=2,#圈线粗度
                    lty=1, #圈线类型
                    fill=c('red',"blue","yellow","#7ad2f6","green"), #填充颜色
                    col=c('red',"blue","yellow","#7ad2f6","green"), #圈线颜色
                    cat.col=c('red',"blue","yellow","#7ad2f6","green"),#A和B的颜色
                    cat.cex = 4,# A和B的大小
                    rotation.degree = 0,#旋转角度
                    main = "",#主标题内容
                    main.cex = 2,#主标题大小
                    sub = "",#亚标题内容
                    sub.cex = 1,#亚标题字大小
                    cex=3,#里面交集字的大小
                    alpha = 0.5,#透明度
                    reverse=TRUE,
                    scaled     = FALSE)
    grid::grid.draw(T)
    dev.off()
    # filename33 = paste(path,"ven",".jpg",sep = "",collapse="_")
    # jpeg(file=filename33)
    grid::grid.draw(T)
    dev.off()
    grid::grid.draw(T)
  }else if (length(names(ven3)) == 6) {

    print("ven not use for more than 6")
  }
 return(list(T,ven2))
}


