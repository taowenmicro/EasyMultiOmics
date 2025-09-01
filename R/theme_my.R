
#' @title Custom Theme and Color Settings for the analyse
#' @description
#' The `theme_my` function provides two custom ggplot2 themes (`mytheme1` and `mytheme2`) and dynamically assigns color palettes based on the number of unique groups (`gnum`) in a phyloseq object. This function is specifically designed for creating consistent and aesthetically pleasing visualizations of microbial community data.
#'
#' @param ps A `phyloseq` object containing the microbial community data and metadata.
#' @param gnum Calculate the Number of Unique Groups in Sample Metadata
#'
#' @author Tao Wen \email{2018203048@@njau.edu.cn}, Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' res = theme_my(ps.16s)
#' mytheme1 = res[[1]]
#' mytheme2 = res[[2]]
#' colset1 = res[[3]]
#' @export
theme_my<- function(ps= ps,gnum=NULL){

#--设定颜色
    if (is.null(gnum)) {
      gnum <- phyloseq::sample_data(ps)$Group %>% unique() %>% length()
    }

#-包括几个常用的主题

# 设置主题#-----
mytheme1 = theme_bw() + theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  legend.position="right",

  legend.title = element_blank(),
  legend.background=element_blank(),
  legend.key=element_blank(),
  # legend.text= element_text(size=7),
  # text=element_text(),
  # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
  plot.title = element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
  axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
  axis.text = element_text(size = 20,face = "bold"),
  axis.text.x = element_text(colour = "black",size = 14),
  axis.text.y = element_text(colour = "black",size = 14),
  legend.text = element_text(size = 15,face = "bold")
)

mytheme2 = theme_bw() + theme(
  panel.background=element_blank(),
  panel.grid=element_blank(),
  legend.position="right",

  legend.title = element_blank(),
  legend.background=element_blank(),
  legend.key=element_blank(),
  # legend.text= element_text(size=7),
  # text=element_text(),
  # axis.text.x=element_text(angle=45,vjust=1, hjust=1)
  plot.title = element_text(vjust = -8.5,hjust = 0.1),
  axis.title.y =element_text(size = 24,face = "bold",colour = "black"),
  axis.title.x =element_text(size = 24,face = "bold",colour = "black"),
  axis.text = element_text(size = 20,face = "bold"),
  axis.text.x = element_text(colour = "black",size = 14,angle = 90),
  axis.text.y = element_text(colour = "black",size = 14),
  legend.text = element_text(size = 15,face = "bold")
)

# #--设定颜色
# gnum <- unique(sample_data(ps)$Group) %>% length()
# gnum

if (gnum <= 9) {
  #设定颜色#------------
  #调用所有这个包中的调色板
  RColorBrewer::display.brewer.all()
  #提取特定个数的调色板颜色，会出图显示
  # RColorBrewer::display.brewer.pal(9,"Set1")
  colset1 <- c(brewer.pal(9,"Set1"))
  colset2 <- brewer.pal(12,"Paired")
  colset3 <- c(brewer.pal(12,"Set1"),brewer.pal(9,"Pastel1"))
  colset4 = colset3
}


if (gnum > 9) {
  #设定颜色#------------
  #调用所有这个包中的调色板
  RColorBrewer::display.brewer.all()
  #提取特定个数的调色板颜色，会出图显示
  # RColorBrewer::display.brewer.pal(9,"Set1")
  colset1 <- colorRampPalette(RColorBrewer::brewer.pal(9,"Set1"))(gnum)
  colset2 <- brewer.pal(12,"Paired")
  colset3 <- c(brewer.pal(12,"Set1"),brewer.pal(9,"Pastel1"))
  colset4 = colset3
}



return(list(
  mytheme1 = mytheme1,
  mytheme2 = mytheme2,
  colset1 = colset1,
  colset2 = colset2,
  colset3 = colset3,
  colset4 = colset4
))

}



