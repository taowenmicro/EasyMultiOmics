#' @title Perform Orthogonal Partial least squares Discriminant Analysis (OPLS-DA) on metabolite data
#' @description
#' The oplsda.ms function conducts Orthogonal Partial least squares Discriminant Analysis (OPLS-DA) on metabolite data using the ropls package.
#' It visualizes the results of OPLS-DA including score plots, vip plots, and Permutation plots.

#' @param ps A phyloseq format file used as an alternative for the input containing metabolite composition table,
#' metabolite classification table, and sample metadata.
#' @param Group Column name for groupID in map table(sample metadata).
#' @param orthoI_num The number of orthogonal components in the OPLS-DA.
#' @param oplspath Path to save the output files.
#' @param vipplot_top Top number of VIP features to consider in plots.
#' @param ncol Number of columns for arranging plots.
#' @param nrow Number of rows for arranging plots.
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @export
#' @examples
#' {library(pacman)
#' library(ggsci)
#' library(ropls)
#' res = oplsda.ms(ps = ps.ms, ncol=3,nrow = 1)
#' p1 = res[[1]]
#' p1
#' p2 = res[[2]]
#' p2
#' p3 = res[[3]]
#' p3
#' dat1 = res[[4]]
#' dat1
#' dat2 = res[[4]]
#' dat2
#' dat3 = res[[6]]
#' dat3}
oplsda.ms <- function(ps=ps,Group="Group",orthoI_num=3,oplspath=NULL,vipplot_top=30,
                      ncol=ncol,nrow=nrow)
  {
  Score_list <- list()
  VIP_list <- list()
  Permutation_list <- list()
 sub_design <- as.data.frame(phyloseq::sample_data(ps))
Desep_group <- as.character(levels(as.factor(sub_design$Group)))
aaa = combn(Desep_group,2)
# i=1
for (i in 1:dim(aaa)[2]) {
  Desep_group = aaa[,i]
  group1 = paste(Desep_group[1],Desep_group[2],sep = "__")
  dat = ps %>%
    subset_samples.wt(Group,c(Desep_group[1] ,Desep_group[2])) %>%
    filter_taxa(function(x) sum(x ) > 0 , TRUE) %>%
    vegan_otu()

  group = ps%>%
    scale_micro() %>%
    subset_samples.wt(Group,c(Desep_group[1] ,Desep_group[2])) %>%
    sample_data() %>%
    .$Group %>%
    as.factor()
  #oplsda分析
  set.seed(1234)
  oplsda = opls(dat,group, predI =1,orthoI = orthoI_num,permI = 200,crossvalI = 8)#自动选择要求第一主成分显著
  #提取模型解释度
  oplsda.eig <- data.frame(oplsda@modelDF)
  #提取检验数据
  permutation_data <- data.frame(oplsda@suppLs[["permMN"]])
  permutation_name <- paste("permutation",1:200,sep = "")
  rownames(permutation_data) <- c("model",permutation_name)
  #提取模型中样本得分
  sample.score = oplsda@scoreMN %>%  #得分矩阵
    as.data.frame() %>%
    mutate(Group = group,
           o1=oplsda@orthoScoreMN[,"o1"],
           o2=oplsda@orthoScoreMN[,"o2"],
           o3=oplsda@orthoScoreMN[,"o3"]) #正交矩阵
  head(sample.score)#查看
  #变量重要性，VIP
  #VIP 值帮助寻找重要的代谢物,添加orthol为模型vip值,不添加为第一主成分即组间差异vip值
  ovip <- data.frame( VIP_o1 = getVipVn(oplsda,orthoL = TRUE), #正交模型的变量重要性值。
                      VIP_p1 = getVipVn(oplsda,orthoL = FALSE)) #预测成分的变量的重要性。
  #根据vip筛选差异微生物
  ovip_select <- ovip %>%  dplyr::filter(VIP_p1 >= 1) %>% # 以VIP>=1为阈值进行筛选。
    arrange(VIP_p1) %>% # 根据VIP值排序
    mutate(feature = factor(rownames(.),levels = rownames(.)))
  ovip_select %>% head()
  dim(ovip_select)#vip大于1的物质个数
  #top30
  vip_plot <- ovip_select[(nrow(ovip_select)-vipplot_top+1):nrow(ovip_select),]
  tem = as.data.frame(t(dat))[row.names(vip_plot), ] %>% as.data.frame()
  vip_plot <- cbind(tem, vip_plot)

  #可视化
  #得分图可视化
  p1 <- ggplot(sample.score, aes(p1, o1, color = Group)) +
    geom_hline(yintercept = 0, linetype = 'dashed', size = 0.5) + #横向虚线
    geom_vline(xintercept = 0, linetype = 'dashed', size = 0.5) +
    geom_point(alpha = 0.6, size = 2)+
    xlab(paste("T score[1] (",round(oplsda.eig["p1",1]*100,2), "%)",sep=""))+
    ylab(paste("Orthogonal T score[1] (",round(oplsda.eig["o1",1]*100,2), "%)",sep=""))+
    ggtitle(paste(group1,"         OPLS-DA Scores Plot",sep = ""))+
    stat_ellipse(level = 0.95, linetype = 'solid',
                 size = 1,alpha=0.7) + #添加置信区间
    theme_bw() +
    theme(panel.border=element_rect(colour= "black",fill=NA,size=0.75),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          plot.background=element_blank(),
          plot.title =element_text(face="bold",color="black",size=12,vjust=0.7,hjust = 0.5),
          axis.title.x = element_text(size = 11,colour = "black",face = "bold",vjust = 1),
          axis.title.y = element_text(size = 11,colour = "black",face = "bold",hjust = 0.5,vjust = 1.6),
          axis.text=element_text(color="black",size=10,vjust = 0.5,hjust = 0.5),
          legend.background=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(face="bold",color="black",size=9.2),
          legend.spacing.x=unit(0.2,"cm"),
          legend.key = element_blank(),
          legend.key.size =unit(0.45,"cm"),
          legend.position=c(0.5,0.99),
          legend.direction = "horizontal",
          legend.justification = c(0.5,0.84))

  p1
  #vip棒棒糖
  p2 <- ggplot(vip_plot,aes(x= VIP_p1,y=feature))+
    geom_segment(aes(x=1,xend=VIP_p1,y=feature,yend=feature),linewidth=0.6)+
    geom_point(color="#377EB8",shape = 19,size=2.8)+  theme_bw()+
    xlab("VIP scores")+ylab("Metabolites")+
    labs(title=paste(group1," Top",vipplot_top," different features",sep = ""))+
    scale_x_continuous(limits=c(1,1.08*max(vip_plot$VIP_p1)-0.08),expand = c(0,0))+
    theme(panel.border=element_rect(colour= "black",fill=NA,size=0.75),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          plot.background=element_blank(),
          plot.title =element_text(face="bold",color="black",size=12,vjust=0.7,hjust = 0.5),
          axis.title.x = element_text(size = 11,colour = "black",face = "bold",vjust = 1),
          axis.title.y = element_text(size = 11,colour = "black",face = "bold",hjust = 0.5,vjust = 1.6),
          axis.text.y = element_text(color="black",size=7.5,vjust = 0.5,hjust = 1)
    )

  p2
  # FileName <- paste(oplspath,group1,"_TOP30_VIP.plot", ".pdf", sep = "")
  # ggsave(FileName, p2,width = 12.4/2.54, height =13/2.54)

  #过拟合检验
  permutation_plot <- reshape2::melt(permutation_data[,c("R2Y.cum.","Q2.cum.","sim")],id="sim")
  #R2,Q2截距,斜率，图例标签
  IR2 <- summary(lm(value~sim,data = permutation_plot[which(permutation_plot$variable=="R2Y.cum."),]))$coefficients["(Intercept)","Estimate"]
  IQ2 <- summary(lm(value~sim,data = permutation_plot[which(permutation_plot$variable=="Q2.cum."),]))$coefficients["(Intercept)","Estimate"]

  labelR2 <- expression(paste("R"^2, "Y(cum)",sep = ""))
  labelQ2 <- expression(paste("Q"^2,"(cum)",sep = ""))
  #出图
  p3 <- ggplot(permutation_plot,aes(x= sim,y=value,color=variable,shape = variable))+
    theme_bw()+geom_point(alpha=0.8,size=1.6)+
    scale_color_manual(values=c("Q2.cum."="#377EB8","R2Y.cum."="#E41A1C"),
                       labels=c("R2Y.cum."=labelR2,"Q2.cum."=labelQ2))+
    scale_shape_manual(values=c(19,17),labels=c(labelR2,labelQ2))+
    geom_segment(x=0,y=IR2,xend=1,yend=permutation_data["model","R2Y.cum."],linetype="dashed",color=rgb(95/255,95/255,95/255),linewidth=0.35)+
    geom_segment(x=0,y=IQ2,xend=1,yend=permutation_data["model","Q2.cum."],linetype="dashed",color=rgb(95/255,95/255,95/255),linewidth=0.35)+
    scale_y_continuous(limits = c(min(permutation_plot$value)-0.1,max(permutation_plot$value)+0.1),expand = c(0,0))+
    ggtitle(paste(group1,"      OPLS-DA  permutation",sep = ""))+
    labs(x="Correlation Coefficient",y = expression(paste("R"^2, "Y(cum) and ","Q"^2,"(cum)",sep = "")),
         subtitle =  bquote(Intercepts ~ ":" ~ R^2 ~ Y(cum) ~ "=" ~ .(round(IR2, 2))~"," ~ Q^2 ~ "(cum)"~"=" ~ .(round(IQ2, 2))))+
    theme(panel.border=element_rect(colour= "black",fill=NA,size=0.75),
          panel.grid.minor=element_blank(),
          panel.background=element_blank(),
          plot.background=element_blank(),
          plot.title =element_text(face="bold",color="black",size=12,vjust=0.7,hjust = 0.5),
          plot.subtitle = element_text(face="bold",color="black",size=11,vjust=0.7,hjust = 0.5),
          axis.title.x = element_text(size = 11,colour = "black",face = "bold",vjust = 1),
          axis.title.y = element_text(size = 11,colour = "black",face = "bold",hjust = 0.5,vjust = 1),
          axis.text = element_text(color="black",size=10,vjust = 0.5,hjust = 0.5),
          legend.background=element_blank(),
          legend.title=element_blank(),
          legend.text=element_text(color="black",size=9.5,face="bold",vjust = 0.5,hjust = 0.5),
          legend.spacing.x=unit(0,"cm"),
          legend.key = element_blank(),
          legend.key.size =unit(0.6,"cm"),
          legend.position=c(0.99,0),
          legend.direction = "vertical",
          legend.justification = c(0.95,0))+
    geom_hline(yintercept=0,linetype=1,color = 'black',size = 0.5) +
    geom_vline(xintercept=0,linetype=1,color = 'black',size = 0.5)

  p3
  Score_list[[i]] <- p1
  VIP_list[[i]] <- p2
  Permutation_list[[i]] <- p3
  if (!is.null(oplspath)) {
    fs::dir_create(oplspath,recursive = TRUE)
    file <- paste(oplspath,group1,"_model_statistic.csv",sep="")
    write.csv(oplsda.eig,file,quote =T)
    file <- paste(oplspath,group1,"_model_permutation.csv",sep="")
    write.csv(permutation_data,file,quote =T)
    file <- paste(oplspath,group1,"_sample_score.csv",sep="")
    write.csv(sample.score,file,quote =T)
    file <- paste(oplspath,group1,"_VIP_score.csv",sep="")
    write.csv(ovip,file,quote =T)
    file <- paste(oplspath,group1,"_VIP_difference_top",vipplot_top,".csv",sep="")
    write.csv(vip_plot,file,quote =T)
    FileName <- paste(oplspath,group1,"_score.plot", ".pdf", sep = "")
    ggsave(FileName, p1,width = 14/2.54, height =11/2.54)
    FileName <- paste(oplspath,group1,"_TOP30_VIP.plot", ".pdf", sep = "")
    ggsave(FileName, p2,width = 12.4/2.54, height =13/2.54)
    FileName <- paste(oplspath,group1,"_permutation.plot", ".pdf", sep = "")
    ggsave(FileName, p3,width = 14/2.54, height =12/2.54)
  }
  permutation_plot$desep_group <- group1
  vip_plot$desep_group <- group1
  sample.score$desep_group <- group1
  if (i==1) {
    permutation_plotdata <- permutation_plot
    vip_plotdata <- vip_plot[,c("VIP_o1","VIP_p1" ,"feature","desep_group")]
    sample.scoredata <- sample.score
  }else{
    permutation_plotdata <- rbind(permutation_plotdata,permutation_plot)
    vip_plotdata <- rbind(vip_plotdata,vip_plot[,c("VIP_o1","VIP_p1" ,"feature","desep_group")])
    sample.scoredata <- rbind(sample.scoredata,sample.score)
  }
}
Score_plot = ggpubr::ggarrange(plotlist = Score_list,common.legend = F,ncol = ncol,nrow = nrow)
VIP_plot = ggpubr::ggarrange(plotlist = VIP_list,common.legend = F,ncol = ncol,nrow = nrow)
Permutation_plot = ggpubr::ggarrange(plotlist = Permutation_list,common.legend = F,ncol = ncol,nrow = nrow)
return(list(Score_plot,VIP_plot,Permutation_plot,sample.scoredata,vip_plotdata,permutation_plotdata,
       Score_list,VIP_list, Permutation_list))
}
