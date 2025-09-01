#' @title Random forest modeling for transcriptome functional composition data
#' @description
#' Random forest modeling for transcriptome functional composition data was used to screen for characteristic genes.
#' @param otu Transcriptome functional composition table.
#' @param map Sample metadata.
#' @param tax transcriptome functional classification table.
#' @param tree Phylogenetic tree.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition data, taxonomy table, and sample metadata.
#' @param optimal Select the number of important genes to be displayed.
#' @param rfcv TURE or FALSE,whether need to do cross-validation.
#' @param nrfcvnum Number of cross-validation.
#' @param min Circle diagram inner diameter adjustment.
#' @param group Column name for groupID in map table(sample metadata).
#' @param max Circle diagram outer diameter adjustment.
#' @return list contain ggplot object and table:
#' \item{p1}{A lollipop diagram of the selected number of important genes.}
#' \item{p2}{A circle diagram of the selected number of important genes.}
#' \item{a3}{Taxonomic annotation information for the selected number of important genes and their importance.}
#' \item{pn}{Combined table of model accuracy rates and confusion matrix.}
#' \item{a2}{Taxonomic annotation information for all genes and their importance.}
#' @author
#' Tao Wen \email{2018203048@njau.edu.cn},
#' Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' library(dplyr)
#' library(ggClusterNet)
#' library(caret)
#' library(randomForest)
#' # Example with a phyloseq object
#' data(ps.trans)
#' ps = ps.trans%>% filter_taxa(function(x) sum(x ) > 5 , TRUE)
#' res <- randomforest.trans( ps= ps, group  = "Group", optimal = 50)
#' p1 = res[[1]];p1
#' p2=res[[2]];p2
#' p3=res[[4]];p3
#' dat1 =res[[3]];dat1
#' dat2 =res[[5]];dat2
#'@export



randomforest.trans <- function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
                              ps = NULL,
                              group  = "Group",
                              optimal = 20,
                              rfcv = FALSE,
                              nrfcvnum = 5,
                              min = -1,
                              max = 5,
                              fill = "Level1",
                              lab = "id"
){


  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group) %>% ggClusterNet::scale_micro()
  map = as.data.frame(phyloseq::sample_data(ps))
  #-scaleing relative abundancce#----
  mapping = as.data.frame(phyloseq::sample_data(ps))
  otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
  colnames(otutab) = paste("wentao",colnames(otutab),sep = "")


  # Set classification info.
  otutab$group = factor(mapping$Group)
  colnames(otutab) <- gsub("-","_",colnames(otutab))
  colnames(otutab) <- gsub("[/]","_",colnames(otutab))
  colnames(otutab) <- gsub("[(]","_",colnames(otutab))
  colnames(otutab) <- gsub("[)]","_",colnames(otutab))
  colnames(otutab) <- gsub("[:]","_",colnames(otutab))
  colnames(otutab) <- gsub("[[]","_",colnames(otutab))
  colnames(otutab) <- gsub("[]]","_",colnames(otutab))
  colnames(otutab) <- gsub("[#]","_",colnames(otutab))
  colnames(otutab) <- gsub("[+]","_",colnames(otutab))
  colnames(otutab) <- gsub(" ","_",colnames(otutab))
  colnames(otutab) <- gsub("[,]","_",colnames(otutab))


  model_rf= randomForest::randomForest(group ~ ., data=otutab,
                                       importance=TRUE, proximity=TRUE)
  print(model_rf)
  #--提取混淆矩阵
  Confusion_matrix <- as.data.frame(model_rf$confusion)
  Confusion_matrix$class.error <- round(Confusion_matrix$class.error,3)
  Confusion_matrix$Group = row.names(Confusion_matrix)

  Confusion_matrix <- dplyr::select(Confusion_matrix , Group, everything())
  #-提取正确率
  model_Accuracy_rates <- paste(round(100-tail(model_rf$err.rate[,1],1)*100,2),"%",sep = "")
  model_Accuracy_rates = data.frame(ID = "model Accuracy rates",model_Accuracy_rates = model_Accuracy_rates)
  colnames(model_Accuracy_rates) = c("Random foreest","model")
  tab2 <- ggpubr::ggtexttable(Confusion_matrix, rows = NULL)
  tab1 <- ggpubr::ggtexttable(model_Accuracy_rates, rows = NULL)
  library(patchwork)
  pn <- tab1/tab2

  # if (rfcv) {
  #   result = Micro.rfcv(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = ps,group  = "Group",optimal = 20,nrfcvnum = 6)
  #
  #   prfcv = result[[1]]# plot rfcv
  #   # result[[2]]# plotdata
  #   rfcvtable = result[[3]]# table rfcv
  # } else{
  #   prfcv = NULL
  #   rfcvtable = NULL
  # }
  #

  #------------k可视化
  ###### inportant OTU picked out and plot
  a=as.data.frame(round(randomForest::importance(model_rf), 2))
  a$id=row.names(a)

  row.names(a)  = gsub("wentao","",row.names(a))
  a$id = gsub("wentao","",a$id)
  a2<- dplyr::arrange(a, desc(MeanDecreaseAccuracy)) %>% as.data.frame()
  row.names(a2)=a2$id

  a3=head(a2,n=optimal)
  head(a3)
  if (is.null(ps@tax_table)) {
    tax = data.frame(row.names =  gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))]),
                     ID =  gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))]),
                     class =  gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))]))

  } else{
    tax <- as.data.frame( ggClusterNet::vegan_tax(ps))
  }


  # head(tax)
  tax$org.id = row.names(tax)
  row.names(tax) = gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))])
  tax = tax[rownames(a3),]
  dim(a3)
  # head(tax)
  tax$id = NULL
  a3 = merge(a3,tax,by = "row.names",all = F)
  row.names(a3) = a3$Row.names
  a3$Row.names = NULL
  # dim(a3)
  # head(a3)
  if (is.null(ps@tax_table)) {
    tax = data.frame(row.names =  gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))]),
                     ID =  gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))]),
                     class =  gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))]))

  } else{
    tax <- as.data.frame( ggClusterNet::vegan_tax(ps))
  }
  # head(tax)
  tax$org.id = row.names(tax)
  row.names(tax) = gsub("wentao","",colnames(otutab)[-(length(colnames(otutab)))])
  tax = tax[rownames(a2),]
  dim(a2)
  # head(a2)
  tax$id = NULL
  a2 = merge(a2,tax,by = "row.names",all = F)
  row.names(a2) = a3$Row.names
  a2$Row.names = NULL


  OTU =  ggClusterNet::vegan_otu(ps)
  ### pice mapping
  design = as.data.frame(phyloseq::sample_data(ps))

  #mean abundance by groups
  iris.split <- split(as.data.frame(OTU),as.factor(design$Group))
  iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
  norm2 <- do.call(rbind,iris.apply)%>% # combine result
    t()
  colnames(norm2) = paste(colnames(norm2),"mean",sep = "")



  ind_fal = merge(a3,norm2,by = "row.names",all = F)
  head(ind_fal)



  # plot fire bar 1
  if (is.null(lab)) {
    lab = "id"
  }

  # head(a3)
  if (is.null(fill)) {
    a3$fill = "class"
    fill = "fill"
  }
  tem = a3 %>% arrange(desc(MeanDecreaseAccuracy))
  head(tem)
  a3$id = factor(a3$id,level =tem$id[length(tem$id):1] )
  p1 <- ggplot(a3, aes(x = MeanDecreaseAccuracy, y = id,
                       fill = !!sym(fill),color = !!sym(fill))) +
    geom_point(size=6,pch=21)+
    geom_segment(aes(yend=id),xend=0,size=3)+
    geom_text(aes(x =MeanDecreaseAccuracy*1.1,  label = !!sym(lab)),size = 3,color = "black")



  # plot 2
  a3<- dplyr::arrange(a3, desc(MeanDecreaseAccuracy))
  a3$iid = paste(1:length(a3$id))
  angle1 = 90 - 360 * ( as.numeric(a3$iid) - 0.5) /length(a3$id)
  a3$id = factor(a3$id,levels = a3$id)
  p2 = a3  %>%
    ggplot(aes(x = factor(id), y = MeanDecreaseAccuracy ,label = id)) +
    geom_bar(stat = 'identity', position = 'dodge',fill = "blue") +
    # scale_fill_manual(values = mi)+
    geom_text(hjust = 0, angle = angle1, alpha = 1) +
    coord_polar() +
    ggtitle('') +
    ylim(c(min,max))+
    theme_void()
  p2

  return(list(p1,p2,a3,pn,a2))
}


## rfcv function
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title For cross-validation of microbiome data
#' @description For cross-validation of microbiome data
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param tax taxonomy table
#' @param ps phyloseq object of microbiome
#' @param Group column name for groupID in map table.
#' @param optimal important OTU number which selected
#' @param rfcv TURE or FELSE,whether need to do cross-validation
#' @param nrfcvnum Number of cross-validation
#' @details
#' @return list contain ggplot object and table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @examples
#' # data form github
#' result = Micro.rfcv(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = ps_rela,group  = "Group",optimal = 20,nrfcvnum = 6)
#' prfcv = result[[1]]# plot rfcv
# result[[2]]# plotdata
#' rfcvtable = result[[3]]# table rfcv
#'@export
# Micro.rfcv = function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
#                       ps = NULL,group  = "Group",optimal = 20,nrfcvnum = 5){
#
#
#   ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
#   otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
#   # Set classification info.
#   otutab$group = factor(mapping$Group)
#   colnames(otutab) <- gsub("-","_",colnames(otutab))
#   # rfcv for select···
#   n = ncol(otutab)-1
#   myotutab_t= otutab[1:n]
#   set.seed(315)
#   result= rfcv(myotutab_t, otutab$group, cv.fold=5, scale = "log", step = 0.9)
#   # with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
#   result1 = result
#   error.cv = data.frame(num = result$n.var, error.1 =  result$error.cv)
#   for (i in 316:(314+ nrfcvnum)){
#     print(i)
#     set.seed(i)
#     result= rfcv(myotutab_t, otutab$group, cv.fold=5, scale = "log", step = 0.9)
#     error.cv = cbind(error.cv, result$error.cv)
#   }
#   n.var = error.cv$num
#   error.cv = error.cv[,2:6]
#   colnames(error.cv) = paste('err',1:5,sep='.')
#   err.mean = apply(error.cv,1,mean)
#   allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
#   head(allerr)
#   data <- gather(allerr, key = "group", value = "value",-num)
#   head(data)
#
#   p <- ggplot() +
#     geom_line(data = data,aes(x = num, y = value,group = group), colour = 'grey') +
#     geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') +
#     coord_trans(x = "log2") +
#     scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
#     labs(title=paste('Training set (n = ', dim(otutab)[1],')', sep = ''),
#          x='Number of families ', y='Cross-validation error rate') +
#     annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep=""))
#   return(list(plot = p,plotdata = data,origdata = allerr))
# }


# Roc function
# You can learn more about package at:
#
#   https://github.com/microbiota/amplicon

#' @title Comparison of three machine methods (randomforest,SVM,GLM).
#' @description Comparison of three machine methods (randomforest,SVM,GLM).
#' @param otu OTU/ASV table;
#' @param map Sample metadata;
#' @param tax taxonomy table
#' @param ps phyloseq object of microbiome
#' @param Group column name for groupID in map table.
#' @param repnum Modeling times
#' @details
#' @return list contain ggplot object and table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn}, Yong-Xin Liu \email{yxliu@@genetics.ac.cn}
#' @references
#'
#' Jingying Zhang, Yong-Xin Liu, Na Zhang, Bin Hu, Tao Jin, Haoran Xu, Yuan Qin, Pengxu Yan, Xiaoning Zhang, Xiaoxuan Guo, Jing Hui, Shouyun Cao, Xin Wang, Chao Wang, Hui Wang, Baoyuan Qu, Guangyi Fan, Lixing Yuan, Ruben Garrido-Oter, Chengcai Chu & Yang Bai.
#' NRT1.1B is associated with root microbiota composition and nitrogen use in field-grown rice.
#' Nature Biotechnology, 2019(37), 6:676-684, DOI: \url{https://doi.org/10.1038/s41587-019-0104-4}
#'
#' @examples
#' result = MicroRoc( ps = ps,group  = "Group")
#' #--提取roc曲线
#' result[[1]]
#' #提取AUC值
#' result[[2]]
#'@export
# MicroRoc <- function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
#                      ps = NULL,group  = "Group",repnum = 5){
#
#
#   ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
#   mapping = as.data.frame(phyloseq::sample_data(ps))
#   otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
#   colnames(otutab) <- gsub("-","_",colnames(otutab))
#   colnames(otutab) <- gsub("-","_",colnames(otutab))
#   colnames(otutab) <- gsub("[/]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[(]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[)]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[:]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[[]","",colnames(otutab))
#   colnames(otutab) <- gsub("[]]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[#]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[+]","_",colnames(otutab))
#   colnames(otutab) <- gsub(" ","_",colnames(otutab))
#   colnames(otutab) <- gsub("[;]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[,]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[?]","_",colnames(otutab))
#   test = as.data.frame((otutab))
#   mapping$Group = as.factor(mapping$Group)
#   test$group = factor(mapping$Group)
#   colnames(test) = paste("OTU",colnames(test),sep = "")
#
#   # random forest
#
#   test = dplyr::select(test,OTUgroup,everything())
#   train = test
#   folds <- createFolds(y=test[,1],k=repnum)
#   AUC =c()
#
#   max=0
#   num=0
#   fc<-as.numeric()#fc 为测试数据真实分组信息
#   mod_pre<-as.numeric()
#   i= 1
#   for(i in 1:repnum){
#     fold_test<-train[folds[[i]],]
#     fold_train<-train[-folds[[i]],]
#
#
#     colnames(fold_test) <- gsub("-","_",colnames(fold_test))
#     colnames(fold_train) <- gsub("-","_",colnames(fold_train))
#
#     model<-randomForest(OTUgroup~.,data=fold_train, importance=TRUE, proximity=TRUE)
#     print(model)
#     model_pre<-predict(model,newdata = fold_test,type="prob")
#     fc<-append(fc,as.factor(fold_test$OTUgroup))
#     mod_pre<-append(mod_pre,model_pre[,2])
#   }
#
#
#   #- pick data and plot
#   pred <- prediction(mod_pre, fc)
#   perf <- performance(pred,"tpr","fpr")
#   x <- unlist(perf@x.values)  ##提取x值
#   y <- unlist(perf@y.values)
#   plotdata <- data.frame(x,y)
#   names(plotdata) <- c("x", "y")
#   AUC[1] = paste("rf AUC:",round(performance(pred,'auc')@y.values[[1]],3),sep = " ")
#   head(plotdata)
#   rf.rocdata = plotdata
#   g0 <- ggplot(plotdata) +
#     geom_path(aes(x = x, y = y, colour = x), size=1,color = "red") +
#     labs(x = "False positive rate", y = "Ture positive rate") + # , title ="Random Forest"
#     annotate("text", x=0.75, y=0.5, label=paste("Red: ",AUC[1],sep = ""))
#   df<-cbind(fc,as.numeric(mod_pre))
#   #-svm
#   max=0
#   num=0
#   fc<-as.numeric()
#   mod_pre<-as.numeric()
#
#   for(i in 1:repnum){
#     fold_test<-train[folds[[i]],]
#     # head(fold_test)
#     fold_train<-train[-folds[[i]],]
#     model<-svm(OTUgroup~.,data=fold_train,probability=TRUE)
#     model
#     model_pre<-predict(model,newdata = fold_test,decision.values = TRUE, probability = TRUE)
#     fc<-append(fc,(fold_test$OTUgroup))
#     mod_pre<-append(mod_pre,as.numeric(attr(model_pre, "probabilities")[,2]))
#
#   }
#
#   pred <- prediction(mod_pre, fc)
#   perf <- performance(pred,"tpr","fpr")
#   x <- unlist(perf@x.values)
#   y <- unlist(perf@y.values)
#   plotdata <- data.frame(x,y)
#   names(plotdata) <- c("x", "y")
#   svm.rocdata = plotdata
#   AUC[2] = paste("svm AUC:",round(performance(pred,'auc')@y.values[[1]],3),sep = " ")
#
#
#   g1 <- g0 +
#     geom_path(data = plotdata,aes(x = x, y = y, colour = x), size=1,color = "blue") +
#     # labs(x = "False positive rate", y = "Ture positive rate") +
#     annotate("text", x=0.75, y=0.4, label=paste("Blue: ",AUC[2],sep = ""))
#
#   df<-cbind(df,cbind(fc,mod_pre))
#
#   #GLM
#   max=0
#   num=0
#   fc<-as.numeric()
#   mod_pre<-as.numeric()
#   for(i in 1:repnum){
#     fold_test<-train[folds[[i]],]
#     fold_train<-train[-folds[[i]],]
#     model<-glm(OTUgroup~.,family='binomial',data=fold_train)
#     model
#     model_pre<-predict(model,type='response',newdata=fold_test)
#     model_pre
#
#     fc<-append(fc,fold_test$OTUgroup)
#     mod_pre<-append(mod_pre,as.numeric(model_pre))
#   }
#
#   pred <- prediction(mod_pre, fc)
#   perf <- performance(pred,"tpr","fpr")
#   x <- unlist(perf@x.values)  ##提取x值
#   y <- unlist(perf@y.values)
#   plotdata <- data.frame(x,y)
#   names(plotdata) <- c("x", "y")
#   glm.rocdata = plotdata
#   AUC[3] = paste("GLM AUC:",round(performance(pred,'auc')@y.values[[1]],3),sep = " ")
#
#   g2 <- g1 +
#     geom_path(data = plotdata,aes(x = x, y = y, colour = x), size=1,color = "black") +
#     labs(x = "False positive rate", y = "Ture positive rate") +
#     annotate("text", x=0.75, y=0.3, label=paste("Black: ",AUC[3],sep = ""))
#
#   g2
#
#   df<-cbind(df,cbind(fc,mod_pre))
#
#
#   return(list(g2,AUC,df,list(rf.rocdata=rf.rocdata,svm.rocdata =svm.rocdata,glm.rocdata = glm.rocdata )))
# }
#
#
#
#
# MicroRF_GC <- function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
#                        ps = NULL,group  = "Group",optimal = 20,rfcv = FALSE,nrfcvnum = 5,min = -1,max = 5
# ){
#
#
#   ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group) %>% ggClusterNet::scale_micro()
#   map = as.data.frame(phyloseq::sample_data(ps))
#   #-scaleing relative abundancce#----
#   mapping = as.data.frame(phyloseq::sample_data(ps))
#   otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
#   tem = colnames(otutab)
#   colnames(otutab) = paste("RE",1:length(colnames(otutab)),sep = "")
#
#   tem2 = data.frame(ID = colnames(otutab),name = tem)
#   # colnames(otutab) = paste("wentao",colnames(otutab),sep = "")
#
#
#   # Set classification info.
#   otutab$group = factor(mapping$Group)
#   # colnames(otutab) <- gsub("-","_",colnames(otutab))
#   model_rf= randomForest::randomForest(group ~ ., data=otutab, importance=TRUE, proximity=TRUE)
#   print(model_rf)
#   #--提取混淆矩阵
#   Confusion_matrix <- as.data.frame(model_rf$confusion)
#   Confusion_matrix$class.error <- round(Confusion_matrix$class.error,3)
#   Confusion_matrix$Group = row.names(Confusion_matrix)
#
#   Confusion_matrix <- dplyr::select(Confusion_matrix , Group, everything())
#   #-提取正确率
#   model_Accuracy_rates <- paste(round(100-tail(model_rf$err.rate[,1],1)*100,2),"%",sep = "")
#   model_Accuracy_rates = data.frame(ID = "model Accuracy rates",model_Accuracy_rates = model_Accuracy_rates)
#   colnames(model_Accuracy_rates) = c("Random foreest","Fu wilt model")
#   tab2 <- ggpubr::ggtexttable(Confusion_matrix, rows = NULL)
#   tab1 <- ggpubr::ggtexttable(model_Accuracy_rates, rows = NULL)
#   library(patchwork)
#   pn <- tab1/tab2
#
#   if (rfcv) {
#     result = Micro.rfcv(otu = NULL,tax = NULL,map = NULL,tree = NULL ,ps = ps_rela,group  = "Group",optimal = 20,nrfcvnum = 6)
#
#     prfcv = result[[1]]# plot rfcv
#     # result[[2]]# plotdata
#     rfcvtable = result[[3]]# table rfcv
#   } else{
#     prfcv = NULL
#     rfcvtable = NULL
#   }
#
#
#   #------------k可视化
#   ###### inportant OTU picked out and plot
#   a=as.data.frame(round(randomForest::importance(model_rf), 2))
#   a$id=row.names(a)
#   head(a)
#
#   row.names(a)  = tem
#   a$id = tem
#   a2<- dplyr::arrange(a, desc(MeanDecreaseAccuracy)) %>% as.data.frame()
#   row.names(a2)=a2$id
#   # optimal = 40
#   a3=head(a2,n=optimal)
#
#   OTU =  ggClusterNet::vegan_otu(ps)
#   ### pice mapping
#   design = as.data.frame(phyloseq::sample_data(ps))
#
#   #mean abundance by groups
#   iris.split <- split(as.data.frame(OTU),as.factor(design$Group))
#   iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
#   norm2 <- do.call(rbind,iris.apply)%>% # combine result
#     t()
#   colnames(norm2) = paste(colnames(norm2),"mean",sep = "")
#
#
#
#   ind_fal = merge(a3,norm2,by = "row.names",all = F)
#   head(ind_fal)
#
#
#
#   # plot fire bar 1
#   head(a3)
#   p1 <- ggplot(a3, aes(x = MeanDecreaseAccuracy, y = reorder(id,MeanDecreaseAccuracy))) +
#     geom_point(size=6,pch=21,fill = "#9ACD32",color = "#9ACD32")+
#     geom_segment(aes(yend=id),xend=0,size=3,color = "#9ACD32")+
#     geom_label(aes(x =MeanDecreaseAccuracy*1.1,  label = id),size = 3) + theme_classic()
#
#
#
#   # plot 2
#   a3<- dplyr::arrange(a3, desc(MeanDecreaseAccuracy))
#   a3$iid = paste(1:length(a3$id))
#   angle1 = 90 - 360 * ( as.numeric(a3$iid) - 0.5) /length(a3$id)
#   a3$id = factor(a3$id,levels = a3$id)
#   p2 = a3  %>%
#     ggplot(aes(x = factor(id), y = MeanDecreaseAccuracy ,label = id)) +
#     geom_bar(stat = 'identity', position = 'dodge',fill = "blue") +
#     # scale_fill_manual(values = mi)+
#     geom_text(hjust = 0, angle = angle1, alpha = 1) +
#     coord_polar() +
#     # ylim(c(min,max))+
#     theme_void()
#   p2
#
#   return(list(p1,p2,prfcv,rfcvtable,a3,pn))
# }
#
#
#
#
#
# # ps = ps01 %>%
# #   scale_micro("sampling") %>%
# #   scale_micro("log") %>%
# #   subset_samples.wt("Group",c("CD_R0" ,"CD_NR0")) %>%
# #   tax_glom_wt(6)
#
# # MicroRoc2 ( ps = ps,group  = "Group",
# #                       repnum = 5,
# #                       method = "bys")
# #
# # MicroRoc2 ( ps = ps,group  = "Group",
# #             repnum = 5,
# #             method = "nnet")
# #
# # MicroRoc2 ( ps = ps,group  = "Group",
# #             repnum = 5,
# #             method = "rf")
# #
# #
# # MicroRoc2 ( ps = ps,group  = "Group",
# #             repnum = 5,
# #             method = "mboost")
#
#
#
# MicroRoc2 <- function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
#                       ps = NULL,
#                       group  = "Group",
#                       repnum = 5,
#                       method = "rf"
# ){
#   ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
#   mapping = as.data.frame(phyloseq::sample_data(ps))
#   otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
#   colnames(otutab) <- gsub("-","_",colnames(otutab))
#   colnames(otutab) <- gsub("-","_",colnames(otutab))
#   colnames(otutab) <- gsub("[/]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[(]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[)]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[:]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[[]","",colnames(otutab))
#   colnames(otutab) <- gsub("[]]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[#]","_",colnames(otutab))
#   colnames(otutab) <- gsub("[+]","_",colnames(otutab))
#   colnames(otutab) <- gsub(" ","_",colnames(otutab))
#   colnames(otutab) <- gsub("[;]","_",colnames(otutab))
#   test = as.data.frame((otutab))
#   mapping$Group = as.factor(mapping$Group)
#   test$group = factor(mapping$Group)
#   colnames(test) = paste("OTU",colnames(test),sep = "")
#
#   # random forest
#
#   test = dplyr::select(test,OTUgroup,everything())
#   train = test
#
#   folds <- createFolds(y=test[,1],k=repnum)
#   AUC =c()
#
#   max=0
#   num=0
#   fc<-as.numeric()#fc 为测试数据真实分组信息
#   mod_pre<-as.numeric()
#   i= 1
#
#   for(i in 1:repnum){
#     fold_test<-train[folds[[i]],]
#     fold_train<-train[-folds[[i]],]
#
#
#     colnames(fold_test) <- gsub("-","_",colnames(fold_test))
#     colnames(fold_train) <- gsub("-","_",colnames(fold_train))
#
#     if (method == "rf" ) {
#       model<-randomForest(OTUgroup~.,data=fold_train, importance=TRUE, proximity=TRUE)
#       print(model)
#       model_pre<-predict(model,newdata = fold_test,type="prob")
#       fc<-append(fc,as.factor(fold_test$OTUgroup))
#       mod_pre<-append(mod_pre,model_pre[,2])
#     }
#
#
#     if (method == "nnet") {
#       library(nnet)
#       set.seed(50)
#       a_nn = nnet(OTUgroup ~ .,
#                   fold_train,size = 2,
#                   rang = 0.3,
#                   decay = 5e-4,
#                   MaxNWts=83581,
#                   maxit = 220)
#       ##type参数为class,因此输出的是预测的类标号而非概率矩阵
#       model_pre = predict(a_nn,
#                           newdata = fold_test[,-1],
#                           type = 'raw')
#       fc<-append(fc,as.factor(fold_test$OTUgroup))
#       mod_pre<-append(mod_pre,model_pre[,1])
#     }
#
#     if (method == "bys") {
#       library(e1071)
#       a_nb <- naiveBayes(OTUgroup ~ .,
#                          fold_train)
#       # a_nb
#       # print(a_nb)
#       # summary(a_nb)
#       #预测结果
#       # predict(a_nb,newdata = fold_test)
#
#       ##type参数为class,因此输出的是预测的类标号而非概率矩阵
#       model_pre = predict(a_nb,
#                           newdata = fold_test[,-1],
#                           type = 'raw')
#       fc<-append(fc,as.factor(fold_test$OTUgroup))
#       mod_pre<-append(mod_pre,model_pre[,2])
#     }
#
#
#     if (method == "mboost") {
#       library(caret)
#       library(mboost)
#
#       a_mboost = blackboost(OTUgroup ~ .,
#                             fold_train, family = Binomial())
#       # a_mboost
#       # predict(a_mboost,fold_test, type = 'class')
#       model_pre = predict(a_mboost,
#                           newdata = fold_test[,-1],
#                           type = 'link')
#       # model_pre[,1] = 1 - model_pre[,1]
#       fc<-append(fc,as.factor(fold_test$OTUgroup))
#       mod_pre<-append(mod_pre,model_pre[,1])
#     }
#
#
#
#   }
#
#
#   #- pick data and plot
#   pred <- prediction(mod_pre, fc)
#   perf <- performance(pred,"tpr","fpr")
#   x <- unlist(perf@x.values)  ##提取x值
#   y <- unlist(perf@y.values)
#   plotdata <- data.frame(x,y)
#   names(plotdata) <- c("x", "y")
#   AUC[1] = paste("rf AUC:",round(performance(pred,'auc')@y.values[[1]],3),sep = " ")
#   head(plotdata)
#   rf.rocdata = plotdata
#   g0 <- ggplot(plotdata) +
#     geom_path(aes(x = x, y = y, colour = x), size=1,color = "red") +
#     labs(x = "False positive rate", y = "Ture positive rate") + # , title ="Random Forest"
#     annotate("text", x=0.75, y=0.5, label=paste("Red: ",AUC[1],sep = ""))
#   g0
#   return(list(g0))
# }
#
#
# randomforest.wt = function(
    #     pst = pst,
#     ROC = F,
#     rfcv = F,
#     optimal = 50
# ){
#
#   # if (ROC ){
#   #   #--三种机器学习方法评测
#   #   result = MicroRoc( ps = pst,group  = "Group",repnum = 5)
#   #   #--提取roc曲线
#   #   p <- result[[1]] +
#   #     mytheme1
#   #
#   #   #提取AUC值
#   #   data <- result[[2]]
#   #   filename = paste(matpath,"/three_method_AUCvalue.csv",sep = "")
#   #   write.csv(data,filename,quote = F)
#   #   data <- result[[3]]
#   #   filename = paste(matpath,"/three_method_AUCdata.csv",sep = "")
#   #   write.csv(data,filename,quote = F)
#   #   filename = paste(matpath,"/three_method_AUC_plot.pdf",sep = "")
#   #   ggsave(filename,p,width = 8,height = 8)
#   #   filename = paste(matpath,"/three_method_AUC_plot.jpg",sep = "")
#   #   ggsave(filename,p,width = 8,height = 8)
#   # }
#   #
#
#   mapping = as.data.frame(phyloseq::sample_data(pst))
#   #--随机森林全套-如果圈图尚未显示前面几个，就设定max大一点
#   result = MicroRF(ps = pst,
#                    group  = "Group",
#                    optimal = optimal,
#                    fill = "Phylum",
#                    rfcv = F,
#                    nrfcvnum = 5,
#                    lab = "Genus",
#                    min = -1,max = 5)
#   #火柴图展示前二十个重要的OTU
#   p1 <- result[[1]] +
#     mytheme1
#   p1
#
#
#   filename = paste(matpath,"/randonforest_loading.pdf",sep = "")
#   ggsave(filename,p1,width = 8,height = optimal/2)
#   filename = paste(matpath,"/randonforest_loading.jpg",sep = "")
#   ggsave(filename,p1,width = 8,height = optimal/2)
#
#   # 圈图展示
#   p <- result[[2]]
#
#   filename = paste(matpath,"/randonforest_loading_circle.pdf",sep = "")
#   ggsave(filename,p,width = 8,height = 10)
#   filename = paste(matpath,"/randonforest_loading_circle.jpg",sep = "")
#   ggsave(filename,p,width = 8,height = 10)
#
#   p <- result[[6]]
#
#   filename = paste(matpath,"/Show_model.pdf",sep = "")
#   ggsave(filename,p,width = 8,height = 4)
#   filename = paste(matpath,"/Show_model.jpg",sep = "")
#   ggsave(filename,p,width = 8,height = 4)
#
#
#   if (rfcv){
#     # 展示交叉验证结果
#     p <- result[[3]]
#     filename = paste(matpath,"/randonforest_cross_check.pdf",sep = "")
#     ggsave(filename,p,width = 8,height = 12)
#     data <- result[[4]]
#     filename = paste(matpath,"/randomforest_cross_data.csv",sep = "")
#     write.csv(data,filename,quote = F)
#   }
#
#   data <- result[[5]]
#   filename = paste(matpath,"/randomforest_data.csv",sep = "")
#   write.csv(data,filename,quote = F)
#
#
#
#   ps1 <- pst %>%
#     subset_taxa.wt("OTU", as.character(data$org.id))
#   ps1
#
#
#
#   otu_table = as.data.frame(t(vegan_otu(ps1)))
#   head(otu_table)
#
#   design = as.data.frame(sample_data(ps1))
#   ## 计算相对丰度，计算每个物种丰度均值，按照均值排序
#   OTU = as.matrix(otu_table)
#   norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
#   norma = norm %>%
#     t() %>% as.data.frame()
#   #数据分组计算平均值
#   iris.split <- split(norma,as.factor(design$Group))
#
#   iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
#   # 组合结果
#   norm2 <- do.call(rbind,iris.apply)%>%
#     t()
#   norm2 = as.data.frame(norm2)
#   norm2$mean=apply(norm2,1,mean)
#   norm2$ID = row.names(norm2)
#   colnames(norm2)
#   abun = norm2
#   head(abun)
#   abun$mean = NULL
#   library(reshape2)
#   abun_a = melt(abun,
#                 id.var = c("ID"),
#                 variable.name = "id",
#                 value.name = "count")
#
#   head(abun_a)
#   abun_a$iid = rep(paste(1:(length(abun_a$id)/2)),2)
#
#   library(plyr)
#   abun_a1 = ddply(abun_a,"iid",transform,percount = count/sum(count)*100)
#   head(abun_a1)
#   mi = c( "firebrick3","navy")
#
#   head(abun_a1)
#   head(data)
#   tem = data.frame(org.id = data$org.id,rf.id = data$id)
#   abun_a1$ID = factor(abun_a1$ID,levels = data$org.id[length(data$org.id):1])
#   abun3 = abun_a1 %>% left_join(tem,by  = c("ID" = "org.id"))
#   head(abun3)
#
#   p2 = abun3  %>%
#     ggplot(aes(y = rf.id, x = percount, fill =id, group =id)) +
#     geom_bar(stat = 'identity') +
#     scale_fill_manual(values = mi)+
#     # theme_classic()+
#     theme_classic()+
#     theme(axis.text.x=element_text(angle=90,vjust=0.5, hjust=1,size = 4))
#   p2
#
#   plotname <- paste(matpath,"/a4_random_forst_loaing_abun",optimal,".pdf",sep = "")
#   ggsave(plotname, p, width = 8, height = 4)
#
#
#   library(aplot)
#
#
#
#   p3 = p1  %>%
#     aplot::insert_right(p2, width=.5)
#
#   plotname <- paste(matpath,"/a4_random_forst_loaing_abun.important",optimal,".pdf",sep = "")
#   ggsave(plotname, p3, width = 15, height = 10)
#
#
#
#   otu_table = as.data.frame(t(vegan_otu(pst)))
#   head(otu_table)
#
#   design = as.data.frame(sample_data(pst))
#   ## 计算相对丰度，计算每个物种丰度均值，按照均值排序
#   OTU = as.matrix(otu_table)
#   norm = t(t(OTU)/colSums(OTU,na=TRUE)) #* 100 # normalization to total 100
#   norma = norm %>%
#     t() %>% as.data.frame()
#   #数据分组计算平均值
#   iris.split <- split(norma,as.factor(design$Group))
#
#   iris.apply <- lapply(iris.split,function(x)colMeans(x,na.rm = TRUE))
#   # 组合结果
#   norm2 <- do.call(rbind,iris.apply)%>%
#     t()
#   norm2 = as.data.frame(norm2)
#   norm2$mean=apply(norm2,1,mean)
#   norm2$ID = row.names(norm2)
#   colnames(norm2)
#   abun = norm2
#   head(abun)
#   colnames(abun) = paste0(colnames(abun),"abundance.tab")
#   data = result[[7]]
#   head(data)
#   abun2 = abun %>% left_join(data,by = c("IDabundance.tab"="org.id")) %>% dplyr::arrange(desc(MeanDecreaseAccuracy))
#   filename = paste(matpath,"/randomforest_data.r.abundance.order.csv",sep = "")
#   write.csv(abun2,filename,quote = F)
# }
#





