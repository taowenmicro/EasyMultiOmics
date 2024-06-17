
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
#'
Roc.micro <- function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
                     ps = NULL,group  = "Group",repnum = 5){


  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  mapping = as.data.frame(phyloseq::sample_data(ps))
  otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
  colnames(otutab) <- gsub("-","_",colnames(otutab))
  colnames(otutab) <- gsub("-","_",colnames(otutab))
  colnames(otutab) <- gsub("[/]","_",colnames(otutab))
  colnames(otutab) <- gsub("[(]","_",colnames(otutab))
  colnames(otutab) <- gsub("[)]","_",colnames(otutab))
  colnames(otutab) <- gsub("[:]","_",colnames(otutab))
  colnames(otutab) <- gsub("[[]","",colnames(otutab))
  colnames(otutab) <- gsub("[]]","_",colnames(otutab))
  colnames(otutab) <- gsub("[#]","_",colnames(otutab))
  colnames(otutab) <- gsub("[+]","_",colnames(otutab))
  colnames(otutab) <- gsub(" ","_",colnames(otutab))
  colnames(otutab) <- gsub("[;]","_",colnames(otutab))
  colnames(otutab) <- gsub("[,]","_",colnames(otutab))
  colnames(otutab) <- gsub("[?]","_",colnames(otutab))
  test = as.data.frame((otutab))
  mapping$Group = as.factor(mapping$Group)
  test$group = factor(mapping$Group)
  colnames(test) = paste("OTU",colnames(test),sep = "")

  # random forest

  test = dplyr::select(test,OTUgroup,everything())
  train = test
  folds <- createFolds(y=test[,1],k=repnum)
  AUC =c()

  max=0
  num=0
  fc<-as.numeric()#fc 为测试数据真实分组信息
  mod_pre<-as.numeric()
  i= 1
  for(i in 1:repnum){
    fold_test<-train[folds[[i]],]
    fold_train<-train[-folds[[i]],]


    colnames(fold_test) <- gsub("-","_",colnames(fold_test))
    colnames(fold_train) <- gsub("-","_",colnames(fold_train))

    model<-randomForest(OTUgroup~.,data=fold_train, importance=TRUE, proximity=TRUE)
    print(model)
    model_pre<-predict(model,newdata = fold_test,type="prob")
    fc<-append(fc,as.factor(fold_test$OTUgroup))
    mod_pre<-append(mod_pre,model_pre[,2])
  }


  #- pick data and plot
  pred <- prediction(mod_pre, fc)
  perf <- performance(pred,"tpr","fpr")
  x <- unlist(perf@x.values)  ##提取x值
  y <- unlist(perf@y.values)
  plotdata <- data.frame(x,y)
  names(plotdata) <- c("x", "y")
  AUC[1] = paste("rf AUC:",round(performance(pred,'auc')@y.values[[1]],3),sep = " ")
  head(plotdata)
  rf.rocdata = plotdata
  g0 <- ggplot(plotdata) +
    geom_path(aes(x = x, y = y, colour = x), size=1,color = "red") +
    labs(x = "False positive rate", y = "Ture positive rate") + # , title ="Random Forest"
    annotate("text", x=0.75, y=0.5, label=paste("Red: ",AUC[1],sep = ""))
  df<-cbind(fc,as.numeric(mod_pre))
  #-svm
  max=0
  num=0
  fc<-as.numeric()
  mod_pre<-as.numeric()

  for(i in 1:repnum){
    fold_test<-train[folds[[i]],]
    # head(fold_test)
    fold_train<-train[-folds[[i]],]
    model<-svm(OTUgroup~.,data=fold_train,probability=TRUE)
    model
    model_pre<-predict(model,newdata = fold_test,decision.values = TRUE, probability = TRUE)
    fc<-append(fc,(fold_test$OTUgroup))
    mod_pre<-append(mod_pre,as.numeric(attr(model_pre, "probabilities")[,2]))

  }

  pred <- prediction(mod_pre, fc)
  perf <- performance(pred,"tpr","fpr")
  x <- unlist(perf@x.values)
  y <- unlist(perf@y.values)
  plotdata <- data.frame(x,y)
  names(plotdata) <- c("x", "y")
  svm.rocdata = plotdata
  AUC[2] = paste("svm AUC:",round(performance(pred,'auc')@y.values[[1]],3),sep = " ")


  g1 <- g0 +
    geom_path(data = plotdata,aes(x = x, y = y, colour = x), size=1,color = "blue") +
    # labs(x = "False positive rate", y = "Ture positive rate") +
    annotate("text", x=0.75, y=0.4, label=paste("Blue: ",AUC[2],sep = ""))

  df<-cbind(df,cbind(fc,mod_pre))

  #GLM
  max=0
  num=0
  fc<-as.numeric()
  mod_pre<-as.numeric()
  for(i in 1:repnum){
    fold_test<-train[folds[[i]],]
    fold_train<-train[-folds[[i]],]
    model<-glm(OTUgroup~.,family='binomial',data=fold_train)
    model
    model_pre<-predict(model,type='response',newdata=fold_test)
    model_pre

    fc<-append(fc,fold_test$OTUgroup)
    mod_pre<-append(mod_pre,as.numeric(model_pre))
  }

  pred <- prediction(mod_pre, fc)
  perf <- performance(pred,"tpr","fpr")
  x <- unlist(perf@x.values)  ##提取x值
  y <- unlist(perf@y.values)
  plotdata <- data.frame(x,y)
  names(plotdata) <- c("x", "y")
  glm.rocdata = plotdata
  AUC[3] = paste("GLM AUC:",round(performance(pred,'auc')@y.values[[1]],3),sep = " ")

  g2 <- g1 +
    geom_path(data = plotdata,aes(x = x, y = y, colour = x), size=1,color = "black") +
    labs(x = "False positive rate", y = "Ture positive rate") +
    annotate("text", x=0.75, y=0.3, label=paste("Black: ",AUC[3],sep = ""))

  g2

  df<-cbind(df,cbind(fc,mod_pre))


  return(list(g2,AUC,df,list(rf.rocdata=rf.rocdata,svm.rocdata =svm.rocdata,glm.rocdata = glm.rocdata )))
}

