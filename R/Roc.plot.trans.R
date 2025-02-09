#' @title Comparison of three machine methods (randomforest,SVM,GLM)
#' @description Comparison of three machine methods (randomforest,SVM,GLM).
#' @param otu transcriptome functional composition table.
#' @param map Sample metadata.
#' @param tax taxonomy table.
#' @param ps A phyloseq format file used as an alternative for the input containing transcriptome functional composition table, tax, and sample metadata.
#' @param Group Column name for groupID in map table.
#' @param repnum Modeling times.
#' @return list contain ggplot object and table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn},  Peng-Hao Xie \email{2019103106@njau.edu.cn}
#' @examples
#' data(ps.trans)
#' ps=ps.trans%>%filter_OTU_ps(Top=100)
#' id = sample_data(ps)$Group %>% unique()
#' aaa = combn(id,2)
#' i= 1
#' group = c(aaa[1,i],aaa[2,i])
#' pst = ps %>% subset_samples.wt("Group",group) %>%filter_taxa(function(x) sum(x ) > 10, TRUE)
#' res = Roc.plot.trans( ps = pst,group  = "Group",repnum = 5)
#' p33.1 =  res[[1]]
#' p33.1
#' p33.2 =  res[[2]]
#' p33.2
#' dat =  res[[3]]
#' dat
#'@export
Roc.plot.trans <- function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
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
    #print(model)
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
  print("finsh rf")
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
  print("finsh svm")
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

  print("finsh glm")
  return(list(g2,AUC,df,list(rf.rocdata=rf.rocdata,svm.rocdata =svm.rocdata,glm.rocdata = glm.rocdata )))
}

