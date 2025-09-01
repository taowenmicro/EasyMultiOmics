
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
#' @return list contain ggplot object and table.
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn},  Peng-Hao Xie \email{2019103106@njqu.edu.cn}
#' @examples
#' # data form github
#' result = rfcv.Micro(ps = ps_rela,group  = "Group",optimal = 20,nrfcvnum = 6)
#' prfcv = result[[1]]# plot rfcv
# result[[2]]# plotdata
#' rfcvtable = result[[3]]# table rfcv
#'@export
rfcv.Micro = function(otu = NULL,tax = NULL,map = NULL,tree = NULL,
                      ps = NULL,group  = "Group",optimal = 20,nrfcvnum = 5){


  ps = ggClusterNet::inputMicro(otu,tax,map,tree,ps,group  = group)
  otutab = as.data.frame((ggClusterNet::vegan_otu(ps)))
  mapping =  sample_data(ps)
  # Set classification info.
  otutab$group = factor(mapping$Group)
  colnames(otutab) <- gsub("-","_",colnames(otutab))
  # rfcv for select···
  n = ncol(otutab)-1
  myotutab_t= otutab[1:n]
  set.seed(315)
  result= rfcv(myotutab_t, otutab$group, cv.fold=5, scale = "log", step = 0.9)
  # with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))
  result1 = result
  error.cv = data.frame(num = result$n.var, error.1 =  result$error.cv)
  for (i in 316:(314+ nrfcvnum)){
    print(i)
    set.seed(i)
    result= rfcv(myotutab_t, otutab$group, cv.fold=5, scale = "log", step = 0.9)
    error.cv = cbind(error.cv, result$error.cv)
  }
  n.var = error.cv$num
  error.cv = error.cv[,2:6]
  colnames(error.cv) = paste('err',1:5,sep='.')
  err.mean = apply(error.cv,1,mean)
  allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
  head(allerr)
  data <- gather(allerr, key = "group", value = "value",-num)
  head(data)

  p <- ggplot() +
    geom_line(data = data,aes(x = num, y = value,group = group), colour = 'grey') +
    geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') +
    coord_trans(x = "log2") +
    scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
    labs(title=paste('Training set (n = ', dim(otutab)[1],')', sep = ''),
         x='Number of families ', y='Cross-validation error rate') +
    annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep=""))
  return(list(plot = p,plotdata = data,origdata = allerr))
}
