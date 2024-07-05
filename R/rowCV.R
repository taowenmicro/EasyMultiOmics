rowSD = function(x){
  apply(x,1, sd)
}


rowCV = function(x){
  rowSD(x)/rowMeans(x)
}
