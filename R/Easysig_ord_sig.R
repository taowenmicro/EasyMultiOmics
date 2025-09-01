#' Reverse the significance level of letter marks
#'
#' @param data a data.frame contain the input data
#' @param  ID col index wtich need to reverse significance level,if NULL was marked, would transform to each col.
#' @examples
#' # data(data_wt)
#' result = KwWlx(data = data_wt, i= 4)
#' datord <- ord_sig(data = result[[1]],ID = NULL)
#' @return data.frame
#' @author Contact: Tao Wen \email{2018203048@@njau.edu.cn} Jun Yuan \email{junyuan@@njau.edu.cn}
#' @references
#'
#' Yuan J, Zhao J, Wen T, Zhao M, Li R, Goossens P, Huang Q, Bai Y, Vivanco JM, Kowalchuk GA, Berendsen RL, Shen Q
#' Root exudates drive the soil-borne legacy of aboveground pathogen infection
#' Microbiome 2018,DOI: \url{doi: 10.1186/s40168-018-0537-x}
#' @export



ord_sigold <- function(data = da,ID = "groups"){
 if (!is.null(ID)) {
   data[,ID] = gsub(" ","",as.character(data[,ID]))
   aa <- unique(as.character(data[,ID] ))

   # sort(aa,decreasing =TRUE)
   # sort(aa,decreasing =FALSE)

   tmp <- data.frame(ori = sort(aa,decreasing =TRUE),new = sort(aa,decreasing =FALSE))

   for (i in 1:nrow(data)) {
     data[,ID][i] <- as.character(tmp$new[match(data[,ID][i],tmp$ori)])

   }

 }

   if (is.null(ID)) {
     for (i in 1:ncol(data)) {
       aa <- unique(as.character(data[,i] ))
       tmp <- data.frame(ori = sort(aa,decreasing =TRUE),new = sort(aa,decreasing =FALSE))
       for (ii in 1:nrow(data)) {
         data[,i][ii] <- as.character(tmp$new[match(data[,i][ii],tmp$ori)])
       }
     }
   }
  return(data)

}





