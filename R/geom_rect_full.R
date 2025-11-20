#' Rectangle path with two rounded corners (A-P1-P2-B)
#'
#' @description
#' Draw a rectangular polyline from A(x,y) to B(xend,yend),
#' going through the two opposite vertices P1 and P2 so that
#' A, P1, P2, B form a rectangle. At P1 and P2, replace sharp
#' right angles with circular/elliptical fillets.
#'
#' @param mapping,data,... ggplot2 arguments. Required: x,y,xend,yend.
#' @param height,width Numeric. Offsets for the rectangle side not defined by A,B.
#'   If A,B share y (horizontal edge), use `height` (rectangle height).
#'   If A,B share x (vertical edge), use `width` (rectangle width).
#' @param rx1,ry1 Fillet radii at P1.
#' @param rx2,ry2 Fillet radii at P2.
#' @param n Points per fillet arc.
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x=0,y=0,xend=4,yend=0)
#'
#' ggplot(df) +
#'   geom_rect_full(aes(x,y,xend=xend,yend=yend),
#'                  height=3,
#'                  rx1=0.4,ry1=0.4,
#'                  rx2=0.4,ry2=0.4,
#'                  linewidth=1) +
#'   coord_equal() + theme_minimal()
#'
#' @export
geom_rect_full <- function(mapping=NULL, data=NULL, ...,
                           height=NULL, width=NULL,
                           rx1=0.3, ry1=0.3,
                           rx2=0.3, ry2=0.3,
                           n=30,
                           na.rm=FALSE, show.legend=NA, inherit.aes=TRUE) {
  ggplot2::layer(
    stat=StatRectFull, data=data, mapping=mapping,
    geom=ggplot2::GeomPath, position="identity",
    show.legend=show.legend, inherit.aes=inherit.aes,
    params=list(height=height,width=width,
                rx1=rx1,ry1=ry1,rx2=rx2,ry2=ry2,
                n=n,na.rm=na.rm,...)
  )
}

StatRectFull <- ggplot2::ggproto(
  "StatRectFull", ggplot2::Stat,
  required_aes=c("x","y","xend","yend"),
  compute_panel=function(data, scales,
                         height=NULL,width=NULL,
                         rx1=0.3, ry1=0.3, rx2=0.3, ry2=0.3,
                         n=30){

    fillet_arc <- function(Pin,P,Pout,rx,ry,n){
      u1 <- Pin-P; L1 <- sqrt(sum(u1^2)); if(L1==0) return(NULL); u1<-u1/L1
      u2 <- Pout-P; L2 <- sqrt(sum(u2^2)); if(L2==0) return(NULL); u2<-u2/L2
      rx <- min(rx,0.9*L1); ry <- min(ry,0.9*L2)
      Ein <- P+u1*rx; Eout<-P+u2*ry; C<-P+u1*rx+u2*ry
      th<-seq(0,pi/2,length.out=n)
      arc<-t(sapply(th,function(t) C-rx*u1*sin(t)-ry*u2*cos(t)))
      list(xs=c(Ein[1],arc[,1],Eout[1]), ys=c(Ein[2],arc[,2],Eout[2]))
    }

    out <- lapply(seq_len(nrow(data)), function(i){
      x1<-data$x[i]; y1<-data$y[i]
      x2<-data$xend[i]; y2<-data$yend[i]

      # 判断是水平边还是垂直边
      if(abs(y1-y2) < 1e-8){
        # AB 水平，矩形高度给定
        h <- if(is.null(height)) 1 else height
        P1 <- c(x1, y1+h)
        P2 <- c(x2, y2+h)
      } else if(abs(x1-x2) < 1e-8){
        # AB 垂直，矩形宽度给定
        w <- if(is.null(width)) 1 else width
        P1 <- c(x1+w, y1)
        P2 <- c(x2+w, y2)
      } else {
        stop("A and B must be aligned horizontally or vertically")
      }

      f1<-fillet_arc(c(x1,y1),P1,P2,rx1,ry1,n)
      f2<-fillet_arc(P1,P2,c(x2,y2),rx2,ry2,n)

      xs<-c(x1,f1$xs,f2$xs,x2)
      ys<-c(y1,f1$ys,f2$ys,y2)
      df<-data.frame(x=xs,y=ys,PANEL=data$PANEL[i],group=i)
      df
    })
    do.call(rbind,out)
  }
)
