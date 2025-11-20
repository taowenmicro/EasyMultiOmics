#' Offset-rectangle path with two rounded corners (A-P1-P2-B)
#'
#' @description
#' Connect A(x,y) and B(xend,yend) which are located on opposite rectangle sides
#' (not necessarily at corners, but somewhere along the sides). The path goes
#' from A up (or down) to side edge P1, across to opposite edge P2, and then
#' down (or up) to B. At P1 and P2 the 90° corners are replaced by quarter-circle/
#' ellipse arcs.
#'
#' @param mapping,data,... ggplot2 args. Required: x,y,xend,yend.
#' @param side "top" or "bottom" if A,B are on left/right sides; "left"/"right"
#'   if A,B are on top/bottom sides. This chooses which adjacent side contains
#'   the imaginary points P1,P2.
#' @param rx1,ry1 Radii for fillet at P1.
#' @param rx2,ry2 Radii for fillet at P2.
#' @param n Number of points per fillet arc.
#'
#' @examples
#' library(ggplot2)
#' df <- data.frame(x=0,y=1.2,xend=4,yend=2.3)
#'
#' # A在左边，B在右边，经过顶边
#' ggplot(df) +
#'   geom_offset_rect(aes(x,y,xend=xend,yend=yend),
#'                    side="top", rx1=0.4, ry1=0.4, rx2=0.4, ry2=0.4,
#'                    linewidth=1) +
#'   coord_equal() + theme_minimal()
#'
#' # A在下边，B在上边，经过右边
#' df2 <- data.frame(x=1,y=0,xend=3,yend=3)
#' ggplot(df2) +
#'   geom_offset_rect(aes(x,y,xend=xend,yend=yend),
#'                    side="right", rx1=0.3, ry1=0.3, rx2=0.3, ry2=0.3) +
#'   coord_equal() + theme_minimal()
#'
#' @export
geom_offset_rect <- function(mapping=NULL, data=NULL, ...,
                             side=c("top","bottom","left","right"),
                             rx1=0.3, ry1=0.3,
                             rx2=0.3, ry2=0.3,
                             n=30,
                             na.rm=FALSE, show.legend=NA, inherit.aes=TRUE) {
  side <- match.arg(side)
  ggplot2::layer(
    stat=StatOffsetRect, data=data, mapping=mapping,
    geom=ggplot2::GeomPath, position="identity",
    show.legend=show.legend, inherit.aes=inherit.aes,
    params=list(side=side,rx1=rx1,ry1=ry1,rx2=rx2,ry2=ry2,n=n,na.rm=na.rm,...)
  )
}

StatOffsetRect <- ggplot2::ggproto(
  "StatOffsetRect", ggplot2::Stat,
  required_aes=c("x","y","xend","yend"),
  compute_panel=function(data, scales,
                         side="top", rx1=0.3,ry1=0.3,rx2=0.3,ry2=0.3,n=30){
    fillet_arc <- function(Pin,P,Pout,rx,ry,n){
      u1<-Pin-P; L1<-sqrt(sum(u1^2)); if(L1==0) return(NULL); u1<-u1/L1
      u2<-Pout-P; L2<-sqrt(sum(u2^2)); if(L2==0) return(NULL); u2<-u2/L2
      rx<-min(rx,0.9*L1); ry<-min(ry,0.9*L2)
      Ein<-P+u1*rx; Eout<-P+u2*ry; C<-P+u1*rx+u2*ry
      th<-seq(0,pi/2,length.out=n)
      arc<-t(sapply(th,function(t) C-rx*u1*sin(t)-ry*u2*cos(t)))
      list(xs=c(Ein[1],arc[,1],Eout[1]), ys=c(Ein[2],arc[,2],Eout[2]))
    }

    out <- lapply(seq_len(nrow(data)), function(i){
      x1<-data$x[i]; y1<-data$y[i]
      x2<-data$xend[i]; y2<-data$yend[i]

      if(side %in% c("top","bottom")){
        y_side <- if(side=="top") max(y1,y2)+abs(y2-y1)+1e-6 else min(y1,y2)-abs(y2-y1)-1e-6
        P1<-c(x1,y_side); P2<-c(x2,y_side)
      } else {
        x_side <- if(side=="right") max(x1,x2)+abs(x2-x1)+1e-6 else min(x1,x2)-abs(x2-x1)-1e-6
        P1<-c(x_side,y1); P2<-c(x_side,y2)
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
