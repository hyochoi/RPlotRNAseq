
#' @export
boxplot.hy=function(value,robust=FALSE,
                    title=NULL,title.cex=2,
                    ylab=NULL,
                    value.labels=NULL,tick.at=NULL,
                    ylim=NULL,
                    indlist=NULL,indlist.col=NULL,
                    indcex=2,text=FALSE,textwhat=NULL,
                    print.outliers=FALSE,outliers.range=NULL,outliers.text=FALSE,
                    colmat=NULL,outliers.col=NULL,
                    box.col="black",print.plot=TRUE) {
  require(wesanderson); require(RColorBrewer)
  n = length(value)
  qvalue = quantile(value)
  iqr = qvalue[4]-qvalue[2]
  if (!robust) {
    if (is.null(outliers.range)) {
      fence=c(qvalue[2]-iqr*1.5,qvalue[4]+iqr*1.5)
    } else {
      fence=c(qvalue[2]-iqr*outliers.range,qvalue[4]+iqr*outliers.range)
    }
  } else {
    iqr.u=qvalue[4]-qvalue[3] # iqr for large values
    iqr.b=qvalue[3]-qvalue[2] # iqr for low values
    if (is.null(outliers.range)) {
      fence=c(qvalue[2]-2*iqr.b*1.5,qvalue[4]+2*iqr.u*1.5)
    } else {
      fence=c(qvalue[2]-2*iqr.b*outliers.range,qvalue[4]+2*iqr.u*outliers.range)
    }
  }
  box.outliers.below=which(value<fence[1])
  box.outliers.above=which(value>fence[2])
  box.outliers = c(which(value<fence[1]),which(value>fence[2]))

  if (print.outliers) {
    # indlist=box.outliers
    if (outliers.text) {
      text=TRUE
    }
  }
  if (print.plot) {
    if ((!is.null(indlist)) & (length(indlist)==0)) {indlist=NULL}
    # Define colors for points
    if (is.null(colmat)) {
      colmat=rep("grey",n)
      if (!is.null(box.outliers)) {
        if (is.null(outliers.col)) {
          outliers.col = "red"
        }
        colmat[box.outliers] = outliers.col;
      }

      if (!is.null(indlist)) {
        if (is.null(indlist.col)) {
          indlist.col=wes_palette(n=5,"Darjeeling1")[2]
        }
        colmat[indlist] = indlist.col;
      }
    }

    # Define x values
    x = rnorm(n=length(value),mean=0.5,sd=0.03)

    # Define ylim if it is null
    if (is.null(ylim)) {
      ylim=yaxis.hy(value)
      if (ylim[1]==ylim[2]) {
        if (ylim[1]==0) {
          ylim[1]=-0.5; ylim[1]=0.5
        } else {
          ylim[1]=ylim[1]*0.5; ylim[2]=ylim[2]*1.5
        }
      }
    }
    if (is.null(indlist)) {
      plot(x,value,axes=F,col=colmat,ylim=ylim,
           xlab=NA,ylab=NA,xlim=c(min(x)-0.05,max(x)+0.05),pch=19,cex=1,lwd=3)
    } else {
      plot(x[-indlist],value[-indlist],axes=F,col=colmat[-indlist],ylim=ylim,
           xlab=NA,ylab=NA,xlim=c(min(x)-0.05,max(x)+0.05),pch=19,cex=1,lwd=3)
    }

    if (is.null(value.labels)) {
      tick.at=seq(ylim[1],ylim[2],length=5)
      if ((qvalue[5]-qvalue[1]>0) & (qvalue[5]<1e-2)) {
        digits=abs(round(log10(qvalue[5])))+1
        value.labels=round(tick.at,digits=digits)
      } else {
        value.labels=round(tick.at,digits=2)
      }
    } else {
      if (is.null(tick.at)) {
        cat("tick.at should be specified if value.labels are defined.","\n")
      }
    }
    axis(side=2,lwd=0.8,cex.axis=1,at=tick.at,labels=value.labels)
    abline(h=ylim[1])
    mtext(side=3, title, line=1, cex=title.cex)
    if (!is.null(ylab)) {
      mtext(side=2,ylab,line=3,cex=2)
    }

    if (!is.null(indlist)) {
      if (text) {
        if (is.null(textwhat)) {
          text(x[indlist], value[indlist], indlist, col = colmat[indlist],
               cex = indcex)
        } else {
          value[indlist] = seq(min(value[indlist]), max(value[indlist]),
                               length.out = length(indlist))
          text(x[indlist], value[indlist], textwhat, col = colmat[indlist],
               cex = indcex)
        }
      } else {
        points(x[indlist],value[indlist],col=colmat[indlist],pch=19,cex=indcex)
      }
    }
    segments(x0=min(x),x1=max(x),y0=qvalue[3],y1=qvalue[3],col=box.col,lwd=5)
    segments(x0=min(x),x1=max(x),y0=qvalue[2],y1=qvalue[2],col=box.col,lwd=2)
    segments(x0=min(x),x1=max(x),y0=qvalue[4],y1=qvalue[4],col=box.col,lwd=2)
    segments(x0=min(x),x1=min(x),y0=qvalue[2],y1=qvalue[4],col=box.col,lwd=2)
    segments(x0=max(x),x1=max(x),y0=qvalue[2],y1=qvalue[4],col=box.col,lwd=2)

    if (print.outliers) {
      seg.range=max(x)-min(x)
      if (fence[1]>=ylim[1]) {
        segments(x0=min(x)+0.2*seg.range,x1=max(x)-0.2*seg.range,
                 y0=fence[1],y1=fence[1],col=box.col,lwd=2)
      }
      segments(x0=0.5,x1=0.5,y0=qvalue[2],y1=max(ylim[1],fence[1]),col=box.col,lwd=2)
      if (fence[2]<=ylim[2]) {
        segments(x0=min(x)+0.2*seg.range,x1=max(x)-0.2*seg.range,
                 y0=fence[2],y1=fence[2],col=box.col,lwd=2)
      }
      segments(x0=0.5,x1=0.5,y0=qvalue[4],y1=min(ylim[2],fence[2]),col=box.col,lwd=2)
    }
  }
  if (print.plot) {
    return(list(outliers=box.outliers,
                outliers.above=box.outliers.above,
                outliers.below=box.outliers.below,
                iqr=iqr,fence=fence))
  } else {
    return(list(outliers=box.outliers,
                outliers.above=box.outliers.above,
                outliers.below=box.outliers.below,
                iqr=iqr,fence=fence))
  }
}

#' @export
boxplot.hy.v0=function(value,title=NULL,ylab=NULL,
                       value.labels=NULL,tick.at=NULL,
                       indlist=NULL,outliers=NULL,
                       colmat=NULL,outlier.col=NULL,
                       box.col="black") {
  require(wesanderson); require(RColorBrewer)
  n = length(value)

  # Define colors for points
  if (is.null(colmat)) {
    colmat=rep("grey",n)
    if (!is.null(outliers)) {
      if (is.null(outlier.col)) {
        outlier.col = "red"
      }
      colmat[outliers] = outlier.col;
    }
    if (!is.null(indlist)) {
      if (is.null(indlist.col)) {
        indlist.col=wes_palette(n=5,"Darjeeling")[2]
      }
      colmat[indlist] = indlist.col;
    }
  }

  x = rnorm(n=length(value),mean=0.5,sd=0.03)

  ylim=yaxis.hy(value)
  # par(mar=c(1,2.5,3,0))
  plot(x,value,axes=F,col=colmat,ylim=ylim,
       xlab=NA,ylab=NA,xlim=c(min(x)-0.05,max(x)+0.05),pch=19,cex=1,lwd=3)

  if (is.null(value.labels)) {
    tick.at=seq(ylim[1],ylim[2],length=5)
    value.labels=round(tick.at,digits=2)
  } else {
    if (is.null(tick.at)) {
      cat("tick.at should be specified if value.labels are defined.","\n")
    }
  }
  axis(side=2,lwd=0.8,cex.axis=1.5,at=tick.at,labels=value.labels)
  abline(h=ylim[1])
  mtext(side=3, title, line=1, cex=1)
  if (!is.null(ylab)) {
    mtext(side=2,ylab,line=3,cex=2)
  }

  qvalue = quantile(value)
  segments(x0=min(x),x1=max(x),y0=qvalue[3],y1=qvalue[3],col=box.col,lwd=5)
  segments(x0=min(x),x1=max(x),y0=qvalue[2],y1=qvalue[2],col=box.col,lwd=2)
  segments(x0=min(x),x1=max(x),y0=qvalue[4],y1=qvalue[4],col=box.col,lwd=2)
  segments(x0=min(x),x1=min(x),y0=qvalue[2],y1=qvalue[4],col=box.col,lwd=2)
  segments(x0=max(x),x1=max(x),y0=qvalue[2],y1=qvalue[4],col=box.col,lwd=2)
  points(x[indlist],value[indlist],col=colmat[indlist],pch=19,cex=2)

}
