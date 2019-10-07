#' @export
plot.hy = function(x,y,indlist=NULL,text=F,
                   colmat=NULL,indcol="red",
                   cex=1,indcex=1.2,...) {
  n=length(x)
  if (length(indlist)==0) {indlist=NULL}
  if (is.null(indlist)) {
    if (is.null(colmat)) {
      colmat=rep("black",length(x))
    }
    plot(x=x,y=y,col=colmat,...)
  } else {
    if (is.null(colmat)) {
      colmat=rep("grey",n)
      colmat[indlist]=rep(indcol,length(indlist))
    }
    plot(x[-indlist],y=y[-indlist],col=colmat[-indlist],cex=cex,xlim=yaxis.hy(x),ylim=yaxis.hy(y),...)
    if (!text) {
      points(x[indlist],y[indlist],col=colmat[indlist],cex=indcex,...)
    } else {
      text(x[indlist],y[indlist],indlist,col=colmat[indlist],cex=indcex,...)
    }
  }
}

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

#' @export
vec2text.hy=function(x,quote.mark=T){
  n = length(x);
  fr = x[1]; # first entry

  if (n>1) {
    if (is.numeric(fr)) {
      for (i in 2:n){
        fr=paste0(fr,",",x[i])
      }# for
    } else {
      for (i in 2:n){
        if (quote.mark) {
          fr=cat(fr,'","',x[i],sep="");
        } else {
          fr=cat(fr,',',x[i],sep="");
        }
      }# for
    }# if
  }# if
  return(fr);
}

#' @export
dendcluster.hy = function(get.hclust,row.col="row",cutoff=NULL,numCluster=NULL,
                          drawplot=TRUE,
                          title="dendrogram",colorBranch=NULL) {
  require(dendextend);

  ## The number of clusters issue is fixed
  if (row.col=="row") {
    hclust.out = get.hclust$row_out
  }
  if (row.col=="col") {
    hclust.out = get.hclust$col_out
  }
  dend = as.dendrogram(hclust.out);

  if (is.null(numCluster)) {
    if (is.null(h)) {
      cutoff = 10 # default
    }
    hc = cutree(tree=hclust.out,h=cutoff,k=NULL)
    dendcluster = t(t(hc))
    k = max(dendcluster)
    numCluster = k
  } else {
    hc = cutree(tree=hclust.out,k=numCluster,h=NULL)
    dendcluster = t(t(hc))
    k = max(hc)
    if (k==(numCluster+1)) {
      warning("The number of detected clusters is one more than numCluster.")
      # For some reason, this happens. I don't know why.
      # If this is the case, I set k as numCluster by setting numCluster <- numCluster-1
      numCluster = numCluster - 1
      hc = cutree(tree=hclust.out,k=numCluster,h=NULL)
      dendcluster = t(t(hc))
      k = max(hc)
    }
  }

  # Draw dendrogram
  if (drawplot) {
    if (is.null(colorBranch)) {
      colorBranch = brewer.pal(8, "Dark2")[1:k]
    }
    dend = color_branches(dend,k=numCluster,col=colorBranch)
    # o.dend = order.dendrogram(dend)
    # labels(dend) = dendcluster[o.dend]
    # labels_colors(dend) = as.integer(dendcluster[o.dend])
    plot(dend,leaflab="none")
    if (!is.null(numCluster)) {
      abline(h=cutoff,col="red",lwd=2)
    }
    title(title)
  }

  # Match clusters and dendrogram
  if (is.null(numCluster)) {
    cluster.order = cutree(hclust.out,h=cutoff,order_clusters_as_data=FALSE)
  } else {
    cluster.order = cutree(hclust.out,k=numCluster,order_clusters_as_data=FALSE)
  }

  return(list(cluster.data=dendcluster,cluster.dend=cluster.order))
}

#' @export
textplot.hy = function(textonplot="text", ...) {
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, textonplot, ...)
}

#' @export
palette.hy = function(set="normal") {
  require(RColorBrewer); require(wesanderson)
  if (set=="normal") {
    palette(c(brewer.pal(8,"Set1"),brewer.pal(8,"Dark2"),
              wes_palette(n=4,"GrandBudapest1"),wes_palette(n=5,"Darjeeling1"),
              wes_palette(n=5,"Cavalcanti1")))
  } else if (set=="group") {
    palette(c(brewer.pal(5,"Set2"),brewer.pal(7,"Dark2"),brewer.pal(8,"Set3")[c(1,3:7)],
              wes_palette(n=4,"GrandBudapest1"),wes_palette(n=5,"Darjeeling1"),
              wes_palette(n=5,"Cavalcanti1")))
  }
}

#' @export
get.hclust.hy = function(value,main="Heatmap",method="complete",
                         col.order=NULL,row.order=NULL,
                         col.order.mat=NULL,row.order.mat=NULL,
                         col.dend=TRUE,row.dend=TRUE) {

  require(gplots)
  X = value;
  col_out = row_out = NULL;
  # Get orders in the column and the row
  if (is.null(col.order)) {
    if (is.null(col.order.mat)) {
      col_out = hclust(dist(t(X)),method=method)
      col_out_order = col_out$order
    } else {
      col_out = hclust(dist(t(col.order.mat)),method=method)
      col_out_order = col_out$order
    }
  } else {
    col_out_order = col.order
  }

  if (is.null(row.order)) {
    if (is.null(row.order.mat)) {
      row_out = hclust(dist(X),method=method)
      row_out_order = rev(row_out$order)
    } else {
      row_out = hclust(dist(row.order.mat),method=method)
      row_out_order = rev(row_out$order)
    }
  } else {
    row_out_order = row.order
  }

  X_sorted = X[row_out_order,col_out_order]
  return(list(X_sorted=X_sorted,
              row_out_order=row_out_order,col_out_order=col_out_order,
              col_out=col_out,row_out=row_out))
}


#' @export
heatmap.hy = function(get.hclust,main="Heatmap",only.heatmap=FALSE,
                      col=NULL,breaks=NULL,color.quantile=NULL,color.key=TRUE,
                      col.dend=TRUE,row.dend=TRUE,
                      col.cutoff=NULL,col.numCluster=NULL,col.colorBranch=NULL,
                      row.cutoff=NULL,row.numCluster=NULL,row.colorBranch=NULL,
                      col.trait=NULL,row.trait=NULL,
                      col.trait.col=NULL,row.trait.col=NULL) {
  # Updated on 10/11/2017
  # Hyo Young Choi

  require(RColorBrewer); require(gplots);
  require(dendextend);

  X_sorted = get.hclust$X_sorted;
  row_out_order = get.hclust$row_out_order;
  col_out_order = get.hclust$col_out_order;
  row_out = get.hclust$row_out;
  col_out = get.hclust$col_out;

  if (is.null(col_out) & is.null(row_out)) {
    col.dend = FALSE; row.dend = FALSE;
  } else {
    if (is.null(col_out)) {
      col.dend = FALSE;
    } else {
      if (is.null(row_out)) {
        row.dend = FALSE;
      }
    }
  }

  if (is.null(col)) {
    if (is.null(color.quantile)) {
      col = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
    } else {
      col = rev(colorRampPalette(brewer.pal(9, "RdBu"))(length(color.quantile)-1))
    }
  }

  if (only.heatmap) {
    if (is.null(color.quantile)) {
      if (is.null(breaks)) {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,axes=FALSE)
      } else {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
      }
    } else {
      breaks = quantile(as.vector(X_sorted),probs=color.quantile)
      image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
    }
    if (!is.null(main)) {
      title(main)
    }
  } else {
    if ((is.null(col.trait)) & (is.null(row.trait))) {
      m=matrix(c(4,2,3,1),ncol=2,byrow=TRUE)
      layout(m,widths=c(0.2,0.8),heights = c(0.2,0.8))
    } else if ((!is.null(col.trait)) & (is.null(row.trait))) {
      m=matrix(c(4,2,3,1,6,5),ncol=2,byrow=TRUE)
      layout(m,widths=c(0.2,0.8),heights=c(0.15,0.8,0.05))
    } else if ((is.null(col.trait)) & (!is.null(row.trait))) {
      m=matrix(c(4,2,6,3,1,5),ncol=3,byrow=TRUE)
      layout(m,widths=c(0.15,0.8,0.05),heights=c(0.2,0.8))
    } else {
      m=matrix(c(4,2,9,3,1,5,7,6,8),ncol=3,byrow=TRUE)
      layout(m,widths=c(0.15,0.8,0.05),heights=c(0.15,0.8,0.05))
    }
    # layout.show(max(m))

    # Figure 1
    par(mar=c(0,0,0,0))
    if (is.null(color.quantile)) {
      if (is.null(breaks)) {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,axes=FALSE)
      } else {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
      }
    } else {
      breaks = quantile(as.vector(X_sorted),probs=color.quantile)
      image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
    }
    par(mar=c(5,4,4,1)+0.1)


    # Figure 2
    par(mar=c(0,0,2,0))
    if (col.dend) {
      coldend = as.dendrogram(col_out)

      if ((is.null(col.numCluster)) & (is.null(col.cutoff))) {
        plot(coldend,leaflab = "none",axes=FALSE,xaxs = "i")
      } else {
        if (is.null(col.numCluster)) {
          colcluster = t(t(cutree(col_out,h=col.cutoff)))
          col.numCluster = length(table(colcluster))
        }
        if (is.null(col.colorBranch)) {
          palette(rainbow(10))
          col.colorBranch = 1:col.numCluster
        }
        coldend = color_branches(coldend,k=col.numCluster,col=col.colorBranch)

        par(mar=c(0,0,2,0))
        plot(coldend,leaflab="none",axes=FALSE,xaxs="i")
        # abline(h=col.cutoff,lwd=2,col="red")
      }
    } else {
      plot.new()
    }
    title(main)
    par(mar=c(5,4,4,1)+0.1)

    # Figure 3
    par(mar=c(0,0,0,0))
    if (row.dend) {
      rowdend = as.dendrogram(row_out);
      if ((is.null(row.numCluster)) & (is.null(row.cutoff))) {
        plot(rowdend,leaflab = "none",axes=FALSE,yaxs = "i",horiz =TRUE)
      } else {
        if (is.null(row.numCluster)) {
          hc = cutree(rowdend,h=row.cutoff)
          k = max(hc)
          row.numCluster = k
        } else {
          hc = cutree(rowdend,k=row.numCluster,h=NULL)
          k = max(hc)
          if (k==(row.numCluster+1)) {
            warning("The number of detected row clusters is one more than numCluster.")
            # For some reason, this happens. I don't know why.
            # If this is the case, I set k as numCluster by setting numCluster <- numCluster-1
            row.numCluster = row.numCluster - 1
            hc = cutree(tree=rowdend,k=row.numCluster,h=NULL)
            k = max(hc)
          }
        }
        if (is.null(row.colorBranch)) {
          # row.colorBranch = brewer.pal(8, "Dark2")
          row.colorBranch = brewer.pal(8, "Dark2")[1:k]
        }

        rowdend = color_branches(rowdend,k=row.numCluster,col=row.colorBranch)
        par(mar=c(0,0,0,0))
        plot(rowdend,leaflab="none",axes=FALSE,yaxs="i",horiz=TRUE)
        # abline(h=col.cutoff,lwd=2,col="red")
      }
    } else {
      plot.new()
    }
    par(mar=c(5,4,4,1)+0.1)

    # Figure 4
    par(mar=c(0,0,0,0))
    if (color.key) {
      scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
      }
      par(mar=c(3,2,3,2),cex=0.75,mgp=c(2,1,0));
      min.raw = min(X_sorted); max.raw = max(X_sorted); interval=c(min.raw,max.raw)
      z = seq(min.raw,max.raw,length.out=1000);
      breaks = quantile(as.vector(X_sorted),probs=color.quantile)
      if (!is.null(color.quantile)) {
        image(matrix(z,ncol=1),col=col,breaks=breaks,xaxt = "n", yaxt = "n")
      } else {
        image(matrix(z,ncol=1),col=col,xaxt = "n", yaxt = "n")
      }

      par(usr = c(0, 1, 0, 1))
      lv <- pretty(c(min.raw,max.raw))
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(side=1,at=xv,labels=lv,cex=0.2,lwd=0,line=-0.5,tck=-0.015)
      mtext(side = 1, "Value", line=1.5,cex=0.8)

      # kernel density
      dens <- density(X_sorted, adjust = 0.25, na.rm = TRUE,
                      from = min.raw, to = max.raw)
      omit <- dens$x < min(interval) | dens$x > max(interval)
      dens$x <- dens$x[!omit]
      dens$y <- dens$y[!omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = "cyan",
            lwd = 1)
      axis(side=2,at=pretty(dens$y)/max(dens$y)*0.95,labels=pretty(dens$y),
           cex=0.2,lwd=0,line=-0.5,tck=-0.015)
      # title
      key.title <- "Color Key"
      title(key.title)
      par(mar=c(5,4,4,1)+0.1)
    } else {
      plot.new()
    }
    par(mar=c(5,4,4,1)+0.1)

    # Figure 5
    if (!is.null(row.trait)) {
      if (is.null(row.trait.col)) {
        row.trait.col = colorRampPalette(brewer.pal(9, "Blues"))(length(seq(0,1,0.01))-1)
        # row.trait.col = colorRampPalette(brewer.pal(9, "BuGn"))(100)
      }
      par(mar=c(0,0,0,0))
      breaks = quantile(as.vector(row.trait[row_out_order]),probs=seq(0,1,0.01))
      image(t(rev(row.trait[row_out_order])),col=row.trait.col,breaks=breaks,axes=FALSE)
      # image(t(matrix(row.trait[row_out_order],ncol=1)),
      #       col=row.trait.col,axes=FALSE)
      par(mar=c(5,4,4,1)+0.1)
    }

    if (!is.null(col.trait)) {
      if (is.null(col.trait.col)) {
        col.trait.col = colorRampPalette(brewer.pal(9, "OrRd"))(100)
      }
      par(mar=c(0,0,0,0))
      image(matrix(col.trait[col_out_order],ncol=1),
            col=col.trait.col,axes=FALSE)
    }
  }
  par(mfrow=c(1,1),mar=c(5,4,4,1)+0.1)
}

#' @export
heatmap3.hy = function(get.hclust,
                       rowlabels=NULL,rowlabel.size=0.5,
                       collabels=NULL,collabel.size=0.5,
                       main="Heatmap",only.heatmap=FALSE,
                       col=NULL,breaks=NULL,color.quantile=NULL,color.key=TRUE,
                       col.cutoff=NULL,col.numCluster=NULL,col.colorBranch=NULL,
                       col.trait=NULL,row.trait=NULL,
                       col.trait.col=NULL,row.trait.col=NULL) {

  # Instead of row dendrogram, display rownames
  # Created on 10/19/2017
  # Hyo Young Choi

  require(RColorBrewer); require(gplots)
  X_sorted = get.hclust$X_sorted;
  row_out_order = get.hclust$row_out_order;
  col_out_order = get.hclust$col_out_order;
  row_out = get.hclust$row_out;
  col_out = get.hclust$col_out;

  col.dend = TRUE; row.dend = TRUE;
  if (is.null(col_out) & is.null(row_out)) {
    col.dend = FALSE; row.dend = FALSE;
  } else {
    if (is.null(col_out)) {
      col.dend = FALSE;
    } else {
      if (is.null(row_out)) {
        row.dend = FALSE;
      }
    }
  }

  if (is.null(col)) {
    if (is.null(color.quantile)) {
      col = rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))
    } else {
      col = rev(colorRampPalette(brewer.pal(9, "RdBu"))(length(color.quantile)-1))
    }
  }

  if (only.heatmap) {
    if (is.null(color.quantile)) {
      if (is.null(breaks)) {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,axes=FALSE)
      } else {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
      }
    } else {
      breaks = quantile(as.vector(X_sorted),probs=color.quantile)
      image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
    }
    if (!is.null(main)) {
      title(main)
    }
    axis(1,at=seq(0,1,length.out=ncol(X_sorted)),labels=collabels[col_out_order],
         tck=-0.015,lwd=0,line=-0.5,cex.axis=collabel.size,las=2,cex=0.1)
    axis(2,at=seq(0,1,length.out=nrow(X_sorted)),labels=rev(rowlabels[row_out_order]),
         tck=-0.015,lwd=0,line=-0.5,cex.axis=rowlabel.size,las=2,cex=0.1)
  } else {
    if ((is.null(col.trait)) & (is.null(row.trait))) {
      m=matrix(c(4,2,3,1),ncol=2,byrow=TRUE)
      layout(m,widths=c(0.1,0.9),heights = c(0.15,0.85))
    } else if ((!is.null(col.trait)) & (is.null(row.trait))) {
      m=matrix(c(4,2,3,1,6,5),ncol=2,byrow=TRUE)
      layout(m,widths=c(0.1,0.9),heights=c(0.15,0.82,0.03))
    } else if ((is.null(col.trait)) & (!is.null(row.trait))) {
      m=matrix(c(4,2,6,3,1,5),ncol=3,byrow=TRUE)
      layout(m,widths=c(0.1,0.85,0.05),heights=c(0.15,0.85))
    } else {
      m=matrix(c(4,2,9,3,1,5,7,6,8),ncol=3,byrow=TRUE)
      layout(m,widths=c(0.1,0.85,0.05),heights=c(0.15,0.82,0.03))
    }
    # layout.show(max(m))

    # Figure 1
    par(mar=c(0,0,0,0))
    if (is.null(color.quantile)) {
      if (is.null(breaks)) {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,axes=FALSE)
      } else {
        image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
      }
    } else {
      breaks = quantile(as.vector(X_sorted),probs=color.quantile)
      image(t(X_sorted[seq(nrow(X_sorted),1),]),col=col,breaks=breaks,axes=FALSE)
    }
    axis(1,at=seq(0,1,length.out=ncol(X_sorted)),labels=collabels[col_out_order],
         tck=-0.015,lwd=0,line=-0.5,cex.axis=collabel.size,las=2,cex=0.1)
    axis(2,at=seq(0,1,length.out=nrow(X_sorted)),labels=rev(rowlabels[row_out_order]),
         tck=-0.015,lwd=0,line=-0.5,cex.axis=rowlabel.size,las=2,cex=0.1)
    par(mar=c(5,4,4,1)+0.1)


    # Figure 2
    library(dendextend); library(RColorBrewer);
    par(mar=c(0,0,2,0))
    if (col.dend) {
      coldend = as.dendrogram(col_out)

      if ((is.null(col.numCluster)) & (is.null(col.cutoff))) {
        plot(coldend,leaflab = "none",axes=FALSE,xaxs = "i")
      } else {
        if (is.null(col.numCluster)) {
          colcluster = t(t(cutree(col_out,h=col.cutoff)))
          col.numCluster = length(table(colcluster))
        }
        if (is.null(col.colorBranch)) {
          palette(rainbow(10))
          col.colorBranch = 1:col.numCluster
        }
        coldend = color_branches(coldend,k=col.numCluster,col=col.colorBranch)

        par(mar=c(0,0,2,0))
        plot(coldend,leaflab="none",axes=FALSE,xaxs="i")
        # abline(h=col.cutoff,lwd=2,col="red")
      }
    } else {
      plot.new()
    }
    title(main)
    par(mar=c(5,4,4,1)+0.1)

    # Figure 3
    par(mar=c(0,0,0,0))
    plot.new()
    par(mar=c(5,4,4,1)+0.1)

    # Figure 4
    if (color.key) {
      scale01 <- function(x, low = min(x), high = max(x)) {
        x <- (x - low)/(high - low)
        x
      }
      # par(mar=c(3,2,3,2),cex=0.75,mgp=c(2,1,0),usr=c(0,1,0,1));
      # par(mar=c(2,1,0,0),cex=0.75,mgp=c(2,1,0),usr=c(0,1,0,1));
      par(mar=c(2,2,1,0))
      min.raw = min(X_sorted); max.raw = max(X_sorted); interval=c(min.raw,max.raw)
      z = seq(min.raw,max.raw,length.out=1000);
      if (!is.null(color.quantile)) {
        breaks = quantile(as.vector(X_sorted),probs=color.quantile)
        image(matrix(z,ncol=1),col=col,breaks=breaks,xaxt = "n", yaxt = "n")
      } else {
        image(matrix(z,ncol=1),col=col,xaxt = "n", yaxt = "n")
      }

      par(usr = c(0, 1, 0, 1))
      lv <- pretty(c(min.raw,max.raw))
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      axis(side=1,at=xv,labels=lv,cex.axis=0.8,lwd=0,line=-0.5,tck=-0.015)
      # mtext(side = 1, "Value", line=1.5,cex=0.8)

      # kernel density
      dens <- density(X_sorted, adjust = 0.25, na.rm = TRUE,
                      from = min.raw, to = max.raw)
      omit <- dens$x < min(interval) | dens$x > max(interval)
      dens$x <- dens$x[!omit]
      dens$y <- dens$y[!omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y/max(dens$y) * 0.95, col = "cyan",
            lwd = 1)
      axis(side=2,at=pretty(dens$y)/max(dens$y)*0.95,labels=pretty(dens$y),
           cex.axis=0.8,lwd=0,line=-0.5,tck=-0.015)
      # title
      # key.title <- "Color Key"
      # title(key.title)
    } else {
      par(mar=c(0,0,0,0))
      plot.new()
    }
    par(mar=c(5,4,4,1)+0.1)

    # Figure 5
    if (!is.null(row.trait)) {
      if (is.null(row.trait.col)) {
        row.trait.col = colorRampPalette(brewer.pal(9, "Blues"))(length(seq(0,1,0.01))-1)
        # row.trait.col = colorRampPalette(brewer.pal(9, "BuGn"))(100)
      }
      par(mar=c(0,0,0,0))
      breaks = quantile(as.vector(row.trait[row_out_order]),probs=seq(0,1,0.01))
      image(t(rev(row.trait[row_out_order])),col=row.trait.col,breaks=breaks,axes=FALSE)
      # image(t(matrix(row.trait[row_out_order],ncol=1)),
      #       col=row.trait.col,axes=FALSE)
      par(mar=c(5,4,4,1)+0.1)
    }

    if (!is.null(col.trait)) {
      if (is.null(col.trait.col)) {
        col.trait.col = colorRampPalette(brewer.pal(9, "OrRd"))(100)
      }
      par(mar=c(0,0,0,0))
      image(matrix(col.trait[col_out_order],ncol=1),
            col=col.trait.col,axes=FALSE)
    }
    par(mfrow=c(1,1),mar=c(5,4,4,1)+0.1)
  }
}


#' @export
split.char <- function(barcode){
  tmp.split <- strsplit(barcode,"-")[[1]] ;
  return(paste0(tmp.split[1],"-",tmp.split[2],"-",tmp.split[3],"-",tmp.split[4]))
}

#' @export
pca.hy <- function(data, subt.mean=TRUE){
  ##  PCA with the covariance matrix (X%*%t(X)/n-1) where X=data
  ##  data: d by r (d is the number of variables, n is the number of subjects)
  ##  If subt.mean=TRUE, subtract mean from the data.
  if (subt.mean){
    x <- data-apply(data,1,mean) ;
    div.fac <- ncol(x)-1;
  } else {
    x <- data ;
    div.fac <- ncol(x);
  }

  x.mda <- svd(x) ;
  dirmat <- x.mda$u ;
  if (dirmat[1,1]<0){  dirmat[,1] <- -dirmat[,1]}
  projmat <- t(dirmat)%*%x ;
  eigenval <- ((x.mda$d)^2)/div.fac ;
  stdprojmat <- diag(1/sqrt(eigenval))%*%projmat		#  standardized
  return(list(dirmat=dirmat, projmat=projmat, eigenval=eigenval, stdprojmat=stdprojmat))
}

#' @export
yaxis.hy <- function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}


#' @export
split.exon = function(exon) {
  a <- exon ;
  b <- strsplit(a,":")[[1]][2] ;
  c <- gsub("-",",",b) ;
  d <- as.numeric(strsplit(c,",")[[1]]) ;
  epl <- d[seq(1,length(d),by=2)] ;
  epr <- d[seq(2,length(d),by=2)] ;
  ep <- cbind(epl,epr) ;

  seqdir = strsplit(exon,":")[[1]][3] ;
  if (seqdir=="+") {
    return(ep)
  } else {
    return(ep[c(nrow(ep):1),c(2:1)])
  }
}

