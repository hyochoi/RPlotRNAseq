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
heatmap2.hy = function(get.hclust,main="Heatmap",only.heatmap=FALSE,
                       col=NULL,breaks=NULL,color.quantile=NULL,color.key=TRUE,
                       col.dend=TRUE,row.dend=TRUE,
                       col.cutoff=NULL,col.numCluster=NULL,col.colorBranch=NULL,
                       row.cutoff=NULL,row.numCluster=NULL,row.colorBranch=NULL,
                       col.trait=NULL,row.trait=NULL,
                       col.label=NULL,row.label=NULL,
                       col.label.size=0.5,row.label.size=0.5,
                       col.trait.col=NULL,row.trait.col=NULL) {
  # Updated on 6/9/2020
  # column and row labels added
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
    if (!is.null(col.label)) {
      axis(1,at=seq(0,1,length.out=ncol(X_sorted)),labels=col.label[col_out_order],
           tck=-0.015,lwd=0,line=-0.5,cex.axis=col.label.size,las=2,cex=0.1,srt=45)
    }
    if (!is.null(row.label)) {
      axis(2,at=seq(0,1,length.out=nrow(X_sorted)),labels=rev(row.label[row_out_order]),
           tck=-0.015,lwd=0,line=-0.5,cex.axis=row.label.size,las=2,cex=0.1)
    }
  } else {
    if ((is.null(col.trait) & is.null(col.label)) & (is.null(row.trait) & is.null(row.label))) {
      m=matrix(c(4,2,3,1),ncol=2,byrow=TRUE)
      layout(m,widths=c(0.2,0.8),heights = c(0.2,0.8))
    } else if ((!(is.null(col.trait) & is.null(col.label))) & (is.null(row.trait) & is.null(row.label))) {
      m=matrix(c(4,2,3,1,6,5),ncol=2,byrow=TRUE)
      layout(m,widths=c(0.2,0.8),heights=c(0.15,0.8,0.05))
    } else if ((is.null(col.trait) & is.null(col.label)) & (!(is.null(row.trait) & is.null(row.label)))) {
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
    if (!is.null(col.label)) {
      axis(1,at=seq(0,1,length.out=ncol(X_sorted)),labels=col.label[col_out_order],
           tck=-0.015,lwd=0,line=-0.5,cex.axis=col.label.size,las=2,cex=0.1,srt=45)
    }
    if (!is.null(row.label)) {
      axis(2,at=seq(0,1,length.out=nrow(X_sorted)),labels=rev(row.label[row_out_order]),
           tck=-0.015,lwd=0,line=-0.5,cex.axis=row.label.size,las=2,cex=0.1)
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
