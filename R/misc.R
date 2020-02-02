#' @export
plot.hy = function(x,y,indlist=NULL,text=F,
                   colmat=NULL,indcol="red",
                   cex=1,indcex=1.2,xlim=NULL,ylim=NULL,
                   xlab=NULL,ylab=NULL,...) {
  n=length(x)
  if (length(indlist)==0) {indlist=NULL}
  if (is.null(xlim)) {
    xlim=yaxis.hy(x)
  }
  if (is.null(ylim)) {
    ylim=yaxis.hy(y)
  }
  if (is.null(indlist)) {
    if (is.null(colmat)) {
      colmat=rep("black",length(x))
    }
    plot(x=x,y=y,col=colmat,xlim=xlim,ylim=ylim,axes=F,ylab=NA,xlab=NA,cex=cex,...)
    box()
    abline(0, 0, lty = 2, col = rgb(0.5, 0.5, 0.5))
    abline(v = 0, lty = 2, col = rgb(0.5, 0.5, 0.5))
    axis(side = 1, tck = -0.015, labels = NA)
    axis(side = 1, lwd = 0, line = -1, cex = 0.2, cex.axis = 0.9)
    axis(side = 2, tck = -0.015, lwd = 0, line = -1, cex.axis = 0.9)
    mtext(side = 1, xlab, line = 1, cex = 1)
    mtext(side = 2, ylab, line = 1, cex = 1)
  } else {
    if (is.null(colmat)) {
      colmat=rep("grey",n)
      colmat[indlist]=rep(indcol,length(indlist))
    }
    plot(x[-indlist],y=y[-indlist],col=colmat[-indlist],cex=cex,
         xlim=xlim,ylim=ylim,
         axes=F,ylab=NA,xlab=NA,...)
    box()
    abline(0, 0, lty = 2, col = rgb(0.5, 0.5, 0.5))
    abline(v = 0, lty = 2, col = rgb(0.5, 0.5, 0.5))
    axis(side = 1, tck = -0.015, labels = NA)
    axis(side = 1, lwd = 0, line = -1, cex = 0.2, cex.axis = 0.9)
    axis(side = 2, tck = -0.015, lwd = 0, line = -1, cex.axis = 0.9)
    mtext(side = 1, xlab, line = 1, cex = 1)
    mtext(side = 2, ylab, line = 1, cex = 1)
    if (!text) {
      points(x[indlist],y[indlist],col=colmat[indlist],cex=indcex,...)
    } else {
      text(x[indlist],y[indlist],indlist,col=colmat[indlist],cex=indcex,...)
    }
  }
}

#' @export
KDEscatter = function(X,high=0.7,low=0.3,indlist=NULL,
                      plot.diag="upper",
                      colmat=NULL,ylim.list=NULL,
                      print.corr=TRUE, ...) {
  n = dim(X)[1]
  m = dim(X)[2]

  par(mfrow=c(m,m),mar=c(1.5,1.5,1.5,1.5))
  for (i in 1:m) {
    for (j in i:m) {
      if (i==j) {
        par(mfg=c(i,j))
        kdeplot.hy(X[,i],high=0.6,low=0.2,indlist=indlist,
                   colmat=colmat,xlim=ylim.list[[i]],...)
        title(main=colnames(X)[i])
      } else {
        if ((plot.diag=="both") | (plot.diag=="upper")) {
          par(mfg=c(i,j))
          plot.hy(x=X[,j],y=X[,i],indlist=indlist,
                  colmat=colmat,
                  xlim=ylim.list[[j]],ylim=ylim.list[[i]],...)
          if (print.corr) {
            legend("topleft",bty="n",
                   legend=paste("corr =",round(cor(X[,j],X[,i]),digits=3)))
          }
        }
        if ((plot.diag=="both") | (plot.diag=="lower")) {
          par(mfg=c(j,i))
          plot.hy(x=X[,i],y=X[,j],indlist=indlist,
                  colmat=colmat,
                  xlim=ylim.list[[i]],ylim=ylim.list[[j]],...)
        }
        # par(mfg=c(i,j))
        # plot.hy(x=X[,i],y=X[,j],colmat=colmat,pch=19)
      }
    }
  }
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
yaxis.hy <- function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}


#' @export
split_exon = function(exon) {
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

