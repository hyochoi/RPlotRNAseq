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

