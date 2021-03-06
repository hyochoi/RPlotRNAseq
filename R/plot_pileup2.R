#' @export
plot_pileup4str = function(Pileup,PileupPS=NULL,Ranges,cases=NULL,logcount=NULL,
                       plot.meanpileup=TRUE,plot.meanpileupPS=TRUE,
                       col.pileup=NULL,col.meanpileup="grey",
                       main=NULL,cex.main=1.2,
                       print.ranges=TRUE,
                       xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,...) {

  ##  % needed variables
  if (missing(Pileup)) {
    stop("Pileup is missing")
  }
  if (missing(Ranges)) {
    stop("Genomic ranges should be needed (See get_Ranges)")
  }
  if (is.null(PileupPS)) {
    samplist = colnames(Pileup)
  } else {
    samplist = colnames(PileupPS)
  }
  n = length(samplist)
  exons = Ranges$lRanges[,c(2,3)]
  if (is.null(cases)) cases = 1:n
  if (!is.numeric(cases)) {
    caseIDs=cases
    cases=which(samplist %in% caseIDs)
    if (length(cases)==0) {
      stop("No sample matches with the given case IDs.")
    }
  }

  require("RColorBrewer")
  candicol1 = c(brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
                brewer.pal(8,"Pastel2")[6],
                brewer.pal(9,"YlOrBr")[1],
                brewer.pal(9,"YlOrRd")[1],
                brewer.pal(9,"YlOrRd")[2],
                brewer.pal(9,"Reds")[1],
                brewer.pal(9,"RdPu")[1],
                brewer.pal(9,"OrRd")[1],
                brewer.pal(9,"OrRd")[2],
                brewer.pal(9,"Oranges")[1],
                brewer.pal(9,"Oranges")[2],
                "aliceblue");
  candicol2 = brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
  candicol3 = brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
  candicol = c(candicol2,candicol3);
  exon.col = candicol1[9]

  # Set plot parameters
  if (is.null(xlim)) xlim = c(0,nrow(Pileup))
  if (is.null(ylim)) ylim = yaxis.hy(Pileup)
  if (is.null(xlab)) xlab = "genomic positions"
  if (is.null(ylab)) ylab = "read depth"
  if (is.null(main)) {
    if (length(cases)>1) {
      main = paste0(Ranges$Gene)
    } else if (length(cases)==1) {
      main = paste0(Ranges$Gene," | sample #",cases," (ID:",samplist[cases],")")
    }
  }
  # pileup colors
  if (is.null(col.pileup)) {
    if (length(cases)>10) {
      x.mda = svd(Pileup) ;
      projmat = diag(x.mda$d)%*%t(x.mda$v) ;
      projmat[1,] = -projmat[1,] ;
      # colset = rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
      colset = colorRampPalette(brewer.pal(10, "Spectral"))(n)[rank(-projmat[1,])]
    } else {
      colset = rep("black",n)
    }
  } else {
    if (length(col.pileup)==n) {
      colset=col.pileup
    } else {
      colset = rep("black",n)
      colset[cases] = rep(col.pileup,length=length(cases))
    }
  }

  # Start plotting
  if (print.ranges) {
    par(mar=c(3.2,3,3,1.5))
  } else {
    par(mar=c(3,3,3,1.5))
  }
  meanPileup = apply(Pileup, 1, median) ;
  plot(meanPileup, type='l', lty=2, lwd=0.5, ylim=c(min(0,ylim[1]),ylim[2]),
       xlim=xlim, axes=F, ylab=NA, xlab=NA, xaxs="i",yaxs="i", col="white") ;
  for (i in 1:nrow(exons)){
    polygon(x=c(rep(exons[i,1],2),rep(exons[i,2],2)),y=c(-10000,(max(Pileup)+10000),(max(Pileup)+10000),-10000),col=exon.col,border=NA) ;
  }
  abline(v=exons[,1],lty=1,col="lightyellow3",lwd=0.1) ;
  abline(v=exons[,2],lty=1,col="lightyellow3",lwd=0.1) ;
  title(main, cex.main=cex.main,font.main=1,line=0.5);

  if (print.ranges) {
    x.tick.at = c(1,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
    x.labels.l = c(Ranges$Gene,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
    x.labels.c = c(Ranges$Gene,Ranges$cRanges[2:nrow(Ranges$lRanges),1],max(Ranges$cRanges))
    x.labels.g = c(Ranges$chr,Ranges$gRanges[2:nrow(Ranges$gRanges),2],max(Ranges$gRanges))
    axis(side=1, tck=-0.01, at=x.tick.at, labels=NA, col.ticks="darkgrey") ;
    axis(side=1, lwd=0, line=-1, cex.axis=0.8,col.axis="darkgrey",
         at=x.tick.at,labels=x.labels.c) ;
    axis(side=1, lwd=0, line=-0.1, cex.axis=0.8,col.axis="darkgrey",
         at=x.tick.at,labels=x.labels.g) ;
    mtext(side=1, xlab, line=2, cex=1) ;
  } else {
    mtext(side=1, xlab, line=1, cex=1)
  }
  if (!is.null(logcount)) {
    if (logcount==1) {
      labels = c(1,5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    } else {
      labels = c(5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    }
    tick.at = log10(labels+logcount)-log10(logcount);
  } else {
    tick.at = NULL;
    labels = TRUE;
  }
  axis(side=2, tck=-0.02, at=tick.at, col.ticks="darkgrey",las=1,
       labels=labels,lwd=0,line=-0.8,cex.axis=0.8,col.axis="darkgrey")
  mtext(side=2, ylab, line=1.5, cex=1) ;

  box(lwd=1.5)
  if (plot.meanpileup){
    points(meanPileup, type='l', lty=1, lwd=2, col=col.meanpileup) ;
  }
  for (case in cases) {
    points(Pileup[,case], type='l', lty=1, lwd=1.5, col=colset[case], ...) ;
  }
  if (length(cases)==1) {
    d = nrow(Pileup)
    if (!is.null(PileupPS)) {
      polygon(cbind(c(1,c(1:d),d),c(0,Pileup[,case],0)),col="pink")
      if (plot.meanpileupPS){
        polygon(cbind(c(1,c(1:d),d),c(0,PileupPS[,case],0)),col="skyblue",border=FALSE)
        meanPileup2 = apply(PileupPS, 1, mean) ;
        points(meanPileup2, type='l', lty=1, lwd=1, col="black") ;
      } else {
        polygon(cbind(c(1,c(1:d),d),c(0,PileupPS[,case],0)),col="skyblue")
      }
    }
  }
}

yaxis.hy <- function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}
