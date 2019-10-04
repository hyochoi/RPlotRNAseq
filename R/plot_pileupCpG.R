#' @import SCISSOR
#' @export
plot_pileupCpG = function(Pileup,Ranges,cpgData,cases=NULL, ...) {
  require(SCISSOR)
  if (missing(cpgData)) {
    stop("cpgData is missing")
  }
  chrnum = Ranges$chr
  cpg.island0 = cpgData[which(cpgData$chr==chrnum),]
  cpg.island = cpg.island0[which((cpg.island0$end>min(Ranges$gRanges)) & (cpg.island0$start<max(Ranges$gRanges))),]

  d = max(Ranges$lRanges)
  if (nrow(cpg.island)>0) {
    cgi.pos=cbind(find_region(gpos=cpg.island$start,Ranges=Ranges),find_region(gpos=cpg.island$end,Ranges=Ranges))

    if (Ranges$strand == "-") {
      cgi.pos=matrix(cgi.pos[,c(2,1)],ncol=2)
    }
    cgi.pos[which(cgi.pos[,1]<=0),1] = 1
    cgi.pos[which(cgi.pos[,2]>d),2] = d
  }

  cgi.bar=rep(0,d)
  for (i in 1:nrow(cgi.pos)) {
    if (cgi.pos[i,1]<cgi.pos[i,2]) {
      cgi.bar[c(cgi.pos[i,1]:cgi.pos[i,2])]=cpg.island$perCpG[i]/100
    }
  }

  m=matrix(c(c(0,1,0.2,1),c(0,1,0,0.2)),byrow=T,ncol=4)
  split.screen(m)

  palette_SCISSOR()
  screen(2,new=TRUE)
  par(mar=c(3,3.5,0.5,2))
  image(matrix(cgi.bar,ncol=1),col=brewer.pal(8,"Greens"),axes=F)
  axis(1,at=apply(cgi.pos,1,mean)/d,labels=cpg.island$perCpG,line=-0.5,lwd=0,cex.axis=1.5)
  par(mar=c(5,4,4,1)+0.1)

  screen(1,new=TRUE)
  par(mar=c(1.5,3.5,2,2))
  plot_pileup(Pileup=Pileup,Ranges=Ranges,cases=cases, ...)
  par(mar=c(5,4,4,1)+0.1)
}
