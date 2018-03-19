mutpos0.hy <- function(x,exon,is.intron=FALSE,num.intron=NULL){
  ## %  First version: mutpos.hy.2
  ## %  Last updated: 11/26/2015
  ## %  Find mutation position in coverage file.
  a <- exon ;
  b <- strsplit(a,":")[[1]][2] ;
  c <- gsub("-",",",b) ;
  d <- as.numeric(strsplit(c,",")[[1]]) ;
  epl <- d[seq(1,length(d),by=2)] ;
  epr <- d[seq(2,length(d),by=2)] ;
  ep <- cbind(epl,epr) ;
  #  exon.count <- length(epl) ;
  #%#%  Find which exon contains the mutation
  which.exon.hy <- function(x,ep){
    if (is.null(dim(ep))){
      ep = matrix(ep,nrow=1)
    }
    int <- which((x-ep[,1]>=0) & (x-ep[,2]<=0));
    if (length(int)>0){
      return(list(int=int,in.exon=TRUE));
    } else {
      int <- max(which(x-ep[,1]>0));
      return(list(int=int,in.exon=FALSE));
    }
  }
  if (!is.intron){		#  Exons only
    mut.exon <- which.exon.hy(x,ep);
    ep.tmp <- as.matrix(ep[1:mut.exon$int,]);
    if (mut.exon$int==1){ ep.tmp <- t(ep.tmp)}
    if (mut.exon$in.exon){
      ep.tmp[mut.exon$int,2] <- x;
    }
    ep2 <- ep.tmp - ep.tmp[,1] + 1;
    if (nrow(ep2)==1){ y <- ep2[1,2] } else {
      new.ep <- ep2 + c(0,cumsum(ep2[1:(nrow(ep2)-1),2]));
      y <- new.ep[mut.exon$int,2];
    }
  } else {				#  Exons + Introns
    ep2 <- ep - ep[1,1] + 1;
    xp <- x - ep[1,1] + 1;
    if (is.null(num.intron)){	#  All introns
      y <- xp;
    } else {					#  Parts of introns
      ep2 <- cbind(ep2[,1]-num.intron,ep2,ep2[,2]+num.intron); ep2[1,1] <- 1;
      if (nrow(ep2)>1) {
        overlap.id <- which(ep2[1:(nrow(ep2)-1),4]>ep2[2:nrow(ep2),1]) ;
        ep2[overlap.id,4] <- ep2[overlap.id,3];
        ep2[(overlap.id+1),1] <- ep2[overlap.id,3] + 1;
        ep2[nrow(ep2),4] <- ep2[nrow(ep2),3];
      }
      #%#%
      y = rep(0,length(xp));
      for (i in 1:length(xp)){
        ixp = xp[i];
        mut.exon <- which.exon.hy(ixp,ep2[,c(1,4)]);
        ep2.tmp <- as.matrix(ep2[1:mut.exon$int,c(1,4)]);
        if (mut.exon$int==1){ ep2.tmp <- t(ep2.tmp) }
        if (mut.exon$in.exon){
          ep2.tmp[mut.exon$int,2] <- ixp;
        }
        ep3 <- ep2.tmp - ep2.tmp[,1] + 1;
        if (nrow(ep3)==1){ y[i] <- ep3[1,2] } else {
          new.ep <- ep3 + c(0,cumsum(ep3[1:(nrow(ep3)-1),2]));
          y[i] <- new.ep[mut.exon$int,2];
        }
      }
    }
  }
  return(y);
}
