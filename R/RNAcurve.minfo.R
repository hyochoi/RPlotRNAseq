RNAcurve.minfo = function(datmat,indlist=NULL,dai,barcode,
                          mpar=NULL,one.plot=FALSE,same.yaxis=TRUE,
                          minfo=NULL,mlegend=TRUE,mcol=NULL,
                          mzoom=FALSE,mtype.zoom="Splice_Site",
                          colmat=NULL,
                          title.id=TRUE, title="Sample", title.cex=1, title.barcode=TRUE,
                          mean.plot=TRUE, mean.col="grey", fullpar=FALSE,ylim=NULL,xlim=NULL,
                          ep.col=NULL,
                          ex.num=NULL,int.num=NULL,ex.num.col=NULL,int.num.col=NULL,
                          ylab=NULL, xlab=NULL,yaxis.logcount=NULL,...) {
  ## % Draw individual curves with the mean curve of datmat and vertical lines on exon edges.
  ## % datmat        : d by n  (columns are each sample)
  ## % indlist       : a vector of lists of sample id to be drawn
  ## % mutation.info : mutation information for one gene to be plotted.
  ##                    - all mutations will be indicated by vertical lines
  ##                    - If NULL, no action. (only curve)
  ## % colmat        : color matrix to be used. (length should be same as n)
  ## % ex.num        : exon list that will be colored. (from the right in figure)
  ## % int.num       : intron list that will be colored. (from the right in figure)
  ## % Updated       : September 15, 2017
  ## % Hyo Young Choi

  #  Create colmat by using SVD
  if (is.null(colmat)) {
    x.mda <- svd(datmat) ;
    projmat <- diag(x.mda$d)%*%t(x.mda$v) ;
    projmat[1,] <- -projmat[1,] ;
    colmat <- rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
  }

  epl = dai$epm[,2]; epr = dai$epm[,3]
  num.intron = dai$intron.len;

  if (!is.null(minfo)) {
    mutation.info = minfo[,c(8,9,6)]; mut_start=NULL; mut_end=NULL;
    if (nrow(mutation.info)>0) {
      mut_start = apply(matrix(as.numeric(minfo[,3]),ncol=1),1,
                        FUN=function(x){mutpos.hy(x,exon=exon,is.intron=TRUE,num.intron=dai$intron.len)})
      mut_end = apply(matrix(as.numeric(minfo[,4]),ncol=1),1,
                      FUN=function(x){mutpos.hy(x,exon=exon,is.intron=TRUE,num.intron=dai$intron.len)})
    }
    mutation.info = cbind(mutation.info,mut_start,mut_end);

    mut_types =  c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
                   "Nonsense_Mutation","Nonstop_Mutation","RNA","Silent","Splice_Site","Translation_Start_Site");
  }

  if (one.plot) {
    if (is.null(indlist)) {
      case=1
      RNAcurveOne.minfo(datmat=datmat,indlist=case,dai,barcode,
                        minfo=NULL,mlegend=FALSE,mcol=mcol,
                        colmat=colmat, same.yaxis=TRUE,
                        title.id=FALSE,title=title, title.cex=title.cex,
                        mean.plot=FALSE, mean.col=,mean.col,ylim=ylim,xlim=xlim,
                        ep.col=ep.col,title.barcode=FALSE,
                        ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                        ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
      for (case in 2:ncol(datmat)) {
        points(datmat[,case],type='l', lty=1, col=colmat[case], ...)
      }
    } else {
      RNAcurveOne.minfo(datmat=datmat,indlist=indlist[1],dai=dai,barcode=barcode,
                        minfo=minfo,mlegend=FALSE,mcol=mcol,
                        colmat=colmat, same.yaxis=TRUE,
                        title.id=FALSE,title=title, title.cex=title.cex,
                        mean.plot=FALSE, mean.col=,mean.col,ylim=ylim,xlim=xlim,
                        ep.col=ep.col,title.barcode=FALSE,
                        ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                        ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
      for (case in indlist[2:length(indlist)]) {
        points(datmat[,case],type='l', lty=1, col=colmat[case], ...)
      }
    }
  } else {
    ncurve = length(indlist);
    if (ncurve==1) {
      if (mzoom) {
        case.barcode = barcode[indlist];
        if (!is.null(minfo)) {
          if (nrow(mutation.info)>0) {
            mutype0 = mutation.info[which(mutation.info$barcode2==case.barcode),3];
            mutpos0 = mutation.info[which(mutation.info$barcode2==case.barcode),4];
            mutpos0 = mutpos0[which(mutype0==mtype.zoom)][1]
            if (length(mutpos0)>0) {
              mutpos1 = mutpos0 ;
              xlim = c(mutpos1-200,mutpos1+200)
            }
          }# (nrow(mutation.info)>0)
        }
      }# if (mzoom)

      RNAcurveOne.minfo(datmat=datmat,indlist=indlist,dai=dai,barcode=barcode,
                        minfo=minfo,mlegend=mlegend,mcol=mcol,
                        colmat=colmat,
                        title.id=title.id,title=title, title.cex=title.cex,
                        mean.plot=mean.plot, mean.col=mean.col,
                        same.yaxis=same.yaxis,ylim=ylim,xlim=xlim,
                        ep.col=ep.col,title.barcode=title.barcode,
                        ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                        ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
    } else {
      if (fullpar){
        par(mfrow=c(6,5), mar=c(2,3,1.5,0.5))
        for (id in indlist){
          if (mzoom) {
            case.barcode = barcode[indlist];
            if (!is.null(minfo)) {
              if (nrow(mutation.info)>0) {
                mutype0 = mutation.info[which(mutation.info$barcode2==case.barcode),3];
                mutpos0 = mutation.info[which(mutation.info$barcode2==case.barcode),4];
                mutpos0 = mutpos0[which(mutype0==mtype.zoom)][1]
                if (length(mutpos0)>0) {
                  mutpos1 = mutpos0 ;
                  xlim = c(mutpos1-200,mutpos1+200)
                }
              }# (nrow(mutation.info)>0)
            }
          }# if (mzoom)

          RNAcurveOne.minfo(datmat=datmat,indlist=id,dai=dai,barcode=barcode,
                            minfo=minfo,mlegend=mlegend,mcol=mcol,
                            colmat=colmat, same.yaxis=same.yaxis,
                            title.id=title.id,title=title, title.cex=title.cex,
                            mean.plot=mean.plot, mean.col=mean.col,ylim=ylim,xlim=xlim,
                            ep.col=ep.col,title.barcode=title.barcode,
                            ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                            ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
        }
      } else {
        if (is.null(mpar)) {
          if ((1 < ncurve) & (ncurve <=4)) {
            par(mfrow=c(1,4), mar=c(2,3,1.5,0.5))
          } else if ((4 < ncurve) & (ncurve <=16)) {
            kk=ceiling(ncurve/4)
            par(mfrow=c(kk,4), mar=c(2,2,1.5,0.5))
          } else if ((16 < ncurve) & (ncurve <=30)) {
            kk=ceiling(ncurve/5)
            par(mfrow=c(kk,5), mar=c(2,3,1.5,0.5))
          } else {
            par(mfrow=c(6,5), mar=c(2,3,1.5,0.5))
          }
        } else {
          par(mfrow=c(6,5), mar=c(2,3,1.5,0.5))
        }
        for (id in indlist){
          if (mzoom) {
            case.barcode = barcode[indlist];
            if (!is.null(minfo)) {
              if (nrow(mutation.info)>0) {
                mutype0 = mutation.info[which(mutation.info$barcode2==case.barcode),3];
                mutpos0 = mutation.info[which(mutation.info$barcode2==case.barcode),4];
                mutpos0 = mutpos0[which(mutype0==mtype.zoom)][1]
                if (length(mutpos0)>0) {
                  mutpos1 = mutpos0 ;
                  xlim = c(mutpos1-200,mutpos1+200)
                }
              }# (nrow(mutation.info)>0)
            }
          }# if (mzoom)

          RNAcurveOne.minfo(datmat=datmat,indlist=id,dai=dai,barcode=barcode,
                            minfo=minfo,mlegend=mlegend,mcol=mcol,
                            colmat=colmat, same.yaxis=same.yaxis,
                            title.id=title.id,title=title, title.cex=title.cex,
                            mean.plot=mean.plot, mean.col=mean.col,ylim=ylim,xlim=xlim,
                            ep.col=ep.col,title.barcode=title.barcode,
                            ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                            ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
        }
      }
    }
  }
}
