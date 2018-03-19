mutpos.hy = function(x,exon,is.intron=FALSE,num.intron=NULL){
  a <- exon ;
  b <- strsplit(a,":")[[1]][2] ;
  c <- gsub("-",",",b) ;
  d <- as.numeric(strsplit(c,",")[[1]]) ;
  seqdir = strsplit(exon,":")[[1]][3] ;

  epl <- d[seq(1,length(d),by=2)] ;
  epr <- d[seq(2,length(d),by=2)] ;
  ep <- cbind(epl,epr) ;

  endpos = epr[length(epr)];
  tlen = mutpos0.hy(x=endpos,exon=exon,is.intron=is.intron,num.intron=num.intron);
  mutpos = mutpos0.hy(x=x,exon=exon,is.intron=is.intron,num.intron=num.intron);

  if (seqdir=="+") {
    return(mutpos);
  } else {
    return(tlen-mutpos+1);
  }
}
