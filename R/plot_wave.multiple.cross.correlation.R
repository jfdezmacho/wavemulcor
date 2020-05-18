plot_wave.multiple.cross.correlation <- #3.1.0.
  function (Lst, lmax, by=3) {
    ##Producing cross-correlation plot
    J <- length(Lst$YmaxR)-1
    xxnames <- names(Lst$data)
    returns.cross.cor <- Lst$xy.mulcor
    returns.lower.ci <- Lst$ci.mulcor$lower
    returns.upper.ci <- Lst$ci.mulcor$upper
    YmaxR <- Lst$YmaxR
    valnames <- c(paste("level",1:J),paste("s",J))
    par(mfrow=c(ceiling((J+1)/2),2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
    ymin <- -0.1
    if (length(Lst$data)<3) ymin <- -1
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    for(j in (J+1):1) {
      matplot((1:(2*lmax+1)),returns.cross.cor[j,], type="l", lty=1, ylim=c(ymin,1), xaxt="n",
              xlab="", ylab="", main=valnames[j])
      if(j<3) {axis(side=1, at=seq(1,2*lmax+1,by=lmax/by), labels=seq(-lmax,lmax,by=lmax/by))}
      #axis(side=2, at=c(-.2, 0, .5, 1))
      abline(h=0,v=lmax+1)              ##Add Straight horiz and vert Lines to a Plot
      lines(returns.lower.ci[j,], lty=2, col=2) ##Add Connected Line Segments to a Plot
      lines(returns.upper.ci[j,], lty=2, col=2)
      mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
      if (length(unique(YmaxR))==1) {
        mtext(xxnames[YmaxR][1], side=3, at=1, line=-1, cex=.8)
      }else {
        mtext(xxnames[YmaxR], side=3, at=seq(1,2*lmax+1,by=lmax/by), line=-1, cex=.5)
      }
    }
    par(las=0)
    mtext('lag', side=1, outer=TRUE, adj=0.5)
    mtext('Wavelet Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5)
    return()
  }
