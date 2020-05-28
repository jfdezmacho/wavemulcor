plot_local.multiple.cross.correlation <- #3.1.0
  function(Lst, lmax, xaxt="s"){
    ##Producing correlation plot
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    val <- Lst$cor$vals
    low.ci <- Lst$cor$lower
    upp.ci <- Lst$cor$upper
    lag.max <- trunc((ncol(val)-1)/2)
    lag0 <- lag.max+1
    YmaxR <- Lst$YmaxR
    N <- length(YmaxR)
    xxnames <- names(Lst$data)
    lag.labs <- c(paste("lead",lag.max:1),paste("lag",0:lag.max))
    ymim <- -0.1
    if (length(Lst$data)<3) ymim <- -1
    par(mfcol=c(lmax+1,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    for(i in c(-lmax:0,lmax:1)+lag0) {
      plot(1:N,val[,i], type="l", lty=1, ylim=c(ymim,1), xaxt=xaxt,
              xlab="", ylab="", main=lag.labs[i])
      abline(h=0)                     ##Add Straight horiz
      lines(low.ci[,i], lty=1, col=2) ##Add Connected Line Segments to a Plot
      lines(upp.ci[,i], lty=1, col=2)
      mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
      # xvar <- seq(1,N,M)
      if (length(unique(YmaxR))==1) {
        mtext(xxnames[YmaxR][1], at=1, side=3, line=-1, cex=.8)
      } else {
        xvaru <- t(t(which(diff(sign(diff(as.matrix(val[,i]))))==-2))+1)
        xvarl <- t(t(which(diff(sign(diff(as.matrix(val[,i]))))==2))+1)
        # xvar <- t(t(which(abs(diff(sign(diff(as.matrix(val[,i])))))==2))+1)
        # xvar2 <- xvar[-length(xvar)]+diff(xvar)/2
        mtext(xxnames[YmaxR][xvaru], at=xvaru, side=3, line=-.5, cex=.5)
        # mtext(xxnames[YmaxR][xvar2], at=xvar2, side=3, line=-1, cex=.5)
        mtext(xxnames[YmaxR][xvarl], at=xvarl, side=3, line=-1, cex=.5)
      }
      if (xaxt!="s") axis(side=1, at=at, labels=label)
    }
    par(las=0)
    mtext('time', side=1, outer=TRUE, adj=0.5)
    mtext('Local Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5)
    return()
  }
