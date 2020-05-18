plot_wave.multiple.correlation <- #3.1.0.
  function(Lst) {
    ##Producing correlation plot
    J <- length(Lst$YmaxR)-1
    xxnames <- names(Lst$data)
    cor <- Lst$xy.mulcor[1:J,]
    YmaxR <- Lst$YmaxR[1:J]
    par(mfrow=c(1,1), las=1, mar=c(5,4,4,2)+.1)
    ymin <- min(cor)
    ymax <- max(cor)
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    matplot(2^(0:(J-1)), cor, ylim=c(ymin-0.1,ymax+0.1),
            type="b", log="x", pch="*LU", cex=c(1.5,.75,.75), xaxt="n", lty=c(1,2,2), col=c(1,4,4),
            xlab="wavelet scale", ylab="Wavelet Multiple Correlation", cex.axis=0.75)
    mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
    abline(h=0)
    axis(side=1, at=2^(0:(J-1)))
    if (length(unique(YmaxR))==1) {
      mtext(xxnames[YmaxR][1], side=3, at=1, line=1, cex=.8)
    }else {
      mtext(xxnames[YmaxR], side=3, at=2^(0:J), line=-0.5, cex=.5)
    }
    return()
  }
