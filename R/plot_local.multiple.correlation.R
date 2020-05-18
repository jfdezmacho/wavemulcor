plot_local.multiple.correlation <- #3.1.0.
  function(Lst, xaxt="s") {
    ##Producing correlation plot
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
        label <- xaxt[[2]]
        xaxt <- "n"
    }
    cor <- as.data.frame(Lst$cor)
    YmaxR <- Lst$YmaxR
    N <- length(YmaxR)
    xxnames <- names(Lst$data)
    par(mfrow=c(1,1), las=1, mar=c(5,4,4,2)+.1)
    ymin <- min(cor)
    ymax <- max(cor)
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    matplot(1:N, cor, type="n", xaxt=xaxt, ylim=c(ymin-0.1,ymax+0.1),
            xlab="", ylab="Local Multiple Correlation")
    abline(h=0)              ##Add Straight horiz and vert Lines to a Plot
    matlines(1:N,cor, lty=c(1,2,2), col=c(1,2,2))
    mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
    if (length(unique(YmaxR))==1) {
      mtext(xxnames[YmaxR][1], at=1, side=3, line=-1, cex=.8)
    } else {
      xvaru <- t(t(which(diff(sign(diff(as.matrix(cor[,"val"]))))==-2))+1)
      xvarl <- t(t(which(diff(sign(diff(as.matrix(cor[,"val"]))))==2))+1)
      xvar <- t(t(which(abs(diff(sign(diff(as.matrix(cor[,"val"])))))==2))+1)
      xvar2 <- xvar[-length(xvar)]+diff(xvar)/2
      mtext(xxnames[YmaxR][xvaru], at=xvaru, side=3, line=-.5, cex=.5)
      mtext(xxnames[YmaxR][xvar2], at=xvar2, side=3, line=-1, cex=.5)
      mtext(xxnames[YmaxR][xvarl], at=xvarl, side=3, line=-1.5, cex=.5)
    }
    if (xaxt!="s") axis(side=1, at=at, labels=label)
    return()
  }
