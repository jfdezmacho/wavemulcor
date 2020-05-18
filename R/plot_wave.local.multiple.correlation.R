plot_wave.local.multiple.correlation <- #3.1.0.
  function(Lst, xaxt="s"){
    ##Producing correlation plot
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    cor <- Lst$cor
    J <- length(Lst$YmaxR)-1
    N <- length(Lst$data[[1]][[1]])
    valnames <- c(paste("level",1:J),paste("s",J))
    par(mfrow=c(ceiling((J+1)/2),2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    for(j in (J+1):1) {
      vj <- as.data.frame(cor[[j]]) #cbind(val[,j],lo[,j],up[,j])
      ymin <- min(vj, na.rm=TRUE)
      ymax <- max(vj, na.rm=TRUE)
      matplot(1:N,vj, ylim=c(ymin-0.1,ymax+0.1), #before:dont remember why vj[,-2] here and lines(lo[,j]...) below
              type="l", lty=c(1,2,2), col=c(1,2,2), xaxt=xaxt,
              xlab="", ylab="", main=valnames[j])
      if(xaxt!="s" & j<=(J+1)) {axis(side=1, at=at, labels=label)}
      mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
      abline(h=0)              ##Add Straight horiz and vert Lines to a Plot
    }
    par(las=0)
    # mtext('time', side=1, outer=TRUE, adj=0.5)
    mtext('Wavelet Local Multiple Correlation', side=2, outer=TRUE, adj=0.5)
  }
