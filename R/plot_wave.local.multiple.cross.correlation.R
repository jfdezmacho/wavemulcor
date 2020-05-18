plot_wave.local.multiple.cross.correlation <- #3.1.0
  function(Lst, lmax, lag.first=FALSE, xaxt="s", pdf.write=NULL) {
    ##Producing cross-correlation plot
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    blocki <- function(lab){
      mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
      plot(1:N,vj$vals[,i], type="l", lty=1, ylim=c(ymin,1), xaxt=xaxt,
           xlab="", ylab="", main=lab)
      abline(h=0)                     ##Add Straight horiz
      lines(vj$lower[,i], lty=1, col=2) ##Add Connected Line Segments to a Plot
      lines(vj$upper[,i], lty=1, col=2)
      mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
      # xvar <- seq(1,N,M)
      if (length(unique(YmaxR[[j]]))==1) {
        mtext(xxnames[YmaxR[[j]]][1], at=1, side=3, line=-1, cex=.5)
      } else {
        xvaru <- t(t(which(diff(sign(diff(as.matrix(vj$vals[,i]))))==-2))+1)
        xvarl <- t(t(which(diff(sign(diff(as.matrix(vj$vals[,i]))))==2))+1)
        # xvar <- t(t(which(abs(diff(sign(diff(as.matrix(vj$vals[,i])))))==2))+1)
        # xvar2 <- xvar[-length(xvar)]+diff(xvar)/2
        mtext(xxnames[YmaxR[[j]]][xvaru], at=xvaru, side=3, line=-.5, cex=.3)
        # mtext(xxnames[YmaxR[[j]]][xvar2], at=xvar2, side=3, line=-1, cex=.3)
        mtext(xxnames[YmaxR[[j]]][xvarl], at=xvarl, side=3, line=-1, cex=.3)
      }
      if (xaxt!="s") axis(side=1, at=at, labels=label)
    }
    cor <- Lst$cor
    YmaxR <- Lst$YmaxR
    J <- length(Lst$YmaxR)-1
    N <- length(Lst$data[[1]][[1]])
    xxnames <- names(Lst$data)
    level.lab <- c(paste("level",1:J),paste("s",J))
    lagnames <- c(paste("Lead",lmax:1),paste("Lag",0:lmax))
    ymin <- -0.1
    if (length(Lst$data)<3) ymin <- -1
    if(lag.first){
      for(i in c(-lmax:0,lmax:1)+lmax+1) {
        if (!is.null(pdf.write))
          cairo_pdf(paste0("plot_",pdf.write,"_WLMCC_",lagnames[i],".pdf"), width=8.27,height=11.69)
        par(mfrow=c(ceiling((J+1)/2),2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,1.2,0))
        for(j in (J+1):1) {
          vj <- cor[[j]]
          blocki(level.lab[j])
        }
        par(las=0)
        mtext('time', side=1, outer=TRUE, adj=0.5)
        mtext('Wavelet Local Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5)
        mtext(lagnames[i], side=3, outer=TRUE, adj=0.5)
        if (!is.null(pdf.write)) dev.off()
      }
    } else{
      for(j in 1:(J+1)) {
        vj <- cor[[j]]
        if (!is.null(pdf.write))
          cairo_pdf(paste0(pdf.write,"_WLMCC_",level.lab[j],".pdf"), width=8.27,height=11.69)
        par(mfcol=c(lmax+1,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,1.2,0))
        for(i in c(-lmax:0,lmax:1)+lmax+1) {
        # for(i in 1:(2*lmax+1)) {
        blocki(lagnames[i])
        }
        par(las=0)
        mtext('time', side=1, outer=TRUE, adj=0.5)
        mtext('Wavelet Local Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5)
        mtext(level.lab[j], side=3, outer=TRUE, adj=0.5)
        if (!is.null(pdf.write)) dev.off()
      }
    }
    return()
  }
