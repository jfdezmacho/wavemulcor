plot_wave.multiple.regression <- #3.1.0.
  function(Lst, nsig=2) {
    ##Producing regression plot
    # requireNamespace(magrittr)
    J <- length(Lst$YmaxR)-1
    vals <- Lst$xy.mulreg$rval[1:J,-1]  #exclude constant
    stdv <- Lst$xy.mulreg$rstd[1:J,-1]
    lows <- Lst$xy.mulreg$rlow[1:J,-1]
    upps <- Lst$xy.mulreg$rupp[1:J,-1]
    pval <- Lst$xy.mulreg$rpva[1:J,-1]
    order <- Lst$xy.mulreg$rord[1:J,-1]-1
    order[order==0] <- NA       #exclude dependent variable
    vars <- t(matrix(colnames(vals),length(Lst$data),J))
    YmaxR <- Lst$YmaxR[1:J]
    xxnames <- names(Lst$data)
    sel <- order<=nsig & pval<=0.05              #select firs nsig 5%signif. predictors
    vals.sig <- vals*sel
    lows.sig <- lows*sel
    upps.sig <- upps*sel
    # pval.sig <- pval*sel
    order.sig <- order*sel
    vals.sig[vals.sig==0] <-
      lows.sig[lows.sig==0] <-
      upps.sig[upps.sig==0] <-
      # pval.sig[pvals.sig==0] <-
      order.sig[order.sig==0] <- NA
    # requireNamespace(RColorBrewer)
    mycolors <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
    par(mfrow=c(1,1), las=1, mar=c(5,4,4,2)+.1)
    ymin <- min(vals.sig,na.rm=TRUE)
    ymax <- max(vals.sig,na.rm=TRUE)
    x <- 2^(0:(J-1))
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    matplot(x,vals, log="x", ylim=c(ymin-0.1,ymax+0.1),
            type="n", xaxt="n", lty=3, col=8, cex.axis=0.75,
            xlab="wavelet scale", ylab="Wavelet Multiple Regression")
    # shade <- 1.96*stdv %>% apply(1,max)
    for (i in ncol(stdv):1){
      shade <- 1.96*stdv[,i]
      polygon(c(x,rev(x)),c(-shade,rev(shade)), col=gray(0.8,alpha=0.2), border=NA)
    }
    matlines(x,vals, log="x", lty=1, col=8)
    # v <- (vals*(pval<=0.05)) %>% replace(.==0,NA) #%>% replace(.==-1.,NA)
    # matlines(x,v, log="x",lty=1, col=mycolors[8])
    # v <- (vals*(rbind(pval[-1,],tail(pval,1))<=0.05)) %>% replace(.==0,NA) #%>% replace(.==-1.,NA)
    # matlines(tail(x,2),v[(J-1):J,], log="x",lty=1, col=mycolors[8])
    if(abs(ymax-ymin)<3) lo<-2 else lo<-4
    mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
    abline(h=seq(floor(ymin),ceiling(ymax),length.out=lo),col=8)
    matlines(x, vals.sig, log="x", type="b", pch="*", lty=1, lwd=2,  col=mycolors)
    matlines(x, lows.sig, log="x", type="b", pch="*",lty=2, col=mycolors)
    matlines(x, upps.sig, log="x", type="b", pch="*",lty=2, col=mycolors)
    col <- (order.sig<=nsig)*1 +(order.sig>=nsig)*8
    # xvar <- seq(1,N,M)
    text(x, vals.sig, labels=vars, col=col, cex=.5)
    text(x, vals.sig, labels=order,pos=1, col=col, cex=.5)
    if (length(unique(YmaxR))==1) {
      mtext(xxnames[YmaxR][1],at=1,  side=3, line=-1, outer=TRUE, cex=.8)
    }else {
      mtext(xxnames[YmaxR], at=2^(0:J),  side=3, line=-0.5, cex=.5)
    }
    axis(side=1, at=x)
    return()
  }
