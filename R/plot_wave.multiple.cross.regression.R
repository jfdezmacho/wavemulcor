
plot_wave.multiple.cross.regression <- #3.1.1
  function(Lst, lmax, nsig=2, by=3) {
    ##Producing cross-regression plot
    # requireNamespace(magrittr)
    J <- length(Lst$YmaxR)-1
    YmaxR <- Lst$YmaxR[1:J]
    xxnames <- names(Lst$data)
    vars <- t(matrix(xxnames,length(Lst$data),2*lmax+1)) #Lst$xy.mulreg$vars[1:J,]
    valnames <- c(paste("level",1:J),paste("s",J))
    # requireNamespace(RColorBrewer)
    mycolors <- RColorBrewer::brewer.pal(n=8, name="Dark2")
    par(mfrow=c(ceiling((J+1)/2),2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    for(j in (J+1):1) {
      vals <- Lst$xy.mulreg$rval[[j]][,-1] #exclude constant
      stdv <- Lst$xy.mulreg$rstd[[j]][,-1]
      lows <- Lst$xy.mulreg$rlow[[j]][,-1]
      pval <- Lst$xy.mulreg$rpva[[j]][,-1]
      upps <- Lst$xy.mulreg$rupp[[j]][,-1]
      order <- Lst$xy.mulreg$rord[[j]][,-1]-1
      order[order==0] <-
        vals[order==0] <-                  #exclude dependent variable
        stdv[order==0] <-
        lows[order==0] <-
        upps[order==0] <-
        pval[order==0] <- NA
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
      # ymin <- min(lows,na.rm=TRUE)
      # ymax <- max(upps,na.rm=TRUE)
      ymin <- min(vals,na.rm=TRUE)
      ymax <- max(vals,na.rm=TRUE)
      matplot(1:(2*lmax+1),vals, ylim=c(ymin-0.1,ymax+0.1), xaxt="n",
              type="n", lty=3,
              xlab="", ylab="", main=valnames[j], col=8)
      # shade <- 1.96*reg.stdv %>% apply(1,max)
      for (i in ncol(stdv):1){
        shade <- 1.96*stdv[,i]
        polygon(c(1:(2*lmax+1),(2*lmax+1):1),c(-shade,rev(shade)), col=gray(0.8,alpha=0.2), border=NA)
      }
      matlines(1:(2*lmax+1),vals, lty=1, col=8)
      if(abs(ymax-ymin)<3) lo<-2 else lo<-4
      abline(h=seq(floor(ymin),ceiling(ymax),length.out=lo),col=8)
      abline(v=lmax+1,col=8)
      matlines(1:(2*lmax+1),vals.sig, lty=1, lwd=2, col=mycolors)
      matlines(1:(2*lmax+1),lows.sig, lty=2, col=mycolors)
      matlines(1:(2*lmax+1),upps.sig, lty=2, col=mycolors)
      mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
      if(j<3) {axis(side=1, at=seq(1, 2*lmax+1, by=lmax/by), labels=seq(-lmax, lmax, by=lmax/by))}
      #axis(side=2, at=c(-.2, 0, .5, 1))
      xvar <- seq(1,2*lmax+1,6)
      col <- (order<=3)*1 +(order>1)*8
      text(xvar,vals[xvar,], labels=vars[xvar,], col=col[xvar,], cex=.5)
      text(xvar,vals[xvar,], labels=order[xvar,],pos=1, col=col[xvar,], cex=.5)
      if (length(unique(YmaxR))==1) {
        mtext(xxnames[YmaxR][1], side=3, at=1, line=-1, cex=.8)
      }else {
        mtext(xxnames[YmaxR], side=3, at=seq(1,2*lmax+1,by=lmax/by), line=-1, cex=.5)
      }
    }
    par(las=0)
    mtext('lag', side=1, outer=TRUE, adj=0.5)
    mtext('Wavelet Multiple Cross-Regression', side=2, outer=TRUE, adj=0.5)
    return()
    }
