plot_local.multiple.cross.regression <- #3.1.0
  function(Lst, lmax, nsig=2, xaxt="s"){
    ##Producing regression plot
    # requireNamespace(magrittr)
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    val <- Lst$cor$vals
    reg.vals <- Lst$reg$rval[,,-1]       #exclude constant
    reg.stdv <- Lst$reg$rstd[,,-1]
    reg.lows <- Lst$reg$rlow[,,-1]
    reg.upps <- Lst$reg$rupp[,,-1]
    reg.pval <- Lst$reg$rpva[,,-1]
    reg.order <- Lst$reg$rord[,,-1]-1
    reg.order[reg.order==0] <-
      reg.vals[reg.order==0] <-            #exclude dependent variable
      reg.stdv[reg.order==0] <-
      reg.lows[reg.order==0] <-
      reg.upps[reg.order==0] <-
      reg.pval[reg.order==0] <- NA
    lag.max <- trunc((ncol(val)-1)/2)
    lag0 <- lag.max+1
    YmaxR <- Lst$YmaxR
    N <- length(YmaxR)
    xxnames <- names(Lst$data)
    lag.labs <- c(paste("lead",lag.max:1),paste("lag",0:lag.max))
    reg.vars <- t(matrix(xxnames,length(Lst$data),N))
    reg.sel <- reg.order<=nsig & reg.pval<=0.05              #select firs nsig 5%signif. predictors
    reg.vals.sig <- reg.vals*reg.sel
    reg.lows.sig <- reg.lows*reg.sel
    reg.upps.sig <- reg.upps*reg.sel
    # reg.pval.sig <- reg.pval*reg.sel
    reg.order.sig <- reg.order*reg.sel
    reg.vals.sig[reg.vals.sig==0] <-
      reg.lows.sig[reg.lows.sig==0] <-
      reg.upps.sig[reg.upps.sig==0] <-
      # reg.pval.sig[reg.pvals.sig==0] <-
      reg.order.sig[reg.order.sig==0] <- NA
    # requireNamespace(RColorBrewer)
    mycolors <- RColorBrewer::brewer.pal(n=8, name="Dark2")
    par(mfcol=c(lmax+1,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
    ymin <- min(reg.vals,na.rm=T) #head(unique(sort(reg.vals[,,-1])))[2]
    ymax <- max(reg.vals,na.rm=T)
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    for(i in c(-lmax:0,lmax:1)+lag0) {
      matplot(1:N,reg.vals[,i,], ylim=c(ymin-0.1,ymax+0.1),
              type="n", xaxt=xaxt, lty=3, col=8,
              xlab="", ylab="", main=lag.labs[i])
      # shade <- 1.96*reg.stdv %>% apply(1,max)
      for (j in dim(reg.stdv)[3]:1){
        shade <- 1.96*reg.stdv[,i,j]
        polygon(c(1:N,rev(1:N)),c(-shade,rev(shade)), col=gray(0.8,alpha=0.2), border=NA)
      }
      matlines(1:N,reg.vals[,i,], lty=1, col=8)
      if(abs(ymax-ymin)<3) lo<-2 else lo<-4
      abline(h=seq(floor(ymin),ceiling(ymax),length.out=lo),col=8)
      matlines(1:N, reg.vals.sig[,i,], lty=1, lwd=2, col=mycolors)
      matlines(1:N, reg.lows.sig[,i,], lty=2, col=mycolors)
      matlines(1:N, reg.upps.sig[,i,], lty=2, col=mycolors)
      mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
      col <- (reg.order[,i,]<=3)*1 +(reg.order[,i,]>3)*8
      # xvar <- seq(1,N,M)
      xvar <- t(t(which(abs(diff(sign(diff(reg.vals[,i,]))))==2,arr.ind=T))+c(1,0))
      text(xvar, reg.vals[xvar,i,], labels=reg.vars[xvar,], col=col,cex=.3)
      text(xvar, reg.vals[xvar,i,], labels=reg.order[xvar,i,],pos=1, col=col,cex=.3)
      if (length(unique(YmaxR))==1) {
        mtext(xxnames[YmaxR][1], at=1, side=3, line=-1, cex=.5)
      } else {
        xvaru <- t(t(which(diff(sign(diff(as.matrix(val[,i]))))==-2))+1)
        xvarl <- t(t(which(diff(sign(diff(as.matrix(val[,i]))))==2))+1)
        # xvar <- t(t(which(abs(diff(sign(diff(as.matrix(val[,i])))))==2))+1)
        # xvar2 <- xvar[-length(xvar)]+diff(xvar)/2
        mtext(xxnames[YmaxR][xvaru], at=xvaru, side=3, line=-.5, cex=.3)
        # mtext(xxnames[YmaxR][xvar2], at=xvar2, side=3, line=-1, cex=.5)
        mtext(xxnames[YmaxR][xvarl], at=xvarl, side=3, line=-1, cex=.3)
      }
      if (xaxt!="s") axis(side=1, at=at, labels=label)
    }
    par(las=0)
    mtext('time', side=1, outer=TRUE, adj=0.5)
    mtext('Local Multiple Cross-Regression', side=2, outer=TRUE, adj=0.5)
    return()
  }
