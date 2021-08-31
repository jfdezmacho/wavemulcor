plot_local.multiple.regression <- #3.1.0.
  function(Lst, nsig=2, xaxt="s") {
    ##Producing regression plot
    # requireNamespace(magrittr)
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    cor <- as.data.frame(Lst$cor)
    reg.vals <- Lst$reg$rval[,-1]       #exclude constant
    reg.stdv <- Lst$reg$rstd[,-1]
    reg.lows <- Lst$reg$rlow[,-1]
    reg.upps <- Lst$reg$rupp[,-1]
    reg.pval <- Lst$reg$rpva[,-1]
    reg.order <- Lst$reg$rord[,-1]-1
    reg.order[reg.order==0] <-
      reg.vals[reg.order==0] <-            #exclude dependent variable
      reg.stdv[reg.order==0] <-
      reg.lows[reg.order==0] <-
      reg.upps[reg.order==0] <-
      reg.pval[reg.order==0] <- NA
    YmaxR <- Lst$YmaxR
    N <- length(YmaxR)
    xxnames <- names(Lst$data)
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
    par(mfrow=c(1,1), las=1, mar=c(5,4,4,2)+.1)
    ymin <- min(reg.vals.sig,na.rm=TRUE)
    ymax <- max(reg.vals.sig,na.rm=TRUE)
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    matplot(1:N,reg.vals, ylim=c(ymin-0.1,ymax+0.1),
            type="n", xaxt=xaxt, lty=3, col=8,
            xlab="", ylab="Local Multiple Regression")
    # shade <- 1.96*reg.stdv %>% apply(1,max)

    # for (i in 1:6){
    #   y1 <- yearcrisis[1,i]; y2 <- yearcrisis[2,i]
    #   ndy <- ifelse (y2<20, 260, 79)          #nof days in year = 5*52=260, except last year 2020=79
    #   days <- c(((y1-1)*260+1),((y2-1)*260+ndy)-10)
    #   polygon(c(days,rev(days)), c(ymin-0.1,ymin-0.1,ymax+0.1,ymax+0.1), col="grey95", border=NA)
    # }

    for (i in ncol(reg.stdv):1){
      shade <- 1.96*reg.stdv[,i]
      polygon(c(1:N,rev(1:N)),c(-shade,rev(shade)), col=gray(0.8,alpha=0.2), border=NA)
    }
    matlines(1:N,reg.vals, lty=1, col=8)
    # v <- (reg.vals*(reg.pval<=0.05) %>% replace(.==0,NA)) #%>% replace(.==-1.,NA)
    # matlines(1:N,v, lty=1, col=mycolors[8])
    if(abs(ymax-ymin)<3) lo<-2 else lo<-4
    abline(h=seq(floor(ymin),ceiling(ymax),length.out=lo),col=8)
    matlines(1:N, reg.vals.sig, lty=1, lwd=2,  col=mycolors)
    matlines(1:N, reg.lows.sig, lty=2, col=mycolors)
    matlines(1:N, reg.upps.sig, lty=2, col=mycolors)
    mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
    col <- (reg.order.sig<=nsig)*1 +(reg.order.sig>=nsig)*8
    # xvar <- seq(1,N,M)
    xvar <- t(t(which(abs(diff(sign(diff(reg.vals.sig))))==2,arr.ind=T))+c(1,0))
    text(xvar[,1], reg.vals.sig[xvar], labels=reg.vars[xvar], col=col[xvar], cex=.8)
    text(xvar[,1], reg.vals.sig[xvar], labels=reg.order[xvar],pos=1, col=col[xvar],cex=.5)
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
