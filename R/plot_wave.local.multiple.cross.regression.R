plot_wave.local.multiple.cross.regression <- #3.1.0
  function(Lst, lmax, nsig=2, xaxt="s", pdf.write=NULL){
    ##Producing cross-regression plot
    # requireNamespace(magrittr)
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    cor <- Lst$cor
    reg <- Lst$reg
    YmaxR <- Lst$YmaxR
    J <- length(Lst$YmaxR)-1
    N <- length(Lst$data[[1]][[1]])
    xxnames <- names(Lst$data)
    level.lab <- c(paste("level",1:J),paste("s",J))
    lagnames <- c(paste("Lead",lmax:1),paste("Lag",0:lmax))
    reg.vars <- t(matrix(xxnames,length(Lst$data),N))
    # requireNamespace(RColorBrewer)
    mycolors <- RColorBrewer::brewer.pal(n=8, name="Dark2")
    reg.vars <- t(matrix(xxnames,length(Lst$data),N)) #xy.mulcor$xy.mulreg$vars[1:J,]
    for(j in 1:(J+1)) {
      vj <- lapply(reg[[j]],function(x){x[,,-1]}) #exclude constant
      vj$rord <- vj$rord-1
      vj$rord[vj$rord==0] <-
        vj$rval[vj$rord==0] <-            #exclude dependent variable
        vj$rstd[vj$rord==0] <-
        vj$rlow[vj$rord==0] <-
        vj$rupp[vj$rord==0] <-
        vj$rpva[vj$rord==0] <- NA
      vj.sel <- vj$rord<=nsig & vj$rpva<=0.05              #select firs nsig 5%signif. predictors
      vj$rval.sig <- vj$rval*vj.sel
      vj$rlow.sig <- vj$rlow*vj.sel
      vj$rupp.sig <- vj$rupp*vj.sel
      # vj$rpva.sig <- vj$rpva*vj.sel
      vj$rord.sig <- vj$rord*vj.sel
      vj$rval.sig[vj$rval.sig==0] <-
        vj$rlow.sig[vj$rlow.sig==0] <-
        vj$rupp.sig[vj$rupp.sig==0] <-
        vj$rord.sig[vj$rord.sig==0] <- NA
      if (!is.null(pdf.write))
        cairo_pdf(paste0("plot_",pdf.write,"_WLMCR_",level.lab[j],".pdf"), width=8.27,height=11.69)
      mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
      par(mfcol=c(lmax+1,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,1.2,0))
      ymin <- min(vj$rval,na.rm=TRUE) #head(unique(sort(vj$rval)))[2]
      ymax <- max(vj$rval,na.rm=TRUE)
      for(i in c(-lmax:0,lmax:1)+lmax+1) {
        matplot(1:N,vj$rval[,i,], ylim=c(ymin-0.1,ymax+0.1),
                type="n", xaxt=xaxt, lty=3, col=8,
                xlab="", ylab="", main=lagnames[i])
        # shade <- 1.96*vj$rstd %>% apply(1,max)
        for (k in dim(vj$rstd)[3]:1){
          shade <- 1.96*vj$rstd[,i,k]
          polygon(c(1:N,rev(1:N)),c(-shade,rev(shade)), col=gray(0.8,alpha=0.2), border=NA)
        }
        matlines(1:N,vj$rval[,i,], lty=1, col=8)
        if(abs(ymax-ymin)<3) lo<-2 else lo<-4
        abline(h=seq(floor(ymin),ceiling(ymax),length.out=lo),col=8)
        matlines(1:N, vj$rval.sig[,i,], lty=1, lwd=2, col=mycolors)
        matlines(1:N, vj$rlow.sig[,i,], lty=2, col=mycolors)
        matlines(1:N, vj$rupp.sig[,i,], lty=2, col=mycolors)
        mtext(mark, side=1, line=-1, adj=1, col=rgb(0,0,0,.1),cex=.2)
        col <- (vj$rord[,i,]<=3)*1 +(vj$rord[,i,]>3)*8
        # xvar <- seq(1,N,M)
        xvar <- t(t(which(abs(diff(sign(diff(vj$rval[,i,]))))==2,arr.ind=TRUE))+c(1,0))
        text(xvar, vj$rval[xvar,i,2], labels=reg.vars[xvar,2], col=col,cex=.3)
        text(xvar, vj$rval[xvar,i,], labels=vj$rord[xvar,i,],pos=1, col=col,cex=.3)
        if (length(unique(YmaxR[[j]]))==1) {
          mtext(xxnames[YmaxR[[j]]][1], at=1, side=3, line=-1, cex=.5)
        } else {
          xvaru <- t(t(which(diff(sign(diff(as.matrix(vj$rval[,i,]))))==-2))+1)
          xvarl <- t(t(which(diff(sign(diff(as.matrix(vj$rval[,i,]))))==2))+1)
          # xvar <- t(t(which(abs(diff(sign(diff(as.matrix(vj$rval[,i])))))==2))+1)
          # xvar2 <- xvar[-length(xvar)]+diff(xvar)/2
          mtext(xxnames[YmaxR[[j]]][xvaru], at=xvaru, side=3, line=-1, cex=.3)
          # mtext(xxnames[YmaxR[[j]]][xvar2], at=xvar2, side=3, line=-1, cex=.3)
          mtext(xxnames[YmaxR[[j]]][xvarl], at=xvarl, side=3, line=-1, cex=.3)
        }
        if (xaxt!="s") axis(side=1, at=at, labels=label)
      }
      par(las=0)
      mtext('time', side=1, outer=TRUE, adj=0.5)
      mtext('Wavelet Local Multiple Cross-Regression', side=2, outer=TRUE, adj=0.5)
      mtext(level.lab[j], side=3, outer=TRUE, adj=0.5)
      if (!is.null(pdf.write)) dev.off()
    }
    return()
  }
