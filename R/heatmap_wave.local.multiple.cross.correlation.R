heatmap_wave.local.multiple.cross.correlation <- #3.1.1
  function(Lst, lmax, lag.first=FALSE, xaxt="s", ci=NULL, pdf.write=NULL) {
    ##Producing cross-correlation plot
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    cor <- Lst$cor
    lag.max <- trunc((ncol(cor$d1$vals)-1)/2)
    lag0 <- lag.max+1
    YmaxR <- Lst$YmaxR
    J <- length(Lst$YmaxR)-1
    level.labs <- c(paste("level",1:J),paste("s",J))
    scale.labs <- c(paste0(2^(1:J),"-",2^((1:J)+1)),"smooth")
    lag.labs <- c(paste("lead",lag.max:1),paste("lag",0:lag.max))
    if (is.null(ci)||ci=="center") {clim <- range(sapply(cor,function(x){x[["vals"]]}))
    } else if (ci=="lower") {clim <- range(sapply(cor,function(x){x[["lower"]]}))
    } else if (ci=="upper"){clim <- range(sapply(cor,function(x){x[["upper"]]}))}
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    if(lag.first){
      par(mfcol=c(lmax+1,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(0.1,1.2,1.2,1.2))
      if (!is.null(pdf.write))
        cairo_pdf(paste0("heat_",pdf.write,"_WLMCC_lags.pdf"), width=8.27,height=11.69)
      for(i in c(-lmax:0,lmax:1)+lag0) {
        # val <- sapply(cor,function(x){x[["vals"]][,i]})
        if (is.null(ci)||ci=="center") {val <- sapply(cor,function(x){x[["vals"]][,i]})
        } else if (ci=="lower") {val <- sapply(cor,function(x){x[["lower"]][,i]})
        } else if (ci=="upper"){val <- sapply(cor,function(x){x[["upper"]][,i]})}
        colkey <- FALSE
        if (i==lag0) colkey <- list(cex.axis=0.75)
        plot3D::image2D(z=val, x=1:nrow(val), y=1:ncol(val),
                        main=lag.labs[i], sub="", xlab="", ylab="",
                        xaxt=xaxt, yaxt="n", cex.axis=0.75,
                        colkey=colkey, clim=clim, clab=expression(varphi),
                        rasterImage=T, contour=list(lwd=2, col=plot3D::jet.col(11)))
        text(x=grconvertX(0.1,from="npc"), y=grconvertY(0.97,from="npc"), labels=mark,
             col=rgb(0,0,0,.1),cex=.2)
        par(las=0)
        if(xaxt!="s") {axis(side=1, at=at, labels=label, cex.axis=0.75)}
        axis(side=2, at = if(i<=lag0) 1:ncol(val) else FALSE,
             labels = if(i<=lag0) scale.labs else FALSE, las=1, cex.axis=0.75)
      }
      title(main='Wavelet Local Multiple Cross-Correlation', outer=TRUE)
      mtext("period", side=2, outer=TRUE, adj=0.5)
      if (!is.null(pdf.write)) dev.off()
    } else{
      par(mfrow=c(ceiling((J+1)/2),2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(0.1,1.2,1.2,1.2))
      if (!is.null(pdf.write))
        cairo_pdf(paste0("heat_",pdf.write,"_WLMCC_levels.pdf"), width=8.27,height=11.69)
      for(j in (J+1):1) {
        # val <- cor[[j]]$vals
        if (is.null(ci)||ci=="center") {val <- cor[[j]]$vals[,lag0+(-lmax:lmax)]
        } else if (ci=="lower") {val <- cor[[j]]$lower[,lag0+(-lmax:lmax)]
        } else if (ci=="upper"){val <- cor[[j]]$upper[,lag0+(-lmax:lmax)]}
        colkey <- FALSE
        if (j==1) colkey <- list(cex.axis=0.75)
        plot3D::image2D(z=val, x=1:nrow(val), y=1:ncol(val),
                        main=level.labs[j], sub="", xlab="", ylab="",
                        xaxt=xaxt, yaxt="n", cex.axis=0.75,
                        colkey=colkey, clim=clim, clab=expression(varphi),
                        rasterImage=T, contour=list(lwd=2, col=plot3D::jet.col(11)))
        text(x=grconvertX(0.1,from="npc"), y=grconvertY(0.97,from="npc"), labels=mark,
             col=rgb(0,0,0,.15),cex=.2)
        par(las=0)
        if(xaxt!="s") {axis(side=1, at=at, labels=label, cex.axis=0.75)}
        axis(side=2, at = if((J%%2==0&j%%2==1)||(J%%2==1&j%%2==0)) 1:ncol(val) else FALSE,
             labels = if((J%%2==0&j%%2==1)||(J%%2==1&j%%2==0)) lag.labs[lag0+(-lmax:lmax)] 
             else FALSE, las=1, cex.axis=0.75)
      }
      title(main='Wavelet Local Multiple Cross-Correlation', outer=TRUE)
      mtext("lead / lag", side=2, outer=TRUE, adj=0.5)
      if (!is.null(pdf.write)) dev.off()
    }
    return()
  }
