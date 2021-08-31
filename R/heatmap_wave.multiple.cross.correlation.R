heatmap_wave.multiple.cross.correlation <- #3.1.0.
  function(Lst, lmax, by=3, ci=NULL, pdf.write=NULL) {
  ##Producing heat map
    if (is.null(ci)||ci=="center") {cor <- Lst$xy.mulcor
    } else if (ci=="lower") {cor <- Lst$ci.mulcor$lower
    } else if (ci=="upper"){cor <- Lst$ci.mulcor$upper}
    J <- length(Lst$YmaxR)-1
    par(mfcol=c(1,1), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(0.1,1.2,1.2,1.2))
    scale.labs <- c(paste0(2^(1:J),"-",2^((1:J)+1)),"smooth")
    if (!is.null(pdf.write))
      cairo_pdf(paste0("heat_",pdf.write,"_WMCC.pdf"), width=11.69,height=8.27)
    val <- t(cor) #sapply(cor,function(x){x[["val"]]})
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    plot3D::image2D(z=val, x=1:nrow(val), y=1:ncol(val),
            main="", sub="", xlab="", ylab="", cex.axis=0.75,
            xaxt="n", yaxt="n", clab = expression(varphi), colkey=list(cex.axis=0.75),
            rasterImage=TRUE, contour=list(lwd=2, col=plot3D::jet.col(11)))
    axis(side=1, at=seq(1,2*lmax+1,by=lmax/by), labels=seq(-lmax,lmax,by=lmax/by), cex.axis=0.75)
    axis(side=2, at=1:ncol(val),labels=scale.labs, las=1, cex.axis=0.75)
    text(x=grconvertX(0.1,from="npc"), y=grconvertY(0.98,from="npc"), labels=mark,
         col=rgb(0,0,0,.1),cex=.2)
    par(las=0)
    title(main='Wavelet Multiple Cross-Correlation', outer=TRUE)
    mtext('lag', side=1, outer=FALSE, adj=0.5)
    mtext("period", side=2, outer=TRUE, adj=0.5)
    if (!is.null(pdf.write)) dev.off()
    return()
  }
