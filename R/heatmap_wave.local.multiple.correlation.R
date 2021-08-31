heatmap_wave.local.multiple.correlation <- #3.1.0.
  function(Lst, xaxt="s", ci=NULL, pdf.write=NULL) {
  ##Producing heat map
    if (xaxt[1]!="s"){
      at <- xaxt[[1]]
      label <- xaxt[[2]]
      xaxt <- "n"
    }
    cor <- Lst$cor
    J <- length(Lst$YmaxR)-1
    par(mfcol=c(1,1), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(0.1,1.2,1.2,1.2))
    scale.labs <- c(paste0(2^(1:J),"-",2^((1:J)+1)),"smooth")
    if (!is.null(pdf.write))
      cairo_pdf(paste0("heat_",pdf.write,"_WLMC.pdf"), width=11.69,height=8.27)
    if (is.null(ci)||ci=="center") {val <- sapply(cor,function(x){x[["val"]]})
    } else if (ci=="lower") {val <- sapply(cor,function(x){x[["lo"]]})
    } else if (ci=="upper"){val <- sapply(cor,function(x){x[["up"]]})}
    clim <- range(sapply(cor,function(x){x[["val"]]}))
    mark <- paste0("\u00A9jfm-wavemulcor3.1.0_",Sys.time()," ")
    plot3D::image2D(z=val, x=1:nrow(val), y=1:ncol(val),
            main="", sub="", xlab="", ylab="", 
            xaxt=xaxt, yaxt="n", cex.axis=0.75,
            colkey=list(cex.axis=0.75), clim=clim, clab=expression(varphi),
            rasterImage=TRUE, contour=list(lwd=2, col=plot3D::jet.col(11)))
    if(xaxt!="s") {axis(side=1, at=at, labels=label, cex.axis=0.75)}
    axis(side=2, at=1:ncol(val),labels=scale.labs, las=1, cex.axis=0.75)
    text(x=grconvertX(0.1,from="npc"), y=grconvertY(0.98,from="npc"), labels=mark,
         col=rgb(0,0,0,.1),cex=.2)
    par(las=0)
    title(main='Wavelet Local Multiple Correlation', outer=TRUE)
    mtext("period", side=2, outer=TRUE, adj=0.5)
    if (!is.null(pdf.write)) dev.off()
    return()
  }
