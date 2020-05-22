test_that("wave.local.multiple.cross.* works", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  M <- 10
  window <- "gauss"
  J <- 3 #trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  ##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  # ---------------------------
  
  Lst <- wave.local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.local.multiple.cross.regression(xx, M, window=window, lag.max=lmax) #, ymaxr=NULL) 
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax, ymaxr=1)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.local.multiple.cross.regression(xx, M, window=window, lag.max=lmax, ymaxr=1) 
  expect_true(is.list(Lst))
})


test_that("plot_WLMCC and plot_WLMCR work with ymax=NULL", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  M <- 10
  window <- "gauss"
  J <- 3 #trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  
  # Note: WLMCR may take more than 10 seconds of CPU time on some systems
  
  Lst <- wave.local.multiple.cross.regression(xx, M, window=window, lag.max=lmax) #, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE) #, xaxt="s", ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE) #, xaxt="s", ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other") #, xaxt="s") #, pdf.write=NULL)
  )
  
  ##Producing cross-correlation plot
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE) #, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE) #, xaxt="s")
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, xaxt="s") #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, xaxt="s") #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper", xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other", xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other", xaxt="s") #, pdf.write=NULL)
  )

  ##Producing cross-correlation plot
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE, xaxt="s")
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, xaxt=xaxt) #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, xaxt=xaxt) #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other", xaxt=xaxt) #, pdf.write=NULL)
  )

  ##Producing cross-correlation plot
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE, xaxt=xaxt)
  )
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE, xaxt=xaxt)
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt=xaxt)
  )
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation heat map
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE,xaxt="blah") #, ci=null) #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE,xaxt="blah") #, ci=null) #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other",xaxt="blah") #, pdf.write=null)
  )

  ##Producing cross-correlation plot
  expect_error(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE, xaxt="blah")
  )
  expect_error(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE, xaxt="blah")
  )
  
  ##Producing cross-regression plot
  expect_error(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="blah")
  )
  
})

test_that("plot_WLMCC and plot_WLMCR work with ymax=1", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  M <- 10
  window <- "gauss"
  J <- 1 #trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  
  # Note: WLMCR may take more than 10 seconds of CPU time on some systems
  
  Lst <- wave.local.multiple.cross.regression(xx, M, window=window, lag.max=lmax, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE) #, xaxt="s", ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE) #, xaxt="s", ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other") #, xaxt="s") #, pdf.write=NULL)
  )
  
  ##Producing cross-correlation plot
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE) #, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE) #, xaxt="s")
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, xaxt="s") #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, xaxt="s") #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper", xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other", xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other", xaxt="s") #, pdf.write=NULL)
  )
  
  ##Producing cross-correlation plot
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE, xaxt="s")
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, xaxt=xaxt) #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, xaxt=xaxt) #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other", xaxt=xaxt) #, pdf.write=NULL)
  )
  
  ##Producing cross-correlation plot
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE, xaxt=xaxt)
  )
  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE, xaxt=xaxt)
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt=xaxt)
  )
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation heat map
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE,xaxt="blah") #, ci=null) #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE,xaxt="blah") #, ci=null) #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="center",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="center",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="lower",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="lower",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="upper",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="upper",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=FALSE, ci="other",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
                                                  lag.first=TRUE, ci="other",xaxt="blah") #, pdf.write=null)
  )
  
  ##Producing cross-correlation plot
  expect_error(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE, xaxt="blah")
  )
  expect_error(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=TRUE, xaxt="blah")
  )
  
  ##Producing cross-regression plot
  expect_error(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="blah")
  )
  
})