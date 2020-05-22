test_that("wave.multiple.cross.* works", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  # M <- 10
  # window <- "gauss"
  J <- trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  ##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  # ---------------------------
  
  Lst <- wave.multiple.cross.correlation(xx, lag.max=lmax) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.multiple.cross.regression(xx, lag.max=lmax) #, ymaxr=NULL) 
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.multiple.cross.correlation(xx, lag.max=lmax, ymaxr=1)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.multiple.cross.regression(xx, lag.max=lmax, ymaxr=1) 
  expect_true(is.list(Lst))
})


test_that("plot_WMCC and plot_WLMCR work with ymax=NULL", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  # M <- 10
  # window <- "gauss"
  J <- trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.multiple.cross.regression(xx, lag.max=lmax) #, ymaxr=1)
  
  # ---------------------------
  
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="center") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="lower") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="upper") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="other") #, pdf.write=NULL)
  )
  
  # ---------------------------
  
  ##Producing cross-correlation plot
  expect_null(
    plot_wave.multiple.cross.correlation(Lst, lmax)
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.multiple.cross.regression(Lst, lmax, nsig=2)
  )
  
  ##Producing cross-correlation plot
  expect_error(
    plot_wave.multiple.cross.correlation(Lst, lmax, xaxt="blah")
  )
  
  ##Producing cross-regression plot
  expect_error(
    plot_wave.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="blah")
  )
  
  # ---------------------------
  
})


test_that("plot_WMCC and plot_WLMCR work with ymax=1", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  # M <- 10
  # window <- "gauss"
  J <- 1 #trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.multiple.cross.regression(xx, lag.max=lmax, ymaxr=1)
  
  # ---------------------------
  
  ##Producing cross-correlation heat map
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="center") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="lower") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="upper") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.multiple.cross.correlation(Lst, lmax, ci="other") #, pdf.write=NULL)
  )
  
  # ---------------------------
  
  ##Producing cross-correlation plot
  expect_null(
    plot_wave.multiple.cross.correlation(Lst, lmax)
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_wave.multiple.cross.regression(Lst, lmax, nsig=2)
  )
  
  ##Producing cross-correlation plot
  expect_error(
    plot_wave.multiple.cross.correlation(Lst, lmax, xaxt="blah")
  )
  
  ##Producing cross-regression plot
  expect_error(
    plot_wave.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="blah")
  )
  
  # ---------------------------
  
})