test_that("wave.multiple.* works", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  # M <- 10
  # window <- "gauss"
  J <- trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  ##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  # ---------------------------
  
  Lst <- wave.multiple.correlation(xx) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.multiple.regression(xx) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.multiple.correlation(xx, ymaxr=1)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.multiple.regression(xx, ymaxr=1)
  expect_true(is.list(Lst))
  
})


test_that("plot_WMC and plot_WMR work with ymax=NULL", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  # M <- 10
  # window <- "gauss"
  J <- trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.multiple.regression(xx) #, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing correlation plot
  expect_null(
    plot_wave.multiple.correlation(Lst)
  )
  expect_null(
    plot_wave.multiple.correlation(Lst)
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.multiple.regression(Lst, nsig=2)
  )
  
})

test_that("plot_WMC and plot_WMR work with ymax=1", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  # M <- 10
  # window <- "gauss"
  J <- trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.multiple.regression(xx, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing correlation plot
  expect_null(
    plot_wave.multiple.correlation(Lst)
  )
  expect_null(
    plot_wave.multiple.correlation(Lst)
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.multiple.regression(Lst, nsig=2)
  )
  
})