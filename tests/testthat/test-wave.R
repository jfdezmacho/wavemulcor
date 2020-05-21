test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("wave.local.multiple.cross.correlation works", {
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
  
  ##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax)
  expect_true(is.list(Lst))
})

test_that("wave.local.multiple.cross.regression works", {
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
  
  ##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.local.multiple.cross.regression(xx, M, window=window, lag.max=lmax) 
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  ##Producing cross-correlation plot

  expect_null(
    plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE) #, xaxt="s")
  )
  
  ##Producing cross-regression plot
  
  expect_null(
    plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2) #, xaxt="s")
  )
})
