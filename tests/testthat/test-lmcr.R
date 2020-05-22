test_that("local.multiple.cross.* works", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  # wf <- "d4"
  M <- 10
  window <- "gauss"
  # J <- trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd <- returns[,"DEM.USD"]
  jpyusd <- returns[,"JPY.USD"]
  rand   <- rnorm(length(returns[,"DEM.USD"]))
  
  ##xx <- list(demusd, jpyusd)
  xx <- data.frame(demusd, jpyusd, rand)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  # ---------------------------
  
  Lst <- local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- local.multiple.cross.regression(xx, M, window=window, lag.max=lmax) #, ymaxr=NULL) 
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax, ymaxr=1)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- local.multiple.cross.regression(xx, M, window=window, lag.max=lmax, ymaxr=1) 
  expect_true(is.list(Lst))
})


test_that("plot_LMCC and plot_LMCR work with ymaxr=NULL", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  # wf <- "d4"
  M <- 10
  window <- "gauss"
  # J <- trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd <- returns[,"DEM.USD"]
  jpyusd <- returns[,"JPY.USD"]
  rand   <- rnorm(length(returns[,"DEM.USD"]))
  
  ##xx <- list(demusd, jpyusd)
  xx <- data.frame(demusd, jpyusd, rand)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  
  # Note: WLMCR may take more than 10 seconds of CPU time on some systems
  
  Lst <- local.multiple.cross.regression(xx, M, window=window, lag.max=lmax) #, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 

  ##Producing cross-correlation plot
  expect_null(
    plot_local.multiple.cross.correlation(Lst, lmax) #, xaxt="s")
  )

  ##Producing cross-regression plot
  expect_null(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"

  ##Producing cross-correlation plot
  expect_null(
    plot_local.multiple.cross.correlation(Lst, lmax, xaxt="s")
  )

  ##Producing cross-regression plot
  expect_null(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation plot
  expect_null(
    plot_local.multiple.cross.correlation(Lst, lmax, xaxt=xaxt)
  )

  ##Producing cross-regression plot
  expect_null(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt=xaxt)
  )
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation plot
  expect_error(
    plot_local.multiple.cross.correlation(Lst, lmax, xaxt="blah")
  )

  ##Producing cross-regression plot
  expect_error(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="blah")
  )
  
})

test_that("plot_LMCC and plot_LMCR work with ymaxr=1", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  # wf <- "d4"
  M <- 10
  window <- "gauss"
  # J <- 1 #trunc(log2(N))-3
  lmax <- 1
  
  set.seed(140859)
  
  demusd <- returns[,"DEM.USD"]
  jpyusd <- returns[,"JPY.USD"]
  rand   <- rnorm(length(returns[,"DEM.USD"]))
  
  ##xx <- list(demusd, jpyusd)
  xx <- data.frame(demusd, jpyusd, rand)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  
  # Note: WLMCR may take more than 10 seconds of CPU time on some systems
  
  Lst <- local.multiple.cross.regression(xx, M, window=window, lag.max=lmax, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  
  ##Producing cross-correlation plot
  expect_null(
    plot_local.multiple.cross.correlation(Lst, lmax) #, xaxt="s")
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"
  
  ##Producing cross-correlation plot
  expect_null(
    plot_local.multiple.cross.correlation(Lst, lmax, xaxt="s")
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation plot
  expect_null(
    plot_local.multiple.cross.correlation(Lst, lmax, xaxt=xaxt)
  )
  
  ##Producing cross-regression plot
  expect_null(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt=xaxt)
  )
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing cross-correlation plot
  expect_error(
    plot_local.multiple.cross.correlation(Lst, lmax, xaxt="blah")
  )
  
  ##Producing cross-regression plot
  expect_error(
    plot_local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="blah")
  )
  
})