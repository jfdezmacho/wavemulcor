test_that("local.multiple.* works", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  # wf <- "d4"
  M <- 10
  window <- "gauss"
  # J <- trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd <- returns[,"DEM.USD"]
  jpyusd <- returns[,"JPY.USD"]
  rand   <- rnorm(length(returns[,"DEM.USD"]))
  
  ##xx <- list(demusd, jpyusd)
  xx <- data.frame(demusd, jpyusd, rand)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  # ---------------------------
  
  Lst <- local.multiple.correlation(xx, M, window=window) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  Lst <- local.multiple.correlation(xx[,1:2], M, window=window) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- local.multiple.regression(xx, M, window=window) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  Lst <- local.multiple.regression(xx[,1:2], M, window=window) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- local.multiple.correlation(xx, M, window=window, ymaxr=1)
  expect_true(is.list(Lst))
  Lst <- local.multiple.correlation(xx[,1:2], M, window=window, ymaxr=1)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- local.multiple.regression(xx, M, window=window, ymaxr=1)
  expect_true(is.list(Lst))
  Lst <- local.multiple.regression(xx[,1:2], M, window=window, ymaxr=1)
  expect_true(is.list(Lst))
  
})


test_that("plot_LMC and plot_LMR work with ymax=NULL", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  M <- 10
  window <- "gauss"
  J <- trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd <- returns[,"DEM.USD"]
  jpyusd <- returns[,"JPY.USD"]
  rand   <- rnorm(length(returns[,"DEM.USD"]))
  
  ##xx <- list(demusd, jpyusd)
  xx <- data.frame(demusd, jpyusd, rand)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- local.multiple.regression(xx, M, window=window) #, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing correlation plot
  expect_null(
    plot_local.multiple.correlation(Lst) #, xaxt="s")
  )
  expect_null(
    plot_local.multiple.correlation(Lst) #, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_local.multiple.regression(Lst, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"
  ##Producing correlation plot
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt="s")
  )
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_local.multiple.regression(Lst, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 

  ##Producing correlation plot
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt=xaxt)
  )
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt=xaxt)
  )
  
  ##Producing regression plot
  expect_null(
    plot_local.multiple.regression(Lst, nsig=2, xaxt=xaxt)
  )
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing correlation plot
  expect_error(
    plot_local.multiple.correlation(Lst, xaxt="blah")
  )
  expect_error(
    plot_local.multiple.correlation(Lst, xaxt="blah")
  )
  
  ##Producing regression plot
  expect_error(
    plot_local.multiple.regression(Lst, nsig=2, xaxt="blah")
  )
  
})

test_that("plot_LMC and plot_LMR work with ymax=1", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  M <- 10
  window <- "gauss"
  J <- 1 #trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd <- returns[,"DEM.USD"]
  jpyusd <- returns[,"JPY.USD"]
  rand   <- rnorm(length(returns[,"DEM.USD"]))
  
  ##xx <- list(demusd, jpyusd)
  xx <- data.frame(demusd, jpyusd, rand)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- local.multiple.regression(xx, M, window=window, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing correlation plot
  expect_null(
    plot_local.multiple.correlation(Lst) #, xaxt="s")
  )
  expect_null(
    plot_local.multiple.correlation(Lst) #, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_local.multiple.regression(Lst, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"
  ##Producing correlation plot
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt="s")
  )
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_local.multiple.regression(Lst, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing correlation plot
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt=xaxt)
  )
  expect_null(
    plot_local.multiple.correlation(Lst, xaxt=xaxt)
  )
  
  ##Producing regression plot
  expect_null(
    plot_local.multiple.regression(Lst, nsig=2, xaxt=xaxt)
  )
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing correlation plot
  expect_error(
    plot_local.multiple.correlation(Lst, xaxt="blah")
  )
  expect_error(
    plot_local.multiple.correlation(Lst, xaxt="blah")
  )
  
  ##Producing regression plot
  expect_error(
    plot_local.multiple.regression(Lst, nsig=2, xaxt="blah")
  )
  
})