test_that("wave.local.multiple.* works", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  M <- 10
  window <- "gauss"
  J <- 3 #trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  ##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  # ---------------------------
  
  Lst <- wave.local.multiple.correlation(xx, M, window=window) #, ymaxr=NULL)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.local.multiple.regression(xx, M, window=window) #, ymaxr=NULL)
  expect_true(is.list(Lst))

    # ---------------------------
  
  Lst <- wave.local.multiple.correlation(xx, M, window=window, ymaxr=1)
  expect_true(is.list(Lst))
  
  # ---------------------------
  
  Lst <- wave.local.multiple.regression(xx, M, window=window, ymaxr=1)
  expect_true(is.list(Lst))
  
})


test_that("plot_WLMC and plot_WLMR work with ymaxr=NULL", {
  data(exchange)
  returns <- diff(log(exchange))
  returns <- ts(returns, start=1970, freq=12)
  N <- dim(returns)[1]
  wf <- "d4"
  M <- 10
  window <- "gauss"
  J <- 3 #trunc(log2(N))-3
  # lmax <- 1
  
  set.seed(140859)
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.local.multiple.regression(xx, M, window=window) #, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing correlation heat map
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst) #, xaxt="s", ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="center") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other") #, xaxt="s") #, pdf.write=NULL)
  )

  ##Producing correlation plot
  expect_null(
    plot_wave.local.multiple.correlation(Lst) #, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.correlation(Lst) #, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.local.multiple.regression(Lst, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"
  ##Producing correlation heat map
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, xaxt="s") #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="center", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper", xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other", xaxt="s") #, pdf.write=NULL)
  )

  ##Producing correlation plot
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.local.multiple.regression(Lst, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing correlation heat map
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, xaxt=xaxt) #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="center", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other", xaxt=xaxt) #, pdf.write=NULL)
  )

  ##Producing correlation plot
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt=xaxt)
  )
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt=xaxt)
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.local.multiple.regression(Lst, nsig=2, xaxt=xaxt)
  )
  
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing correlation heat map
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst,xaxt="blah") #, ci=null) #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="center",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other",xaxt="blah") #, pdf.write=null)
  )

  ##Producing correlation plot
  expect_error(
    plot_wave.local.multiple.correlation(Lst, xaxt="blah")
  )
  expect_error(
    plot_wave.local.multiple.correlation(Lst, xaxt="blah")
  )
  
  ##Producing regression plot
  expect_error(
    plot_wave.local.multiple.regression(Lst, nsig=2, xaxt="blah")
  )
  
})

test_that("plot_WLMC and plot_WLMR work with ymaxr=1", {
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
  
  demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
  jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
  rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)
  
  xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
  names(xx) <- c("DEM.USD","JPY.USD","rand")
  
  Lst <- wave.local.multiple.regression(xx, M, window=window, ymaxr=1)
  
  # ---------------------------
  
  #xaxt NULL 
  ##Producing correlation heat map
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst) #, xaxt="s", ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="center") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper") #, xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other") #, xaxt="s") #, pdf.write=NULL)
  )
  
  ##Producing correlation plot
  expect_null(
    plot_wave.local.multiple.correlation(Lst) #, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.correlation(Lst) #, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.local.multiple.regression(Lst, nsig=2) #, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt ="s"
  ##Producing correlation heat map
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, xaxt="s") #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="center", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower", xaxt="s") #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper", xaxt="s") #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other", xaxt="s") #, pdf.write=NULL)
  )
  
  ##Producing correlation plot
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt="s")
  )
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt="s")
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.local.multiple.regression(Lst, nsig=2, xaxt="s")
  )
  
  # ---------------------------
  
  #xaxt a list(at=,label=)
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing correlation heat map
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, xaxt=xaxt) #, ci=NULL) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="center", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_null(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper", xaxt=xaxt) #, pdf.write=NULL)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other", xaxt=xaxt) #, pdf.write=NULL)
  )
  
  ##Producing correlation plot
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt=xaxt)
  )
  expect_null(
    plot_wave.local.multiple.correlation(Lst, xaxt=xaxt)
  )
  
  ##Producing regression plot
  expect_null(
    plot_wave.local.multiple.regression(Lst, nsig=2, xaxt=xaxt)
  )
  
  # ---------------------------
  
  #xaxt a list anything else
  xaxt <- list(at=seq(1,5),label=paste0(seq(1,5))) 
  
  ##Producing correlation heat map
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst,xaxt="blah") #, ci=null) #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="center",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="lower",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="upper",xaxt="blah") #, pdf.write=null)
  )
  expect_error(
    heatmap_wave.local.multiple.correlation(Lst, ci="other",xaxt="blah") #, pdf.write=null)
  )
  
  ##Producing correlation plot
  expect_error(
    plot_wave.local.multiple.correlation(Lst, xaxt="blah")
  )
  expect_error(
    plot_wave.local.multiple.correlation(Lst, xaxt="blah")
  )
  
  ##Producing regression plot
  expect_error(
    plot_wave.local.multiple.regression(Lst, nsig=2, xaxt="blah")
  )
  
})