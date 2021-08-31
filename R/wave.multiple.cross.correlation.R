wave.multiple.cross.correlation <- #2.2.3
  function(xx, lag.max=NULL, p=.975, ymaxr=NULL) {
    sum.of.squares <- function(x) { sum(x^2, na.rm=TRUE) / sum(!is.na(x)) }
    sum.of.not.squares <- function(x) { sum(x, na.rm=TRUE) / sum(!is.na(x)) }
    d <- length(xx)       #number of series
    dd <- d*(d-1)/2       #number of correlations
    l <- length(xx[[1]])  #number of scales J+1 (wavelet coefficients at levels 1 to J plus the scaling coeffs at level J+1)
    N <- length(xx[[1]][[1]]) #number of observations
    if(is.null(lag.max)) {lag.max <- trunc(sqrt(length(xx[[1]][[l]]))/2)}
    lm <- min(length(xx[[1]][[l]])-1, lag.max, na.rm=TRUE)
    x.var <- vector("list", d)
    for(j in 1:d) {
      x.var[[j]] <- unlist(lapply(xx[[j]], sum.of.squares))
    }
    xy.cor <- vector("list", dd)
    xy <- vector("list", l)
    jk <- 0
    for(k in 1:(d-1)) {
      for(j in (k+1):d) {
        jk <- jk+1
        for(i in 1:l) {
          xy[[i]] <- as.vector(xx[[j]][[i]] * xx[[k]][[i]])
        }
        xy.cov <- unlist(lapply(xy, sum.of.not.squares))
        xy.cor[[jk]] <- xy.cov / sqrt(x.var[[j]] * x.var[[k]])
      }}
    xy.cor.vec <- matrix(unlist(xy.cor),l,dd)
    xy.mulcor <- matrix(0, l, 2*lm+1)
    ## Note: [    1 2 ... lm-1 lm | lm+1 | lm+2 ... 2*lm 2*lm+1 ]
    ##       [ ...Pimax var leads |   0  | Pimax var lags...    ]
    YmaxR <- vector("numeric",l)
    for(i in 1:l) {
      r <- xy.cor.vec[i,]
      P <- diag(d)/2
      P[lower.tri(P)] <- r
      P <- P+t(P)
      Pidiag <- diag(solve(P))
      if(is.null(ymaxr)) {
        YmaxR[i] <- Pimax <- which.max(Pidiag) ## detect i | x[i] on rest x gives max R2
      } else {YmaxR[i] <- Pimax <- ymaxr}
      sgnr <- 1
      if (dd==1) sgnr <- sign(r)
      xy.mulcor[i,lm+1] <- sgnr*sqrt(1-1/Pidiag[Pimax])
      ## lag=0: this must be same as in wave.multiple.correlation
      if(lm>0) {
        x <- sapply(xx[-Pimax],'[[',i)
        y <- sapply(xx[Pimax],'[[',i)
        z <- y
        vlength <- length(y)
        for(j in 1:lm) {       ## now we obtain R2 of var[Pimax] with lagged values
          y <- c(y[2:vlength], NA)
          z <- c(NA, z[1:(vlength-1)])
          lm_yx <- summary(lm(formula = y ~ x))
          lm_zx <- summary(lm(formula = z ~ x))
          sgnr <- 1
          if (dd==1) sgnr <- sign(cor(y,x,use="complete.obs"))
          xy.mulcor[i,lm+1+j] <- sgnr*sqrt( lm_yx$r.squared )
          ## Note: var[Pimax] lags behind the others: y[t+j]<--x[t]hat
          sgnr <- 1
          if (dd==1) sgnr <- sign(cor(z,x,use="complete.obs"))
          xy.mulcor[i,lm+1-j] <- sgnr*sqrt( lm_zx$r.squared )
          ## Note: var[Pimax] leads the others:       x[t]hat<--z[t-j]
        }}
    } # end of for(i in 1:l) loop
    lags <- length(-lm:lm)
    oldw <- getOption("warn")
    options(warn = -1)
    sqrtn <- sqrt(matrix(trunc(N/2^(1:l)), nrow=l, ncol=lags) - 3)
    options(warn = oldw)
    alow <- atanh(xy.mulcor)-qnorm(p)/sqrtn
    if (dd>1) alow <- pmax(alow,0) ## wavemulcor can only be negative in bivariate case
    aupp <- atanh(xy.mulcor)+qnorm(p)/sqrtn
    ci.mulcor <- list( lower=tanh(alow), upper=tanh(aupp) )
    Lst <- list(xy.mulcor=xy.mulcor,ci.mulcor=ci.mulcor,YmaxR=YmaxR)
    return(Lst)
  }

