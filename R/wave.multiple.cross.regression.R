wave.multiple.cross.regression <- #3.0.0
  function(xx, lag.max=NULL, p = .975, ymaxr=NULL) {
    sum.of.squares <- function(x) { sum(x^2, na.rm=TRUE) / sum(!is.na(x)) }
    sum.of.not.squares <- function(x) { sum(x, na.rm=TRUE) / sum(!is.na(x)) }
    asnum <- function(x){setNames(as.numeric(x),names(x))}

    d <- length(xx)       #number of series
    dd <- d*(d-1)/2       #number of correlations
    l <- length(xx[[1]])  #number of scales J+1 (wavelet coefficients at levels 1 to J plus the scaling coeffs at level J+1)
    N <- length(xx[[1]][[1]]) #number of observations

    if(is.null(lag.max)) {lag.max=trunc(sqrt(length(xx[[1]][[l]]))/2)}
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
    rval <- rstd <- rlow <- rupp <- rtst <- rpva <- rord <- vector("list", l)
    # xy.mulreg <- lapply(vector("list", l), function(x) {vector("list",2*lm+1)})
    ## Note: [    1 2 ... lm-1 lm | lm+1 | lm+2 ... 2*lm 2*lm+1 ]
    ##       [ ...Pimax var leads |   0  | Pimax var lags...    ]
    YmaxR <- vector("numeric",l)
    for(i in 1:l) {
      xy.mulreg <- vector("list",2*lm+1)
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
      x0 <- sapply(xx[-Pimax],'[[',i)
      y0 <- sapply(xx[Pimax],'[[',i)
      depvar <- matrix(c(-1,0,Inf,0),1,4)
      if (is.null(names(xx))) row.names(depvar) <- "Y"
      else row.names(depvar) <- names(xx[Pimax])
      z0 <- summary(lm(formula = y0 ~ x0))$coefficients
      if(Pimax<nrow(z0)) z1 <- z0[(Pimax+1):nrow(z0),,drop=FALSE]
      else z1 <- NULL
      z0 <- rbind(z0[1:Pimax,,drop=FALSE], depvar, z1)
      if (is.null(names(xx))) row.names(z0)[1] <- "b0"
      else row.names(z0) <- c("b0",names(xx))
      xy.mulreg[[lm+1]] <- z0
      # xy.mulreg <- z[order(z[,4]),1:2]  #coefficients (and their stdvs) ordered from most to least significant
      ## lag=0: this must be same as in wave.multiple.regression
      if(lm>0) {
        x <- x0
        z <- y <- y0
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
          zjr <- lm_yx$coefficients
          if(Pimax<nrow(zjr)) z1 <- zjr[(Pimax+1):nrow(zjr),,drop=FALSE]
          else z1 <- NULL
          zjr <- rbind(zjr[1:Pimax,,drop=FALSE], depvar, z1)
          if (is.null(names(xx))) row.names(zjr)[1] <- "b0"
          else row.names(zjr) <- c("b0",names(xx))
          xy.mulreg[[lm+1+j]] <- zjr
          sgnr <- 1
          if (dd==1) sgnr <- sign(cor(z,x,use="complete.obs"))
          xy.mulcor[i,lm+1-j] <- sgnr*sqrt( lm_zx$r.squared )
          ## Note: var[Pimax] leads the others:       x[t]hat<--z[t-j]
          zjl <- lm_zx$coefficients
          if(Pimax<nrow(zjl)) z1 <- zjl[(Pimax+1):nrow(zjl),,drop=FALSE]
          else z1 <- NULL
          zjl <- rbind(zjl[1:Pimax,,drop=FALSE], depvar, z1)
          if (is.null(names(xx))) row.names(zjl)[1] <- "b0"
          else row.names(zjl) <- c("b0",names(xx))
          xy.mulreg[[lm+1-j]] <- zjl
        }}
      rval[[i]] <- lapply(xy.mulreg, function(x){x[,'Estimate']})
      rstd[[i]] <- lapply(xy.mulreg, function(x){x[,'Std. Error']})
      rlow[[i]] <- lapply(xy.mulreg, function(x){x[,'Estimate']-qt(p,N-d)*x[,'Std. Error']})
      rupp[[i]] <- lapply(xy.mulreg, function(x){x[,'Estimate']+qt(p,N-d)*x[,'Std. Error']})
      rtst[[i]] <- lapply(xy.mulreg, function(x){x[,'t value']})
      rpva[[i]] <- lapply(xy.mulreg, function(x){x[,'Pr(>|t|)']})
      rord[[i]] <- lapply(xy.mulreg, function(x){match(abs(x[,'t value']),sort(abs(x[,'t value']),decreasing=TRUE))})
      rval[[i]] <- t(sapply(rval[[i]],asnum))
      rstd[[i]] <- t(sapply(rstd[[i]],asnum))
      rlow[[i]] <- t(sapply(rlow[[i]],asnum))
      rupp[[i]] <- t(sapply(rupp[[i]],asnum))
      rtst[[i]] <- t(sapply(rtst[[i]],asnum))
      rpva[[i]] <- t(sapply(rpva[[i]],asnum))
      rord[[i]] <- t(sapply(rord[[i]],asnum))
    } # end for(i in 1:l) loop
    lags <- length(-lm:lm)
    oldw <- getOption("warn")
    options(warn = -1)
    sqrtn <- sqrt(matrix(trunc(N/2^(1:l)), nrow=l, ncol=lags) - 3)
    options(warn = oldw)
    alow <- atanh(xy.mulcor)-qnorm(p)/sqrtn
    if (dd>1) alow <- pmax(alow,0) ## wavemulcor can only be negative in bivariate case
    aupp <- atanh(xy.mulcor)+qnorm(p)/sqrtn
    ci.mulcor <- list( lower=tanh(alow), upper=tanh(aupp) )
    xy.mulreg <- list( rval=rval, rstd=rstd, rlow=rlow, rupp=rupp,
                    rtst=rtst, rord=rord, rpva=rpva )
    xy.mulreg <- lapply(xy.mulreg,setNames,names(xx[[1]]))
    #) #vars=t(sapply(rval,names)) --> names of variables: useful if somehow reordered
    Lst <- list(xy.mulcor=xy.mulcor,ci.mulcor=ci.mulcor,xy.mulreg=xy.mulreg,YmaxR=YmaxR,data=xx)
    return(Lst)
  }

