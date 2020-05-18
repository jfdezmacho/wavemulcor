wave.multiple.correlation <- #2.2.3
function(xx, N=length(xx[[1]][[1]]), p = .975, ymaxr=NULL) {
  sum.of.squares <- function(x) { sum(x^2, na.rm=TRUE) / sum(!is.na(x)) }
  sum.of.not.squares <- function(x) { sum(x, na.rm=TRUE) / sum(!is.na(x)) }

  d <- length(xx)       #number of series
  dd <- d*(d-1)/2       #number of correlations
  l <- length(xx[[1]])  #number of scales J+1 (wavelet coefficients at levels 1 to J plus the scaling coeffs at level J+1)
  # N <- length(xx[[1]][[1]]) #number of observations
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

  xy.mulcor <- vector("numeric", l) ##xy.mulcor <- matrix(NA,l,1,dimnames=list(names(xy.cor[[1]])))
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
    xy.mulcor[[i]] <- sgnr*sqrt(1-1/Pidiag[Pimax]) ## max(sqrt(1-1/diag(solve(P))))
  }

  n <- trunc(N/2^(1:l))
  sqrtn <- sapply(n,function(x){if (x>=3) sqrt(x-3) else NaN})
  alow <- atanh(xy.mulcor)-qnorm(p)/sqrtn
  if (dd>1) alow <- pmax(alow,0) ## wavemulcor can only be negative in bivariate case 
  aupp <- atanh(xy.mulcor)+qnorm(p)/sqrtn
  out <- data.frame( wavemulcor=xy.mulcor, lower=tanh(alow), upper=tanh(aupp) )
  Lst <- list(xy.mulcor=as.matrix(out),YmaxR=YmaxR)
  return(Lst)
}

