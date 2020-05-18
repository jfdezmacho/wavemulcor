#----------------------------------------------------------
local.multiple.correlation <- #2.2.2
  function(xx, M, window="gauss", p=.975, ymaxr=NULL) {
    window <- substr(tolower(window),1,4)
    if (window=="unif"){
      #uniform window:
      weights <- function(z,M) {w<-z<=M; w/sum(w)}
    } else if (window=="clev"||window=="tric"){
      #Cleveland tricube window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M)^3)^3; w/sum(w)}
    } else if (window=="epan"||window=="para"){
      #Epanechnikov parabolic window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M)^2); w/sum(w)}
    } else if (window=="bart"||window=="tria"){
      #Bartlett triangular window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M)); w/sum(w)}
    } else if (window=="wend"){
      #Wendland 1995 window:
      weights <- function(z,M) {w<-z<=M; w<-w*(1-abs(z/M))^4*(1+4*abs(z/M)); w/sum(w)}
    } else if (window=="gaus"){
      #gauss window:
      weights <- function(z,M) {w<-1/exp((z/M)^2); w/sum(w)}
    } else {stop("wrong window type")}

    weighted <- function(x,s,M) {z<-abs( seq_along(x)-s ); w<-weights(z,M); w*x} #s=observation number; z=distance |t-s|
    weightedsq <- function(x,s,M) {z<-abs( seq_along(x)-s ); w<-weights(z,M); w^2*x}
    sum.of.squares <- function(x) { sum(x^2, na.rm=TRUE) / sum(!is.na(x)) }
    sum.of.not.squares <- function(x) { sum(x, na.rm=TRUE) / sum(!is.na(x)) }
    d <- length(xx)       #number of series
    dd <- d*(d-1)/2       #number of correlations
    N <- nrow(xx) #number of observations
    val <- lo <- up <- YmaxR <- matrix(nrow=N)
    for (s in 1:N) {
      x.w <- x.var <- vector("list", d)
      for(j in 1:d) {
        x.w[[j]] <- weighted( xx[[j]] ,s,M)  #xx.w[[j]][[i]][[s]]
        x.var[[j]] <- sum.of.squares(x.w[[j]])
      }
      xy.cor <- vector("list", dd)
      jk <- 0
      for(k in 1:(d-1)) {
        for(j in (k+1):d) {
          jk <- jk+1
          if (is.null(attr(xx,"wave"))) {
            xy.cor[[jk]] <- cor(x.w[[j]],x.w[[k]],use="complete.obs")
          } else {
            xy.w <- as.vector( x.w[[j]] * x.w[[k]] )
            xy.cov <- sum.of.not.squares(xy.w)
            xy.cor[[jk]] <- xy.cov / sqrt(x.var[[j]] * x.var[[k]])
          }
          if(sum(is.infinite(xy.cor[[jk]]))>0) browser()
        }}
      r <- unlist(xy.cor)
      if (is.na(sum(r))||is.nan(sum(r))){ Pimax <- xy.mulcor <- NA
      } else{
        P <- diag(d)/2
        P[lower.tri(P)] <- r
        P <- P+t(P)
        if (qr(P)$rank < d) {xy.mulcor <- 1; Pimax <- NA}
        else {
          Pidiag <- diag(solve(P))
          if(is.null(ymaxr)) {
            Pimax <- which.max(Pidiag) ## detect i | x[i] on rest x gives max R2
          } else {Pimax <- ymaxr}
          sgnr <- 1
          if (dd==1) sgnr <- sign(r)
          xy.mulcor <- sgnr*sqrt(1-1/Pidiag[Pimax]) ## max(sqrt(1-1/diag(solve(P))))
          if (abs(xy.mulcor)>1) browser()
        }}
      #}
      swsq <- sum( weightedsq( !is.na(xx[[1]]) ,s,M) ,na.rm=TRUE)
      Nhat <- 1/swsq
      if (Nhat>=3) sqrtN <- sqrt(Nhat-3) else sqrtN <- NaN
      val[s] <- xy.mulcor
      lo[s]  <- tanh(atanh(xy.mulcor)-qnorm(p)/sqrtN)
      if (dd>1) lo[s] <- pmax(lo[s],0) ## wavemulcor can only be negative in bivariate case
      up[s]  <- tanh(atanh(xy.mulcor)+qnorm(p)/sqrtN)
      YmaxR[s] <- Pimax
    } # end for (s in 1:N) loop
    # val <- as.data.frame(val)
    # lo <- as.data.frame(lo)
    # up <- as.data.frame(up)
    # YmaxR <- as.data.frame(YmaxR)
    # names(val) <- names(lo) <- names(up) <- names(YmaxR) <- names(xx[[1]])
    Lst <- list(val=val,lo=lo,up=up,YmaxR=YmaxR)
    return(Lst)
  }
#--------------------------------------------------------------
