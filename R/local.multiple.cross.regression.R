#----------------------------------------------------------
local.multiple.cross.regression <- #3.0.0
  function(xx, M, window="gauss", lag.max=NULL, p=.975, ymaxr=NULL) {
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
    N <- nrow(xx)         #number of observations
    if(is.null(lag.max)) {lag.max=trunc(sqrt(length(xx[[1]]))/2)}
    lm <- min(length(xx[[1]])-1, lag.max, na.rm=TRUE)
    val <- lo <- up <- matrix(nrow=N,ncol=2*lm+1)
    YmaxR <- matrix(nrow=N)
    rval<- rstd<- rlow<- rupp<- rtst<- rpva<- rord<- array(dim=c(N,2*lm+1,d+1))
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
      xy.mulcor <- matrix(nrow=2*lm+1)
      if (is.na(sum(r))||is.nan(sum(r))){ Pimax <- xy.mulcor[lm+1] <- NA
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
          xy.mulcor[lm+1] <- sgnr*sqrt(1-1/Pidiag[Pimax]) ## max(sqrt(1-1/diag(solve(P))))
          if (abs(xy.mulcor[lm+1])>1) browser()
          ## lag=0: this must be same as in local.multiple.correlation
        }}
      #}
      depvar <- matrix(c(-1,0,Inf,0),1,4)
      if (is.null(names(xx))) row.names(depvar) <- "Y"
      else row.names(depvar) <- names(xx[Pimax])
      x <- as.matrix(as.data.frame(x.w[-Pimax]))
      y <- x.w[[Pimax]]
      ## Note: var[Pimax] contemporaneous with the others: y[t]<--x[t]hat
      lm_yx <- summary(lm(formula = y ~ x))
      xy.mulreg <- lm_yx$coefficients
      if(Pimax<nrow(xy.mulreg)) xy.mulreg.r <- xy.mulreg[(Pimax+1):nrow(xy.mulreg),,drop=FALSE]
      else xy.mulreg.r<-NULL
      xy.mulreg <- rbind(xy.mulreg[1:Pimax,,drop=FALSE], depvar, xy.mulreg.r)
      if (is.null(names(xx))) row.names(z)[1] <- "b0"
      else row.names(xy.mulreg) <- c("b0",names(xx))
      rval[s,lm+1,] <- xy.mulreg[,'Estimate']
      rstd[s,lm+1,] <- xy.mulreg[,'Std. Error']
      rlow[s,lm+1,] <- xy.mulreg[,'Estimate']-qt(p,N-d)*xy.mulreg[,'Std. Error']
      rupp[s,lm+1,] <- xy.mulreg[,'Estimate']+qt(p,N-d)*xy.mulreg[,'Std. Error']
      rtst[s,lm+1,] <- xy.mulreg[,'t value']
      rpva[s,lm+1,] <- xy.mulreg[,'Pr(>|t|)']
      rord[s,lm+1,] <- match(abs(xy.mulreg[,'t value']),sort(abs(xy.mulreg[,'t value']),decreasing=TRUE))
      ## next 3 lines to check xy.mulcor above
      # sgnr <- 1
      # if (dd==1) sgnr <- sign(cor(y,x,use="complete.obs"))
      # xy.mulcor2[lm+1] <- sgnr*sqrt( lm_yx$r.squared )
      if(lm>0) {
        z <- y
        vlength <- length(y)
        for(j in 1:lm) {       ## now we obtain R2 of var[Pimax] with lagged values
          y <- c(y[2:vlength], NA)
          ## Note: var[Pimax] lags behind the others: y[t+j]<--x[t]hat
          lm_yx <- summary(lm(formula = y ~ x))
          sgnr <- 1
          if (dd==1) sgnr <- sign(cor(y,x,use="complete.obs"))
          xy.mulcor[lm+1+j] <- sgnr*sqrt( lm_yx$r.squared )
          xy.mulreg <- lm_yx$coefficients
          if(Pimax<nrow(xy.mulreg)) xy.mulreg.r <- xy.mulreg[(Pimax+1):nrow(xy.mulreg),,drop=FALSE]
          else xy.mulreg.r<-NULL
          xy.mulreg <- rbind(xy.mulreg[1:Pimax,,drop=FALSE], depvar, xy.mulreg.r)
          if (is.null(names(xx))) row.names(z)[1] <- "b0"
          else row.names(xy.mulreg) <- c("b0",names(xx))
          rval[s,lm+1+j,] <- xy.mulreg[,'Estimate']
          rstd[s,lm+1+j,] <- xy.mulreg[,'Std. Error']
          rlow[s,lm+1+j,] <- xy.mulreg[,'Estimate']-qt(p,N-d)*xy.mulreg[,'Std. Error']
          rupp[s,lm+1+j,] <- xy.mulreg[,'Estimate']+qt(p,N-d)*xy.mulreg[,'Std. Error']
          rtst[s,lm+1+j,] <- xy.mulreg[,'t value']
          rpva[s,lm+1+j,] <- xy.mulreg[,'Pr(>|t|)']
          rord[s,lm+1+j,] <- match(abs(xy.mulreg[,'t value']),sort(abs(xy.mulreg[,'t value']),decreasing=TRUE))
          z <- c(NA, z[1:(vlength-1)])
          ## Note: var[Pimax] leads the others:       x[t]hat<--z[t-j]
          lm_zx <- summary(lm(formula = z ~ x))
          sgnr <- 1
          if (dd==1) sgnr <- sign(cor(z,x,use="complete.obs"))
          xy.mulcor[lm+1-j] <- sgnr*sqrt( lm_zx$r.squared )
          xy.mulreg <- lm_zx$coefficients
          if(Pimax<nrow(xy.mulreg)) xy.mulreg.r <- xy.mulreg[(Pimax+1):nrow(xy.mulreg),,drop=FALSE]
          else xy.mulreg.r<-NULL
          xy.mulreg <- rbind(xy.mulreg[1:Pimax,,drop=FALSE], depvar, xy.mulreg.r)
          if (is.null(names(xx))) row.names(z)[1] <- "b0"
          else row.names(xy.mulreg) <- c("b0",names(xx))
          rval[s,lm+1-j,] <- xy.mulreg[,'Estimate']
          rstd[s,lm+1-j,] <- xy.mulreg[,'Std. Error']
          rlow[s,lm+1-j,] <- xy.mulreg[,'Estimate']-qt(p,N-d)*xy.mulreg[,'Std. Error']
          rupp[s,lm+1-j,] <- xy.mulreg[,'Estimate']+qt(p,N-d)*xy.mulreg[,'Std. Error']
          rtst[s,lm+1-j,] <- xy.mulreg[,'t value']
          rpva[s,lm+1-j,] <- xy.mulreg[,'Pr(>|t|)']
          rord[s,lm+1-j,] <- match(abs(xy.mulreg[,'t value']),sort(abs(xy.mulreg[,'t value']),decreasing=TRUE))
        }}
      swsq <- sum( weightedsq( !is.na(xx[[1]]) ,s,M) ,na.rm=TRUE)
      Nhat <- 1/swsq
      if (Nhat>=3) sqrtN <- sqrt(Nhat-3) else sqrtN <- NaN
      val[s,] <- xy.mulcor
      lo[s,]  <- tanh(atanh(xy.mulcor)-qnorm(p)/sqrtN)
      if (dd>1) lo[s,] <- pmax(lo[s,],0) ## wavemulcor can only be negative in bivariate case
      up[s,]  <- tanh(atanh(xy.mulcor)+qnorm(p)/sqrtN)
      YmaxR[s] <- Pimax
    } # end for (s in 1:N) loop
    dimnames(rval)[[3]] <- dimnames(rlow)[[3]] <- dimnames(rupp)[[3]] <- row.names(xy.mulreg)
    outcor <- list(vals=val,lower=lo,upper=up)
    outreg <- list( rval=rval, rstd=rstd,
                    rlow=rlow, rupp=rupp,
                    rtst=rtst, rord=rord, rpva=rpva )
    Lst <- list(cor=outcor,reg=outreg,YmaxR=YmaxR,data=xx)
    return(Lst)
  }
#--------------------------------------------------------------
