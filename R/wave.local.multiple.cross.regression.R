wave.local.multiple.cross.regression <- #3.0.0
  function(xx, M, window="gauss", lag.max=NULL, p=.975, ymaxr=NULL) {
    df.swap.list <- function(xx){
      yy <- list()
      for (i in seq_along(xx[[1]])){
        yy[[i]] <- as.data.frame( lapply(xx,'[[',i) )
        names(yy[[i]]) <- names(xx)
        attr(yy[[i]],"wave") <- attr(xx[[1]],"class")
      }
      names(yy) <- names(xx[[1]])
      return(yy)
    }
    d <- length(xx)       #number of series
    dd <- d*(d-1)/2       #number of correlations
    l <- length(xx[[1]])  #number of scales J+1 (wavelet coefficients at levels 1 to J plus the scaling coeffs at level J+1)
    # N <- length(xx[[1]][[1]]) #number of observations
    wav <- attr(xx[[1]],"wavelet")
    if(substr(wav,1,1)=="h") {L<-2
    } else if(substr(wav,1,1)=="d") {L<-as.numeric(substring(wav,2))
    } else if(substr(wav,1,1)=="l") {L<-as.numeric(substring(wav,3))} else L<-8

    ###IMPORTANT: must shift wave coeffs to left by L_i/2 at each level to be in phase with time series!!! (Getal, p.145)
    for(j in 1:d) {
      xx[[j]] <- phase.shift(xx[[j]], wav)
    }

    cor <- reg <- YmaxR <- list()

    x <- df.swap.list(xx)

    cat("\nlev:")
    for(i in 1:l) { cat(sprintf("%s",i))
      out <- local.multiple.cross.regression(x[[i]], M, window=window, lag.max=lag.max, p=p, ymaxr=ymaxr)
      cor[[i]] <- out$cor
      reg[[i]] <- out$reg
      YmaxR[[i]] <- out$YmaxR
    }
    YmaxR <- as.data.frame(YmaxR)
    names(cor) <- names(reg) <- names(YmaxR) <- names(xx[[1]])
    Lst <- list(cor=cor,reg=reg,YmaxR=YmaxR,data=xx)
    return(Lst)
  }
#--------------------------------------------------------------
