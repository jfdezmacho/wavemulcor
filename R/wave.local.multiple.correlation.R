wave.local.multiple.correlation <- #2.2.2
  function(xx, M, window="gauss", p=.975, ymaxr=NULL) {
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

    val <- lo <- up <- YmaxR <- list()

    x <- df.swap.list(xx)

    cat("\nlev:")
    for(i in 1:l) { cat(sprintf("%s",i))
      out <- local.multiple.correlation(x[[i]], M, window=window, p=p, ymaxr=ymaxr)
      val[[i]] <- out$val
      lo[[i]] <-  out$lo
      up[[i]] <-  out$up
      YmaxR[[i]]<- out$YmaxR
    }
    val <- as.data.frame(val)
    lo <- as.data.frame(lo)
    up <- as.data.frame(up)
    YmaxR <- as.data.frame(YmaxR)
    names(val) <- names(lo) <- names(up) <- names(YmaxR) <- names(xx[[1]])
    Lst <- list(val=val,lo=lo,up=up,YmaxR=YmaxR)
    return(Lst)
  }
#--------------------------------------------------------------
