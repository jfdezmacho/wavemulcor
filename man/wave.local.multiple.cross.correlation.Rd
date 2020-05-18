\name{wave.local.multiple.cross.correlation}
\alias{wave.local.multiple.cross.correlation}
\alias{wlmcc}
\title{
Wavelet routine for local multiple cross-correlation
}
\description{
Produces an estimate of the multiscale local multiple cross-correlation
(as defined below) along with approximate confidence intervals.
}
\usage{
wave.local.multiple.cross.correlation(xx, M,
     window="gauss", lag.max=NULL, p=.975, ymaxr=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{A list of \eqn{n} (multiscaled) time series, usually the outcomes of dwt or modwt, \emph{i.e.}
 xx <- list(v1.modwt.bw, v2.modwt.bw, v3.modwt.bw)}
  \item{M}{length of the weight function or rolling window.}
  \item{window}{type of weight function or rolling window.
Six types are allowed, namely
the uniform window,
Cleveland or tricube window,
Epanechnikov or parabolic window,
Bartlett or triangular window,
Wendland window and
the gaussian window.
The letter case and length of the argument are not relevant
as long as at least the first four characters are entered.}
\item{lag.max}{
 maximum lag (and lead). If not set, it defaults to half the square root of the length of the
 original series.}
  \item{p}{one minus the two-sided p-value for the confidence interval, \emph{i.e.} the cdf value.}
  \item{ymaxr}{index number of the variable whose correlation is calculated against a linear combination
 of the rest, otherwise at each wavelet level wlmc chooses the one maximizing the multiple correlation.}
}
\details{
The routine calculates \eqn{J+1} sets of wavelet multiple cross-correlations,one per wavelet level, out of \eqn{n} variables, that can be plotted each as lags and leads time series plots.}
\value{
List of four elements:
\item{val: }{ list of \eqn{J+1} dataframes, each (rows = #observations, columns = #levels)
providing the point estimates for the wavelet local multiple correlation.}
\item{lo: }{ list of \eqn{J+1} dataframes, each (rows = #observations, columns = #lags and leads) providing the lower bounds from the confidence interval.}
\item{up: }{ list of \eqn{J+1} dataframes, each (rows = #observations, columns = #lags and leads) providing the upper bounds from the confidence interval.}
\item{YmaxR: }{ numeric matrix (rows = #observations, columns = #levels) giving, at each wavelet level and time, the index number of the variable whose correlation is calculated against
a linear combination of the rest.
By default, \emph{wlmc} chooses at each wavelet level and value in time the variable maximizing the multiple correlation.}
}
\references{
Fernández-Macho, J., 2018. Time-localized wavelet multiple regression and correlation, Physica A:
Statistical Mechanics, vol. 490, p. 1226--1238. <DOI:10.1016/j.physa.2017.11.050>
}
\author{
Javier Fernández-Macho,
Dpt. of Quantitative Methods,
University of the Basque Country, Agirre Lehendakari etorb. 83, E48015 BILBAO, Spain.
(email: javier.fernandezmacho at ehu.eus).
}
\note{
 Needs \emph{waveslim} package to calculate \emph{dwt} or \emph{modwt} coefficients as inputs to the routine (also for data in the example).
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

## Based on data from Figure 7.9 in Gencay, Selcuk and Whitcher (2001)
## plus one random series.


library(wavemulcor)

data(exchange)
returns <- diff(log(exchange))
returns <- ts(returns, start=1970, freq=12)
N <- dim(returns)[1]
wf <- "d4"
M <- 30
window <- "gauss"
J <- 3 #trunc(log2(N))-3
lmax <- 2

set.seed(140859)

demusd.modwt <- modwt(returns[,"DEM.USD"], wf, J)
demusd.modwt.bw <- brick.wall(demusd.modwt, wf)
jpyusd.modwt <- modwt(returns[,"JPY.USD"], wf, J)
jpyusd.modwt.bw <- brick.wall(jpyusd.modwt, wf)
rand.modwt <- modwt(rnorm(length(returns[,"DEM.USD"])), wf, J)
rand.modwt.bw <- brick.wall(rand.modwt, wf)

##xx <- list(demusd.modwt.bw, jpyusd.modwt.bw)
xx <- list(demusd.modwt.bw, jpyusd.modwt.bw, rand.modwt.bw)
names(xx) <- c("DEM.USD","JPY.USD","rand")

\dontrun{
# Note: WLMCR may take more than 10 seconds of CPU time on some systems

Lst <- wave.local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax)
val <- Lst$val
low.ci <- Lst$lo
upp.ci <- Lst$up
YmaxR <- Lst$YmaxR

# ---------------------------

##Producing cross-correlation plot

xvar <- seq(1,N,M)
level.lab <- c(paste("Level",1:J),paste("Smooth",J))
ymin <- -0.1
if (length(xx)<3) ymin <- -1
for(j in 1:(J+1)) {
  par(mfcol=c(lmax+1,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,1.2,0))
  # xaxt <- c(rep("n",lmax),"s",rep("n",lmax))
  for(i in c(-lmax:0,lmax:1)+lmax+1) {
    matplot(1:N,val[[j]][,i], type="l", lty=1, ylim=c(ymin,1), #xaxt=xaxt[i],
            xlab="", ylab="", main=paste("Lag",i-lmax-1))
    abline(h=0)              ##Add Straight horiz
    lines(low.ci[[j]][,i], lty=1, col=2) ##Add Connected Line Segments to a Plot
    lines(upp.ci[[j]][,i], lty=1, col=2)
    text(xvar,1, labels=names(xx)[YmaxR[[j]]][xvar], adj=0.25, cex=.8)
  }
  par(las=0)
  mtext('time', side=1, outer=TRUE, adj=0.5)
  mtext('Local Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5)
  mtext(level.lab[j], side=3, outer=TRUE, adj=0.5)
}

}
}
% Add one or more standard concepts, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ multivariate wavelet }
\concept{ local correlation }
\concept{ local regression }
\concept{ moving regression }
\concept{ wavelet cross-correlation }
