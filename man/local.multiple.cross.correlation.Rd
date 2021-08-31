\name{local.multiple.cross.correlation}
\alias{local.multiple.cross.correlation}
\alias{lmcc}
\title{
Routine for local multiple cross-correlation
}
\description{
Produces an estimate of local multiple cross-correlations
(as defined below) along with approximate confidence intervals.
}
\usage{
 local.multiple.cross.correlation(xx, M, window="gauss", lag.max=NULL, p=.975, ymaxr=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xx}{A list of \eqn{n} time series, \emph{e.g.}
 xx <- list(v1, v2, v3)}
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
\item{lag.max}{ maximum lag (and lead). If not set, it defaults to half the square root of the length of the original series.}
  \item{p}{one minus the two-sided p-value for the confidence interval, \emph{i.e.} the cdf value.}
  \item{ymaxr}{index number of the variable whose correlation is calculated against a linear combination of the rest, otherwise at each wavelet level lmc chooses the one maximizing the multiple correlation.}
}
\details{
The routine calculates a set of  time series of  multiple cross-correlations, one per lag and lead) out of \eqn{n} variables.
}
\value{
List of four elements:
\item{vals: }{ numeric matrix (rows = #observations, cols = #lags and leads) providing the point estimates for the local multiple cross-correlation.}
\item{lower: }{ numeric vmatrix (rows = #observations, cols = #lags and leads) providing the lower bounds from the confidence interval.}
\item{upper: }{ numeric matrix (rows = #observations, cols = #lags and leads) providing the upper bounds from the confidence interval.}
\item{YmaxR: }{ numeric matrix (rows = #observations, cols = #lags and leads) giving, at each value in time, the index number of the variable
whose correlation is calculated against a linear combination of the rest.
By default, \emph{lmcc} chooses at each value in time the variable maximizing the multiple correlation.}
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

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{

## Based on  Figure 4 showing correlation structural breaks in Fernandez-Macho (2018).

library(wavemulcor)

data(exchange)
returns <- diff(log(as.matrix(exchange)))
returns <- ts(returns, start=1970, freq=12)
N <- dim(returns)[1]

M <- 30
window <- "gauss"
lmax <- 1

demusd <- returns[,"DEM.USD"]
jpyusd <- returns[,"JPY.USD"]

set.seed(140859)

xrand <- rnorm(N)

xx <- data.frame(demusd, jpyusd, xrand)
##exchange.names <- c(colnames(returns), "RAND")

Lst <- local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax)
val <- Lst$vals
low.ci <- Lst$lower
upp.ci <- Lst$upper
YmaxR <- Lst$YmaxR

# ---------------------------

##Producing correlation plot

colnames(val) <- paste("Lag",-lmax:lmax)
xvar <- seq(1,N,M)
par(mfcol=c(lmax+1,2), las=1, pty="m", mar=c(2,3,1,0)+.1, oma=c(1.2,1.2,0,0))
ymin <- -0.1
if (length(xx)<3) ymin <- -1
for(i in c(-lmax:0,lmax:1)+lmax+1) {
  matplot(1:N,val[,i], type="l", lty=1, ylim=c(ymin,1), #xaxt="n",
          xlab="", ylab="", main=colnames(val)[i])
  # if(i==lmax+1) {axis(side=1, at=seq(0,N+50,50))}
  #axis(side=2, at=c(-.2, 0, .5, 1))
  abline(h=0)              ##Add Straight horiz
  lines(low.ci[,i], lty=1, col=2) ##Add Connected Line Segments to a Plot
  lines(upp.ci[,i], lty=1, col=2)
  text(xvar,1, labels=names(xx)[YmaxR][xvar], adj=0.25, cex=.8)
}
par(las=0)
mtext('time', side=1, outer=TRUE, adj=0.5)
mtext('Local Multiple Cross-Correlation', side=2, outer=TRUE, adj=0.5)

}
% Add one or more standard concepts, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ local correlation }
\concept{ local regression }
\concept{ moving regression }
\concept{ wavelet cross-correlation }
