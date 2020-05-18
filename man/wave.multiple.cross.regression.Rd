\name{wave.multiple.cross.regression}
\alias{wave.multiple.cross.regression}
\alias{wmcr}
\title{
 Wavelet routine for multiple cross-regression
}
\description{
 Produces an estimate of the multiscale multiple cross-regression
 (as defined below).
}
\usage{
 wave.multiple.cross.regression(xx, lag.max=NULL, p = .975, ymaxr=NULL)
}
\arguments{
  \item{xx}{ A list of \eqn{n} (multiscaled) time series, usually the outcomes of dwt or modwt, \emph{i.e.} xx <- list(v1.modwt.bw, v2.modwt.bw, v3.modwt.bw)}
  \item{lag.max}{ maximum lag (and lead). If not set, it defaults to half the square root of the length of the original series.}
  \item{p}{one minus the two-sided p-value for the confidence interval, \emph{i.e.} the cdf value.}
  \item{ymaxr}{ index number of the variable whose correlation is calculated against a linear combination of the rest, otherwise at each wavelet level wmc chooses the one maximizing the multiple correlation.}
}
\details{
 The routine calculates one single set of wavelet multiple cross-regressions out of \eqn{n} variables that can be plotted as one single set of graphs (one per wavelet level).
}
\value{
 List of four elements:
 \item{xy.mulcor}{ List of three elements:}
 \itemize{
 \item{wavemulcor:}{ numeric matrix (rows = #levels, #cols = #lags and leads) with as many rows as levels in the wavelet transform object.
 The columns provide the point estimates for the wavelet multiple cross-correlations at different lags (and leads).
 The central column (lag=0) replicates the wavelet multiple correlations.
 Columns to the right (lag>0) give wavelet multiple cross-correlations with positive lag,
 \emph{i.e.} with y=var[Pimax] lagging behind a linear combination of the rest: x[t]hat --> y[t+j].
 Columns to the left (lag<0) give wavelet multiple cross-correlations with negative lag,
 \emph{i.e.} with y=var[Pimax] leading a linear combination of the rest: y[t-j] --> x[t]hat.}
 \item{lower: }{ numeric matrix (rows = #levels, #cols = #lags and leads) of lower bounds of the confidence interval.}
 \item{upper: }{ numeric matrix (rows = #levels, #cols = #lags and leads) of upper bounds of the confidence interval.}
 }
 \item{xy.mulreg:}{ List of seven elements:}
\itemize{
\item{rval: }{ numeric array (1stdim = #levels, 2nddim = #lags and leads, 3rddim = #regressors+1) of regression estimates.}
\item{rstd: }{ numeric array (1stdim = #levels, 2nddim = #lags and leads, 3rddim = #regressors+1) of their standard deviations.}
\item{rlow: }{ numeric array (1stdim = #levels, 2nddim = #lags and leads, 3rddim = #regressors+1) of their lower bounds.}
\item{rupp: }{ numeric array (1stdim = #levels, 2nddim = #lags and leads, 3rddim = #regressors+1) of their upper bounds.}
\item{rtst: }{ numeric array (1stdim = #levels, 2nddim = #lags and leads, 3rddim = #regressors+1) of their t statistic values.}
\item{rord: }{ numeric array (1stdim = #levels, 2nddim = #lags and leads, 3rddim = #regressors+1) of their index order when sorted by significance.}
\item{rpva: }{ numeric array (1stdim = #levels, 2nddim = #lags and leads, 3rddim = #regressors+1) of their p values.}
}
 \item{YmaxR:}{ numeric vector giving, at each wavelet level, the index number of the variable
 whose correlation is calculated against a linear combination of the rest. By default,
 \emph{wmcr} chooses at each wavelet level the variable maximizing the multiple correlation.}
 \item{data: }{ dataframe (rows = #levels, cols = #regressors) of original data.}
 }
\references{
 Fernández-Macho, J., 2018. Time-localized wavelet multiple regression and correlation, Physica A:
Statistical Mechanics, vol. 490, p. 1226--1238. <DOI:10.1016/j.physa.2017.11.050>
}

\author{
 Javier Fernández-Macho,
 Dpt. of Quantitative Methods,
 University of the Basque Country,
 Agirre Lehendakari etorb. 83, E48015 BILBAO, Spain. (email: javier.fernandezmacho at ehu.eus).
}
\note{
 Needs \emph{waveslim} package to calculate \emph{dwt} or \emph{modwt} coefficients as inputs to the routine (also for data in the example).
}
\examples{
## Based on data from Figure 7.9 in Gencay, Selcuk and Whitcher (2001)
## plus one random series.

library(wavemulcor)

data(exchange)
returns <- diff(log(exchange))
returns <- ts(returns, start=1970, freq=12)
N <- dim(returns)[1]

wf <- "d4"
J <- trunc(log2(N))-3
lmax <- 36

set.seed(140859)

demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)

# ---------------------------

xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
names(xx) <- c("DEM.USD","JPY.USD","rand")

Lst <- wave.multiple.cross.regression(xx, lmax)

# ---------------------------

##Producing correlation plot

plot_wave.multiple.cross.correlation(Lst, lmax) #, by=2)

##Producing regression plot

plot_wave.multiple.cross.regression(Lst, lmax) #, by=2)

}
\concept{ multivariate wavelet}
\concept{ wavelet cross-correlation }
\concept{ wavelet cross-regression }


