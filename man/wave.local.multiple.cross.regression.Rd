\name{wave.local.multiple.cross.regression}
\alias{wave.local.multiple.cross.regression}
\alias{wlmcr}
\title{
Wavelet routine for local multiple cross-regression
}
\description{
Produces an estimate of the multiscale local multiple cross-regression
(as defined below) along with approximate confidence intervals.
}
\usage{
wave.local.multiple.cross.regression(xx, M,
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
The routine calculates \eqn{J+1} sets of wavelet multiple cross-regressions, one per wavelet level, out of \eqn{n} variables, that can be plotted each as lags and leads time series plots.}
\value{
List of four elements:
\item{cor: }{ List of three elements:}
\itemize{
\item{val: }{ list of \eqn{J+1} dataframes, each (rows = #observations, columns = #levels)
providing the point estimates for the wavelet local multiple correlation.}
\item{lo: }{ list of \eqn{J+1} dataframes, each (rows = #observations, columns = #lags and leads) providing the lower bounds from the confidence interval.}
\item{up: }{ list of \eqn{J+1} dataframes, each (rows = #observations, columns = #lags and leads) providing the upper bounds from the confidence interval.}
}
\item{reg: }{ List of \eqn{J+1} elements, one per wavelet level, each with:}
\itemize{
\item{rval: }{ numeric array (1st_dim = #observations, 2nd-dim = #lags and leads, 3rd_dim = #regressors+1) of local regression estimates.}
\item{rstd: }{ numeric array (1st_dim = #observations, 2nd-dim = #lags and leads, 3rd_dim = #regressors+1) of their standard deviations.}
\item{rlow: }{ numeric array (1st_dim = #observations, 2nd-dim = #lags and leads, 3rd_dim = #regressors+1) of their lower bounds.}
\item{rupp: }{ numeric array (1st_dim = #observations, 2nd-dim = #lags and leads, 3rd_dim = #regressors+1) of their upper bounds.}
\item{rtst: }{ numeric array (1st_dim = #observations, 2nd-dim = #lags and leads, 3rd_dim = #regressors+1) of their t statistic values.}
\item{rord: }{ numeric array (1st_dim = #observations, 2nd-dim = #lags and leads, 3rd_dim = #regressors+1) of their index order when sorted by significance.}
\item{rpva: }{ numeric array (1st_dim = #observations, 2nd-dim = #lags and leads, 3rd_dim = #regressors+1) of their p values.}
}
\item{YmaxR: }{ dataframe (rows = #observations, columns = #levels) giving, at each wavelet level and time, the index number of the variable whose correlation is calculated against
a linear combination of the rest.
By default, \emph{wlmcr} chooses at each wavelet level and value in time the variable maximizing the multiple correlation.}
\item{data: }{ dataframe (rows = #observations, cols = #regressors) of original data.}
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

demusd.modwt <- brick.wall(modwt(returns[,"DEM.USD"], wf, J), wf)
jpyusd.modwt <- brick.wall(modwt(returns[,"JPY.USD"], wf, J), wf)
rand.modwt   <- brick.wall(modwt(rnorm(length(returns[,"DEM.USD"])), wf, J), wf)

xx <- list(demusd.modwt, jpyusd.modwt, rand.modwt)
names(xx) <- c("DEM.USD","JPY.USD","rand")

\dontrun{
# Note: WLMCR may take more than 10 seconds of CPU time on some systems

Lst <- wave.local.multiple.cross.regression(xx, M, window=window, lag.max=lmax) #, ymaxr=1)

# ---------------------------

##Producing cross-correlation plot

plot_wave.local.multiple.cross.correlation(Lst, lmax, lag.first=FALSE) #, xaxt="s")

##Producing cross-regression plot

plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2) #, xaxt="s")
}
}
% Add one or more standard concepts, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ multivariate wavelet}
\concept{ local correlation }
\concept{ local regression }
\concept{ moving regression }
\concept{ wavelet cross-correlation }
\concept{ wavelet cross-regression }
