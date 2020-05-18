\name{plot_local.multiple.cross.regression}
\alias{plot_local.multiple.cross.regression}
\alias{plot_lmcr}
\title{
Auxiliary routine for plotting local multiple cross-regressions
}
\description{
Produces a plot of local multiple cross-regressions.
}
\usage{
plot_local.multiple.cross.regression(Lst, lmax, nsig=2, xaxt="s")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from local.multiple.cross.regression.}
    \item{lmax}{ maximum lag (and lead). }
  \item{nsig}{An optional value for the number of significant variables to plot_ Default is 2.}
  \item{xaxt}{An optional vector of labels for the "x" axis. Default is 1:n.}
}
\details{
The routine produces time series plots of local multiple cross-regressions
with their confidence interval for every lag and lead.
Also, at every upturn and downturn of the corresponding local multiple cross-correlation,
the name of the variable that maximizes that multiple correlation against the rest is shown on top.
The others are named ordered by significance when they are relevant.}
\value{Plot.}
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

% Add one or more standard concepts, see file 'KEYWORDS' in the
% R documentation directory.
\concept{ local correlation }
\concept{ local regression }
\concept{ moving regression }
\concept{ wavelet cross-correlation }
\concept{ wavelet cross-regression }
\concept{ plot }

