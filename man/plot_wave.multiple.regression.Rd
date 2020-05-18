\name{plot_wave.multiple.regression}
\alias{plot_wave.multiple.regression}
\alias{plot_wmr}
\title{
Auxiliary routine for plotting wave multiple regressions
}
\description{
Produces a plot of wave multiple regressions.
}
\usage{
plot_wave.multiple.regression(Lst, nsig=2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from wave.multiple.regression.}
  \item{nsig}{An optional value for the number of significant variables to plot_ Default is 2.}
}
\details{
The routine produces a plot of wave multiple regressions with their confidence interval.
Also, the name of the variable that maximizes that multiple correlation against the rest is shown on top. 
The others are named with their order of significance when they are relevant.}
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
\concept{ multivariate wavelet}
\concept{ local correlation }
\concept{ local regression }
\concept{ moving regression }
\concept{ wavelet correlation }
\concept{ wavelet regression }
\concept{ plot }
