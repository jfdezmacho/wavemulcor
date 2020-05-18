\name{plot_wave.multiple.correlation}
\alias{plot_wave.multiple.correlation}
\alias{plot_wmc}
\title{
Auxiliary routine for plotting wave multiple correlations
}
\description{
Produces a plot of wave multiple correlations.
}
\usage{
plot_wave.multiple.correlation(Lst)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from wave.multiple.regression or wave.multiple.correlation.}
}
\details{
The routine produces a plot of wave multiple correlations, at each wavelet level, with its confidence interval.
Also, at each wavelet level, the name of the variable that maximizes its multiple correlation against the rest is shown.}
\value{Plot.}
\references{
Fernández-Macho, J., 2012. Wavelet multiple correlation and cross-correlation: A multiscale analysis of Eurozone stock markets. Physica A: Statistical Mechanics and its Applications 391, 1097--1104. <DOI:10.1016/j.physa.2011.11.002>

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
\concept{ wave correlation }
\concept{ wave regression }
\concept{ moving regression }
\concept{ wavelet correlation }
\concept{ plot }