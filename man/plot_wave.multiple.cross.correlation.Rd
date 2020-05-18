\name{plot_wave.multiple.cross.correlation}
\alias{plot_wave.multiple.cross.correlation}
\alias{plot_wmcc}
\title{
Auxiliary routine for plotting wave multiple cross-correlations
}
\description{
Produces a plot of wave multiple cross-correlations.
}
\usage{
plot_wave.multiple.cross.correlation(Lst, lmax, by=3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from wave.multiple.cross.regression or wave.multiple.cross.correlation.}
  \item{lmax}{ maximum lag (and lead). }
  \item{by}{ labels are printed every lmax/by. Default is 3. }
}
\details{
The routine produces a set of plots of wave multiple cross-correlations, one per wavelet level, 
with their confidence intervals.
Also, at each wavelet level, the name of the variable that maximizes its multiple 
cross-correlation against the rest is shown.}
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
\concept{ wavelet cross-correlation }
\concept{ plot }