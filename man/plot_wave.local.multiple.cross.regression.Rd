\name{plot_wave.local.multiple.cross.regression}
\alias{plot_wave.local.multiple.cross.regression}
\alias{plot_wlmcr}
\title{
Auxiliary routine for plotting wave local multiple cross-regressions
}
\description{
Produces a set of plots of wave local multiple cross-regressions.
}
\usage{
plot_wave.local.multiple.cross.regression(Lst, lmax, nsig=2, 
     xaxt="s", pdf.write=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from wave.multiplecross.regression.}
  \item{lmax}{ maximum lag (and lead). }
  \item{nsig}{ An optional value for the number of significant variables to plot_ Default is 2.}
  \item{xaxt}{ An optional vector of labels for the "x" axis. Default is 1:n.}
  \item{pdf.write}{ Optional name leader to save files to pdf format.
  The actual filename of each page  
  "plot_<pdf.write>_WLMCC_<j>.pdf", where "j" is the wavelet level \eqn{j=1...(J+1)}.}  
}
\details{
The routine produces \eqn{J+1} pages, one per wavelet level, each with time series plots of wave multiple cross-regressions at different lags and leads, each with their confidence interval.
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
\concept{ local correlation }
\concept{ local regression }
\concept{ moving regression }
\concept{ wavelet cross-correlation }
\concept{ wavelet cross-regression }
\concept{ plot }