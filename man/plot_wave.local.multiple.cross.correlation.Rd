\name{plot_wave.local.multiple.cross.correlation}
\alias{plot_wave.local.multiple.cross.correlation}
\alias{plot_wlmcc}
\title{
Auxiliary routine for plotting wave local multiple cross-correlations
}
\description{
Produces a plot of wave local multiple cross-correlations.
}
\usage{
plot_wave.local.multiple.cross.correlation(Lst, lmax,
     lag.first=FALSE, xaxt="s", pdf.write=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from local.multiple.cross.regression.}
  \item{lmax}{ maximum lag (and lead). }
  \item{lag.first}{ if TRUE, it produces lag-lead pages with \eqn{J+1} wavelet plots each.
  Otherwise (default) it gives wavelet pages with \eqn{2*lmax+1} lag-lead plots each. }
  \item{xaxt}{An optional vector of labels for the "x" axis. Default is 1:n.}
  \item{pdf.write}{ Optional name leader to save files to pdf format.
  The actual filename of each page is either
  "plot_<pdf.write>_WLMCC_<i>.pdf",  where "i" is the lag/lead \eqn{i=-lmax...+lmax}, or,
  "plot_<pdf.write>_WLMCC_<j>.pdf", where "j" is the wavelet level \eqn{j=1...(J+1)}.}
}
\details{
The routine produces time series plots of wave local multiple cross-correlations with
  their confidence intervals.
  Also, at every upturn and downturn, the name of the variable that maximizes its
  multiple cross-correlation against the rest is shown.}
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
\concept{ plot }
