\name{heatmap_wave.local.multiple.cross.correlation}
\alias{heatmap_wave.local.multiple.cross.correlation}
\alias{heatmap_wlmcc}
\title{
Auxiliary routine for heatmaping wave local multiple cross-correlations
}
\description{
Produces heatmaps of wave local multiple cross-correlations.
}
\usage{
heatmap_wave.local.multiple.cross.correlation(Lst, lmax,
     lag.first=FALSE, xaxt="s", ci=NULL, pdf.write=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from wave.local.multiple.cross.regression.}
  \item{lmax}{ maximum lag (and lead). }
  \item{lag.first}{ if TRUE, it produces lag-lead pages with \eqn{J+1} wavelet heatmaps each.
  Otherwise (default) it gives wavelet pages with \eqn{2*lmax+1} lag-lead heatmaps each. }
  \item{xaxt}{An optional vector of labels for the "x" axis. Default is 1:n.}
  \item{ci}{ value to plot: "center" value of confidence interval (\emph{i.e.} the estimated cross-correlation),
  the  "lower" bound, or the "upper" bound. Default is "center".}
  \item{pdf.write}{ Optional name leader to save files to pdf format.
  The actual name of the file is either
  "heat_<pdf.write>_WLMCC_lags.pdf" or, "heat_<pdf.write>_WLMCC_levels.pdf".}
}
\details{
The routine produces a set of time series vs. wavelet periods heatmaps of wave local multiple
  cross-correlations at different lags and leads.}
\value{Heat map.}
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
\concept{ heat map }
