\name{heatmap_wave.multiple.cross.correlation}
\alias{heatmap_wave.multiple.cross.correlation}
\alias{heatmap_wmcc}
\title{
Auxiliary routine for heatmaping wave multiple cross-correlations
}
\description{
Produces heatmaps of wave multiple cross-correlations.
}
\usage{
heatmap_wave.multiple.cross.correlation(Lst, lmax, by=3, ci=NULL, pdf.write=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from wave.multiple.cross.regression or wave.multiple.cross.correlation.}
  \item{lmax}{ maximum lag (and lead). }
  \item{by}{ labels are printed every lmax/by. Default is 3. }
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
\concept{ heat map }
