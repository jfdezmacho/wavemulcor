\name{wavemulcor-package}
\alias{wavemulcor-package}
\alias{wavemulcor}
\docType{package}
\title{
\packageTitle{wavemulcor}
}
\description{
\packageDescription{wavemulcor}
}
\details{
Wavelet routines that calculate single sets of wavelet multiple
regressions and correlations (WMR and WMC), and cross-regressions
and cross-correlations (WMCR and WMCC) from a multivariate time series.
Dynamic versions of the routines allow the wavelet local multiple
(cross-)regressions (WLMR and WLMCR) and (cross-)correlations
(WLMC and WLMCC) to evolve over time.
The output from these Wavelet statistics can later be plotted in single graphs,
as an alternative to trying to make sense out of several sets of wavelet
correlations or wavelet cross-correlations.
The code is based on the calculation, at each wavelet scale,
of the square root of the coefficient of determination in a
linear combination of variables for which such coefficient of
determination is a maximum. The code provided here is based on
the wave.correlation routine in Brandon Whitcher's waveslim R
package Version: 1.6.4, which in turn is based on wavelet
methodology developed in Percival and Walden (2000),
Gençay, Selçuk and Whitcher (2002) and others.
Version 2 incorporates wavelet local multiple correlations
(WLMC). These are like the previous global WMC but consisting
in one single set of multiscale correlations along time. That
is, at each time t, they are calculated by letting a window of
weighted wavelet coefficients around t move along time. Six
weight functions are provided. Namely, the uniform window,
Cleveland's tricube window, Epanechnikov's parabolic window,
Bartlett's triangular window and Wendland's truncated power
window and the Gaussian window. Version 2.2 incorporates
auxiliary functions that calculate local multiple correlations
and cross-correlations (LMC, LMCC). They are calculated by
letting move along time a window of weighted time series values
around t. Any of the six weight functions mentioned above can
be used. They also feed a new routine to compute wavelet local
multiple cross-correlation (WLMCC).
Version 3 extends all the previous correlation routines
(WMC, WMCC, LMC, WLMC, WLMCC) to handle wavelet regressions
(WMR, WMCR, LMR, WLMR, WLMCR) that provide regression coefficients
and statistics across wavelet scales.
Auxiliary plot_ and heatmap_ routines are also provided
to visualize the wavmulcor statistics.

}
\author{
Javier Fernández-Macho,
Dpt. of Quantitative Methods,
University of the Basque Country, Agirre Lehendakari etorb. 83, E48015 BILBAO, Spain.
(email: javier.fernandezmacho at ehu.eus).
}
\references{
Fernández-Macho, J., 2012. Wavelet multiple correlation and cross-correlation: A multiscale analysis of Eurozone stock markets. Physica A: Statistical Mechanics and its Applications 391, 1097--1104. <DOI:10.1016/j.physa.2011.11.002>

Fernández-Macho, J., 2018. Time-localized wavelet multiple regression and correlation, Physica A:
Statistical Mechanics, vol. 490, p. 1226--1238. <DOI:10.1016/j.physa.2017.11.050>
}
%%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\keyword{ ts }
\keyword{ multivariate }
\keyword{ htest }
\keyword{ regression }
