\name{xrand}
\alias{xrand}
\docType{data}
\title{
Correlation structural breaks data
}
\description{
Simulated data showing correlation structural breaks in Figure 4  of Fernández-Macho (2018).
}
\usage{data("xrand")}
\format{
  A data frame with 512 observations on the following 2 variables.
  \describe{
    \item{\code{xrand1}}{a numeric vector}
    \item{\code{xrand2}}{a numeric vector}
  }
}
\details{
\eqn{xrand1[t]} and \eqn{xrand2[t]}
are highly correlated at low frequencies (long timescales) but uncorrelated at high frequencies (short timescales).
However, during a period of time spanning the second third of the sample (\eqn{T/3<t<2T/3}) that behavior is reversed so that data become highly correlated at short timescales but uncorrelated at low frequencies.
}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Fernández-Macho, J., 2018. Time-localized wavelet multiple regression and correlation, Physica A, 492, https://doi.org/10.1016/j.physa.2017.11.050
}
\examples{
data(xrand)
## maybe str(xrand) ; plot(xrand) ...
}
\keyword{datasets}
