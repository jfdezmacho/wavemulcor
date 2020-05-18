\name{plot_local.multiple.correlation}
\alias{plot_local.multiple.correlation}
\alias{plot_lmc}
\title{
Auxiliary routine for plotting local multiple correlations
}
\description{
Produces a plot of local multiple correlations.
}
\usage{
plot_local.multiple.correlation(Lst, xaxt="s")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Lst}{A list from local.multiple.regression or local.multiple.correlation.}
  \item{xaxt}{An optional vector of labels for the "x" axis. Default is 1:n.}
}
\details{
The routine produces a time series plot of local multiple correlations with its confidence interval.
Also, at every upturn and downturn, the name of the variable that maximizes its multiple correlation against the rest is shown.
Note that the routine is optimize for \emph{local.multiple.regression}.
If you want to use output from the legacy \emph{local.multiple.correlation} function then you must create an empty list and put that output into a list element named \emph{cor} like this:
Lst <- list();
Lst$cor <- local.multiple.cross.correlation(xx, M, window=window, lag.max=lmax);
Lst$YmaxR <- Lst2$cor$YmaxR;
Lst$cor$YmaxR <- NULL.
}
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
\concept{ wavelet correlation }
\concept{ plot }
