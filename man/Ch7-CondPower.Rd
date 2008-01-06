\name{CondPower}
\alias{CondPower}
\title{Conditional type I and type II error probabilities given current value of the test statistic}

\description{
  Computes conditional type I and type II error probabilities given current value of the test statistic for monitoring based
  upon stochastic curtailment. This is now obsolete and included in the functionality of \dQuote{GrpSeqBnds} and is here for
  instructional purposes only.
}
\usage{
CondPower(Z, frac, drift, drift.end, err.I, sided = 1)
}
\arguments{
  \item{Z}{Current value of test statistic standardized to unit variance.}
  \item{frac}{Current value of the information fraction (variance fraction).}
  \item{drift}{Current value of the drift, i.e. the expected value of the test statistic normalized to have variance equal to
  the information fraction. Required if you want to compute conditional type II error, otherwise enter 0.}
  \item{drift.end}{Projected value of the drift at the end of the trial.}
  \item{err.I}{Overall (total) type I error probability}
  \item{sided}{Enter 1 or 2 for sided-ness of the test.}
}
\value{
  A named numeric vector containing the two components \dQuote{Pr.cond.typeIerr} and \dQuote{Pr.cond.typeIIerr}
}
\references{ 
A General Theory on Stochastic Curtailment for Censored Survival Data
D. Y. Lin, Q. Yao, Zhiliang Ying
Journal of the American Statistical Association, Vol. 94, No. 446 (Jun., 1999), pp. 510-521
}
\author{Grant Izmirlian <izmirlian@nih.gov>}

\seealso{\code{\link{GrpSeqBnds}}}
\examples{
## None as yet
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{design}
\keyword{htest}
