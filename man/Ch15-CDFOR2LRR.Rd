\name{CDFOR2LRR}
\alias{CDFOR2LRR}
\title{Convert CDF Odds Ratio to Logged Relative Risks}
\description{
  Given the values of the baseline hazard and odds ratio of the CDF at
  a grid of time points find the corresponding logged risk ratio.
}
\usage{
  CDFOR2LRR(tcut, tmax, h0, CDFOR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tcut}{Grid of time points (left endpoints)}
  \item{tmax}{The right endpoint of the last interval}
  \item{h0}{Values of the baseline hazard function on given intervals}
  \item{CDFOR}{Values of the odds ratio of the CDF's on the given intervals}
}
\value{
  An m by 2 matrix, where m=length(tcut), having columns 'tcut' and
  logged RR.
}
\author{Grant Izmirlian <izmirlian@nih.gov}
\keyword{survival}
