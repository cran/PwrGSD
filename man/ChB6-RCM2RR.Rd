\name{RCM2RR}
\alias{RCM2RR}
\title{Relative cumulative mortality to Relative Risk}
\description{
  Given the relative cumulative mortality (ratio of CDFs),
  the baseline hazard and censoring hazard at a grid of time points,
  calculates the corresponding risk ratio at a second specified
  grid of time points.
}
\usage{
  RCM2RR(tlook, tcut.i, h.i, hOth, accru, rcm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tlook}{Second grid of time points at which you desire risk ratios}
  \item{tcut.i}{First grid of time points at which baseline hazard,
    censoring hazard and relative cumulative mortality are specified
    (left hand endpoints of intervals)}
  \item{h.i}{Values of baseline hazard on intervals given by tcut.i}
  \item{hOth}{Values of censoring hazard on intervals given by tcut.i}
  \item{accru}{Time at which uniform accrual is completed (starts at 0)}
  \item{rcm}{Values of relative cumulative mortality (ratio of CDFs) on
    interals given by tcut.i}
}

\value{
  Values of risk ratio on intervals given by tlook
}
\author{Grant Izmirlian <izmirlian@nih.gov>}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
