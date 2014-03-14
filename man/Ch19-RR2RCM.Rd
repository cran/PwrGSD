\name{RR2RCM}
\alias{RR2RCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Relative risk to Relative Cumulative Mortality}
\description{
  Relative risk to Relative Cumulative Mortality
}
\usage{
RR2RCM(tlook, tcut.i, tcut.ii, h, rr, hOth, accru)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tlook}{Grid of time points at which you desire cumlative
    relative mortality}
  \item{tcut.i}{Grid of time points at which baseline hazard,
    censoring hazard and relative cumulative mortality are specified
    (left hand endpoints of intervals)}
  \item{tcut.ii}{Grid of time points at which study arm hazard
    is specified (left hand endpoints of intervals)}
  \item{h}{Values of baseline hazard on intervals given by tcut.i}
  \item{rr}{Values of risk ratio on intervals given by tcut.i}
  \item{hOth}{Values of censoring hazard on intervals given by tcut.i}
  \item{accru}{Time at which uniform accrual is completed (starts at 0)}
}

\value{
  Values of relative cumulative mortality (ratio of CDFs) on interals
  given by tlook
}
\author{Grant Izmirlian <izmirlian@nih.gov}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
