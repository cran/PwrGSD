\name{mysurvfit}
\alias{mysurvfit}

\title{My Survfit}
\description{
  Computes numbers at risk, numbers of events at each unique event
  time within levels of a blocking factor
}
\usage{
  mysurvfit(formula = formula(data), data = parent.frame(), subset, na.action = na.fail)
}

\arguments{
  \item{formula}{Should be a formula of the form \code{Surv(ti, ev) ~ block}
  where \code{block} is the blocking factor. It need not be a 
  factor per se but should have relatively few discrete levels.  Sorry,
  no staggered entry allowed at present}
  \item{data}{a dataframe}
  \item{subset}{you can subset the analysis via logical expression
    in variables in the dataframe}
  \item{na.action}{pass a method for handling \code{NA} values in
    \code{block} such as \code{na.omit}, or \code{na.fail}}
}
\value{
  A dataframe of \code{2*NLEV + 1} columns where \code{NLEV} is the
  number of levels of the factor \code{block}.
  \item{time}{The sorted vector of unique event times from all blocks}
  \item{nrisk1}{The number at risk in block level 1 at each event time}
  \item{nevent1}{The number of events in block level 1 at each event time}
  \item{...}{}
  \item{nriskNLEV}{The number at risk in block level NLEV at each event time}
  \item{neventNLEV}{The number of events in block level NLEV at each event time}
}
\author{Grant Izmirlian <izmirlian@nih.gov>}

\examples{
  library(PwrGSD)
  data(lung)

  fit.msf <- mysurvfit(Surv(time, I(status==2)) ~ sex, data=lung)

  fit.msf
  \dontrun{
  plot(fit.msf)
  }
}
\keyword{survival}
