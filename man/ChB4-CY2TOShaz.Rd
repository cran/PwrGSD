\name{CY2TOShaz}
\alias{CY2TOShaz}
\title{Calender year rates to Study Year Rates}
\description{
  Given the cutpoints at which the hazard is to be constant, the 
  values taken by the calender year rates and the calender time
  offset from the start of the trial at which randomization ended,
  this function converts to time on study rates, assuming uniform
  accrual.
}
\usage{
  CY2TOShaz(tcut, t.eor, m, verbose = FALSE)
}
\arguments{
  \item{tcut}{Left hand endpoints of intervals on which time on study
    hazard is taken to be constant}
  \item{t.eor}{Time offsest from the beginning of the trial at which
    randomization ended}
  \item{m}{Annual calender time rates}
  \item{verbose}{do you want to see alot of debugging info--defaults to
    \code{FALSE}}
}
\value{
  hazard = h, table = attr(obj., "tbl")
  \item{hazard}{time on study hazard values taken on intervals specified
    by the argument \code{tcut}}
  \item{table}{a table containg the observed and fitted values}
}
\author{Grant Izmirlian <izmirlian@nih.gov>}

\examples{
  ## none as yet
}
\keyword{survival}
