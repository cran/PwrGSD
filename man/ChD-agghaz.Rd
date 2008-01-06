\name{agghaz}
\alias{agghaz}
\title{Aggregated Hazard}
\description{
  Computes the MLE for the model that assumes piecewise constant hazards
  on intervals defined by a grid of points. One applications for example 
  is to calculate monthly hazard rates given numbers of events, numbers at 
  risk and event times reported to the day. Can also handle time to event
  data stratified on a blocking factor.
}
\usage{
agghaz(t.agg, time, nrisk, nevent)
}
\arguments{
  \item{t.agg}{Vector defining intervals upon which the user wants
    constant hazard rates.}
  \item{time}{Event times, possibly stratified on a blocking factor into
    multiple columns, in units that occur in enough numbers per interval
    specified above. If there is just a single column then it must be in
    column form (see example below).} 
  \item{nrisk}{Numbers at risk at specified event times}
  \item{nevent}{Numbers of events at specified event times}
}
\value{
   \item{time.a}{User supplied left-hand endpoints of intervals of
     hazard constancy}
  \item{nrisk.a}{Numbers at risk on specified intervals}
  \item{nevent.a}{Numbers of events on specified intervals}
}
\author{Grant Izmirlian <izmirlian@nih.gov>}
\examples{
  library(PwrGSD)
  data(lung)
  fit.msf <- mysurvfit(Surv(time, I(status==2)) ~ sex, data=lung)

  ## A single stratum:
  with(fit.msf$Table, agghaz(30, time, cbind(nrisk1), cbind(nevent1)))

  ## Multiple strata--pooled and group 1:
  with(fit.msf$Table, agghaz(30, time, cbind(nrisk1+nrisk2,nrisk1), cbind(nevent1+nevent2,nevent1)))
}
\keyword{survival}
\keyword{manip}
