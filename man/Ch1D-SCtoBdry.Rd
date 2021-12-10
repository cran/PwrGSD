\name{SCtoBdry}
\alias{SCtoBdry}
\title{Converts a stochastic curtailment boundary (conditional type I or
  II error probability) into a (efficacy or futility) boundary on the
  standardized Z scale}
\description{Converts a stochastic curtailment boundary (conditional type I or
  II error probability) into a (efficacy or futility) boundary on the
  standardized Z scale}
\usage{
SCtoBdry(prob, frac, be.end, drift = NULL, drift.end = NULL)
}
\arguments{
  \item{prob}{The stochastic curtailment thresh-hold probability, which
    is the complement of the type I (efficacy) or II (futility) error.
    We typically use 0.90 which will stop for efficacy if the
    probability under the null that the final analysis results in an
    efficacious decision given the data so far exceeds 0.90, and stops
    for futility of the probability under the alternative corresponding
    to the drift arguments, that the final analysis results in a
    futility decision given the data so far, exceeds 0.90. } 
  \item{frac}{The variance ratio. See the \code{\link{GrpSeqBnds}}
    documentation for details.}
  \item{be.end}{Value of efficacy (futility) boundary at the final analysis}
  \item{drift}{The drift function. See the \code{\link{GrpSeqBnds}}
    documentation for details.}
  \item{drift.end}{Required if using a futility boundary. This is the
    value of the drift function at the final analysis. Must be projected
    using the trial design.}
}

\value{A efficacy or futility boundary on the standard normal scale}
\author{Grant Izmirlian}
\examples{

  ## Here we show how to convert a stochastic curtailment procedure for
  ## futility into a futility boundary on the standard normal scale
  library(PwrGSD)

  ## Values of the information fraction at interim analyses --
  ## the sequence does not have to include the last analysis
  frac <- c(0.16, 0.32, 0.54, 0.83, 1.0)

  ## values drift at interim analyses corresponding to values of
  ## frac given above
  drift <- c(0.69, 1.09, 1.54, 2.08, 2.35)

  ## value of the drift at the final analysis (from the design or
  ## projected
  drift.end <- drift[5]

  ## value of the efficacy boundary at the final analysis 
  be.end <- 1.69

  ## stochastic curtailment threshhold probability -- if the probability of rejecting the
  ## null hypothesis by the scheduled end of the trial, under the alternative hypothesis,
  ## and conditional upon the current value of the statistic, is not greater than
  ## prob.thresh, then stop for futility.
  prob.thresh <- 0.90

  ## computes equivalent futility boundary points on the standard normal scale
  SCtoBdry(prob.thresh, frac=frac, be.end=be.end, drift=drift, drift.end=drift.end)
}

\keyword{htest}
\keyword{design}
