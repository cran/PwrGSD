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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (prob, frac, be.end, drift = NULL, drift.end = NULL) 
{
    ef <- !missing(drift)
    if (ef && missing(drift.end)) 
        stop("You must specify both arguments 'drift' and 'drift.end' for a futility boundary " \%,\% 
            "using stochastic curtailment\n")
    nlook <- length(frac)
    if (ef) {
        if (length(drift) > 0 && length(drift) < nlook) 
            stop("Argument 'drift' must be either unspecified (efficacy) or " \%,\% 
                "of the same length as 'frac': length(drift)=" \%,\% 
                length(drift) \%,\% ", " \%,\% "length(frac)=" \%,\% 
                length(frac) \%,\% ".\n")
    }
    if (!ef) {
        drift <- rep(0, nlook)
        drift.end <- 0
    }
    if (length(prob) != 1 && length(prob) < nlook) 
        stop("Argument 'prob' must be either of length 1 or of the same length as " \%,\% 
            "'frac': length(prob)=" \%,\% length(prob) \%,\% ", length(frac)=" \%,\% 
            length(frac) \%,\% ".\n")
    if (length(prob) == 1) 
        prob <- rep(prob, nlook)
    b <- NULL
    frac.old <- 0
    for (l in 1:nlook) {
        frac.new <- frac[l]
        b.l <- .C("StCu2Bnds", pmu = as.double(c(drift[l], drift.end)), 
            pfrac = as.double(frac.new), pzcrit = as.double(be.end), 
            prho = as.double(prob[l]), pef = as.integer(ef), 
            b = double(1), PACKAGE = "PwrGSD")$b
        b <- c(b, b.l)
        frac.old <- frac.new
    }
    b
  }
}
