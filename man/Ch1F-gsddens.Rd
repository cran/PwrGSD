\name{gsd.dens}
\alias{gsd.dens}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A function for computing the probability density for the group
sequentially monitored test statistic.
}
\description{
  A function for computing the probability density for a sequentially monitored test. This
  is the joint density, in the rejection region, of (X_K, K), where X_K is the observed
  value of the test statistic upon efficacy boundary crossing, and K is the analysis
  number at which the efficacy boundary was crossed.
}
\usage{
gsd.dens(x, frac = NULL, scale="Standard")
}
\arguments{
  \item{x}{
    The main argument, \code{x}, is either a object of class \dQuote{boundaries} or a
    numeric vector. If it is of class \dQuote{boundaries} then no other arguments are
    required. If it is a numeric vector then the \code{frac} argument must be
    specified. See below. In this case, \code{x} will be the observed values of the
    statistic at the current and all prior analyses, either on the standard normal scale
    (the default) or on the \dQuote{Brownian} scale. For \dQuote{Brownian} scale, set
    argument \code{scale} to \dQuote{Brownian}.
  }
  \item{frac}{
    Required only when the main argument, \code{x}, is a numeric vector, and must be a
    vector of the same length. In this case, \code{frac} will be the information at the
    current and all prior interim analyses.
  }
  \item{scale}{
    Required only when the main argument, \code{x}, is a numeric vector. A switch
    indicating whether the elements of the numeric vector, \code{x}, are specified on the
    standard normal scale, \code{x}=\dQuote{Standard}, or on the Brownian scale,
    \code{x}=\dQuote{Brownian}.
  }
}
\value{
  A list with elements \code{x}, \code{dF}, \code{x1c}, and
  \code{dF1c}:
  \item{x}{
    Node points used in Gaussian quadrature. See examples below.
  }
  \item{dF}{
    Probability mass at each node point. See examples below.
  }
  \item{x1c}{
    Node points in the continuation region at the first analysis.
  }
  \item{dF1c}{
    Probability mass at each node point in the continuation region
    at the first analysis.
  }
}

\references{
Emerson, S. S. (1993). Computation of the uniform minimum variance unibiased estimator of
a normal mean following a group sequential trialdiscrete sequential boundaries for
clinical trials. Computers and Biomedical Research 26 68--73.

Izmirlian, G. (2014). Estimation of the relative risk following group sequential procedure
based upon the weighted log-rank statistic. Statistics and its Interface 00 00--00
}
\author{
Grant Izmirlian <izmirlig@mail.nih.gov>
}
\note{
  Also used in computation of Rao-Blackwell-ized bias adjusted point
  estimate for statistic observed to cross the efficacy boundary.
}
\seealso{
  \code{\link{EX1gXK}}
}
\examples{
  # Information fraction
  frac <- c(0.15, 0.37, 0.64, 0.76)

  # Efficacy Boundary 
  gsb <- GrpSeqBnds(frac=frac, EfficacyBoundary=LanDemets(spending=ObrienFleming, alpha=0.05))

  # To compute the p-value under the stagewise ordering, for an observed
  #  value of the monitoring statistic 2.1, crossing the efficacy
  #  boundary at the 4th analysis, we do the following

  be <- gsb$table[,"b.e"] 
  be[4] <- 2.1

  sum(gsd.dens(be, frac, scale="Standard")$dF)

}
\keyword{Stage-wise P-value}
\keyword{sequential statistic pdf}
