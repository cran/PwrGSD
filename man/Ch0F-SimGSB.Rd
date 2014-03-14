\name{SimGSB}
\alias{SimGSB}
\title{Verifies the results of "GrpSeqBnds" via simulation}
\description{
  Verifies the results of \code{GrpSeqBnds} via simulation
}
\usage{
SimGSB(object, nsim = 1e+05, ...)
}
\arguments{
  \item{object}{an object of class either \code{boundaries} or \code{PwrGSD}}
  \item{nsim}{number of simulations to do}
  \item{\dots}{if \code{object} is of class \code{PwrGSD} and there are
    more than one statistic under investigation, then you may specify an
    argument \code{stat}. The default value is 1, meaning the first one.}
}
\value{
  A tabulation of the results
}
\author{Grant Izmirlian <izmirlian@nih.gov}

\seealso{\code{\link{GrpSeqBnds}}}
\examples{
  ## none as yet
}
\keyword{misc}
