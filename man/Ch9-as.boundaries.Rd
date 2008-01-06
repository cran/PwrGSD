\name{as.boundaries}
\alias{as.boundaries}
\title{Convert a "PwrGSD" object to a "boundaries" object}
\description{
  Convert a \code{PwrGSD} object to a \code{boundaries} object
}
\usage{
  as.boundaries(object, ...)
}
\arguments{
  \item{object}{an object of class \code{PwrGSD}}
  \item{\dots}{if \code{object} is of class \code{PwrGSD} and there are
    more than one statistic under investigation, then you may specify an
    argument \code{stat}. The default value is 1, meaning the first one.}
}
\value{
  an object of class \code{boundaries}.  See the documentation for
  \code{\link{GrpSeqBnds}}
}
\author{Grant Izmirlian <izmirlian@nih.gov}
\seealso{\code{\link{GrpSeqBnds}}}
\examples{
  ## none as yet
}
\keyword{design}
\keyword{htest}
