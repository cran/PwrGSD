\name{CRRtoRR}
\alias{CRRtoRR}
\title{Cumulative-risk ratios to risk ratios}
\description{
  Given a vector of cumulative-risk ratios, computes risk ratios
}
\usage{
 CRRtoRR(CRR, DT, h = NULL)
}
\arguments{
  \item{CRR}{vector of cumulative risk ratios of length \code{m}}
  \item{DT}{vector of time increments upon which the cumulative ratios
    represent. For example if the hazard takes values
    $h_1, h_2, \ldots, h_m$ on the intervals $[t_1, t_2), [t_2, t_3),
    \ldots, [t_m, t_{m+1})$ then \code{DT} will be c($t_2-t_1,t_3 -t_2,
    \ldots, t_{m+1} - t_m )$
  }
  \item{h}{The hazard in the reference arm, of length \code{m}}
}
\value{
  The vector of risk ratios at the \code{m} time points
}
\author{Grant Izmirlian <izmirlian@nih.gov>}
\examples{
  ## none as yet
}
\keyword{survival}
