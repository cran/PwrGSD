\name{Haybittle}
\alias{Haybittle}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Haybittle method of Boundary Construction}
\description{
  The function \code{Haybittle} is used in calls to the functions
  \code{GrpSeqBnds} and \code{PwrGSD} as a possible setting for the
  argument \code{EfficacyBoundary}. NOTE: the Haybittle method is not
  used for construction a futility boundary (does this even make
  sense?). The Haybittle method is one of four currently availiable
  choices, the others being \code{LanDemets}, \code{SC} (stochastic
  curtailment), and user specified.
}
\usage{
Haybittle(alpha, b.Haybittle, from = NULL, to = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{The total probability of type I error.}
  \item{b.Haybittle}{User specified efficacy boundary at all but the
    last analysis.}
  \item{from}{WARNING EXPERIMENTAL: See the documentation under
    \code{LanDemets} or \code{SC}. I'm not quite sure if this works
    or even makes sense. Don't use it, ok?}
  \item{to}{See above.}
}
\details{
  The Haybittle method is the simplest of all methods for efficacy
  boundary construction. You just specify a value such as 3 which is
  so large that multiple testing with standard normal criterion 3
  (upper quantile is 0.00135) means that after 20 analyses you've
  spent only 0.027.  The value of the criterion for the final analysis
  is determined so that the total type I error is obtained.
}
\value{
  An object of class \code{boundary.construction.method} which is really a list
  with the following components. The print method displays the original
  call.
  \item{type}{Gives the boundary construction method type, which is the character
    string "Haybittle"}
  \item{alpha}{The numeric value passed to the argument 'alpha' which is
    the total probability of type I error.}
  \item{b.Haybittle}{The numeric value passed to the argument 'b.Haybittle' which is
    the user specified efficacy boundary at all but the
    last analysis.}
  \item{from}{Description of 'comp2'}
  \item{to}{You're not using this, right? }
  \item{call}{see above.}
}
\references{see references under \code{\link{PwrGSD}}}
\author{Grant Izmirlian}
\note{The print method returns the call by default} 
\seealso{\code{\link{LanDemets}}, \code{\link{SC}}, \code{\link{GrpSeqBnds}}, and \code{\link{PwrGSD}}}
\examples{
## example 1: what is the result of calling a Boundary Construction Method function
    ## A call to 'Haybittle' just returns the call
    Haybittle(alpha=0.05, b.Haybittle=3)
    
    ## It does arguement checking...this results in an error
    \dontrun{
      Haybittle(alpha=0.05)
    }
    
    ## but really its value is a list with the a component containing
    ## the boundary method type, "LanDemts", and components for each
    ## of the arguments.
    names(Haybittle(alpha=0.05, b.Haybittle=3))

    Haybittle(alpha=0.05, b.Haybittle=3)$type
    Haybittle(alpha=0.05, b.Haybittle=3)$alpha
    Haybittle(alpha=0.05, b.Haybittle=3)$b.Haybittle
    Haybittle(alpha=0.05, b.Haybittle=3)$call

## example 2: ...But the intended purpose of the spending functions its
    ## intended purpose is in constructing calls to 'GrpSeqBnds' and to 'PwrGSD':
     

    frac <- c(0.07614902,0.1135391,0.168252,0.2336901,0.3186155,
              0.4164776,0.5352199,0.670739,0.8246061,1)

    test <- GrpSeqBnds(frac=frac, EfficacyBoundary=Haybittle(alpha=0.025, b.Haybittle=3))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{design}
