\name{LanDemets}
\alias{LanDemets}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Lan-Demets method of Boundary Construction}
\description{
  The function \code{LanDemets} is used in calls to the functions
  \code{GrpSeqBnds} and \code{PwrGSD} as a possible setting for the
  arguments \code{EfficacyBoundary} and \code{FutilityBoundary}, in
  specification of the method whereby efficacy and or futility
  boundaries are to be constructed. The Lan-Demets method is one of
  four currently availiable choices, the others being \code{SC}
  (stochastic curtailment), \code{Haybittle} and user specified.
}
\usage{
LanDemets(alpha, spending, from = NULL, to = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{If \code{LanDemets} is used to specify the \code{EfficacyBoundary} 
    then the argument \code{alpha} is the total probability of type I error.  
    If \code{LanDemets} is used to specify the \code{FutilityBoundary} then the 
    argument \code{alpha} is the total probability of type II error.}
  \item{spending}{Specify the alpha spending function. Set this to
    \code{ObrienFleming}, \code{Pow(rho=<x>)}, or \code{Pocock}. See
    help files for these spending functions.}
  \item{from}{WARNING EXPERIMENTAL: you can actually construct boundaries via a
    hybrid of the 3 boundary construction methods, \code{LanDemets}, \code{SC},
    and 'user specified'. When using a hybrid boundry, set the argument
    \code{EfficacyBoundary} or \code{FutilityBoundary} respectively,
    to a list with components \code{LanDemets}, \code{SC}, or user
    specified numbers. In the former two cases, \code{from} and
    \code{to} are used in \code{LanDemets} and also in \code{SC}
    to stipulate how many interim analyses they are in effect.
    See the help for \code{GrpSeqBnds} and \code{PwrGSD}}
  \item{to}{See above.}
}
\details{
    The cornerstone of the Lan-Demets method is that the
    amount of alpha (type I or II error probability) that is "spent" at
    a given interim analysis is determined via a user specified
    "spending function". A spending function is a monotone increasing
    mapping on (0,1) with range (0,alpha).  The 'alpha' spent at a given
    analysis is determined by the increment in the values of the spending
    function at the current and at the most recent information fractions.
}
\value{
  An object of class \code{boundary.construction.method} which is really a list
  with the following components. The print method displays the original
  call.
  \item{type}{Gives the boundary construction method type, which is the character
    string "LanDemets"}
  \item{alpha}{The numeric value passed to the argument 'alpha' which is
    the total probability of type I (efficacy) or type II (futility) error.}
  \item{spending}{The spending function that was passed to the argument
    'spending'. Note that this will be of class 'name'  for
    'ObrienFleming' and 'Pocock', but will be of class 'function' for
    'Pow'}
  \item{from}{The numeric value passed to the argument 'from'. See above.}
  \item{to}{The numeric value passed to the argument 'to'. See above.}
  \item{call}{returns the call}
}
\references{see references under \code{\link{PwrGSD}}}
\author{Grant Izmirlian}
\note{The print method returns the call by default} 
\seealso{\code{\link{SC}}, \code{\link{ObrienFleming}},
  \code{\link{Pow}}, \code{\link{Pocock}}, \code{\link{GrpSeqBnds}}, and
  \code{\link{PwrGSD}}} 
\examples{
## example 1: what is the result of calling a Boundary Construction Method function
    ## A call to 'LanDemets' just returns the call
    LanDemets(alpha=0.05, spending=ObrienFleming)
    
    ## It does arguement checking...this results in an error
    \dontrun{
      LanDemets(alpha=0.05)
    }
    
    ## but really its value is a list with the a component containing
    ## the boundary method type, "LanDemts", and components for each
    ## of the arguments.
    names(LanDemets(alpha=0.05, spending=ObrienFleming))

    LanDemets(alpha=0.05, spending=ObrienFleming)$type
    LanDemets(alpha=0.05, spending=ObrienFleming)$alpha
    LanDemets(alpha=0.05, spending=ObrienFleming)$spending
    class(LanDemets(alpha=0.05, spending=ObrienFleming)$spending)
    LanDemets(alpha=0.05, spending=Pow(2))$spending
    class(LanDemets(alpha=0.05, spending=Pow(2))$spending)
    LanDemets(alpha=0.05, spending=ObrienFleming)$call

## example 2: ...But the intended purpose of the spending functions its
    ## intended purpose is in constructing calls to 'GrpSeqBnds' and to 'PwrGSD':
     

    frac <- c(0.07614902,0.1135391,0.168252,0.2336901,0.3186155,
              0.4164776,0.5352199,0.670739,0.8246061,1)
    drift <- c(0.3836636,0.5117394,0.6918584,0.8657705,1.091984,
               1.311094,1.538582,1.818346,2.081775,2.345386)

    test <- GrpSeqBnds(frac=frac, EfficacyBoundary=LanDemets(alpha=0.05, spending=ObrienFleming),
                       FutilityBoundary=LanDemets(alpha=0.10, spending=Pocock),
                       drift=drift)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{design}
