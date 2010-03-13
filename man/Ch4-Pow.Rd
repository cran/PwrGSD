\name{Pow}
\alias{Pow}

\title{The Wang-Tsiatis Power Alpha Spending Function}
\description{Stipulates alpha spending according to the Wang-Tsiatis
  Power function in the Lan-Demets boundary construction method. Its
  intended purpose is in constructing calls to \code{GrpSeqBnds} and
  \code{PwrGSD}. 
} 
\usage{
Pow(rho)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rho}{The exponent for the Wang-Tsiatis power spending function}
}
\details{Larger \code{rho} results in more conservative
    boundaries.  \code{rho}=3 is roughly equivalent to Obrien-Fleming
    spending. \code{rho}=1 spends \code{alpha} linearly in the
    information fraction
}
\value{
  An object of class \code{spending.function} which is really a list
  with the following components. The print method displays the original
  call.
  \item{type}{Gives the spending function type, which is the character
    string "Pow"}
  \item{rho}{the numeric value passed to the single argument, \code{rho}}
  \item{call}{returns the call}
}
\references{see references under \code{\link{PwrGSD}}}
\author{Grant Izmirlian}
\note{The print method returns the call by default} 
\seealso{\code{\link{LanDemets}}, \code{\link{ObrienFleming}}, \code{\link{Pocock}},
  \code{\link{GrpSeqBnds}}, \code{\link{PwrGSD}}}
\examples{
## example 1: what is the result of calling a spending function
    ## A call to 'Pow' just returns the call
    Pow(rho=2)

    ## It does argument checking...the following results in an error:
    \dontrun{
      Pow()
    }
    
    ## it doesn't matter whether the argument is named or not,
    ## either produces the same result
    Pow(2)

    ## but really its value is a list with a component named
    ## 'type' equal to "Pow", a component named 'rho' equal
    ## to the numeric value passed to the single argument 'rho'
    ## and a component  named 'call' equal to the call.
    names(Pow(rho=2))
    
    names(Pow(2))
    
    Pow(rho=2)$type
    Pow(rho=2)$rho
    Pow(rho=2)$call    

## example 2: ...But the intended purpose of the spending functions its
    ## intended purpose is in constructing calls to 'GrpSeqBnds' and to 'PwrGSD':
     

    frac <- c(0.07614902,0.1135391,0.168252,0.2336901,0.3186155,
              0.4164776,0.5352199,0.670739,0.8246061,1)
    drift <- c(0.3836636,0.5117394,0.6918584,0.8657705,1.091984,
               1.311094,1.538582,1.818346,2.081775,2.345386)

    test <- GrpSeqBnds(frac=frac, EfficacyBoundary=LanDemets(alpha=0.05, spending=Pow(2)),
                       FutilityBoundary=LanDemets(alpha=0.10, spending=ObrienFleming),
                       drift=drift)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{design}
