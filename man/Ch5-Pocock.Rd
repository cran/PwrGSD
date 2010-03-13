\name{Pocock}
\alias{Pocock}

\title{The Pocock Alpha Spending Function}
\description{Stipulates alpha spending according to the Pocock 
  spending function in the Lan-Demets boundary construction method. Its
  intended purpose is in constructing calls to \code{GrpSeqBnds} and
  \code{PwrGSD}. 
}
\usage{
Pocock()
}
%- maybe also 'usage' for other objects documented here.
\value{
  An object of class \code{spending.function}
  \item{type}{Gives the spending function type, which is the character
    string "Pocock"}
  \item{call}{returns the call}
}
\references{see references under \code{\link{PwrGSD}}}
\author{Grant Izmirlian}
\note{The print method returns the call by default} 
\seealso{\code{\link{LanDemets}}, \code{\link{ObrienFleming}}, \code{\link{Pow}},
  \code{\link{GrpSeqBnds}}, \code{\link{PwrGSD}}}
\examples{
## example 1: what is the result of calling a spending function

    ## A call to 'Pocock' just returns the call
    Pocock()

    ## but really its value is a list with a component named
    ## 'type' equal to "Pocock" and a component named
    ## 'call' equal to the call.
    names(Pocock)

    Pocock()$type

    Pocock()$call

## example 2: ...But the intended purpose of the spending functions its
    ## intended purpose is in constructing calls to 'GrpSeqBnds' and to 'PwrGSD':
     

    frac <- c(0.07614902,0.1135391,0.168252,0.2336901,0.3186155,
              0.4164776,0.5352199,0.670739,0.8246061,1)
    drift <- c(0.3836636,0.5117394,0.6918584,0.8657705,1.091984,
               1.311094,1.538582,1.818346,2.081775,2.345386)

    test <- GrpSeqBnds(frac=frac, EfficacyBoundary=LanDemets(alpha=0.05, spending=Pocock),
                       FutilityBoundary=LanDemets(alpha=0.10, spending=ObrienFleming),
                       drift=drift)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{design}
