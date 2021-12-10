\name{detail}
\alias{detail}
\title{
  Function to extract the 'detail' component from a \code{PwrGSD} object
}
\description{
  Extracts the 'detail' component from an object of class \code{PwrGSD}
}
\usage{
detail(obj)
}
\arguments{
  \item{obj}{
     An object of class \code{PwrGSD} returned from the fucntion
     \code{PwrGSD} or a component of the list returned by the funtion 
     \code{cpd.PwrGSD}
  }
}
\value{
  The 'detail' component of the object. For the Asymptotic method,
  this will be most of the quantities involved in the computation,
  the input parameters such as the various incidence rates, cross over
  rates etc, as well as intermediate computations such as the drift
  function variance function as well. For the simulation method, some of
  these are returned an in addition, the simulated event histories.
}
\author{
  Grant Izmirlian <izmirlian at nih dot gov>
}
