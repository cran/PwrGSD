\name{cpd.PwrGSD}
\alias{cpd.PwrGSD}

\title{Create a skeleton compound PwrGSD object}
  
\description{
 Given a user defined indexing dataframe as its only argument,
  creates a skeleton compound PwrGSD object having a component
  \code{Elements}, a list of \code{PwrGSD} objects, of length
  equal to the number of rows in the indexing dataframe
}
\usage{
cpd.PwrGSD(descr)
}

\arguments{
  \item{descr}{A dataframe of a number of rows equal to the length
    of the resulting list, \code{Elements}, of \code{PwrGSD} objects.
    The user defines the mapping between rows of \code{descr}
    and components of \code{Elements} and uses it to set up a loop
    over scenarios.  There are several S3 classes and methods
    for example \code{plot.cpd.PwrGSD}, which exploit this mapping
    between characteristics of a run and the rows of \code{desr} for
    subsetting and constructing conditioned plots. See the example below.}
}
\value{
  An object of class \code{cpd.PwrGSD} containing elements:
  \item{date}{the POSIX date that the object was created--its quite useful}
  \item{Elements}{a list of length equal to the number of rows of \code{descr}
    which will later contain objects of class \code{PwrGSD}}
  \item{descr}{a copy of the indexing dataframe argument for use in
    navigating the compound object in subsequent calls to other
    functions such as the related \code{plot} method, and the subset
    extractor, \code{Elements}}
}

\note{
  A \code{cpd.PwrGSD} object essentially a list of \code{PwrGSD} objects
  that a user may set up in order to investigate the space of possible trial
  scenarios, test statistics, and boundary construction options. One
  could store a list of results without appealing at all to these
  internal indexing capabilities. The advantage of setting up a
  \code{cpd.PwrGSD} object is the nice summarization functionality
  provided, for example the plot method for the \code{cpd.PwrGSD} class.

  The key ingredient to (i) the construction of the empty
  object, (ii) and summarizing the results in tabular or plotted form
  via its manipulation in subsequent function calls, is
  the indexing dataset, \code{descr} (for description).  The correspondence
  between rows of \code{descr} and elements in the list of \code{PwrGSD}
  objects is purposely left very loose. In the example outlined below,
  the user creates a \dQuote{base case} call to \code{PwrGSD} and then decides
  which quantities in this \dQuote{base case} call to vary in order to
  navigate the space of possible trial scenarios, monitoring statistics
  and boundary construction methods.  Next, for each one of these
  settings being varied, a variable with levels that determine each
  possible setting is created.  The dataset \code{descr} is created
  with one line corresponding to each combination of the selection
  variables so created.  In order to ensure that there is 1-1
  correspondence between the order of the rows in \code{descr} and the
  order in the list \code{Elements} of \code{PwrGSD} objects, the user
  carries out the computation in a loop over rows of \code{descr} in
  which the values of the selection variables in each given row of
  \code{descr} are used to create the corresponding component of
  \code{Elements} via an update the \dQuote{base case} call.
}

\author{Grant Izmirlian <izmirlian@nih.gov>}

\seealso{\code{\link{Elements}}, \code{\link{plot.cpd.PwrGSD}} and \code{\link{Power}}}

\examples{
## don't worry--these examples are guaranteed to work,
## its just inconvenient for the package checker
\dontrun{
  library(PwrGSD)

## In order to set up a compound object of class `cpd.PwrGSD'
## we first construct a base case: a two arm trial randomized in just
## under eight years with a maximum of 20 years of follow-up.
## We compute power at a specific alternative, `rhaz', under
## an interim analysis plan with roughly one annual analysis, some
## crossover between intervention and control arms, with Efficacy
## and futility boundaries constructed via the Lan-Demets procedure
## with O'Brien-Fleming spending on the hybrid scale. Investigate
## the behavior of three weighted log-rank statistics.

test.example <-
  PwrGSD(EfficacyBoundary = LanDemets(alpha = 0.05, spending = ObrienFleming),
         FutilityBoundary = LanDemets(alpha = 0.1, spending = ObrienFleming),
	 RR.Futility = 0.82, sided="1<",method="A",accru =7.73, accrat =9818.65,
         tlook =c(7.14, 8.14, 9.14, 10.14, 10.64, 11.15, 12.14, 13.14,
                  14.14, 15.14, 16.14, 17.14, 18.14, 19.14, 20.14),
         tcut0 =0:19, h0 =c(rep(3.73e-04, 2), rep(7.45e-04, 3),
                            rep(1.49e-03, 15)),
         tcut1 =0:19, rhaz =c(1, 0.9125, 0.8688, 0.7814, 0.6941,
                              0.6943, 0.6072, 0.5202, 0.4332, 0.6520,
                              0.6524, 0.6527, 0.6530, 0.6534, 0.6537,
                              0.6541, 0.6544, 0.6547, 0.6551, 0.6554),
         tcutc0 =0:19, hc0 =c(rep(1.05e-02, 2), rep(2.09e-02, 3),
                              rep(4.19e-02, 15)),
         tcutc1 =0:19, hc1 =c(rep(1.05e-02, 2), rep(2.09e-02, 3),
                              rep(4.19e-02, 15)),
         tcutd0B =c(0, 13), hd0B =c(0.04777, 0),
         tcutd1B =0:6, hd1B =c(0.1109, 0.1381, 0.1485, 0.1637, 0.2446,
                               0.2497, 0),
         noncompliance =crossover, gradual =TRUE,
         WtFun =c("FH", "SFH", "Ramp"),
         ppar =c(0, 1, 0, 1, 10, 10))

## we will construct a variety of alternate hypotheses relative to the
## base case specified above

  rhaz <- 
    c(1, 0.9125, 0.8688, 0.7814, 0.6941, 0.6943, 0.6072, 0.5202, 0.4332,
    0.652, 0.6524, 0.6527, 0.653, 0.6534, 0.6537, 0.6541, 0.6544,
    0.6547, 0.6551, 0.6554)

  max.effect <- 0.80 + 0.05*(0:8)
  n.me <- length(max.effect)

## we will also vary extent of censoring relative to the base case
## specified above

  hc <- c(rep(0.0105, 2), rep(0.0209, 3), rep(0.0419, 15))

  cens.amt <- 0.75 + 0.25*(0:2)
  n.ca <- length(cens.amt)

## we may also wish to compare the Lan-Demets/O'Brien-Fleming efficacy
## boundary with a Lan-Demets/linear spending boundary

  Eff.bound.choice <- 1:2
  ebc.nms <- c("LanDemets(alpha=0.05, spending=ObrienFleming)",
               "LanDemets(alpha=0.05, spending=Pow(1))")
  n.ec <- length(Eff.bound.choice)

## The following line creates the indexing dataframe, `descr', with one
## line for each possible combination of the selection variables we've
## created. 


  descr <- as.data.frame(
              cbind(Eff.bound.choice=rep(Eff.bound.choice, each=n.ca*n.me),
                    cens.amt=rep(rep(cens.amt, each=n.me), n.ec),
                    max.effect=rep(max.effect, n.ec*n.ca)))

  descr$Eff.bound.choice <- ebc.nms[descr$Eff.bound.choice]

## Now descr contains one row for each combination of the levels of
## the user defined selection variables, `Eff.bound.choice',
## `max.effect' and `cens.amt'. Keep in mind that the names and number
## of these variables is arbitrary. Next we create a skeleton
## `cpd.PwrGSD' object with a call to the function `cpd.PwrGSD' with
## argument `descr' 

  test.example.set <- cpd.PwrGSD(descr)

## Now, the newly created object, of class `cpd.PwrGSD', contains
## an element `descr', a component `date', the date created 
## and a component `Elements', an empty list of length equal
## to the number of rows in `descr'.  Next we do the computation in
## a loop over the rows of `descr'.

  n.descr <- nrow(descr)

  for(k in 1:n.descr){

    ## First, we copy the original call to the current call,
    ## `Elements[[k]]$call'

    test.example.set$Elements[[k]]$call <- test.example$call

    ## Use the efficacy boundary choice in the kth row of `descr'
    ## to set the efficacy boundary choice in the current call

    test.example.set$Elements[[k]]$call$EfficacyBoundary <- 
    parse(text=as.character(descr[k,"Eff.bound.choice"]))[[1]]

    ## Derive the `rhaz' defined by the selection variable "max.effect"
    ## in the kth row of `descr' and use this to set the `rhaz'
    ## components of the current call

    test.example.set$Elements[[k]]$call$rhaz <-
                            exp(descr[k,"max.effect"] * log(rhaz))

    ## Derive the censoring components from the selection variable
    ## "cens.amt" in the kth row of `descr' and place that result
    ## into the current call
    
    test.example.set$Elements[[k]]$call$hc0 <-
    test.example.set$Elements[[k]]$call$hc1 <- descr[k, "cens.amt"] * hc

    ## Now the current call corresponds exactly to the selection
    ## variable values in row `k' of `descr'. The computation is
    ## done by calling `update'

    test.example.set$Elements[[k]] <- update(test.example.set$Elements[[k]])
    cat(k/n.descr, "\r")
  }

  ## We can create a new `cpd.PwrGSD' object by subsetting on
  ## the selection variables in `descr':

  Elements(test.example.set, 
           subset=(substring(Eff.bound.choice, 32, 34)=="Obr" &
                            max.effect >= 1))


  ## or we can plot the results -- see the help under `plot.cpd.PwrGSD'

  plot(test.example.set, formula = ~ max.effect | stat * cens.amt,
       subset=(substring(Eff.bound.choice, 32, 34)=="Obr"))

  plot(test.example.set, formula = ~ max.effect | stat * cens.amt,
       subset=(substring(Eff.bound.choice, 32, 34)=="Pow"))

  ## Notice the appearance of the selection variable `stat' which was
  ## not defined in the dataset `descr'. 

  ## Recall that each single "PwrGSD" object can contain results
  ## for a list of test statistics, as in the example shown here where
  ## we have results on three statistics per component of `Elements'.
  ## For this reason the variable `stat' can be also be referenced in
  ## the `subset' or `formula' arguments of calls to this `plot' method,
  ## and in the `subset' argument of the function `Power' shown below.

  ## The function `Power' is used to convert the `cpd.PwrGSD' object
  ## into  a dataframe, stacked by rows of `descr' and by `stat'
  ## (there are three statistics being profiled per each component of
  ## `Elements'), for generating tables or performing other 
  ## computations.

  Power(test.example.set,
        subset=(substring(Eff.bound.choice, 32, 34)=="Pow" & stat \%in\% c(1,3)))

}
}
\keyword{data}
