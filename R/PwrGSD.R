# Nodes and weights for Gauss-Legendre quadrature on [-1,1] for n=24

glegx24 <- c(-0.995187219997, -0.974728555971, -0.938274552003, -0.886415527004,
             -0.820001985974, -0.740124191579, -0.648093651937, -0.545421471389,
             -0.433793507626, -0.315042679696, -0.191118867474, -0.0640568928626,
              0.0640568928626, 0.191118867474,  0.315042679696,  0.433793507626,
              0.545421471389,  0.648093651937,  0.740124191579,  0.820001985974,
              0.886415527004,  0.938274552003,  0.974728555971,  0.995187219997)

glegw24 <- c( 0.0123412297185, 0.0285313860439, 0.0442774398676, 0.0592985850337,
              0.0733464816501, 0.0861901617117, 0.0976186522496, 0.107444270284,
              0.115505668234,  0.121670473118,  0.125837456543,  0.127938195546,
              0.127938195546,  0.125837456543,  0.121670473118,  0.115505668234,
              0.107444270284,  0.0976186522496, 0.0861901617117, 0.0733464816501,
              0.0592985850337, 0.0442774398676, 0.0285313860439, 0.0123412297185)

"%,%" <-
function(x,y) paste(x,y,sep="")

DX <- function(x)c(x[1],diff(x))

"PwrGSD" <- 
function(EfficacyBoundary = LanDemets(alpha=0.05,spending=ObrienFleming),
         FutilityBoundary = LanDemets(alpha=0.10,spending=ObrienFleming),
         sided =c("2",">","<"),method=c("S","A"),accru,accrat,tlook,
         tcut0 = NULL,h0 = NULL,s0 = NULL,tcut1 = NULL,rhaz = NULL,
         h1 = NULL,s1 = NULL,tcutc0 = NULL,hc0 = NULL,sc0 = NULL,tcutc1 = NULL,hc1 = NULL,
         sc1 = NULL,tcutd0A = NULL,hd0A = NULL,sd0A = NULL,tcutd0B = NULL,hd0B = NULL,sd0B = NULL,
         tcutd1A = NULL,hd1A = NULL,sd1A = NULL,tcutd1B = NULL,hd1B = NULL,sd1B = NULL,
         tcutx0A = NULL,hx0A = NULL,sx0A = NULL,tcutx0B = NULL,hx0B = NULL,sx0B = NULL,
         tcutx1A = NULL,hx1A = NULL,sx1A = NULL,tcutx1B = NULL,hx1B = NULL,sx1B = NULL,
         noncompliance = c("none","crossover","mixed","user"),gradual = FALSE,
         WtFun = c("FH","SFH","Ramp"),ppar = cbind(c(0,0)),
         Spend.Info=c("Variance","Events","Hybrid(k)","Calendar"),
         RR.Futility = NULL,qProp.one.or.Q = c("one","Q"),Nsim=NULL,
         detail = FALSE,StatType = c("WLR","ISD"))
{
  .call. <- match.call()
  prfx <- c("Sim","Asy")
  if(missing(method)) method <- "A"
  .call.[[1]] <- as.name(prfx[1+(method=="A")] %,% as.character(.call.[[1]]))
  .call.$method <- NULL
  if(method=="A") .call.$Nsim <- .call.$detail <- .call.$StatType <- NULL
  if(method=="S"){
    if(missing(Nsim)) stop("Argument 'Nsim' is required for method==\"S\"")
    if(missing(StatType)) StatType <- "WLR"
  }
  ans <- eval(.call.)
  ans$call[[1]] <- as.name("PwrGSD")
  ans$call$method <- method
  ans
}

"AsyPwrGSD" <- 
function(EfficacyBoundary = LanDemets(alpha=0.05,spending=ObrienFleming),
         FutilityBoundary = LanDemets(alpha=0.10,spending=ObrienFleming),
         sided=c("2",">","<"),accru,accrat,tlook,tcut0=NULL,h0=NULL,s0=NULL,tcut1=NULL,
         rhaz=NULL, h1=NULL, s1=NULL, tcutc0=NULL, hc0 = NULL, sc0=NULL, 
         tcutc1=NULL, hc1 = NULL, sc1 = NULL, tcutd0A=NULL, hd0A=NULL, 
         sd0A=NULL, tcutd0B=NULL, hd0B=NULL, sd0B=NULL, tcutd1A=NULL, hd1A=NULL, 
         sd1A=NULL, tcutd1B=NULL, hd1B=NULL, sd1B=NULL, tcutx0A=NULL, hx0A=NULL, 
         sx0A=NULL, tcutx0B=NULL, hx0B=NULL, sx0B=NULL, tcutx1A=NULL, hx1A=NULL, 
         sx1A=NULL, tcutx1B=NULL, hx1B=NULL, sx1B=NULL, 
         noncompliance=c("none", "crossover", "mixed", "user"), 
         gradual = FALSE, WtFun=c("FH","SFH","Ramp"), ppar=c(0, 0),
         Spend.Info=c("Variance", "Events", "Hybrid(k)", "Calendar"),
         RR.Futility=NULL, V.end = NULL,qProp.one.or.Q = c("one","Q"))
{

#
#          Alpha.Efficacy=0.05,Alpha.Futility=0.90,
#          Boundary.Efficacy = c("Lan-Demets","Haybittle","SC"),
#          Boundary.Futility = c("Lan-Demets", "Haybittle","SC"),
#          Spending.Efficacy = c("Obrien-Fleming", "Pocock", "Power"),
#          Spending.Futility = c("Obrien-Fleming", "Pocock", "Power"), 
#          rho.Efficacy=0, rho.Futility=0, b.Haybittle = 3,
#          rho.Eff.SC = 1.0, rho.Fut.SC = 1.0,
#
        .call. <- match.call()
        normut <- 8.20953615160139
        psimin <- 1.0e-15

        if(missing(WtFun)) {
          WtFun <- "FH"
          wttyp <- 0
          ppar <- c(0,0)
          nstat <- 1
        }
        else{
          if(missing(ppar) && any(WtFun!="FH"))
            stop("You must specify parameters for the chosen weight function(s) in the 'ppar' argument.")
          nstat <- length(WtFun)
          wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1
        }
        nlook <- length(tlook)
	tend <- tlook[nlook]
	if (missing(sided)) sided <- "2"
	if(!(sided%in%c("2",">","<"))) stop("Argument 'sided' must be " %,%
                                            "equal to '2', '>' or '<'")
        sided <- c(2,1,-1)[grep(sided, c("2",">","<"))]
        
        # determine type of Efficacy Boundary Specification
        mode.E <- "WRONG"
        BE.missing <- missing(EfficacyBoundary)
        if(BE.missing) mode.E <- "NULL"
        if(!BE.missing){
          try.BE <- try(EfficacyBoundary, silent=TRUE)
          if(is.null(attr(try.BE, "class"))) mode.E <- mode(EfficacyBoundary)
          else if(attr(try.BE, "class")=="try-error") mode.E <- "call"
        }
        if (!(mode.E %in% c("call", "numeric", "NULL"))) 
          stop("Argument 'EfficacyBoundary' must be of mode 'call', 'numeric' or 'NULL'")
        is.myE <- FALSE
        n.myE <- 0
        if (mode.E == "numeric") {
          n.myE <- length(EfficacyBoundary)
          is.myE <- TRUE
          my.Efficacy <- EfficacyBoundary
          .call.$EfficacyBoundary <- as.call(expression(LanDemets, alpha=0.05, spending=ObrienFleming))
        }
        do.efficacy <- !BE.missing || is.myE
        if (mode.E == "NULL" && do.efficacy) 
          .call.$EfficacyBoundary <- as.call(expression(LanDemets, alpha=0.05, spending=ObrienFleming))
        
        # determine type of Futility Boundary Specification
        mode.F <- "WRONG"
        BF.missing <- missing(FutilityBoundary)
        if(BF.missing) mode.F <- "NULL"
        if(!BF.missing){
          try.BF <- try(FutilityBoundary, silent=TRUE)
          if(is.null(attr(try.BF, "class"))) mode.F <- mode(FutilityBoundary)
          else if(attr(try.BF, "class")=="try-error") mode.F <- "call"
        }
        if (!(mode.F %in% c("call", "numeric", "NULL"))) 
          stop("Argument 'FutilityBoundary' must be of mode 'call', 'numeric' or 'NULL'")
        is.myF <- FALSE
        n.myF <- 0
        if (mode.F == "numeric") {
          n.myF <- length(FutilityBoundary)
          is.myF <- TRUE
          my.Futility <- FutilityBoundary
          .call.$FutilityBoundary <- as.call(expression(LanDemets, alpha=0.10, spending=ObrienFleming))
        }
        do.futility <- !BF.missing || is.myF
        if(!do.futility) Alpha.Futility <- psimin
        if (mode.F == "NULL" && do.futility)
          .call.$FutilityBoundary <- as.call(expression(LanDemets, alpha=0.10, spending=ObrienFleming))
        
        
        # Parse Boundary Methods:
        ..call. <- .call.
        ..call.[[1]] <- as.name("ParseBoundaryMethods")
        nms.call <- names(..call.)[-1]
        ind.del <- which(!(nms.call %in% c("EfficacyBoundary", "FutilityBoundary")))
        for(k in -sort(-ind.del))
          ..call.[[1+k]] <- NULL
        ..call.$nlooks <- nlook
        eval(..call.)
        
        # Result:  Now the following are defined within this scope:
        #
        #   Alpha.Efficacy, Alpha.Futility, nbnd.e, nbnd.f, nsf.e, nsf.f, rho.Efficacy, rho.Futility,
        #   b.Haybittle, drift.end, crit.e, crit.f
        #
        
        if (!do.efficacy && !do.futility) 
        stop("You must specify one or both of the arguments 'EfficacyBoundary' and 'FutilityBoundary', " %,% 
             "see the documentation")

        b.e <- rep(0, nlook)
        if (is.myE){
          if(n.myE != nlook )
            stop("User supplied efficacy boundary in 'EfficacyBoundary' must be of the same length as 'frac'")
          b.e <- my.Efficacy
        }
        
        b.f <- rep(0, nlook)
        if (is.myF){
          if(n.myF != nlook )
            stop("User supplied futility boundary in 'FutilityBoundary' must be of the same length as 'frac'")
          b.f <- my.Futility
        }
    
#	if(missing(Alpha.Efficacy)) stop("Missing argument 'Alpha.Efficacy', the total type I error")	
#	Alpha.Efficacy <- Alpha.Efficacy/2^(sided == 2)
#        
#       if (missing(Boundary.Efficacy))
#          Boundary.Efficacy <- NULL
#       mode.E <- mode(Boundary.Efficacy)
#       if(!(mode.E %in% c("character", "numeric", "NULL")))
#         stop("Argument 'Boundary.Efficacy' must be of character, numeric or NULL mode")
#       b.e <- rep(0,nlook)
#       is.myE <- FALSE
#       if(mode.E=="numeric"){
#         n.myE <- length(Boundary.Efficacy)
#         if(n.myE != nlook)
#           stop("User supplied efficacy boundary in 'Boundary.Efficacy' must be of the same length as 'frac'")
#         is.myE <- TRUE
#         b.e <- Boundary.Efficacy
#         Boundary.Efficacy <- "Lan-Demets"
#       }
#       if(mode.E=="NULL") Boundary.Efficacy <- "Lan-Demets"
#       nbnd.e <- grep(Boundary.Efficacy, c("Lan-Demets", "Haybittle", "SC"))
#       if(nbnd.e==3 && missing(rho.Eff.SC))
#         stop("Argument 'rho.Eff.SC' (total type I error for Stochastic Curtailment criterion) must be specified")

#       if (missing(Spending.Efficacy))
#         Spending.Efficacy <- "Obrien-Fleming"
#	if(Spending.Efficacy=="Power" && missing(rho.Efficacy))
#	  stop("Argument 'rho.Efficacy' must be supplied")
#       nsf.e <- grep(Spending.Efficacy, c("Obrien-Fleming", "Pocock", "Power"))

#       if(missing(Boundary.Futility))
#         Boundary.Futility <- NULL          
#       mode.F <- mode(Boundary.Futility)       
#       if(!(mode.F %in% c("character", "numeric", "NULL")))
#         stop("Argument 'Boundary.Futility' must be of character, numeric or NULL mode")
#       b.f <- rep(0,nlook)
#       is.myF <- FALSE
#       if(mode.F=="numeric"){
#         n.myF <- length(Boundary.Futility)
#         if(n.myF != nlook)
#           stop("User supplied futility boundary in 'Boundary.Futility' must be of the same length as 'frac'")
#         is.myF <- TRUE
#         b.f <- Boundary.Futility
#         Boundary.Futility <- "Lan-Demets"
#       }

#     	dofu <- 1
#       if(missing(Alpha.Futility) && mode.F != "NULL" && mode.F!="numeric")
#         stop("You specified a futility boundary construction method in the argument 'Boundary.Futility' that requires" %,%
#              " a probability of total type II error, 'Alpha.Futility'")
#       if (missing(Alpha.Futility) && mode.F == "NULL" && (is.myF==0)){
#         dofu <- 0
#         Alpha.Futility <- psimin
#       }
        
#       if(mode.F=="NULL") Boundary.Futility <- "Lan-Demets"
#       nbnd.f <- grep(Boundary.Futility, c("Lan-Demets", "Haybittle", "SC"))
#       if(nbnd.f==3 && missing(rho.Fut.SC))
#         stop("Argument 'rho.Fut.SC' (total type II error for Stochastic Curtailment criterion) must be specified")

#	if(nbnd.f == 2) stop("Haybittle futility boundary makes no sense so that " %,%
#                            "understandibly, it is not supported")
        
#	if(missing(Spending.Futility)) Spending.Futility <- "Obrien-Fleming"
#	if(Spending.Futility=="Power" && missing(rho.Futility))
#	  stop("Argument 'rho.Futility' must be supplied")
#       nsf.f <- grep(Spending.Futility, c("Obrien-Fleming", "Pocock", "Power"))
        
        no.SpndInfo <- missing(Spend.Info)
        if(no.SpndInfo) Spend.Info <- "Variance"
        if(!no.SpndInfo){
          Spend.Info <- as.character(.call.$Spend.Info)
          if(!(Spend.Info[1] %in% c("Variance", "Events", "Hybrid", "Calendar")))
            stop("Argument 'Spend.Info' is an expression of the form 'Variance', 'Events', 'Hybrid(k)', or 'Calendar'")
        }
        spend.info <- grep(Spend.Info[1], c("Variance", "Events", "Hybrid", "Calendar")) - 1
        spend.info.k <- 0
        if(spend.info==2 && length(Spend.Info)>1) spend.info.k <- as.integer(as.numeric(Spend.Info[2])) - 1

        user.V.end <- !missing(V.end)
        if(!user.V.end) V.end <- 0

        if(missing(qProp.one.or.Q))
          qProp.one.or.Q <- 0
        else{
          if(!(qProp.one.or.Q %in% c("one","Q")))
            stop("Argument 'qProp.one.or.Q' must be either \"one\" or \"Q\"")
          qProp.one.or.Q <- grep(qProp.one.or.Q, c("one", "Q")) - 1
        }
        
        nunq <- length(unique(sort(c(tcut0, tcut1, tcutc0, tcutc1, tcutd0A, tcutd0B, tcutd1A, tcutd1B, tcutx0A,
                                     tcutx0B, tcutx1A, tcutx1B))))    
        glegx <- glegx24
        glegw <- glegw24

        NGaussQ <- length(glegx)
	stoh <- function(tcut, s)
	{
		ncut <- length(tcut)
		Sold <- c(1, s[ - (ncut - 1)])
		dt <- diff(tcut)
		h <- log.0(s/Sold)/dt
		h
	}

	no.t0 <- missing(tcut0)
	no.s0 <- missing(s0)
	no.h0 <- missing(h0)
	if((no.s0 && no.h0) || no.t0)
		stop("Must specify 'tcut0' and ('h0' or 's0').")
	if(no.h0)
		h0 <- stoh(tcut0, s0)
	ncut0 <- length(tcut0)

	no.t1 <- missing(tcut1)
	no.s1 <- missing(s1)
	no.h1 <- missing(h1)
	no.rhaz <- missing(rhaz)
	if((no.s1 && no.rhaz && no.h1) || no.t1)
		stop("Must specify 'tcut1' and ('rhaz' or 'h1' or 's1').")
	if(!no.s1)
		h1 <- stoh(tcut1, s1)
	if(!no.rhaz){
		tcut.01 <- support(c(tcut0, tcut1))
		h0.01 <- lookup(tcut0, h0, tcut.01)$y
		rhaz.01 <- lookup(tcut1, rhaz, tcut.01)$y
		tcut1 <- tcut.01
		h1 <- rhaz.01 * h0.01
	}
	ncut1 <- length(tcut1)

        use.rhaz.fu <- 0
        if(missing(RR.Futility)){
          RR.Futility <- rep(1,ncut0)
          if(do.futility==1) use.rhaz.fu <- 1
        }
        if(length(RR.Futility)<ncut0) RR.Futility <- rep(RR.Futility[1], ncut0)
        
	no.tc0 <- missing(tcutc0)
	no.sc0 <- missing(sc0)
	no.hc0 <- missing(hc0)
	if((no.sc0 && no.hc0) || no.tc0)
		stop("Must specify 'tcutc0' and ('hc0' or 'sc0').")
	if(no.hc0)
		hc0 <- stoh(tcutc0, sc0)
	ncutc0 <- length(tcutc0)
	no.tc1 <- missing(tcutc1)
	no.sc1 <- missing(sc1)
	no.hc1 <- missing(hc1)
	if((no.sc1 && no.hc1) || no.tc1)
		stop("Must specify 'tcutc1' and ('hc1' or 'sc1').")
	if(no.hc1)
		hc1 <- stoh(tcutc1, sc1)
	ncutc1 <- length(tcutc1)
	noncompliance <- as.character(.call.$noncompliance)
	no.noncomp <- (length(noncompliance) == 0)
	if(no.noncomp)
		noncompliance <- "none"
	switch(noncompliance,
		none = {
			tcutd0A <- 0
			hd0A <- 1e-7
			ncutd0A <- 1

			tcutd0B <- 0
			hd0B <- 1e-7
			ncutd0B <- 1

			tcutd1A <- 0
			hd1A <- 1e-7
			ncutd1A <- 1

			tcutd1B <- 0
			hd1B <- 1e-7
			ncutd1B <- 1

			tcutx0A <- tcut0
			hx0A <- h0
			ncutx0A <- ncut0

			tcutx0B <- tcut0
			hx0B <- h0
			ncutx0B <- ncut0

			tcutx1A <- tcut1
			hx1A <- h1
			ncutx1A <- ncut1

			tcutx1B <- tcut1
			hx1B <- h1
			ncutx1B <- ncut1
		}
		,
		crossover = {
			no.td0B <- missing(tcutd0B)
			no.sd0B <- missing(sd0B)
			no.hd0B <- missing(hd0B)

			no.td1B <- missing(tcutd1B)
			no.sd1B <- missing(sd1B)
			no.hd1B <- missing(hd1B)

			no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) || 
                                  no.td0B || no.td1B
			if(no.dB)
                          stop("crossover option requires specification of \n'tcutd0B', 'tcutd1B', ('hd0B' or 'sd0B') and " %,%
                               "('hd1B' or 'sd1B').\n")

			if(no.hd0B) hd0B <- stoh(tcutd0B, sd0B)
			ncutd0B <- length(tcutd0B)

			if(no.hd1B) hd1B <- stoh(tcutd1B, sd1B)
			ncutd1B <- length(tcutd1B)

			tcutd0A <- 0
			hd0A <- 1e-7
			ncutd0A <- 1

			tcutd1A <- 0
			hd1A <- 1e-7
			ncutd1A <- 1

			tcutx0A <- tcut0
			hx0A <- h0
			ncutx0A <- ncut0

			tcutx1A <- tcut1
			hx1A <- h1
			ncutx1A <- ncut1

			tcutx0B <- tcut1
			hx0B <- h1
			ncutx0B <- ncut1

			tcutx1B <- tcut0
			hx1B <- h0
			ncutx1B <- ncut0
		}
		,
		mixed = {
			no.td0A <- missing(tcutd0A)
			no.sd0A <- missing(sd0A)
			no.hd0A <- missing(hd0A)

			no.td0B <- missing(tcutd0B)
			no.sd0B <- missing(sd0B)
			no.hd0B <- missing(hd0B)

			no.td1A <- missing(tcutd1A)
			no.sd1A <- missing(sd1A)
			no.hd1A <- missing(hd1A)

			no.td1B <- missing(tcutd1B)
			no.sd1B <- missing(sd1B)
			no.hd1B <- missing(hd1B)

			no.tx0A <- missing(tcutx0A)
			no.sx0A <- missing(sx0A)
			no.hx0A <- missing(hx0A)

			no.tx1A <- missing(tcutx1A)
			no.sx1A <- missing(sx1A)
			no.hx1A <- missing(hx1A)

			no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) || 
                                  no.td0A || no.td1A
			no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) || 
                                  no.td0B || no.td1B
			no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) || 
                                  no.tx0A || no.tx1A

			if(no.dA || no.dB || no.xA) 
                          stop("mixed option requires specification of \n'tcutd0A', 'tcutd1A', ('hd0A' or 'sd0A')," %,%
                               " ('hd1A' or 'sd1A'),\n'tcutd0B', 'tcutd1B', ('hd0B' or 'sd0B')," %,%
                               " ('hd1B' or 'sd1B'),\n'tcutx0A', 'tcutx1A', ('hx0A' or 'sx0A')," %,%
                               " and ('hx1A' or 'sx1A').\n")

			if(no.hd0A) hd0A <- stoh(tcutd0A, sd0A)
			ncutd0A <- length(tcutd0A)

			if(no.hd1A) hd1A <- stoh(tcutd1A, sd1A)
			ncutd1A <- length(tcutd1A)

			if(no.hd0B) hd0B <- stoh(tcutd0B, sd0B)
			ncutd0B <- length(tcutd0B)

			if(no.hd1B) hd1B <- stoh(tcutd1B, sd1B)
			ncutd1B <- length(tcutd1B)

			if(no.hx0A) hx0A <- stoh(tcutx0A, sx0A)
			ncutx0A <- length(tcutx0A)

			if(no.hx1A) hx1A <- stoh(tcutx1A, sx1A)
			ncutx1A <- length(tcutx1A)

			tcutx0B <- tcut1
			hx0B <- h1
			ncutx0B <- ncut1

			tcutx1B <- tcut0
			hx1B <- h0
			ncutx1B <- ncut0
		}
		,
		user = {
			no.td0A <- missing(tcutd0A)
			no.sd0A <- missing(sd0A)
			no.hd0A <- missing(hd0A)

			no.td0B <- missing(tcutd0B)
			no.sd0B <- missing(sd0B)
			no.hd0B <- missing(hd0B)

			no.td1A <- missing(tcutd1A)
			no.sd1A <- missing(sd1A)
			no.hd1A <- missing(hd1A)

			no.td1B <- missing(tcutd1B)
			no.sd1B <- missing(sd1B)
			no.hd1B <- missing(hd1B)

			no.tx0A <- missing(tcutx0A)
			no.sx0A <- missing(sx0A)
			no.hx0A <- missing(hx0A)

			no.tx1A <- missing(tcutx1A)
			no.sx1A <- missing(sx1A)
			no.hx1A <- missing(hx1A)

			no.tx0B <- missing(tcutx0B)
			no.sx0B <- missing(sx0B)
			no.hx0B <- missing(hx0B)

			no.tx1B <- missing(tcutx1B)
			no.sx1B <- missing(sx1B)
			no.hx1B <- missing(hx1B)

			no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) || 
                                  no.td0A || no.td1A
			no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) || 
                                  no.td0B || no.td1B
			no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) || 
                                  no.tx0A || no.tx1A
			no.xB <- (no.sx0B && no.hx0B) || (no.sx1B && no.hx1B) || 
                                  no.tx0B || no.tx1B

			if((no.dA || no.xA) && (no.dB || no.xB)) 
                          stop("user option requires specification of \n('tcutd0A', 'tcutd1A', ('hd0A' or 'sd0A')," %,%
                               " ('hd1A' or 'sd1A') and\n'tcutx0A', 'tcutx1A', ('hx0A' or 'sx0A')," %,%
                               " ('hx1A' or 'sx1A')) or \n('tcutd0B', 'tcutd1B', ('hd0B' or 'sd0B')," %,%
                               " and ('hd1B' or 'sd1B') and\n 'tcutx0B', 'tcutx1B', ('hx0B' or 'sx0B')," %,%
                               " and ('hx1B' or 'sx1B')).\n")

			if(!no.dA) {
				if(no.hd0A) hd0A <- stoh(tcutd0A, sd0A)
				ncutd0A <- length(tcutd0A)

				if(no.hd1A) hd1A <- stoh(tcutd1A, sd1A)
				ncutd1A <- length(tcutd1A)

				if(no.hx0A) hx0A <- stoh(tcutx0A, sx0A)
				ncutx0A <- length(tcutx0A)

				if(no.hx1A) hx1A <- stoh(tcutx1A, sx1A)
				ncutx1A <- length(tcutx1A)

				if(no.dB){
					tcutd0B <- 0
					hd0B <- 1e-7
					ncutd0B <- 1

					tcutd1B <- 0
					hd1B <- 1e-7
					ncutd1B <- 1

					tcutx0B <- tcut0
					hx0B <- h0
					ncutx0B <- ncut0

					tcutx1B <- tcut1
					hx1B <- h1
					ncutx1B <- ncut1
				}
			}
			if(!no.dB){
				if(no.hd0B) hd0B <- stoh(tcutd0B, sd0B)
				ncutd0B <- length(tcutd0B)

				if(no.hd1B) hd1B <- stoh(tcutd1B, sd1B)
				ncutd1B <- length(tcutd1B)

				if(no.hx0B) hx0B <- stoh(tcutx0B, sx0B)
				ncutx0B <- length(tcutx0B)

				if(no.hx1B) hx1B <- stoh(tcutx1B, sx1B)
				ncutx1B <- length(tcutx1B)

				if(no.dA){
					tcutd0A <- 0
					hd0A <- 1e-7
					ncutd0A <- 1

					tcutd1A <- 0
					hd1A <- 1e-7
					ncutd1A <- 1

					tcutx0A <- tcut0
					hx0A <- h0
					ncutx0A <- ncut0

					tcutx1A <- tcut1
					hx1A <- h1
					ncutx1A <- ncut1
				}
			}
		}
		)
        nmax <- max(ncut0,ncut1,nlook)
        cumsum.nppar <- 0
        stat.nms <- NULL
        for(j in 1:nstat){
          prfx <- c("FH-", "SFH-", "R-")[1+wttyp[j]]
          nppar <- c(2,3,1)[1+wttyp[j]]
          ppar.j <- ppar[(cumsum.nppar+1):(cumsum.nppar+nppar)]
          par.string <- glue.2.string(ppar.j, sep="|")
          stat.nms <- c(stat.nms, prfx %,% par.string)
          cumsum.nppar <- cumsum.nppar + nppar
        }
        nbetyp <- length(nbnd.e)
        nbftyp <- length(nbnd.f)
        nsum <- ncut0 + ncut1 + nlook - 2
	ints <- c(nlook,nstat,NGaussQ,ncut0,ncut1,ncutc0,ncutc1,ncutd0A,ncutd0B,ncutd1A,
                  ncutd1B,ncutx0A,ncutx0B,ncutx1A,ncutx1B,gradual,nbnd.e,nbnd.f,
                  nsf.e,nsf.f,do.futility,use.rhaz.fu,spend.info,user.V.end,is.myE,is.myF,spend.info.k,
                  qProp.one.or.Q, sided)        
	ints.nms <- c("nlook","nstat","NGaussQ","ncut0","ncut1","ncutc0",
                      "ncutc1","ncutd0A","ncutd0B","ncutd1A","ncutd1B","ncutx0A","ncutx0B",
                      "ncutx1A","ncutx1B","gradual","nbnd.e." %,% (1:nlook),"nbnd.f." %,% (1:nlook),
                      "nsf.e." %,% (1:nlook),"nsf.f." %,% (1:nlook),"do.futility","use.rhaz.fu",
                      "spend.info","user.V.end","is.myE","is.myF","spend.info.k","qis1orQ","sided")

        dbls <- c(accru, accrat, rho.Efficacy, rho.Futility, crit.e, crit.f)
        ntrial <- floor(accru*accrat)
        dbls.nms <- c("accru", "accrat", "rho.Efficacy." %,%(1:nlook), "rho.Futility." %,% (1:nlook),
                      "crit.e." %,%(1:nlook), "crit.f." %,%(1:nlook))
        
	names(ints) <- ints.nms
	ans <- .C(name = "AsyPwrGSD",
		ints = as.integer(ints),
                dbls = as.double(dbls),
                pttlook = as.double(tlook),
		palphatot = as.double(c(Alpha.Efficacy,Alpha.Futility)),
		lrrf = as.double(log(RR.Futility)),
                bHay = as.double(b.Haybittle),
		ppar = as.double(ppar),
		pgqxw = as.double(c(glegx,glegw)),
		tcut0 = as.double(tcut0),
		h0 = as.double(h0),
		tcut1 = as.double(tcut1),
		h1 = as.double(h1),
		tcutc0 = as.double(tcutc0),
		hc0 = as.double(hc0),
		tcutc1 = as.double(tcutc1),
		hc1 = as.double(hc1),
		tcutd0A = as.double(tcutd0A),
		hd0A = as.double(hd0A),
		tcutd0B = as.double(tcutd0B),
		hd0B = as.double(hd0B),
		tcutd1A = as.double(tcutd1A),
		hd1A = as.double(hd1A),
		tcutd1B = as.double(tcutd1B),
		hd1B = as.double(hd1B),
		tcutx0A = as.double(tcutx0A),
		hx0A = as.double(hx0A),
		tcutx0B = as.double(tcutx0B),
		hx0B = as.double(hx0B),
		tcutx1A = as.double(tcutx1A),
		hx1A = as.double(hx1A),
		tcutx1B = as.double(tcutx1B),
		hx1B = as.double(hx1B),
                wttyp = as.integer(wttyp),
                V.end = as.double(V.end),
                pinffrac = double(nstat*nlook),
                pinffrac.ii= double(nstat*nlook),
		pbounds = as.double(c(rep(b.e, nstat), rep(b.f, nstat))),  
		mufu = double(nstat*nlook),
		mu = double(nstat*nlook),
		palpha0vec = double(2*nstat*nlook),
		palpha1vec = double(2*nstat*nlook),
		RR = double(3*nsum),
                pnnjmp = integer(1),
                E.NT = double(nstat*nsum),
                Var = double(nstat*nsum),
                Eta = double(nstat*nsum),
                betabdry = double(2*nstat*nlook),
		PACKAGE = "PwrGSD")
	detail <- ans
        detail$ntrial <- ntrial
	Pwr <- t(matrix(detail$palpha1vec,nlook,2*nstat))
	dErrII <- Pwr[nstat + (1:nstat),,drop=FALSE]
	Pwr <- Pwr[1:nstat,,drop=FALSE]
	dimnames(dErrII) <- dimnames(Pwr) <- list(stat.nms,tlook)
        njmp <- detail$pnnjmp
	detail$pinffrac <- matrix(detail$pinffrac, nlook, nstat)
        detail$pinffrac.ii <- matrix(detail$pinffrac.ii, nlook, nstat)
	detail$pbounds <- matrix(detail$pbounds,nlook,2*nstat)
        detail$betabdry <- matrix(detail$betabdry, nlook, 2*nstat)
	detail$mu <- matrix(detail$mu, nlook, nstat)
	detail$mufu <- matrix(detail$mufu, nlook, nstat)
        detail$Var <- matrix(detail$Var[1:(njmp*nstat)], njmp, nstat)
        detail$Eta <- matrix(detail$Eta[1:(njmp*nstat)], njmp, nstat)
	detail$palpha0vec <- matrix(detail$palpha0vec, nlook, 2*nstat)
	detail$palpha1vec <- matrix(detail$palpha1vec, nlook, 2*nstat)
	detail$RR <- matrix(detail$RR[1:(3*njmp)],njmp,3)
        detail$E.NT <- 4*matrix(detail$E.NT[1:(njmp*nstat)], njmp, nstat)
	dimnames(detail$pinffrac) <-
 	dimnames(detail$pinffrac.ii) <-
	dimnames(detail$mu) <- list(tlook, stat.nms)
	dimnames(detail$mufu) <- list(tlook, stat.nms)
	dimnames(detail$Var) <-
      	dimnames(detail$Eta) <-
        dimnames(detail$E.NT) <- list(round(detail$RR[,1],7), stat.nms)
	dimnames(detail$pbounds) <- 
        dimnames(detail$betabdry) <- 
	dimnames(detail$palpha0vec) <- 
	dimnames(detail$palpha1vec) <- list(tlook, c(outer(stat.nms,c("-e","-f"),
					    FUN="%,%")))
	dimnames(detail$RR) <- list(rep("",njmp), c("t","htlde0","RR"))
	names(detail$ints) <- ints.nms
        names(detail$dbls) <- dbls.nms
	out <- list(dPower = Pwr, dErrorII = dErrII, detail = detail,call=.call.)
	class(out) <- "PwrGSD"
	out
}

"AsyPwrGSDcall" <-
function(accru,accrat,tlook, tcut0=NULL,h0=NULL,s0=NULL,tcut1=NULL, 
         rhaz=NULL, h1=NULL, s1=NULL, tcutc0=NULL, hc0 = NULL, sc0=NULL, 
         tcutc1=NULL, hc1 = NULL, sc1 = NULL, tcutd0A=NULL, hd0A=NULL, 
         sd0A=NULL, tcutd0B=NULL, hd0B=NULL, sd0B=NULL, tcutd1A=NULL, hd1A=NULL, 
         sd1A=NULL, tcutd1B=NULL, hd1B=NULL, sd1B=NULL, tcutx0A=NULL, hx0A=NULL, 
         sx0A=NULL, tcutx0B=NULL, hx0B=NULL, sx0B=NULL, tcutx1A=NULL, hx1A=NULL, 
         sx1A=NULL, tcutx1B=NULL, hx1B=NULL, sx1B=NULL, 
         noncompliance=c("none", "crossover", "mixed", "user"), 
         gradual = FALSE, WtFun=c("FH","SFH","Ramp"), ppar=c(0, 0),
         Boundary.Efficacy = c("Lan-Demets", "Haybittle"), 
         Boundary.Futility = c("Lan-Demets", "Haybittle"),
         Alpha.Efficacy=0,  Alpha.Futility=0, sided=c("2",">","<"),
         RR.Futility=NULL, Spend.Info=c("Variance", "Events", "Hybrid(k)", "Calendar"),
         Spending.Efficacy = c("Obrien-Fleming", "Pocock", "Power"), 
         Spending.Futility = c("Obrien-Fleming", "Pocock", "Power"), 
         rho.Efficacy=0, rho.Futility=0, b.Haybittle = 3.0,
         rho.Eff.SC = 1.0, rho.Fut.SC = 1.0, V.end = NULL)  
{
    match.call()
}

"print.PwrGSD" <- function (x, ...)
{
    print(x$call)
    nlook <- x$detail$ints["nlook"]
    nstat <- x$detail$ints["nstat"]
    pwr <- c(x$dPower %*% rep(1, nlook))
    dErrII <- 0*x$dPower
    do.futility <- (x$detail$ints["do.futility"]==1)
    if(do.futility) dErrII <- x$dErrorII
    ed <- c(apply(x$dPower + dErrII, 1, 
	  FUN = function(a, b, do.futility) 
		{
        	     n <- length(a)
		     y <- a
        	     if(!do.futility) y <- c(a[-n], 1 - sum(a[-n]))
        	     sum(y * b)
    		}, b = eval(x$call$tlook), do.futility=do.futility))
    eII <- c(dErrII %*% rep(1, nlook))
    stat.nms <- dimnames(x$dPower)[[1]]
    ans <- cbind(pwr, eII, ed)
    dimnames(ans) <- list(stat.nms, c("Power", "Type.II.Err", "ExpDur"))
    print(ans)
    invisible(x)
}

"summary.PwrGSD" <- 
function(object, ...)
{
    nlook <- object$detail$ints[1]
    nstat <- object$detail$ints[2]
    pwr <- c(object$dPower %*% rep(1,nlook))
    dErrII <- 0*object$dPower
    do.futility <- (x$detail$ints["do.futility"]==1)
    if(do.futility) dErrII <- object$dErrorII
    ed <- c(apply(object$dPower + dErrII, 1, 
	  FUN=function(x, v, do.futility)
	      { 
	  	  n <- length(x)
		  y <- x
	  	  if(!do.futility) y <- c(x[-n],1-sum(x[-n]))
	  	  sum(y*v)
	      }, v = object$call$tlook, do.futility=do.futility))
	eII <- c(dErrII %*% rep(1,nlook))
	stat.nms <- dimnames(object$dPower)[[1]]
	ans <- cbind(pwr,eII,ed)
	dimnames(ans) <- list(stat.nms, c("Power","Type.II.Err","ExpDur"))
	out <- object
	out$Tbl <- ans
	out
}

"print.cpd.PwrGSD" <- 
function(x, ...)
{
        print(x$call)
        sapply(x$Elements, FUN = print.PwrGSD)
}

"summary.cpd.PwrGSD" <- 
function(object, ...)
{
        sapply(object$Elements, FUN = summary.PwrGSD)
}

"cpd.PwrGSD" <-
function(descr)
{
  res <- list()
  d.descr <- dim(descr)
  n <- d.descr[1]
  p <- d.descr[2]
  length(res) <- n
  if(!("index" %in% names(descr)))
    descr$index <- 1:n
  
  ans <- list(date=format(Sys.time(), "%Y.%m.%d.%H:%M:%S"), Elements=res, descr=descr)
  class(ans) <- "cpd.PwrGSD"
  ans
}

"Elements" <-
function(object, subset, na.action=na.pass)
{
  m <- .call. <- match.call()
  descr <- object$descr
  n <- dim(descr)[1]
  index <- rep(TRUE, n)
  if(!missing(subset))
    index <- with(descr, eval(.call.$subset))
  
  if(length(index>1)){
    ans <- list(date=object$date, Elements=object$Elements[index], descr=descr[index, ])
    class(ans) <- "cpd.PwrGSD"
  }
  else{
    ans <- object$Elements[index]
    class(ans) <- "PwrGSD"
  }
  ans
}

"as.data.frame.cpd.PwrGSD" <- 
function(x, row.names, optional, ...) 
{
    N <- length(x$Elements)
    idx <- x$descr$index
    det.1 <- x$Elements[[1]]$detail
    det.a <- c(det.1$ntrial, det.1$ints["nlook"], det.1$ints["nstat"])
    names(det.a) <- c("ntrial", "nlook", "nstat")

    nstat <- det.a["nstat"]
    nlook <- det.a["nlook"]
    nms <- outer(1:nstat, 1:nlook, FUN = function(x, y) {
        "nstat[" %,% x %,% "].nlook[" %,% y %,% "]"
    })
    ans <- NULL
    for (k in 1:nstat) {
      tmp.b <- sapply(x$Elements, FUN = function(x, stat) as.boundaries(x, 
                      stat)$table, stat = k)
      d.tmp.b <- dim(tmp.b)
      col.gsbt <- d.tmp.b[1]/nlook
      tmp.b <- matrix(aperm(array(tmp.b, c(nlook, col.gsbt, N))[,c(1,2,5),],c(3,1,2)), N, 3*nlook)
      nstatis <- "stat[" %,% k %,% "]"
      col.nms <- c(t(outer(c("frac.","bFut.","bEff."), nstatis %,% ".nlook[" %,% 1:nlook %,% "]",
                           FUN = "%,%")))
      dimnames(tmp.b) <- list(idx, col.nms)
      tmp.p <- t(sapply(x$Elements, FUN=function(x, stat)c(x$dE[stat,],x$dP[stat,]), stat=k))
      col.nms <- c(t(outer(c("dErrII.","dP."),nstatis %,% ".nlook[" %,% 1:nlook %,% "]",FUN="%,%")))
      dimnames(tmp.p) <- list(idx, col.nms)
      ans <- cbind(ans, tmp.b, tmp.p)
    }

    ans <- as.data.frame(ans)
    ans$index <- x$descr$index
    ans <- merge(x$descr, ans, by = "index")
    ans <- mystack(ans, fu.vars=c("frac","bEff","bFut","dP","dE"))
    ans$stat <- floor(ans$fu/nlook)+1
    ans$nlook <- (ans$fu %% nlook) + 1
    ans$fu <- NULL
    p <- dim(x$descr)[2]
    pp <- dim(ans)[2]
    ans <- ans[,c(1:p, pp-1, pp, (p+1):(pp-2))]
    attr(ans, "detail") <- det.a
    ans
}

"Power" <-
function(object, subset){
  m <- .call. <- match.call()
  
  m1 <- .call.
  m1[[1]] <- as.name("Elements")
  m1 <- eval(m1, sys.parent())
  descr <- m1$descr
  n <- length(m1$Elements)
  nlook <- object$Elements[[1]]$detail$ints["nlook"]
  Pow <- matrix(0, n, 2)
  m1 <- as.data.frame(m1)
  for(k in 1:n)
    Pow[k,] <- rep(1,nlook)%*%as.matrix(m1[nlook*(k-1) + 1:nlook,c("dP", "dE")])

  dimnames(Pow)[[2]] <- c("Power","Type II error")
  
  ans <- list(table=cbind(descr, Pow), call=.call.)
  ans
}

"GrpSeqBnds"<-
function (EfficacyBoundary = LanDemets(alpha=0.05, spending=ObrienFleming),
          FutilityBoundary = LanDemets(alpha=0.10, spending=ObrienFleming),
          frac, frac.ii = NULL, drift = rep(0, length(frac)))
{
    .call. <- match.call()
    is.frac.ii <- TRUE
    if (missing(frac.ii)) {
        frac.ii <- frac
        is.frac.ii <- FALSE
    }
    if (is.frac.ii && length(frac.ii) != length(frac)) 
        stop("Lengths of 'frac' and 'frac.ii' must agree")
    normut <- 8.20953615160139
    psimin <- 1e-15
    nlooks <- length(frac)
    
# determine type of Efficacy Boundary Specification
    mode.E <- "WRONG"
    BE.missing <- missing(EfficacyBoundary)
    if(BE.missing) mode.E <- "NULL"
    if(!BE.missing){
      try.BE <- try(EfficacyBoundary, silent=TRUE)
      if(is.null(attr(try.BE, "class"))) mode.E <- mode(EfficacyBoundary)
      else if(attr(try.BE, "class")=="try-error") mode.E <- "call"
    }
    if (!(mode.E %in% c("call", "numeric", "NULL"))) 
      stop("Argument 'EfficacyBoundary' must be of mode 'call', 'numeric' or 'NULL'")
    is.myE <- FALSE
    n.myE <- 0
    if (mode.E == "numeric") {
      n.myE <- length(EfficacyBoundary)
      is.myE <- TRUE
      my.Efficacy <- EfficacyBoundary
      .call.$EfficacyBoundary <- as.call(expression(LanDemets, alpha=0.05, spending=ObrienFleming))
    }
    do.efficacy <- !BE.missing || is.myE
    if (mode.E == "NULL" && do.efficacy) 
      .call.$EfficacyBoundary <- as.call(expression(LanDemets, alpha=0.05, spending=ObrienFleming))
    
 # determine type of Futility Boundary Specification
    mode.F <- "WRONG"
    BF.missing <- missing(FutilityBoundary)
    if(BF.missing) mode.F <- "NULL"
    if(!BF.missing){
      try.BF <- try(FutilityBoundary, silent=TRUE)
      if(is.null(attr(try.BF, "class"))) mode.F <- mode(FutilityBoundary)
      else if(attr(try.BF, "class")=="try-error") mode.F <- "call"
    }
    if (!(mode.F %in% c("call", "numeric", "NULL"))) 
      stop("Argument 'FutilityBoundary' must be of mode 'call', 'numeric' or 'NULL'")
    is.myF <- FALSE
    n.myF <- 0
    if (mode.F == "numeric") {
      n.myF <- length(FutilityBoundary)
      is.myF <- TRUE
      my.Futility <- FutilityBoundary
      .call.$FutilityBoundary <- as.call(expression(LanDemets, alpha=0.10, spending=ObrienFleming))
    }
    do.futility <- !BF.missing || is.myF
    if(!do.futility) Alpha.Futility <- psimin
    if (mode.F == "NULL")
      .call.$FutilityBoundary <- as.call(expression(LanDemets, alpha=0.10, spending=ObrienFleming))
    
# Parse Boundary Methods:
    ..call. <- .call.
    ..call.[[1]] <- as.name("ParseBoundaryMethods")
    nms.call <- names(..call.)[-1]
    ind.del <- which(!(nms.call %in% c("EfficacyBoundary", "FutilityBoundary")))
    for(k in -sort(-ind.del))
      ..call.[[1+k]] <- NULL
    ..call.$nlooks <- nlooks
    eval(..call.)

# Result:  Now the following are defined within this scope:
#
#   Alpha.Efficacy, Alpha.Futility, nbnd.e, nbnd.f, nsf.e, nsf.f, rho.Efficacy, rho.Futility,
#   b.Haybittle, drift.end, crit.e, crit.f
#

    # for the time being
#   nbnd.e <- nbnd.e[1]
#   nbnd.f <- nbnd.f[1]
#   nsf.e <- nsf.e[1]
#   nsf.f <- nsf.f[1]
#   rho.Efficacy <- rho.Efficacy[1]
#   rho.Futility <- rho.Futility[1]
#   b.Haybittle <- b.Haybittle[1]
#   drift.end <- drift.end[1]
#   crit.e <- crit.e[1]
#   crit.f <- crit.f[1]
    # until I implement changing boundary types
    
    if (!do.efficacy && !do.futility) 
        stop("You must specify one or both of the arguments 'EfficacyBoundary' and 'FutilityBoundary', " %,% 
             "see the documentation")
    
    if (is.myE && (n.myE != nlooks)) 
        stop("User supplied efficacy boundary in 'EfficacyBoundary' must be of the same length as 'frac'")
    if (is.myF && (n.myF != nlooks)) 
        stop("User supplied futility boundary in 'FutilityBoundary' must be of the same length as 'frac'")
    nunq <- nlooks
    glegx <- glegx24
    glegw <- glegw24
    ngqnodes <- length(glegx)
    nbnd.e.sv <- nbnd.e
    nbnd.f.sv <- nbnd.f
    b.e <- alpha.e <- rep(0, nlooks)
    b.f <- alpha.f <- rep(0, nlooks)
    l <- 1
    l.act.e <- 1
    l.act.f <- 1
    fracold <- rep(0, 2)
    fracold.ii <- rep(0, 2)
    fracnew <- frac[1]
    fracnew.ii <- frac.ii[1]
#    sc.drift.factor <- max(frac.ii)^0.5
    sc.drift.factor <- 1
    mu <- drift[1]
    mu.end <- drift.end[1]
    if (nbnd.e[1] == 1 || nbnd.e.sv[1] == 3)
        bold.e <- normut
    if (nbnd.e[1] == 2)
        bold.e <- b.Haybittle[l]
    if (nbnd.f[1] == 1 || nbnd.f.sv[1] == 3)
        bold.f <- -normut
    if (nbnd.f[1] == 2) 
        bold.f <- b.Haybittle[l]
    y.e <- tmp.e <- rep(0, ngqnodes)
    y.f <- tmp.f <- rep(0, ngqnodes)
    x.e <- intgrndx.e <- rep(0, ngqnodes)
    x.f <- intgrndx.f <- rep(0, ngqnodes)
    if (nbnd.e.sv[1] == 3) 
        my.Efficacy <- rep(0, nlooks)
    if (nbnd.f.sv[1] == 3) 
        my.Futility <- rep(0, nlooks)
    while (l <= nlooks) {
        if (nbnd.e.sv[l] == 3) {
            my.Efficacy[l] <-
              .C("StCu2Bnds",
                 pmu = double(2),
                 pfrac = as.double(fracnew.ii),
                 palpha = as.double(c(Alpha.Efficacy,Alpha.Futility)),
                 psided = as.integer(1),
                 prho = as.double(crit.e[l]),
                 pef = as.integer(0),
                 b = double(1),
                 PACKAGE = "PwrGSD")$b
            is.myE <- TRUE
            nbnd.e[l] <- 1
            nsf.e[l] <- 1
        }
        if (nbnd.f.sv[l] == 3) {
            my.Futility[l] <-
              .C("StCu2Bnds",
                 pmu = as.double(c(mu,mu.end)*sc.drift.factor),
                 pfrac = as.double(fracnew.ii),
                 palpha = as.double(c(Alpha.Efficacy, Alpha.Futility)),
                 psided = as.integer(1),
                 prho = as.double(crit.f[l]), 
                 pef = as.integer(1),
                 b = double(1),
                 PACKAGE = "PwrGSD")$b
            is.myF <- TRUE
            nbnd.f[l] <- 1
            nsf.f[l] <- 1
        }
        bnew <- c(normut, -normut)
        if (is.myE || nbnd.e.sv[l] == 3) 
            bnew[1] <- my.Efficacy[l]
        if ((is.myF || nbnd.f.sv[l] == 3) && (1 - fracnew.ii >= 
            1e-06)) 
            bnew[2] <- my.Futility[l]

        ans <- .C("grpseqbnds",
                  do.futility = as.integer(do.futility), 
                  nbnd = as.integer(c(nbnd.e[l], nbnd.f[l])),
                  nsf = as.integer(c(nsf.e[l],nsf.f[l])),
                  rho = as.double(c(rho.Efficacy[l], rho.Futility[l])), 
                  pnlook = as.integer(c(l.act.e, l.act.f)),
                  palphtot = as.double(c(Alpha.Efficacy, Alpha.Futility)),
                  palpha = double(2),
                  psimin = as.double(psimin), 
                  dlact = integer(2),
                  pfracold = as.double(fracold), 
                  pfracnew = as.double(fracnew),
                  pfracold.ii = as.double(fracold.ii), 
                  pfracnew.ii = as.double(fracnew.ii),
                  x = as.double(c(x.e, x.f)),
                  y = as.double(c(y.e, y.f)),
                  tmp = as.double(c(tmp.e, tmp.f)),
                  intgrndx = as.double(c(intgrndx.e, intgrndx.f)), 
                  gqxw = as.double(c(glegx, glegw)),
                  pngqnodes = as.integer(ngqnodes), 
                  mu = as.double(mu),
                  bold = as.double(c(bold.e, bold.f)), 
                  bnew = as.double(bnew),
                  mybounds = as.integer(c(is.myE, is.myF)),
                  PACKAGE = "PwrGSD")
        dlact <- ans$dlact
        if (dlact[1] == 1) {
            if (nbnd.e[l] == 1) 
                b.e[l] <- bold.e <- ifelse(is.myE, bnew[1], ans$bnew[1])
            if (nbnd.e[l] == 2) {
                if (1 - fracnew.ii >= 1e-06) {
                  b.e[l] <- bold.e <- b.Haybittle[l]
                  Alpha.Efficacy <- Alpha.Efficacy - ans$palpha[1]
                }
                else b.e[l] <- ans$bnew[1]
            }
            alpha.e[l] <- ans$palpha[1]
            x.e <- ans$x[1:ngqnodes]
            intgrndx.e <- ans$intgrndx[1:ngqnodes]
            fracold[1] <- ans$pfracnew
            fracold.ii[1] <- ans$pfracnew.ii
            l.act.e <- l.act.e + 1
        }
        else {
            b.e[l] <- bold.e <- ifelse(is.myE, bnew[1], normut)
            alpha.e[l] <- ifelse(is.myE, ans$palpha[1], psimin)
            fracold.ii[1] <- ans$pfracnew.ii
        }
        if (do.futility && dlact[2] == 1) {
            if (nbnd.f[l] == 1) 
                b.f[l] <- bold.f <- ifelse(is.myF && (1 - fracnew.ii >= 1e-06),
                                           bnew[2], ans$bnew[2])
            if (nbnd.f[l] == 2) {
                if (l < nlooks) {
                  b.f[l] <- bold.f <- b.Haybittle[l]
                  Alpha.Futility <- Alpha.Futility - ans$palpha[2]
                }
                else b.f[l] <- ans$bnew[2]
            }
            alpha.f[l] <- ans$palpha[2]
            x.f <- ans$x[ngqnodes + (1:ngqnodes)]
            intgrndx.f <- ans$intgrndx[ngqnodes + (1:ngqnodes)]
            fracold[2] <- ans$pfracnew
            fracold.ii[2] <- ans$pfracnew.ii
            l.act.f <- l.act.f + 1
        }
        if (do.futility && dlact[2] == 0) {
            b.f[l] <- bold.f <- ifelse(is.myF || (l<nlooks), bnew[2], 
                -normut)
            alpha.f[l] <- ifelse(is.myF, ans$palpha[2], psimin)
            fracold.ii[2] <- ans$pfracnew.ii
        }
        l <- l + 1
        fracnew <- frac[l]
        fracnew.ii <- frac.ii[l]
        mu <- drift[l]
    }
    out <- frac.ii
    nms <- "frac"
    if (do.futility) {
        out <- cbind(out, b.f, alpha.f, cumsum(alpha.f))
        nms <- c(nms, "b.f", "alpha.f", "cum-alpha.f")
    }
    out <- cbind(out, b.e, alpha.e, cumsum(alpha.e))
    nms <- c(nms, "b.e", "alpha.e", "cum-alpha.e")
    dimnames(out) <- list(1:nlooks, nms)
    ans <- list(table = out, frac=frac, frac.ii=frac.ii, drift=drift, call = .call.)
    class(ans) <- "boundaries"
    ans
}

"ParseBoundaryMethods" <- 
function(EfficacyBoundary, FutilityBoundary, nlooks)
{
  .call. <- match.call()
  nbnd.e <- nbnd.f <- nsf.e <- nsf.f <- rho.Efficacy <- rho.Futility <-
  b.Haybittle <- drift.end <- crit.e <- crit.f <- from.e <- to.e <- 
  from.f <- to.f <- -9999

  do.efficacy <- TRUE
  do.futility <- TRUE
  if(missing(EfficacyBoundary)){
    do.efficacy <- FALSE
    .call.$EfficacyBoundary <- as.call(expression(LanDemets,alpha=0.05, spending=ObrienFleming))
  }
  if(missing(FutilityBoundary)){
    do.futility <- FALSE
    .call.$FutilityBoundary <- as.call(expression(LanDemets,alpha=0.10, spending=ObrienFleming))
  }

  if(as.character(.call.$EfficacyBoundary)[1]=="c"){
    n.BE.types <- length(.call.$EfficacyBoundary) - 1
    nbnd.e <- nsf.e <- rho.Efficacy <- b.Haybittle <- crit.e <- from.e <- 
    to.e <- rep(-9999, n.BE.types)
    for(k in 1:n.BE.types){
      .call.k <- .call.$EfficacyBoundary[[k+1]]
      BE.type.nm <- as.character(.call.k)[1]
      if(!(BE.type.nm %in% c("LanDemets", "Haybittle", "SC")))
        stop("Boundary type must be one of the following: 'LanDemets()', 'Haybittle()' or 'SC()'")
      nbnd.e[k] <- grep(as.character(.call.k)[1], c("LanDemets", "Haybittle", "SC"))
      if(nbnd.e[k]==1){
        BE.type.class <- class(.call.k)
        if(k==1 && BE.type.class!="call")
          stop("You must set a value for 'alpha' total probability of type I error in the " %,%
	       "call to the 'LanDemets()' method")
        if(BE.type.class=="name") nsf.e[k] <- 1
        if(BE.type.class=="call"){
          BEarg.nms <- names(.call.k)
          if(is.null(BEarg.nms))
            stop("Arguments to 'LanDemets()' must be named")
          if(!("spending" %in% BEarg.nms))
            stop("The 'LanDemets()' method requires you to set the argument 'spending' within " %,%
                 "the call to that method")
          SE.type.nm <- as.character(.call.k$spending)[1]
          if(!(SE.type.nm %in% c("ObrienFleming", "Pocock", "Power")))
            stop("Argument 'spending' must be one of the following: ObrienFleming, Pocock, or Power(rho)")
          nsf.e[k] <- grep(SE.type.nm, c("ObrienFleming", "Pocock", "Power"))
          if(nsf.e[k]==3){
            PE.type.class <- class(.call.k$spending)
            if(PE.type.class!="call")
              stop("The Power spending function requires you to set the argument 'rho' within " %,%
                   "the call")
            rho.Efficacy[k] <- eval(.call.k$spending[[2]])
          }
          if(k==1 && !("alpha" %in% BEarg.nms))
            stop("You must set a value for 'alpha' total probability of type I error in the " %,%
	         "first call to the 'LanDemets()' method")
          if(k==1) Alpha.Efficacy <- eval(.call.k$alpha)
        }
      }
      if(nbnd.e[k]==2){
        BE.type.class <- class(.call.k)
        if(BE.type.class!="call")
          stop("The 'Haybittle()' method requires you to set the argument(s) 'b.Haybittle' and " %,%
               " 'alpha' if it is the first listed method, within the call")
        else{
          BEarg.nms <- names(.call.k)
          if(is.null(BEarg.nms))
            stop("Arguments to 'Haybittle()' must be named")
          if(!("b.Haybittle" %in% BEarg.nms))
            stop("The 'Haybittle()' method requires you to set the argument 'b.Haybittle' within " %,%
	         "the call to that method")
          b.Haybittle[k] <- eval(.call.k$b.Haybittle)
          if(k==1 && !("alpha" %in% BEarg.nms))
            stop("You must set a value for 'alpha' total probability of type I error in the " %,%
	         "first call to the 'Haybittle()' method")
          if(k==1) Alpha.Efficacy <- eval(.call.k$alpha)
        }
      }
      if(nbnd.e[k]==3){
        BE.type.class <- class(.call.k)
        if(BE.type.class!="call")
          stop("The 'SC()' method requires you to set the argument 'crit' within the call to " %,%
               "that method")
        else{
          BEarg.nms <- names(.call.k)
          if(is.null(BEarg.nms))
            stop("Argument to 'SC()' must be named")
          if(!("crit" %in% BEarg.nms))
	    stop("The 'SC()' method requires you to set the argument 'crit' within the call to " %,%
                 "that method")
          crit.e[k] <- eval(.call.k$crit)
          if(k==1 && !("alpha" %in% BEarg.nms))
            stop("You must set a value for 'alpha' total probability of type I error in the " %,%
	         "first call to the 'SC()' method")
          if(k==1) Alpha.Efficacy <- eval(.call.k$alpha)
        }
      }
      if(BE.type.class=="call"){
        is.from <- ("from" %in% BEarg.nms)
        is.to <- ("to" %in% BEarg.nms)
        if(!is.from || !is.to) 
          stop("You must specify both 'from' and 'to' within each efficacy boundary specification")
        from.e[k] <- eval(.call.k$from)
	to.e[k] <- eval(.call.k$to)
      }
    }
    covered <- apply(t(matrix(c(from.e,to.e),n.BE.types,2)), 2, FUN=function(x)x[1]:x[2])
    unl.covered <- unlist(covered)
    covered.all <- all(diff(unl.covered)==1) 
    first.is.1 <- from.e[1]==1
    right.length <- length(unl.covered)==nlooks
    covered.all <- covered.all && first.is.1 && right.length
    if(!covered.all) 
      stop("Your efficacy boundary specification does not cover the full schedule of analyses thus far")
    tmp <- matrix(NA, 5, nlooks)
    for(k in 1:n.BE.types){
      tmp[1,covered[[k]]] <- nbnd.e[k]
      tmp[2,covered[[k]]] <- nsf.e[k]
      tmp[3,covered[[k]]] <- rho.Efficacy[k]
      tmp[4,covered[[k]]] <- b.Haybittle[k]
      tmp[5,covered[[k]]] <- crit.e[k]
    }
    nbnd.e <- tmp[1,]
    nsf.e <- tmp[2,]
    rho.Efficacy <- tmp[3,]
    b.Haybittle <- tmp[4,]
    crit.e <- tmp[5,]
  }
  else{
    BE.type.nm <- as.character(.call.$EfficacyBoundary)[1]
    if(!(BE.type.nm %in% c("LanDemets", "Haybittle", "SC")))
      stop("Boundary type must be one of the following: 'LanDemets()', 'Haybittle()' or 'SC()'")
    nbnd.e <- grep(as.character(.call.$EfficacyBoundary)[1], c("LanDemets", "Haybittle", "SC"))
    if(nbnd.e==1){
      BE.type.class <- class(.call.$EfficacyBoundary)
      if(BE.type.class!="call")
        stop("You must set the argument 'alpha', total probability of type I error, and the " %,%
             "argument 'spending' in the call to the 'LanDemets()' method")
      BEarg.nms <- names(.call.$EfficacyBoundary)
      if(is.null(BEarg.nms))
        stop("Arguments to 'LanDemets()' must be named")
      if(!("spending" %in% BEarg.nms))
        stop("The 'LanDemets()' method requires you to set the argument 'spending' argument within " %,%
             "the call")
      SE.type.nm <- as.character(.call.$EfficacyBoundary$spending)[1]
      if(!(SE.type.nm %in% c("ObrienFleming", "Pocock", "Power")))
        stop("Argument 'spending' must be one of the following: ObrienFleming, Pocock, or Power(rho)")
      nsf.e <- grep(SE.type.nm, c("ObrienFleming", "Pocock", "Power"))
      if(nsf.e==3){
        PE.type.class <- class(.call.$EfficacyBoundary$spending)
        if(PE.type.class!="call")
          stop("The Power spending function requires the argument rho be set to a power")
        rho.Efficacy <- eval(.call.$EfficacyBoundary$spending[[2]])
      }
      if(!("alpha" %in% BEarg.nms))
        stop("You must set a value for 'alpha' total probability of type I error in the " %,%
             "call to the 'LanDemets()' method")
      Alpha.Efficacy <- eval(.call.$EfficacyBoundary$alpha)
    }
    if(nbnd.e==2){
      BE.type.class <- class(.call.$EfficacyBoundary)
      if(BE.type.class!="call")
        stop("The 'Haybittle()' method requires you to set the arguments 'b.Haybittle' and " %,%
               "'alpha' within the call")
      else{
        BEarg.nms <- names(.call.$EfficacyBoundary)
        if(is.null(BEarg.nms))
          stop("Arguments to 'Haybittle()' must be named")
        if(!("b.Haybittle" %in% BEarg.nms))
          stop("The 'Haybittle()' method requires you to set the argument 'b.Haybittle' within " %,%
	       "the call to that method")
	b.Haybittle <- eval(.call.$EfficacyBoundary$b.Haybittle)
        if(!("alpha" %in% BEarg.nms))
          stop("You must set a value for 'alpha' total probability of type I error in the " %,%
               "call to the 'Haybittle()' method")
        Alpha.Efficacy <- eval(.call.$EfficacyBoundary$alpha)
      }
    }
    if(nbnd.e==3){
      BE.type.class <- class(.call.$EfficacyBoundary)
      if(BE.type.class!="call")
        stop("The 'SC()' method requires you to set the argument 'crit' within the call")
      else{
        BEarg.nms <- names(.call.$EfficacyBoundary)
        if(is.null(BEarg.nms))
          stop("Argument to 'SC()' must be named")
        if(!("crit" %in% BEarg.nms))
	  stop("The 'SC()' method requires you to set the argument 'crit' within the call")
        crit.e <- eval(.call.$EfficacyBoundary$crit)
        if(!("alpha" %in% BEarg.nms))
          stop("You must set a value for 'alpha' total probability of type I error in the " %,%
	       "call to the 'SC()' method")
        Alpha.Efficacy <- eval(.call.$EfficacyBoundary$alpha)
      }
    }
    from.e <- 1
    to.e <- nlooks
    nbnd.e <- rep(nbnd.e, nlooks)
    nsf.e <- rep(nsf.e, nlooks)
    rho.Efficacy <- rep(rho.Efficacy, nlooks)
    b.Haybittle <- rep(b.Haybittle, nlooks)
    crit.e <- rep(crit.e, nlooks)
  }

  if(as.character(.call.$FutilityBoundary)[1]=="c"){
    n.BF.types <- length(.call.$FutilityBoundary) - 1
    nbnd.f <- nsf.f <- rho.Futility <- b.Haybittle <- drift.end <- crit.f <- 
    from.f <- to.f <- rep(-9999, n.BF.types)
    for(k in 1:n.BF.types){
      .call.k <- .call.$FutilityBoundary[[k+1]]
      BF.type.nm <- as.character(.call.k)[1]
      if(!(BF.type.nm %in% c("LanDemets", "Haybittle", "SC")))
        stop("Boundary type must be LanDemets, Haybittle or SC")
      nbnd.f[k] <- grep(as.character(.call.k)[1], c("LanDemets", "Haybittle", "SC"))
      if(nbnd.f[k]==1){
        BF.type.class <- class(.call.k)
        if(k==1 && BF.type.class!="call")
          stop("You must set a value for 'alpha' total probability of type II error in the " %,%
	       "specification of the first futility boundary method")
        if(BF.type.class=="name") nsf.f[k] <- 1
        if(BF.type.class=="call"){
          BFarg.nms <- names(.call.k)
          if(is.null(BFarg.nms))
            stop("Arguments to 'LanDemets()' must be named")
          if(!("spending" %in% BFarg.nms))
            stop("The 'LanDemets()' method requires you to set the argument 'spending' within " %,%
                 "the call to that method")
          SF.type.nm <- as.character(.call.k$spending)[1]
          if(!(SF.type.nm %in% c("ObrienFleming", "Pocock", "Power")))
            stop("Argument 'spending' must be one of the following: ObrienFleming, Pocock, or Power(rho)")
          nsf.f[k] <- grep(SF.type.nm, c("ObrienFleming", "Pocock", "Power"))
          if(nsf.f[k]==3){
            PF.type.class <- class(.call.k$spending)
            if(PF.type.class!="call")
              stop("The Power spending function requires you to set the argument 'rho' within " %,%
                   "the call")
            rho.Futility[k] <- eval(.call.k$spending[[2]])
          }
          if(k==1 && !("alpha" %in% BFarg.nms))
            stop("You must set a value for 'alpha' total probability of type II error in the " %,%
	         "specification of the first futility boundary method")
          if(k==1) Alpha.Futility <- eval(.call.k$alpha)
        }
      }
      if(nbnd.f[k]==2){
        BF.type.class <- class(.call.k)
        if(BF.type.class!="call")
          stop("The 'Haybittle()' method requires you to set the argument(s) 'b.Haybittle' and " %,%
               " 'alpha' if it is the first listed method, within the call")          
        else{
          BFarg.nms <- names(.call.k)
          rename <- is.null(BFarg.nms)
          if(rename) names(.call.k)[2] <- "b.Haybittle"
          if(is.null(.call.k$b.Haybittle))
            stop("Method 'Haybittle()' requires the argument 'b.Haybittle'")
          b.Haybittle[k] <- eval(.call.k$b.Haybittle)
          if(k==1 && !("alpha" %in% BFarg.nms))
            stop("You must set a value for 'alpha' total probability of type II error in the " %,%
	         "specification of the first futility boundary method")
          if(k==1) Alpha.Futility <- eval(.call.k$alpha)
        }
      }
      if(nbnd.f[k]==3){
        BF.type.class <- class(.call.k)
        if(BF.type.class!="call")
          stop("Stochastic Curtailment efficacy boundary requires specification of " %,% 
               "the argument 'crit'")
        else{
          BFarg.nms <- names(.call.k)
          if(is.null(BFarg.nms)) 
            stop("Argument to 'SC()' must be named")
          if(is.null(.call.k$crit)||is.null(.call.k$drift.end))
            stop("'SC()' requires the arguments 'crit' and 'drift.end'")
          crit.f[k] <- eval(.call.k$crit)
          drift.end[k] <- eval(.call.k$drift.end)
          if(k==1 && !("alpha" %in% BFarg.nms))
            stop("You must set a value for 'alpha' total probability of type II error in the " %,%
	         "specification of the first futility boundary method")
          if(k==1) Alpha.Futility <- eval(.call.k$alpha)
        }
      }
      if(BF.type.class=="call"){
        is.from <- ("from" %in% BFarg.nms)
        is.to <- ("to" %in% BFarg.nms)
        if(!is.from || !is.to) 
	  stop("You must specify both 'from' and 'to' in each futility boundary specification")
        from.f[k] <- eval(.call.k$from)
        to.f[k] <- eval(.call.k$to)
      }
    }
    covered <- apply(t(matrix(c(from.f,to.f),n.BF.types,2)), 2, FUN=function(x)x[1]:x[2])
    unl.covered <- unlist(covered)
    covered.all <- all(diff(unl.covered)==1) 
    first.is.1 <- from.f[1]==1
    right.length <- length(unl.covered)==nlooks
    covered.all <- covered.all && first.is.1 && right.length
    if(!covered.all) 
      stop("Your futility boundary specification does not cover the full schedule of analyses thus far")

    tmp <- matrix(NA, 5, nlooks)
    for(k in 1:n.BF.types){
      tmp[1,covered[[k]]] <- nbnd.f[k]
      tmp[2,covered[[k]]] <- nsf.f[k]
      tmp[3,covered[[k]]] <- rho.Futility[k]
      tmp[4,covered[[k]]] <- b.Haybittle[k]
      tmp[5,covered[[k]]] <- crit.f[k]
    }
    nbnd.f <- tmp[1,]
    nsf.f <- tmp[2,]
    rho.Futility <- tmp[3,]
    b.Haybittle <- tmp[4,]
    crit.f <- tmp[5,]
  }
  else{
    BF.type.nm <- as.character(.call.$FutilityBoundary)[1]
    if(!(BF.type.nm %in% c("LanDemets", "Haybittle", "SC")))
      stop("Boundary type must be LanDemets, Haybittle or SC")
    nbnd.f <- grep(as.character(.call.$FutilityBoundary)[1], c("LanDemets", "Haybittle", "SC"))
    if(nbnd.f==1){
      BF.type.class <- class(.call.$FutilityBoundary)
      if(BF.type.class!="call")
        stop("You must set the argument 'alpha', total probability of type II error, and the " %,%
             "argument 'spending' in the call to the 'LanDemets()' method")
      BFarg.nms <- names(.call.$FutilityBoundary)
      if(is.null(BFarg.nms))
        stop("Arguments to 'LanDemets()' must be named")
      if(!("spending" %in% BFarg.nms))
        stop("The 'LanDemets()' method requires you to set the argument 'spending' argument within " %,%
             "the call")
      SF.type.nm <- as.character(.call.$FutilityBoundary$spending)[1]
      if(!(SF.type.nm %in% c("ObrienFleming", "Pocock", "Power")))
        stop("Argument 'spending' must be one of the following: ObrienFleming, Pocock, or Power(rho)")
      nsf.f <- grep(SF.type.nm, c("ObrienFleming", "Pocock", "Power"))
      if(nsf.f==3){
        PF.type.class <- class(.call.$FutilityBoundary$spending)
        if(PF.type.class!="call")
          stop("The Power spending function requires the argument rho be set to a power")
        rho.Futility <- eval(.call.$FutilityBoundary$spending[[2]])
      }
      if(!("alpha" %in% BFarg.nms))
        stop("You must set a value for 'alpha' total probability of type II error in the " %,%
             "call to the 'LanDemets()' method")
      Alpha.Futility <- eval(.call.$FutilityBoundary$alpha)
    }
    if(nbnd.f==2){
      BF.type.class <- class(.call.$FutilityBoundary)
      if(BF.type.class!="call")
        stop("The 'Haybittle()' method requires you to set the arguments 'b.Haybittle' and " %,%
               "'alpha' within the call")
      else{
        BFarg.nms <- names(.call.$FutilityBoundary)
        if(is.null(BFarg.nms))
          stop("Arguments to 'Haybittle()' must be named")
        if(!("b.Haybittle" %in% BFarg.nms))
          stop("The 'Haybittle()' method requires you to set the argument 'b.Haybittle' within " %,%
	       "the call to that method")
	b.Haybittle <- eval(.call.$FutilityBoundary$b.Haybittle)
        if(!("alpha" %in% BFarg.nms))
          stop("You must set a value for 'alpha' total probability of type II error in the " %,%
               "call to the 'Haybittle()' method")
        Alpha.Futility <- eval(.call.$FutilityBoundary$alpha)
      }
    }
    if(nbnd.f==3){
      BF.type.class <- class(.call.$FutilityBoundary)
      if(BF.type.class!="call")
        stop("The 'SC()' method requires you to set the arguments 'crit' and 'drift.end' within " %,%
             "the call") 
      else{
        BFarg.nms <- names(.call.$FutilityBoundary)
        if(is.null(BFarg.nms))
          stop("Argument to 'SC()' must be named")
        if(!("crit" %in% BFarg.nms))
	  stop("The 'SC()' method requires you to set the argument 'crit' within the call")
        crit.f <- eval(.call.$FutilityBoundary$crit)
        if(!("drift.end" %in% BFarg.nms))
	  stop("The 'SC()' method requires you to set the argument 'drift.end' within the call")
	drift.end <- eval(.call.$FutilityBoundary$drift.end)
        if(!("alpha" %in% BFarg.nms))
          stop("You must set a value for 'alpha' total probability of type II error in the " %,%
	       "call to the 'SC()' method")
        Alpha.Futility <- eval(.call.$FutilityBoundary$alpha)
      }
    }
    from.f <- 1
    to.f <- nlooks
    nbnd.f <- rep(nbnd.f, nlooks)
    nsf.f <- rep(nsf.f, nlooks)
    rho.Futility <- rep(rho.Futility, nlooks)
    b.Haybittle <- rep(b.Haybittle, nlooks)
    crit.f <- rep(crit.f, nlooks)
  }
  ans <- 
  list(Alpha.Efficacy=Alpha.Efficacy, Alpha.Futility=Alpha.Futility, nbnd.e=nbnd.e, nbnd.f=nbnd.f,
    nsf.e=nsf.e, nsf.f=nsf.f, rho.Efficacy=rho.Efficacy, rho.Futility=rho.Futility,
    b.Haybittle=b.Haybittle, drift.end=drift.end, crit.e=crit.e, crit.f=crit.f, from.e=from.e,
    to.e=to.e, from.f=from.f, to.f=to.f) 
  nms.ans <- names(ans)
  n.ans <- length(nms.ans)
  for(k in 1:n.ans)  assign(nms.ans[k], ans[[k]], parent.frame())
  #
  # note that the intended result of a call to this function is to define and set values 
  # into the objects named above (in the list 'ans'), and to conduct this evaluation within
  # the function that called this one.  Hence, since no arguments need to be returned, we 
  # return zero invisibly, to signify that we have 'finished successfully' just for S&G 
  # if for nothing else.
  #
  invisible(0) 
}

"print.boundaries" <- function(x, ...)
{
  dots <- c(...)
  nms.dots <- names(dots)
  nocall <- FALSE
  if("nocall" %in% nms.dots)
    nocall <- dots["nocall"]  
  ans <- x$table
  if(!nocall){
    cat("call:\n")
    print(x$call)
  }
  print(ans)
  invisible(ans)
}

"plot.boundaries" <- 
function (x, ...) 
{
    .call. <- match.call(expand=TRUE)
    m <- length(.call.) - 2
    nms <- names(.call.)[-(1:2)]
    xtra <- list()
    if(m>0) for(k in 1:m) xtra[[names(.call.)[k]]] <- .call.[[2+k]]
    names(xtra) <- nms
    tbl <- x$table
    d.tbl <- dim(tbl)
    is.fut <- (d.tbl[2] == 7)
    x <- tbl[, "frac"]
    eff <- tbl[, "b.e"]
    fut <- NULL
    if (is.fut) 
        fut <- tbl[, "b.f"]
    yrng <- c(min(c(eff, fut)), max(c(eff, fut)))
    plot.cmd <- as.call(expression(plot))
    plot.cmd$x <- c(0, 1)
    plot.cmd$y <- yrng
    plot.cmd$type <- "n"
    plot.cmd$xlab <- "Information Fraction"
    plot.cmd$ylab <- "z-score"
    if(m>0){
      for(k in 1:m) plot.cmd[[6+k]] <- xtra[[k]]
      names(plot.cmd)[6+(1:m)] <- nms
    }
    eval(plot.cmd)
    lines(x, eff)
    points(x, eff)
    if (is.fut) {
        lines(x, fut)
        points(x, fut)
    }
    invisible(x)
}

"SimGSB" <-
function(object, nsim=1e+05, ...)
{
  UseMethod("SimGSB")
}

"SimGSB.boundaries" <-
function (object, nsim = 1e+05, ...)
{
    .call. <- match.call()
    tab <- object$table
    nms <- dimnames(tab)[[2]]
    .fr. <- eval(object$call$frac)
    .fr.ii <- eval(object$call$frac.ii)
    if(is.null(.fr.ii)) .fr.ii <- .fr.
    n <- length(.fr.)
    do.fu <- length(grep("b.f", nms)) > 0
    a <- rep(-Inf, n)
    al.f <- rep(0, n)
    if (do.fu) {
        .dr. <- eval(object$call$drift)
        a <- .fr.^0.5 * tab[, "b.f"]
        al.f <- tab[, "cum-alpha.f"]
    }
    b <- .fr.^0.5 * tab[, "b.e"]
    al.e <- tab[, "cum-alpha.e"]
    sided <- eval(object$call$sided)
    if (is.null(sided))
        sided <- 1
    dW <- matrix(rnorm(nsim * n), n, nsim)
    W0 <- NULL
    W.i <- 0 * dW[1, ]
    f.old <- 0
    for (i in 1:n) {
        f.new <- .fr.[i]
        df <- f.new - f.old
        W.i <- W.i + df^0.5 * dW[i, ]
        W0 <- rbind(W0, W.i)
        f.old <- f.new
    }
    if (do.fu)
        W1 <- W0 + .dr.
    acc.H0.yet <- rej.H0.yet <- rep(FALSE, nsim)    
    if (sided == 1) {
        cont0.old <- cont1.old <- rep(TRUE, nsim)
        rej.H0 <- acc.H0 <- NULL
        for (i in 1:n) {
            cont0.new <- cont0.old & (W0[i, ] > a[i] & W0[i,
                ] < b[i])
            rej.H0.i <- cont0.old & (!acc.H0.yet) & (W0[i, ] >= b[i])
            rej.H0 <- rbind(rej.H0, rej.H0.i)
            cont0.old <- cont0.new
            if (do.fu) {
                cont1.new <- cont1.old & (W1[i, ] > a[i] & W1[i, ] < b[i])
                acc.H0.i <- cont1.old & (!rej.H0.yet) & (W1[i, ] <= a[i])
                acc.H0 <- rbind(acc.H0, acc.H0.i)
                acc.H0.yet <- acc.H0.yet | acc.H0.i
                cont1.old <- cont1.new
            }
            rej.H0.yet <- rej.H0.yet | rej.H0.i
        }
    }
    if (sided == 2) {
        cont0.old <- cont1.old <- rep(TRUE, nsim)
        rej.H0 <- acc.H0 <- NULL
        for (i in 1:n) {
            cont0.new <- cont0.old & (((W0[i, ] > a[i]) & (W0[i, ] < b[i])) | ((W0[i, ] > -b[i]) & (W0[i, ] < -a[i])))
            rej.H0.i <- cont0.old & (!acc.H0.yet) & ((W0[i, ] >= b[i]) | (W0[i, ] <= -b[i]))
            rej.H0 <- rbind(rej.H0, rej.H0.i)
            cont0.old <- cont0.new
            if (do.fu && a[i] > 0) {
                cont1.new <- cont1.old & (((W1[i, ] > a[i]) & (W1[i, ] < b[i])) | ((W1[i, ] > -b[i]) & (W1[i, ] < -a[i])))
                acc.H0.i <- cont1.old & (!rej.H0.yet) & (W1[i, ] >= -a[i]) & (W1[i, ] <= a[i])
                acc.H0 <- rbind(acc.H0, acc.H0.i)
                acc.H0.yet <- acc.H0.yet | acc.H0.i
                cont1.old <- cont1.new
            }
            rej.H0.yet <- rej.H0.yet | rej.H0.i
        }
    }
    ErrI <- cumsum(c(rej.H0 %*% rep(1/nsim, nsim)))
    if (do.fu)
        ErrII <- cumsum(c(acc.H0 %*% rep(1/nsim, nsim)))
    ans <- cbind(ErrI, al.e)
    nms <- c("eI.est", "eI.act")
    if (do.fu) {
        n.diff <- n - length(ErrII)
        zeros <- rep(0, n.diff)
        ans <- cbind(ans, c(zeros, ErrII), al.f)
        nms <- c(nms, "eII.est", "eII.act")
    }
    dimnames(ans) <- list(signif(.fr.ii, 4), nms)
    ans
}

"SimGSB.PwrGSD" <-
function(object, nsim=1e+05, ...)
{
  object <- as.boundaries(object)
  SimGSB(object)
}

"as.boundaries" <-
function(object, ...){
  UseMethod("as.boundaries")
}

"as.boundaries.boundaries" <-
function(object, ...){
  object
}

"as.boundaries.PwrGSD" <-
function(object, ...)
{
  dots <- c(...)
  nms.dots <- names(dots)
  stat <- 1
  if("stat" %in% nms.dots)
    stat <- dots["stat"]
  nlook <- object$detail$ints["nlook"]
  dofu <- object$detail$ints["do.futility"]
  nstat <- object$detail$ints["nstat"]
  if(dofu==0){
    cl <- as.call(expression(GrpSeqBnds)) 
    cl$frac <- frac <- c(object$detail$pinffrac[,stat])
    if(sum(abs(c(object$detail$pinffrac[,stat])-c(object$detail$pinffrac.ii[,stat]))) > 1e-8)
      cl$frac.ii <- frac.ii <- c(object$detail$pinffrac.ii[,stat])
    else frac.ii <- frac
    cl$EfficacyBoundary <- object$call$EfficacyBoundary
    out <- cbind(object$detail$pinffrac.ii[,stat], object$detail$pbounds[,stat], object$detail$palpha0vec[,stat],
                  cumsum(object$detail$palpha0vec[,stat]))
    dimnames(out) <- list(1:nlook, c("frac","b.e", "alpha.e", "cum-alpha.e"))
    ans <- list(table=out, frac=frac, frac.ii=frac.ii, call=cl)
  }
  if(dofu==1){
    cl <- as.call(expression(GrpSeqBnds))
    cl$frac <- frac <- c(object$detail$pinffrac[,stat])
    if(sum(abs(c(object$detail$pinffrac[,stat])-c(object$detail$pinffrac.ii[,stat]))) > 0)
      cl$frac.ii <- frac.ii <- c(object$detail$pinffrac.ii[,stat])
    else frac.ii <- frac
    cl$EfficacyBoundary <- object$call$EfficacyBoundary
    cl$FutilityBoundary <- object$call$FutilityBoundary
    cl$drift <- drift <- c(object$detail$mufu[,stat])
    if(as.character(cl$FutilityBoundary[[1]])=="SC")
      cl$FutilityBoundary$drift.end <- object$detail$mufu[nlook, stat]
    out <- cbind(object$detail$pinffrac.ii[,stat], object$detail$pbounds[,nstat + stat], object$detail$palpha0vec[,nstat + stat],
                  cumsum(object$detail$palpha0vec[,nstat + stat]), object$detail$pbounds[,stat],
                  object$detail$palpha0vec[,stat], cumsum(object$detail$palpha0vec[,stat]))
    dimnames(out) <- list(1:nlook, c("frac","b.f","alpha.f","cum-alpha.f","b.e", "alpha.e", "cum-alpha.e"))
    ans <- list(table=out, frac=frac, frac.ii=frac.ii, drift=drift, call=cl)
  }
  class(ans) <- "boundaries"
  ans
}

"SimPwrGSD" <-
function(EfficacyBoundary = LanDemets(alpha=0.05,spending=ObrienFleming),
         FutilityBoundary = LanDemets(alpha=0.10,spending=ObrienFleming),
         sided =c("2",">","<"),accru,accrat,tlook,
         tcut0 = NULL,h0 = NULL,s0 = NULL,tcut1 = NULL,rhaz = NULL,
         h1 = NULL,s1 = NULL,tcutc0 = NULL,hc0 = NULL,sc0 = NULL,tcutc1 = NULL,hc1 = NULL,
         sc1 = NULL,tcutd0A = NULL,hd0A = NULL,sd0A = NULL,tcutd0B = NULL,hd0B = NULL,sd0B = NULL,
         tcutd1A = NULL,hd1A = NULL,sd1A = NULL,tcutd1B = NULL,hd1B = NULL,sd1B = NULL,
         tcutx0A = NULL,hx0A = NULL,sx0A = NULL,tcutx0B = NULL,hx0B = NULL,sx0B = NULL,
         tcutx1A = NULL,hx1A = NULL,sx1A = NULL,tcutx1B = NULL,hx1B = NULL,sx1B = NULL,
         noncompliance = c("none","crossover","mixed","user"),gradual = FALSE,
         WtFun = c("FH","SFH","Ramp"),ppar = cbind(c(0,0)),
         Spend.Info=c("Variance","Events","Hybrid(k)","Calendar"),
         RR.Futility = NULL,qProp.one.or.Q = c("one","Q"),Nsim=NULL,
         detail = FALSE,StatType = c("WLR","ISD"))                        
{

    if(!missing(StatType) && any(!(StatType %in% c("WLR","ISD"))))
      stop("Elements of the vector argument 'StatType' be either \"WLR\" or \"ISD\"")
    if(!missing(StatType) && !missing(WtFun) && length(StatType)!=length(WtFun))
      stop("Vector arguments 'StatType' and 'WtFun' must be of the same length")

    if(missing(WtFun)) {
      WtFun <- "FH"
      wttyp <- 0
      ppar <- c(0,0)
      nstat <- 1
    }
    else{
      if(missing(ppar) && any(WtFun!="FH"))
        stop("You must specify parameters for the chosen weight function(s) in the 'ppar' argument.")
      nstat <- length(WtFun)
      wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1
    }
    if(missing(StatType)) StatType <- rep("WLR", nstat)

    stattype <- 1*(StatType == "ISD")

    .call. <- match.call()
    ppar <- t(ppar)
    ndbl <- floor(accru * accrat)
    ndbl <- ndbl + (ndbl%%2)
    n <- ndbl/2
    nlook <- length(tlook)
    if (missing(sided))
        sided <- "2"
    if (!all(sided %in% c("2", ">", "<")))
        stop("Elements of argument 'sided' must be " %,% "equal to '2', '>' or '<'")
    if(length(sided)!=nstat)
        stop("Argument 'sided' must be a vector of length " %,% nstat %,% "containing elements '2', '>' or '<'")
    sided <- c(2, 1, -1)[sapply(sided, c("2", ">", "<"), FUN=grep)]
    
    # determine type of Efficacy Boundary Specification
    mode.E <- "WRONG"
    BE.missing <- missing(EfficacyBoundary)
    if(BE.missing) mode.E <- "NULL"
    if(!BE.missing){
      try.BE <- try(EfficacyBoundary, silent=TRUE)
      if(is.null(attr(try.BE, "class"))) mode.E <- mode(EfficacyBoundary)
      else if(attr(try.BE, "class")=="try-error") mode.E <- "call"
    }
    if (!(mode.E %in% c("call", "numeric", "NULL"))) 
      stop("Argument 'EfficacyBoundary' must be of mode 'call', 'numeric' or 'NULL'")
    is.myE <- FALSE
    n.myE <- 0
    if (mode.E == "numeric") {
      n.myE <- length(EfficacyBoundary)
      is.myE <- TRUE
      my.Efficacy <- EfficacyBoundary
      .call.$EfficacyBoundary <- as.call(expression(LanDemets, alpha=0.05, spending=ObrienFleming))
    }
    do.efficacy <- !BE.missing || is.myE
    if (mode.E == "NULL" && do.efficacy) 
      .call.$EfficacyBoundary <- as.call(expression(LanDemets, alpha=0.05, spending=ObrienFleming))
    
    # determine type of Futility Boundary Specification
    mode.F <- "WRONG"
    BF.missing <- missing(FutilityBoundary)
    if(BF.missing) mode.F <- "NULL"
    if(!BF.missing){
      try.BF <- try(FutilityBoundary, silent=TRUE)
      if(is.null(attr(try.BF, "class"))) mode.F <- mode(FutilityBoundary)
      else if(attr(try.BF, "class")=="try-error") mode.F <- "call"
    }
    if (!(mode.F %in% c("call", "numeric", "NULL"))) 
      stop("Argument 'FutilityBoundary' must be of mode 'call', 'numeric' or 'NULL'")
    is.myF <- FALSE
    n.myF <- 0
    if (mode.F == "numeric") {
      n.myF <- length(FutilityBoundary)
      is.myF <- TRUE
      my.Futility <- FutilityBoundary
      .call.$FutilityBoundary <- as.call(expression(LanDemets, alpha=0.10, spending=ObrienFleming))
    }
    do.futility <- !BF.missing || is.myF
    if(!do.futility) Alpha.Futility <- psimin
    if (mode.F == "NULL" && do.futility)
      .call.$FutilityBoundary <- as.call(expression(LanDemets, alpha=0.10, spending=ObrienFleming))
    
    
    # Parse Boundary Methods:
    ..call. <- .call.
    ..call.[[1]] <- as.name("ParseBoundaryMethods")
    nms.call <- names(..call.)[-1]
    ind.del <- which(!(nms.call %in% c("EfficacyBoundary", "FutilityBoundary")))
    for(k in -sort(-ind.del))
      ..call.[[1+k]] <- NULL
    ..call.$nlooks <- nlook
    eval(..call.)
    
    # Result:  Now the following are defined within this scope:
    #
    #   Alpha.Efficacy, Alpha.Futility, nbnd.e, nbnd.f, nsf.e, nsf.f, rho.Efficacy, rho.Futility,
    #   b.Haybittle, drift.end, crit.e, crit.f
    #
        
    if (!do.efficacy && !do.futility) 
      stop("You must specify one or both of the arguments 'EfficacyBoundary' and 'FutilityBoundary', " %,% 
           "see the documentation")
    
    b.e <- rep(0, nlook)
    if (is.myE){
      if(n.myE != nlook)
        stop("User supplied efficacy boundary in 'EfficacyBoundary' must be of the same length as 'frac'")
      b.e <- my.Efficacy
    }

    b.f <- rep(0, nlook)
    if (is.myF){
      if(n.myF != nlook) 
        stop("User supplied futility boundary in 'FutilityBoundary' must be of the same length as 'frac'")
      b.f <- my.Futility
    }
    
#    if(missing(Alpha.Efficacy)) stop("Missing argument 'Alpha.Efficacy', the total type I error")
#    Alpha.Efficacy <- Alpha.Efficacy/2^(sided == 2)

#    if (missing(Boundary.Efficacy))
#       Boundary.Efficacy <- NULL
#    mode.E <- mode(Boundary.Efficacy)
#    if(!(mode.E %in% c("character", "numeric", "NULL")))
#      stop("Argument 'Boundary.Efficacy' must be of character, numeric or NULL mode")
#    b.e <- rep(0,nlook)
#    is.myE <- FALSE
#    if(mode.E=="numeric"){
#      n.myE <- length(Boundary.Efficacy)
#      if(n.myE != nlook)
#        stop("User supplied efficacy boundary in 'Boundary.Efficacy' must be of the same length as 'frac'")
#      is.myE <- TRUE
#      b.e <- Boundary.Efficacy
#      Boundary.Efficacy <- "Lan-Demets"
#    }
#    if(mode.E=="NULL") Boundary.Efficacy <- "Lan-Demets"    
#    nbnd.e <- grep(Boundary.Efficacy, c("Lan-Demets", "Haybittle", "SC"))
#    if(nbnd.e==3 && missing(rho.Eff.SC))
#      stop("Argument 'rho.Eff.SC' (total type I error for Stochastic Curtailment criterion) must be specified")
#    
#    if (missing(Spending.Efficacy))
#        Spending.Efficacy <- "Obrien-Fleming"
#    if(Spending.Efficacy=="Power" && missing(rho.Efficacy))
#      stop("Argument 'rho.Efficacy' must be supplied")
#    nsf.e <- grep(Spending.Efficacy, c("Obrien-Fleming", "Pocock", "Power"))

    
#    if(missing(Boundary.Futility))
#      Boundary.Futility <- NULL
#    mode.F <- mode(Boundary.Futility)
#    if(!(mode.F %in% c("character", "numeric", "NULL")))
#      stop("Argument 'Boundary.Futility' must be of character, numeric or NULL mode")
#    b.f <- rep(0,nlook)
#    is.myF <- FALSE
#    if(mode.F=="numeric"){
#      n.myF <- length(Boundary.Futility)
#      if(n.myF != nlook)
#        stop("User supplied futility boundary in 'Boundary.Futility' must be of the same length as 'frac'")
#      is.myF <- TRUE
#      b.f <- Boundary.Futility
#      Boundary.Futility <- "Lan-Demets"
#    }

#    dofu <- 1
#    if(missing(Alpha.Futility) && mode.F != "NULL" && mode.F!="numeric")
#      stop("You specified a futility boundary construction method in the argument 'Boundary.Futility' that requires" %,%
#           " a probability of total type II error, 'Alpha.Futility'")
#    if (missing(Alpha.Futility) && mode.F == "NULL" && (is.myF==0)){
#        dofu <- 0
#        Alpha.Futility <- psimin
#    }
    
#    if(dofu==1 && any(stattype==1))
#      stop("Futility Boundaries for Integrated Survival Statistic Not Currently Supported")

#    if(mode.F=="NULL") Boundary.Futility <- "Lan-Demets"
#    nbnd.f <- grep(Boundary.Futility, c("Lan-Demets", "Haybittle", "SC"))
#    if(nbnd.f==3 && missing(rho.Fut.SC))
#      stop("Argument 'rho.Fut.SC' (total type II error for Stochastic Curtailment criterion) must be specified")

#    if(nbnd.f == 2) stop("Haybittle futility boundary makes no sense so that " %,%
#                         "understandibly, it is not supported")
#        
#    if (missing(Spending.Futility)) Spending.Futility <- "Obrien-Fleming"
#    if(Spending.Futility=="Power" && missing(rho.Futility))
#      stop("Argument 'rho.Futility' must be supplied")
#     nsf.f <- grep(Spending.Futility, c("Obrien-Fleming", "Pocock", "Power"))

    no.SpndInfo <- missing(Spend.Info)
    if(no.SpndInfo) Spend.Info <- "Variance"
    if(!no.SpndInfo){
      Spend.Info <- as.character(.call.$Spend.Info)
      if(!(Spend.Info[1] %in% c("Variance", "Events", "Hybrid", "Calendar")))
        stop("Argument 'Spend.Info' is an expression of the form 'Variance', 'Events', 'Hybrid(k)', or 'Calendar'")
    }
    spend.info <- grep(Spend.Info[1], c("Variance", "Events", "Hybrid", "Calendar")) - 1
    spend.info.k <- 0
    if(spend.info==2 && length(Spend.Info)>1) spend.info.k <- as.integer(as.numeric(Spend.Info[2])) - 1

    if(missing(qProp.one.or.Q))
      qProp.one.or.Q <- 0
    else{
      if(!(qProp.one.or.Q %in% c("one","Q")))
        stop("Argument 'qProp.one.or.Q' must be either \"one\" or \"Q\"")
      qProp.one.or.Q <- grep(qProp.one.or.Q, c("one", "Q")) - 1
    }

    stoh <- function(tcut, s) {
        ncut <- length(tcut)
        Sold <- c(1, s[-(ncut - 1)])
        dt <- diff(tcut)
        h <- log.0(s/Sold)/dt
        h
    }
    no.t0 <- missing(tcut0)
    no.s0 <- missing(s0)
    no.h0 <- missing(h0)
    if ((no.s0 && no.h0) || no.t0)
        stop("Must specify 'tcut0' and ('h0' or 's0').")
    if (no.h0)
        h0 <- stoh(tcut0, s0)
    ncut0 <- length(tcut0)
    
    no.t1 <- missing(tcut1)
    no.s1 <- missing(s1)
    no.h1 <- missing(h1)
    no.rhaz <- missing(rhaz)
    if ((no.s1 && no.rhaz && no.h1) || no.t1)
        stop("Must specify 'tcut1' and ('rhaz' or 'h1' or 's1').")
    if (!no.s1)
        h1 <- stoh(tcut1, s1)
    if (!no.rhaz) {
        tcut.01 <- support(c(tcut0, tcut1))
        h0.01 <- lookup(tcut0, h0, tcut.01)$y
        rhaz.01 <- lookup(tcut1, rhaz, tcut.01)$y
        tcut1 <- tcut.01
        h1 <- rhaz.01 * h0.01
    }
    ncut1 <- length(tcut1)

    use.rhaz.fu <- 0
    if(missing(RR.Futility)) {
      RR.Futility <- rep(1,ncut0)
      if(do.futility==1) use.rhaz.fu <- 1
    }
    if(length(RR.Futility)<ncut0) RR.Futility <- rep(RR.Futility[1], ncut0)
    
    no.tc0 <- missing(tcutc0)
    no.sc0 <- missing(sc0)
    no.hc0 <- missing(hc0)
    if ((no.sc0 && no.hc0) || no.tc0)
        stop("Must specify 'tcutc0' and ('hc0' or 'sc0').")
    if (no.hc0)
        hc0 <- stoh(tcutc0, sc0)
    ncutc0 <- length(tcutc0)
    no.tc1 <- missing(tcutc1)
    no.sc1 <- missing(sc1)
    no.hc1 <- missing(hc1)
    if ((no.sc1 && no.hc1) || no.tc1)
        stop("Must specify 'tcutc1' and ('hc1' or 'sc1').")
    if (no.hc1)
        hc1 <- stoh(tcutc1, sc1)
    ncutc1 <- length(tcutc1)
    noncompliance <- as.character(.call.$noncompliance)
    no.noncomp <- (length(noncompliance) == 0)
    if (no.noncomp)
        noncompliance <- "none"
    switch(noncompliance, none = {
        tcutd0A <- 0
        hd0A <- 0
        ncutd0A <- 1
        tcutd0B <- 0
        hd0B <- 0
        ncutd0B <- 1
        tcutd1A <- 0
        hd1A <- 0
        ncutd1A <- 1
        tcutd1B <- 0
        hd1B <- 0
        ncutd1B <- 1
        tcutx0A <- tcut0
        hx0A <- h0
        ncutx0A <- ncut0
        tcutx0B <- tcut0
        hx0B <- h0
        ncutx0B <- ncut0
        tcutx1A <- tcut1
        hx1A <- h1
        ncutx1A <- ncut1
        tcutx1B <- tcut1
        hx1B <- h1
        ncutx1B <- ncut1
    }, crossover = {
        no.td0B <- missing(tcutd0B)
        no.sd0B <- missing(sd0B)
        no.hd0B <- missing(hd0B)
        no.td1B <- missing(tcutd1B)
        no.sd1B <- missing(sd1B)
        no.hd1B <- missing(hd1B)
        no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) ||
            no.td0B || no.td1B
        if (no.dB)
            stop("crossover option requires specification of\n 'tcutd0B', 'tcutd1B', ('hd0B' or " %,%
                "'sd0B') and ('hd1B' or 'sd1B').\n")
        if (no.hd0B)
            hd0B <- stoh(tcutd0B, sd0B)
        ncutd0B <- length(tcutd0B)
        if (no.hd1B)
            hd1B <- stoh(tcutd1B, sd1B)
        ncutd1B <- length(tcutd1B)
        tcutd0A <- 0
        hd0A <- 0
        ncutd0A <- 1
        tcutd1A <- 0
        hd1A <- 0
        ncutd1A <- 1
        tcutx0A <- tcut0
        hx0A <- h0
        ncutx0A <- ncut0
        tcutx1A <- tcut1
        hx1A <- h1
        ncutx1A <- ncut1
        tcutx0B <- tcut1
        hx0B <- h1
        ncutx0B <- ncut1
        tcutx1B <- tcut0
        hx1B <- h0
        ncutx1B <- ncut0
    }, mixed = {
        no.td0A <- missing(tcutd0A)
        no.sd0A <- missing(sd0A)
        no.hd0A <- missing(hd0A)
        no.td0B <- missing(tcutd0B)
        no.sd0B <- missing(sd0B)
        no.hd0B <- missing(hd0B)
        no.td1A <- missing(tcutd1A)
        no.sd1A <- missing(sd1A)
        no.hd1A <- missing(hd1A)
        no.td1B <- missing(tcutd1B)
        no.sd1B <- missing(sd1B)
        no.hd1B <- missing(hd1B)
        no.tx0A <- missing(tcutx0A)
        no.sx0A <- missing(sx0A)
        no.hx0A <- missing(hx0A)
        no.tx1A <- missing(tcutx1A)
        no.sx1A <- missing(sx1A)
        no.hx1A <- missing(hx1A)
        no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) ||
            no.td0A || no.td1A
        no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) ||
            no.td0B || no.td1B
        no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) ||
            no.tx0A || no.tx1A
        if (no.dA || no.dB || no.xA)
            stop("mixed option requires specification of \n'tcutd0A', 'tcutd1A', ('hd0A' or " %,%
                 "'sd0A'), ('hd1A' or 'sd1A'),\n'tcutd0B', 'tcutd1B', ('hd0B' or " %,%
                 "'sd0B'), ('hd1B' or 'sd1B'),\n'tcutx0A', 'tcutx1A', ('hx0A' or " %,%
                 "'sx0A'), and ('hx1A' or 'sx1A').\n")
        if (no.hd0A)
            hd0A <- stoh(tcutd0A, sd0A)
        ncutd0A <- length(tcutd0A)
        if (no.hd1A)
            hd1A <- stoh(tcutd1A, sd1A)
        ncutd1A <- length(tcutd1A)
        if (no.hd0B)
            hd0B <- stoh(tcutd0B, sd0B)
        ncutd0B <- length(tcutd0B)
        if (no.hd1B)
            hd1B <- stoh(tcutd1B, sd1B)
        ncutd1B <- length(tcutd1B)
        if (no.hx0A)
            hx0A <- stoh(tcutx0A, sx0A)
        ncutx0A <- length(tcutx0A)
        if (no.hx1A)
            hx1A <- stoh(tcutx1A, sx1A)
        ncutx1A <- length(tcutx1A)
        tcutx0B <- tcut1
        hx0B <- h1
        ncutx0B <- ncut1
        tcutx1B <- tcut0
        hx1B <- h0
        ncutx1B <- ncut0
    }, user = {
        no.td0A <- missing(tcutd0A)
        no.sd0A <- missing(sd0A)
        no.hd0A <- missing(hd0A)
        no.td0B <- missing(tcutd0B)
        no.sd0B <- missing(sd0B)
        no.hd0B <- missing(hd0B)
        no.td1A <- missing(tcutd1A)
        no.sd1A <- missing(sd1A)
        no.hd1A <- missing(hd1A)
        no.td1B <- missing(tcutd1B)
        no.sd1B <- missing(sd1B)
        no.hd1B <- missing(hd1B)
        no.tx0A <- missing(tcutx0A)
        no.sx0A <- missing(sx0A)
        no.hx0A <- missing(hx0A)
        no.tx1A <- missing(tcutx1A)
        no.sx1A <- missing(sx1A)
        no.hx1A <- missing(hx1A)
        no.tx0B <- missing(tcutx0B)
        no.sx0B <- missing(sx0B)
        no.hx0B <- missing(hx0B)
        no.tx1B <- missing(tcutx1B)
        no.sx1B <- missing(sx1B)
        no.hx1B <- missing(hx1B)
        no.dA <- (no.sd0A && no.hd0A) || (no.sd1A && no.hd1A) ||
            no.td0A || no.td1A
        no.dB <- (no.sd0B && no.hd0B) || (no.sd1B && no.hd1B) ||
            no.td0B || no.td1B
        no.xA <- (no.sx0A && no.hx0A) || (no.sx1A && no.hx1A) ||
            no.tx0A || no.tx1A
        no.xB <- (no.sx0B && no.hx0B) || (no.sx1B && no.hx1B) ||
            no.tx0B || no.tx1B
        if ((no.dA || no.xA) && (no.dB || no.xB))
            stop("user option requires specification of \n('tcutd0A', 'tcutd1A', ('hd0A' or " %,%
                 "'sd0A'), ('hd1A' or 'sd1A') and\n'tcutx0A', 'tcutx1A', ('hx0A' or " %,%
                 "'sx0A'), ('hx1A' or 'sx1A')) or \n('tcutd0B', 'tcutd1B', ('hd0B' or " %,%
                 "'sd0B'), and ('hd1B' or 'sd1B') and\n 'tcutx0B', 'tcutx1B', ('hx0B' or " %,%
                 "'sx0B'), and ('hx1B' or 'sx1B')).\n")
        if (!no.dA) {
            if (no.hd0A)
                hd0A <- stoh(tcutd0A, sd0A)
            ncutd0A <- length(tcutd0A)
            if (no.hd1A)
                hd1A <- stoh(tcutd1A, sd1A)
            ncutd1A <- length(tcutd1A)
            if (no.hx0A)
                hx0A <- stoh(tcutx0A, sx0A)
            ncutx0A <- length(tcutx0A)
            if (no.hx1A)
                hx1A <- stoh(tcutx1A, sx1A)
            ncutx1A <- length(tcutx1A)
            if (no.dB) {
                tcutd0B <- 0
                hd0B <- 0
                ncutd0B <- 1
                tcutd1B <- 0
                hd1B <- 0
                ncutd1B <- 1
                tcutx0B <- tcut0
                hx0B <- h0
                ncutx0B <- ncut0
                tcutx1B <- tcut1
                hx1B <- h1
                ncutx1B <- ncut1
            }
        }
        if (!no.dB) {
            if (no.hd0B)
                hd0B <- stoh(tcutd0B, sd0B)
            ncutd0B <- length(tcutd0B)
            if (no.hd1B)
                hd1B <- stoh(tcutd1B, sd1B)
            ncutd1B <- length(tcutd1B)
            if (no.hx0B)
                hx0B <- stoh(tcutx0B, sx0B)
            ncutx0B <- length(tcutx0B)
            if (no.hx1B)
                hx1B <- stoh(tcutx1B, sx1B)
            ncutx1B <- length(tcutx1B)
            if (no.dA) {
                tcutd0A <- 0
                hd0A <- 0
                ncutd0A <- 1
                tcutd1A <- 0
                hd1A <- 0
                ncutd1A <- 1
                tcutx0A <- tcut0
                hx0A <- h0
                ncutx0A <- ncut0
                tcutx1A <- tcut1
                hx1A <- h1
                ncutx1A <- ncut1
            }
        }
    })
    nunq <- length(unique(sort(c(tcut0, tcut1, tcutc0, tcutc1, tcutd0A, tcutd0B, tcutd1A, tcutd1B, tcutx0A,
                                 tcutx0B, tcutx1A, tcutx1B))))
    glegx <- glegx24
    glegw <- glegw24
    
    NGaussQ <- length(glegx)
    
    cumsum.nppar <- 0
    stat.nms <- NULL
    for(j in 1:nstat){
      stat.<- c("WLR[","ISD[")[1 + stattype[j]]
      wt. <- c("FH(", "SFH(", "R(")[1+wttyp[j]]
      nppar <- c(2,3,1)[1+wttyp[j]]
      ppar.j <- ppar[(cumsum.nppar+1):(cumsum.nppar+nppar)]
      par.string <- glue.2.string(ppar.j, sep=",")
      stat.nms <- c(stat.nms, stat. %,% wt. %,% par.string %,% ")]")
      cumsum.nppar <- cumsum.nppar + nppar
    }
    
    nbetyp <- length(nbnd.e)
    nbftyp <- length(nbnd.f)
    ints <- c(nlook,nstat,NGaussQ,ncut0,ncut1,ncutc0,ncutc1,ncutd0A,ncutd0B,ncutd1A,
              ncutd1B,ncutx0A,ncutx0B,ncutx1A,ncutx1B,gradual,nbnd.e,nbnd.f,
              nsf.e,nsf.f,do.futility,use.rhaz.fu,spend.info,Nsim,is.myE,is.myF,spend.info.k,
              qProp.one.or.Q,sided)

    ints.nms <- c("nlook","nstat","NGaussQ","ncut0","ncut1","ncutc0","ncutc1",
                  "ncutd0A","ncutd0B","ncutd1A","ncutd1B","ncutx0A","ncutx0B","ncutx1A",
                  "ncutx1B","gradual","nbnd.e."%,%(1:nlook),"nbnd.f."%,%(1:nlook),
                  "nsf.e."%,%(1:nlook),"nsf.f."%,%(1:nlook),"do.futility",
                  "use.rhaz.fu","spend.info","Nsim","is.myE","is.myF","spend.info.k",
                  "qis1orQ", "sided-" %,% stat.nms)
    
    dbls <- c(accru,accrat, rho.Efficacy, rho.Futility, crit.e, crit.f)
    
    dbls.nms <- c("accru","accrat", "rho.Efficacy."%,%(1:nlook), "rho.Futility."%,%(1:nlook),
                  "crit.e."%,%(1:nlook), "crit.f."%,%(1:nlook))
    
    logRR.F <- log(RR.Futility)
    
    ans <- .C(name = "SimPwrGSD", 
	ints = as.integer(ints),
        dbls = as.double(dbls),
        pttlook = as.double(tlook), 
	palphatot = as.double(c(Alpha.Efficacy,Alpha.Futility)), 
	lrrf = as.double(logRR.F),
	bHay = as.double(b.Haybittle),
        stattype = as.integer(stattype),
        wttyp = as.integer(wttyp),
	ppar = as.double(ppar), 
	pgqxw = as.double(c(glegx,glegw)),
	tcut0 = as.double(tcut0), 
	h0 = as.double(h0),
        tcut1 = as.double(tcut1), 
	h1 = as.double(h1), 
	tcutc0 = as.double(tcutc0),
        hc0 = as.double(hc0), 
	tcutc1 = as.double(tcutc1), 
	hc1 = as.double(hc1),
        tcutd0A = as.double(tcutd0A), 
	hd0A = as.double(hd0A),
        tcutd0B = as.double(tcutd0B), 
	hd0B = as.double(hd0B),
        tcutd1A = as.double(tcutd1A), 
	hd1A = as.double(hd1A),
        tcutd1B = as.double(tcutd1B), 
	hd1B = as.double(hd1B),
        tcutx0A = as.double(tcutx0A), 
	hx0A = as.double(hx0A),
        tcutx0B = as.double(tcutx0B), 
	hx0B = as.double(hx0B),
        tcutx1A = as.double(tcutx1A), 
	hx1A = as.double(hx1A),
        tcutx1B = as.double(tcutx1B), 
	hx1B = as.double(hx1B),
        t0 = double(n), 
	t1 = double(n),
        tc0 = double(n), 
	tc1 = double(n),
        td0A = double(n), 
	td0B = double(n), 
	td1A = double(n), 
	td1B = double(n), 
	code = integer(2 + ndbl), 
	u = double(2 + ndbl), 
	TT = double(2 + ndbl), 
	delta = integer(2 + ndbl), 
	z = integer(2 + ndbl),
	time = double(2 + ndbl), 
	nrisk = integer(2 + 2*ndbl), 
	nevent = integer(2 + 2*ndbl), 
	pstat = double(nstat), 
	pvar = double(nstat), 
	pndths = integer(1), 
	avg.inffrac = double(nstat * nlook),
        avg.inffrac.ii = double(nstat * nlook),
	avg.bounds = as.double(c(rep(b.e, nstat), rep(b.f, nstat))),
        mufu = double(nstat*nlook),
	palphavec = double(2 * nstat * nlook), 
	pRejAcc = integer(2 * nstat * Nsim),
        kstop = integer(nstat * Nsim),
        duration = double(nstat * Nsim), 
	pstatall = double(nstat * Nsim),
        pvarall = double(nstat * Nsim),
	PACKAGE = "PwrGSD")
    is.detail <- !missing(detail)
    if (is.detail) {
        details <- ans
        ndths <- details$pndths
        names(details$ints) <- ints.nms
        names(details$dbls) <- dbls.nms
        details$pinffrac <- InfFrac
        details$pinffrac.ii <- InfFrac.ii
        details$alphatot <- c(Alpha.Efficacy,Alpha.Futility)
        details$pbounds <- Bounds
        details$palpha0vec <- matrix(ans$palphavec, nlook, 2*nstat)
        details$mufu <- matrix(ans$mufu, nlook, nstat)
	details$RejAcc <- matrix(details$pRejAcc, Nsim, 2*nstat)
        details$time <- details$time[1:ndths]
        details$nrisk <- details$nrisk[1:ndths]
        details$nevent <- details$nevent[1:ndths]
        details$time1 <- details$time1[1:ndths]
        details$nrisk1 <- details$nrisk1[1:ndths]
        details$nevent1 <- details$nevent1[1:ndths]
    }
    names(ints) <- ints.nms
    StatLast <- cbind(ans$pstat, ans$pvar)
    ndeaths <- ans$pndths
    RejAcc <- matrix(ans$pRejAcc, Nsim, 2*nstat)
    RejNull <- RejAcc[,1:nstat,drop=FALSE]
    AccNull <- RejAcc[,nstat + (1:nstat)]
    duration <- matrix(ans$duration, Nsim, nstat)
    StatAll <- matrix(ans$pstatall, Nsim, nstat)
    VarAll <- matrix(ans$pvarall, Nsim, nstat)
    InfFrac <- matrix(ans$avg.inffrac, nlook, nstat)
    InfFrac.ii <- matrix(ans$avg.inffrac.ii, nlook, nstat)
    Bounds <- matrix(ans$avg.bounds, nlook, (1 + do.futility)*nstat)
    kstop <- matrix(ans$kstop, Nsim, nstat)
    sum.kstop <- matrix(0, nlook, nstat)
    for(k in 1:nlook) sum.kstop[k,] <- c(rep(1,Nsim) %*% (kstop==k))
    Nstop <- apply(sum.kstop[nlook:1,,drop=FALSE], 2, FUN=cumsum)[nlook:1,]

    InfFrac <- InfFrac/Nstop
    InfFrac.ii <- InfFrac.ii/Nstop
    for(ef in 1:2)
      Bounds[,nstat*(ef-1) + (1:nstat)] <- Bounds[,nstat*(ef-1) + (1:nstat)]/Nstop

    dimnames(StatLast) <- list(stat.nms, c("stat", "var(stat)"))
    names(ndeaths) <- "ndeaths"
    dimnames(RejNull) <- list(1:Nsim, stat.nms)
    dimnames(duration) <- list(1:Nsim, stat.nms)
    dimnames(StatAll) <- list(1:Nsim, stat.nms)
    dimnames(VarAll) <- list(1:Nsim, stat.nms)
    dimnames(InfFrac) <- list(1:nlook, stat.nms)
    dimnames(InfFrac.ii) <- list(1:nlook, stat.nms)
    bstat.nms <- "Eff." %,% stat.nms
    if(do.futility==1) bstat.nms <- c(bstat.nms, "Fut." %,% stat.nms)
    dimnames(Bounds) <- list(1:nlook, bstat.nms)
    dimnames(sum.kstop) <- list(1:nlook, stat.nms)
    rslts <- cbind(StatAll, VarAll^0.5, duration, RejNull)[,outer(nstat*(0:3), 1:nstat, FUN="+")]
    dimnames(rslts) <- list(1:Nsim, c(outer(c("Stat", "SE", "Duration", "RejNull"), stat.nms, FUN=paste, sep="-")))
    dPower <- sapply(tlook, FUN=function(x,tt,rr)apply(rr*(tt==x),2,FUN=mean),tt=duration,rr=RejNull)
    dErrorII <- sapply(tlook, FUN=function(x,tt,rr)apply(rr*(tt==x),2,FUN=mean),tt=duration,rr=AccNull)
    se.dPower <- ((dPower * (1 - dPower))/Nsim)^0.5
    mu.d <- apply(duration, 2, FUN = mean)
    se.d <- (apply(duration, 2, FUN = var)/Nsim)^0.5
    detail <- list(ints=ints,pinffrac=InfFrac, pinffrac.ii=InfFrac.ii, alphatot=c(Alpha.Efficacy,Alpha.Futility),
                   pbounds=Bounds, palpha0vec=matrix(ans$palphavec, nlook, 2*nstat), mufu=matrix(ans$mufu, nlook, nstat))
    out <- list(dPower=dPower, se.dPower=se.dPower, dErrorII=dErrorII, Exp.Dur=mu.d, se.Exp.Dur=se.d,
                detail=detail, Nsim = Nsim, StatLast = StatLast, ndeaths = ndeaths, Elements = rslts, fail = ans$fail,
                stat.nms=stat.nms, RejNull=RejNull, duration=duration, InfFrac=InfFrac, InfFrac.ii=InfFrac.ii,
                Bounds=Bounds, Nstop=sum.kstop, call = .call.)
    if (is.detail)
        out$detail <- details
    class(out) <- c("SimPwrGSD", "PwrGSD")
    out
}

"print.SimPwrGSD" <-
function (x, ...)
{
    print(x$call)
    stat.nms <- x$stat.nms
    nstat <- x$detail$ints["nstat"]
    nlook <- x$detail$ints["nlook"]
    ans <- cbind(x$dPower %*% rep(1,nlook), (x$se.dPower^2) %*% rep(1, nlook), x$Exp.Dur, x$se.Exp.Dur)
    dimnames(ans) <- list(stat.nms, c("Power","se(Power)", "Exp.Dur", "se(Exp.Dur)"))
    print(ans)
    invisible(x)
}

"summary.SimPwrGSD" <- 
function(object, ...)
{
	out <- object
	stat.nms <- object$stat.nms
	nstat <- object$detail$ints["nstat"]
        ans <- cbind(object$dPower %*% rep(1,nlook), (object$se.dPower^2) %*% rep(1, nlook), object$Exp.Dur, se.Exp.Dur)
        dimnames(ans) <- list(stat.nms, c("Power","se(Power)", "Exp.Dur", "se(Exp.Dur)"))
	out$Tbl <- ans
        out
}

"SimPwrGSDcall" <-
function (Nsim, accru, accrat, tlook,
    Alpha.Efficacy=0, Alpha.Futility=0, sided = c("2", ">",
    "<"), tcut0 = NULL, h0 = NULL, s0 = NULL, tcut1 = NULL, 
    rhaz = NULL, h1 = NULL, s1 = NULL, tcutc0 = NULL, hc0 = NULL, 
    sc0 = NULL, tcutc1 = NULL, hc1 = NULL, sc1 = NULL, tcutd0A = NULL, 
    hd0A = NULL, sd0A = NULL, tcutd0B = NULL, hd0B = NULL, 
    sd0B = NULL, tcutd1A = NULL, hd1A = NULL, sd1A = NULL, 
    tcutd1B = NULL, hd1B = NULL, sd1B = NULL, tcutx0A = NULL, 
    hx0A = NULL, sx0A = NULL, tcutx0B = NULL, hx0B = NULL, 
    sx0B = NULL, tcutx1A = NULL, hx1A = NULL, sx1A = NULL,
    tcutx1B = NULL, hx1B = NULL, sx1B = NULL, 
    noncompliance = c("none", "crossover", "mixed", "user"), 
    gradual = FALSE, detail = FALSE, WtFun =c("FH","SFH","Ramp"), ppar = cbind(c(0, 0)), 
    Boundary.Efficacy = c("Lan-Demets", "Haybittle"), 
    Boundary.Futility = c("Lan-Demets", "Haybittle"),
    RR.Futility = NULL,                          
    Spending.Efficacy = c("Obrien-Fleming", "Pocock", "Power"),
    Spending.Futility = c("Obrien-Fleming", "Pocock", "Power"), 
    rho.Efficacy=0, rho.Futility=0, b.Haybittle = 3)  
{
	match.call()
}

"lookup" <-
function (xgrid, ygrid, x, y0 = 0)
{
        nx <- length(x)
        ngrid <- length(xgrid)
        ans <- .C("lookup",
        xgrid = as.double(xgrid),
        ygrid = as.double(ygrid),
        pngrid = as.integer(ngrid),
        x = as.double(x),
        pnx = as.integer(nx),
        py0 = as.double(y0),
        yatx = as.double(rep(0,nx)),
        index = as.integer(rep(0,nx)),
        PACKAGE = "PwrGSD")
        data.frame(x=ans$x,y=ans$yatx,index=ans$index)
}

"support" <- 
function(x)
{
	sort(unique(x))
}

"glue.2.string" <-
function(x, sep="")
{
	str <- x[1]
	n <- length(x)
	if (n>1)
	for(i in 2.:length(x))
		str <- str %,% sep %,% x[i]
	str
}

"html.table" <- 
function(x, file = "", append = FALSE, main = "", center = FALSE, html.tag = TRUE,
        table.attributes = "BORDER", column.attributes = c(
        row.header.column.attributes, rep(data.column.attributes, length = ncol(
        x))), cell.attributes = "", row.header.attributes = "ALIGN=RIGHT",
        row.header.column.attributes = "", data.column.attributes = "", ...)
{
        # formats data.frames, matrices, and vectors as html tables
        # also accepts a list of data.frames, matrices, and/or vectors
        # does not handle more sophisticated structures
        out <- c()
        # add <HTML> tag
        if(html.tag) out <- c(out, "<HTML>")
        # add <CENTER> tag
        if(center) out <- c(out, "<CENTER>")
        # if x is a list, iterate over elements
        # print list component names for each element
        # pad between elements with empty paragraph
        if(is.list(x) && !inherits(x, "data.frame")) {
                if(!is.null(main) && main != "")
                        out <- c(out, paste("<H2>", main, "</H2>"))
                n <- length(x)
                for(i in seq(length = n)) {
                        out <- c(out, html.table(x[[i]], main = names(x)[i]))
                        if(i < n)
                                out <- c(out, "<P> </P>")
                }
        }
        else {
                # if it's not a list, format it as a character matrix
                if(inherits(x, "data.frame")) {
                        dimnames.x <- dimnames(x)
                        dim.x <- dim(x)
                        x <- sapply(x, format, ...)
                        if(dim.x[1] <= 1) {
                                x <- matrix(x, nrow = dim.x[1], ncol = dim.x[
                                        2], byrow = TRUE)
                                dimnames(x) <- dimnames.x
                        }
                }
                else if(is.matrix(x)) {
                        dimnames.x <- dimnames(x)
                        dim.x <- dim(x)
                        x <- sapply(split(x, col(x)), format, ...)
                        if(dim.x[1] <= 1) {
                                x <- matrix(x, nrow = dim.x[1], ncol = dim.x[
                                        2], byrow = TRUE)
                                dimnames(x) <- dimnames.x
                        }
                }
                else {
                        dimnames.x <- list(names(x), NULL)
                        x <- as.matrix(format(x, ...))
                }
                out <- c(out, paste(collapse = " ", "<TABLE", paste(collapse =
                        " ", table.attributes), ">"))
                if(!is.null(main) && main != "")
                        out <- c(out, paste("<CAPTION> <H3>", main,
                                "</H3> </CAPTION>", collapse = " "))
                out <- c(out, paste("<TR>", paste(paste("<TH",
                        column.attributes, ">"), c(" ", dimnames.x[[2]]),
                        "</TH>", collapse = " "), "</TR>", collapse = " "))
                # data rows
                row.header.attributes <- rep(row.header.attributes, len = nrow(
                        x))
                if(is.matrix(cell.attributes)) {
                        if(!identical(dim(cell.attributes), dim(x)))
                                stop("If cell.attributes is a matrix it must " %,% 
                                     "have same dimensions as x")
                        for(i in seq(length = nrow(x))) {
                                out <- c(out, paste(paste("\n   <TR",
                                        row.header.attributes[i], ">"), "<TH>",
                                        dimnames.x[[1]][i], "</TH>", paste(
                                        sep = "", paste("\n      <TD",
                                        cell.attributes[i,  ], ">", sep = " "),
                                        x[i,  ], "</TD>", collapse = " "),
                                        "</TR>", collapse = " "))
                        }
                }
                else {
                        for(i in seq(length = nrow(x))) {
                                out <- c(out, paste(paste("\n   <TR",
                                        row.header.attributes[i], ">"), "<TH>",
                                        dimnames.x[[1]][i], "</TH>", paste(
                                        sep = "", paste("\n      <TD",
                                        cell.attributes, ">", sep = " "), x[
                                        i,  ], "</TD>", collapse = " "),
                                        "</TR>", collapse = " "))
                        }
                }
                out <- c(out, "</TABLE>")
        }
        if(center)
                out <- c(out, "</CENTER>")
        if(html.tag)
                out <- c(out, "</HTML>")
        # return vector of strings or write to file and return file name
        if(file == "") return(out) else {
                write(out, file = file, append = append)
                invisible(file)
        }
}

"encode" <-
function (x, basis)
{
        n <- length(basis)
        if(length(x)!=length(basis)) stop("lengths of 'x' and 'basis' must agree")
        sum(cumprod(c(1,basis[-n]))*x)
}

"decode" <- 
function (x, basis)
{
        n <- length(basis)
        if(length(x)!=1) stop("'x' must be of length 1")
        if(any(basis<=1)|length(basis)<=1) 
                stop("'basis' must be of length greater than 1, " %,%
        "and consist of elements greater than '1'")
        cpb <- cumprod(c(1,basis[-n]))
        ans <- NULL
        res <- x
        for (i in 1:n){
                new <- floor(res/cpb[n-i+1])
                res <- res - new*cpb[n-i+1]
                ans <- c(ans,new)
        }
        ans[n:1]
}

"Min" <- 
function (x, y)
{
    if (length(y) == 1)
        y <- rep(y, length(x))
    ans <- x
    ans[y < x] <- y[y < x]
    ans
}

"Max" <- 
function (x, y)
{
    if (length(y) == 1)
        y <- rep(y, length(x)) 
    ans <- x
    ans[y > x] <- y[y > x]
    ans
}

"RR2RCM" <- 
function (tlook, tcut.i, tcut.ii, h, rr, hOth, accru)
{
        h.i <- h
        h.ii <- h*rr
        mArmz(tlook, tcut.ii, h.ii, hOth, accru)/
        mArmz(tlook,  tcut.i,  h.i, hOth, accru)
}

"RCM2RR" <- 
function (tlook,tcut.i,h.i,hOth,accru,rcm)
{
        nlook <- length(tlook)
        fit <- list(objective=0)
        h.ii <- attr(fit$objective,"h") <- numeric(0)
        obj <-  function(theta,tlk,tcut.i,h.i,tcut.ii,h.ii,hOth,accru){
                nlk <- length(tlk)
                ans <- ((mArmz(tlk,tcut.ii,c(h.ii,exp(theta)),hOth,accru)/
                        mArmz(tlk,tcut.i,h.i,hOth,accru))[nlk] - rcm[nlk])^2
                attr(ans,"h") <- c(h.ii,exp(theta))
                ans
        }
        ans <- rep(0,nlook)
        for(i in 1:nlook){
                theta <- log(1e-12)
                tlk <- tlook[1:i]
                tcut.ii <- tcut.i[1]
                if(i>1) tcut.ii <- c(tcut.ii,tlook[1:(i-1)])
                h.ii <- attr(fit$objective,"h")
                fit <- optimize(f=obj,interval=c(log(1e-12),log(1e10)),tlk=tlk,
                        tcut.i=tcut.i,h.i=h.i,tcut.ii=tcut.ii,h.ii=h.ii,hOth=hOth,accru=accru)
                ans[i] <- exp(fit$min)
        }
        ans
}

"CRRtoRR" <- 
function(CRR, DT, h = NULL)
{
        RR <- CRR[1]
        m <- length(CRR)
        if(length(DT) != m)
                stop("lengths of 'CRR' and 'DT' must agree")
        if(missing(h))
                h <- rep(1, m)
        for(i in 2:m)
                RR <- c(RR, CRR[i] + sum(h[1:(i - 1)] * DT[1:(i - 1)] * (CRR[
                        i] - RR[1:(i - 1)]))/(h[i] * DT[i]))
        RR
}

"CY2TOShaz" <-
function(tcut, t.eor, m, verbose=FALSE)
{

  ti <- c(0, tcut)
  "Psi.1x" <-
  function(x, ti, h){
    L <- length(ti)
    H <- cumsum(DX(ti)*c(0,h))
    k.x <- max(which(x >= ti))
    ans <- 0
    if(x>ti[k.x])
      ans <- exp(-H[k.x])/h[k.x] * (1 - exp(-h[k.x] * (x - ti[k.x])))
    if(k.x >1)
      for(j in 1:(k.x-1)){
        ans <- ans + exp(-H[j])/h[j] * (1 - exp(-h[j] * (ti[j+1] -
        ti[j])))
    }
    ans
  }
  "Psi" <-
    function(x, ti, h){
      if(length(x)==1)
        ans <- Psi.1x(x, ti, h)
      if(length(x)>1)
        ans <- sapply(x, FUN=Psi.1x, ti=ti, h=h)
    ans
  }
  obj <-
  function(theta, ti, t.eor, m){
    x <- ti[-1]
    h <- exp(theta)
    CDF <- NULL
    if(any(x<=t.eor)){
      x.le.t <- x[x<=t.eor]
      CDF <- c(CDF, (x.le.t - Psi(x.le.t, ti, h))/t.eor)
    }
    if(any(x>t.eor)){
      x.gt.t <- x[x>t.eor]
      CDF <- c(CDF, 1 - (Psi(x.gt.t, ti, h) - Psi(x.gt.t - t.eor, ti,
               h))/t.eor)
    }
    n <- length(CDF)
    m.est <- 1 - (1-CDF)/(1-c(0,CDF[-n]))
    ans <- sum((m - m.est)^2)^0.5
    if(verbose) cat("\n objective: ",ans,"h: ",h,"\n")
    attr(ans,"tbl") <- cbind(x=x, CDF=CDF, m=m, m.est=m.est)
    ans
  }
  th0 <- log(-log(1-m)/DX(ti[-1]))
#  return(Psi.1x(tcut[1], ti, exp(th0)))
  fit <- optim(par=th0, fn=obj, method = "Nelder-Mead", ti=ti, t.eor=t.eor,
               m=m, control=list(maxit=10000))
  h <- exp(fit$par)
  obj. <- obj(fit$par, ti, t.eor, m)
  list(hazard=h, table=attr(obj., "tbl"))
}

"CDFOR2LRR" <-
function(tcut, tmax, h0, CDFOR)
{
  ti <- c(tcut, tmax)
  ind <- which(ti>0)
  ti <- ti[ind]
  m <- length(tcut)
  if(length(CDFOR)!=m)
    stop("lengths of 'tcut' and 'CDFR' must agree")
  DX <- function(x)c(x[1],diff(x))
  H0 <- cumsum(DX(ti) * h0)
  CDF.0 <- 1-exp(-H0)
  Odds.0 <- CDF.0/(1-CDF.0)
  Odds.1 <- CDFOR*Odds.0
  CDF.1 <- Odds.1/(1+Odds.1)
  beta <- log((DX(CDF.1)/(1-CDF.1))/(DX(CDF.0)/(1-CDF.0)))
  ans <- cbind(ti, beta)
  dimnames(ans) <- list(rep("", length(ti)), rep("", 2))
  ans
}

"mArmz" <- 
function (tlook, tcut, h, hOth, accru)
{
        fff <- function(tlk, tcut, h, hOth, accru){
                Min <- function (x, y)
                {
                    if (length(y) == 1)
                        y <- rep(y, length(x))
                    ans <- x
                    ans[y < x] <- y[y < x]
                    ans
                }
                Max <- function (x, y)
                {
                    if (length(y) == 1)
                        y <- rep(y, length(x))
                    ans <- x
                    ans[y > x] <- y[y > x]
                    ans
                }
                tER <- accru
                if(tlk - tER %in% tcut) tER <- tER-1e-10
                mu <- hOth
                ncut <- length(tcut)
                if(tcut[ncut]>tlk-tER) JtlkmtER <- max(min(which(tcut>tlk-tER))-1,0)
                else JtlkmtER <- ncut
                H <- cumsum(c(0,h[-ncut]*diff(tcut)))
                H.tlkmtER <- 0
                if(JtlkmtER>0) {
                        DT <- max(tlk-tER,0)-tcut[JtlkmtER]
                        H.tlkmtER <- H[JtlkmtER] + h[JtlkmtER]*DT
                }
                H. <- Max(H,H.tlkmtER)
                tcut. <- c(tcut,Inf)[-1]
                ans1 <- 0
                if(JtlkmtER>=1){
                        .b0. <- Min(tcut.,tlk-tER)
                        .a0. <- tcut
                        ans1 <- sum((h*exp(-(H+mu*.a0.))/(h+mu)*(1-exp(-(h+mu)*(.b0.-.a0.))))[1:JtlkmtER])
                }
                if(tcut[ncut]>tlk) Jtlk <- max(min(which(tcut>tlk))-1,0)
                else Jtlk <- ncut
                .c. <- Min(tcut.,tlk)
                .b. <- Max(tcut,tlk-tER)
                .a. <- tcut
                factor1 <- h*exp(-(H.+mu*.a.))/tER
                t1 <- (tlk - .b.)*exp(-(h+mu)*(.b.-.a.))/(h+mu)
                t2 <- (tlk - .c.)*exp(-(h+mu)*(.c.-.a.))/(h+mu)
                t3 <- exp(-(h+mu)*(.c.-.a.))/(h+mu)^2
                t4 <- exp(-(h+mu)*(.b.-.a.))/(h+mu)^2
                ans2 <- sum((factor1*(t1-t2+t3-t4))[JtlkmtER:Jtlk])
                ans1+ans2
        }
        sapply(tlook, FUN=fff, tcut=tcut, h=h, hOth=hOth, accru=accru)
}

"EX" <- function(x, h){
	nh <- length(h)
	dx <- diff(x)
	H <- cumsum(c(0,h[-nh]*dx))
	sum(exp(-H) * (1.0 - exp(-h*c(dx,1))*c(rep(1,nh-1),0))/h)
}

"thtilde.const" <- function(x,h,hA,hB,lA,lB)
{

	ans <- .C("htildeConst",
	x = as.double(x),
	nx = as.integer(nx),
	h = as.double(h),
	hA = as.double(hA),
	hB = as.double(hB),
	lA = as.double(lA),
	lB = as.double(lB),
	Stilde = as.double(rep(0,nx)),
	htlde = as.double(rep(0,nx)),
	PACKAGE = "PwrGSD")
}

"thtilde" <- 
function(x,xh,h,xhA,hA,xhB,hB,xlA,lA,xlB,lB,gradual=FALSE)
{
	.call. <- match.call()
	nx <- length(x)
	nh <- length(xh)
	nhA <- length(xhA)
	nhB <- length(xhB)
	nlA <- length(xlA)
	nlB <- length(xlB)

	bad <- FALSE
	bad.h  <- (length( h)!=nh)
	bad.hA <- (length(hA)!=nhA)
	bad.hB <- (length(hB)!=nhB)
	bad.lA <- (length(lA)!=nlA)
	bad.lB <- (length(lB)!=nlB)

	bads <- c(bad.h,bad.hA,bad.hB,bad.lA,bad.lB)
	bad <- any(bads)
	msg <- c("h","hA","hB","lA","lB")[which(bads)]

	glue.2.string <- function (x, sep = "")
	{
	    str <- x[1]
	    n <- length(x)
	    if (n > 1)
	        for (i in 2:length(x)) str <- str %,% sep %,% x[i]
	    str
	}
	if(bad) stop("Check that lengths of time and hazard agree in "%,% glue.2.string(msg[bads], sep=" "))

        nunq <- length(unique(sort(c(tcut0, tcut1, tcutc0, tcutc1, tcutd0A, tcutd0B, tcutd1A, tcutd1B, tcutx0A,
                                     tcutx0B, tcutx1A, tcutx1B))))    
        glegx <- glegx24
        glegw <- glegw24

	ngq <- length(glegx)
        
	ans <- .C("htilde",
	x = as.double(x),
	nx = as.integer(nx),
	gqxw = as.double(c(glegx,glegw)),
	ngq = as.integer(ngq),
	xh = as.double(xh),
	h = as.double(h),
	nh = as.integer(nh),
	xhA = as.double(xhA),
	hA = as.double(hA),
	nhA = as.integer(nhA),
	xhB = as.double(xhB),
	hB = as.double(hB),
	nhB = as.integer(nhB),
	xlA = as.double(xlA),
	lA = as.double(lA),
	nlA = as.integer(nlA),
	xlB = as.double(xlB),
	lB = as.double(lB),
	nlB = as.integer(nlB),
	gradual = as.integer(gradual),
	ftlde = as.double(rep(0,nx)),
	Stlde = as.double(rep(0,nx)),
	htlde = as.double(rep(0,nx)),
	PACKAGE = "PwrGSD")
	list(htilde = ans$htlde, ftilde = ans$ftlde, Stilde = ans$Stlde, call=.call.)
}

"trandfromh" <- function(N, tcut, h)
{
	.C("trandfromh",
	pn = as.integer(N),
	tcut = as.double(tcut),
	h = as.double(h),
	pncut = as.integer(length(tcut)),
	t = as.double(rep(0,N)),
	PACKAGE = "PwrGSD")
}

"trandhcdtl" <- function(N, tcut, h, tcutxA, hxA, tdA, tcutxB, hxB, tdB)
{
	.C("trandhcdtl",
	pn = as.integer(N),
	tcut = as.double(tcut),
	h = as.double(h),
	pncut = as.integer(length(tcut)),
	tcutxA = as.double(tcutxA),
	hxA = as.double(hxA),
	pncutxA = as.integer(length(tcutxA)),
	tdA = as.double(tdA),
	tcutxB = as.double(tcutxB),
	hxB = as.double(hxB),
	pncutXB = as.integer(length(tcutxB)),
	tdB = as.double(tdB),
	code = as.integer(rep(0,N)),
	t = as.double(rep(0,N)),
	PACKAGE = "PwrGSD")
}

"wtdlogrank" <- 
function(formula = formula(data), data = parent.frame(), WtFun = c("FH", "SFH", "Ramp"),
         param = c(0, 0), sided = c(2, 1), subset, na.action, w=FALSE) 
{
    if (missing(sided))
        sided <- 2
    m <- .call. <- match.call()
    m[[1]] <- as.name("model.frame")
    m$param <- m$sided <- m$WtFun <- m$w <- NULL
    m <- eval(m, parent.frame())
    mt <- attr(m, "terms")
    if (is.empty.model(mt)) 
        stop("No treatment indicator specified in model")
    Arm <- model.matrix(mt, m)[, -1]
    nlev <- length(unique(Arm))

    R <- model.extract(m, "response")
    if(nlev>2 || class(R)!="Surv")
      stop("Argument 'formula' is expected to have a response of class \"Surv\" and " %,%
           "a single predictor containing two levels")
    ind.too.small <- (R[,1]<1e-10)
    n.too.small <- sum(ind.too.small)
    if(n.too.small > 0) {
      warning("Deleting " %,% n.too.small %,% "observatioins that are less than 10^-10\n")
      R <- R[-which(ind.too.small),]
      Arm <- Arm[-which(ind.too.small)]
    }    
    TOS <- R[, 1]
    Event <- R[, 2]
    ntimes <- length(unique(TOS[Event!=0]))
    n <- length(Event)
    if(missing(WtFun)){
      param <- c(0,0)
      WtFun<-"FH"
    }
    wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1
    
    ans <- .C("WtdLogRank",
              TOS = as.double(TOS),
              Event = as.integer(Event), 
              Arm = as.integer(Arm),
              pn = as.integer(n),
              wttyp = as.integer(wttyp), 
              par = as.double(param),
              time = as.double(rep(0, ntimes)), 
              nrisk = as.integer(rep(0, 2*ntimes)),
              nevent = as.integer(rep(0, 2*ntimes)),
              wt = double(ntimes),              
              pntimes = as.integer(ntimes),              
              stat = double(1),
              var = double(1),
              UQt = double(ntimes),
              varQt = double(ntimes),
              var1t = double(ntimes),
              PACKAGE = "PwrGSD")

    out <- ans
    out$TOS <- out$Event <- out$Arm <- NULL
    out$time <- ans$time
    out$nrisk <- ans$nrisk[2*(1:ntimes) - 1] + ans$nrisk[2*(1:ntimes)]
    out$nevent <- ans$nevent[2*(1:ntimes)-1] + ans$nevent[2*(1:ntimes)]
    out$nrisk1 <- ans$nrisk[2*(1:ntimes)]
    out$nevent1 <- ans$nevent[2*(1:ntimes)]
    out$wt <- ans$wt
    out$pu0 <- sum(TOS[Arm==0])
    out$pu1 <- sum(TOS[Arm==1])
    sA <- support(Arm)
    out$n0 <- sum(Arm == sA[1])
    out$n1 <- sum(Arm == sA[2])
    out$n <- out$n0 + out$n1
    out$Z <- ans$stat/ans$var^0.5
    out$sided <- sided
    class(out) <- "survtest"
    out$call <- .call.
    out
}

"IntSurvDiff" <- 
function(formula = formula(data), data = parent.frame(), WtFun = c("FH", "SFH", "Ramp"),
         param = c(0, 0), sided = c(2, 1), subset, na.action, w=FALSE) 
{
    if (missing(sided))
        sided <- 2
    m <- .call. <- match.call()
    m[[1]] <- as.name("model.frame")
    m$param <- m$sided <- m$WtFun <- m$w <- NULL
    m <- eval(m, parent.frame())
    mt <- attr(m, "terms")
    if (is.empty.model(mt)) 
        stop("No treatment indicator specified in model")
    Arm <- model.matrix(mt, m)[, -1]
    nlev <- length(unique(Arm))

    R <- model.extract(m, "response")
    if(nlev>2 || class(R)!="Surv")
      stop("Argument 'formula' is expected to have a response of class \"Surv\" and " %,%
           "a single predictor containing two levels")
    ind.too.small <- (R[,1]<1e-10)
    n.too.small <- sum(ind.too.small)
    if(n.too.small > 0) {
      warning("Deleting " %,% n.too.small %,% "observatioins that are less than 10^-10\n")
      R <- R[-which(ind.too.small),]
      Arm <- Arm[-which(ind.too.small)]
    }    
    TOS <- R[, 1]
    Event <- R[, 2]
    ntimes <- length(unique(TOS[Event!=0]))
    n <- length(Event)
    if(missing(WtFun)){
      param <- c(0,0)
      WtFun<-"FH"
    }
    wttyp <- charmatch(WtFun, c("FH", "SFH", "Ramp"))-1

    ans <- .C("IntSurvDiff",
              TOS = as.double(TOS),
              Event = as.integer(Event), 
              Arm = as.integer(Arm),
              pn = as.integer(n),
              wttyp = as.integer(wttyp), 
              par = as.double(param),
              time = as.double(rep(0, ntimes)), 
              nrisk = as.integer(rep(0, 2*ntimes)),
              nevent = as.integer(rep(0, 2*ntimes)),
              pntimes = as.integer(ntimes),              
              stat = as.double(0),
              var = as.double(0),
              wt = as.double(rep(0, ntimes)),
              PACKAGE = "PwrGSD")
    out <- ans
    out$TOS <- out$Event <- out$Arm <- NULL
    out$time <- ans$time
    out$nrisk <- ans$nrisk[2*(1:ntimes) - 1] + ans$nrisk[2*(1:ntimes)]
    out$nevent <- ans$nevent[2*(1:ntimes)-1] + ans$nevent[2*(1:ntimes)]
    out$nrisk1 <- ans$nrisk[2*(1:ntimes)]
    out$nevent1 <- ans$nevent[2*(1:ntimes)]
    out$wt <- ans$wt
    out$pu0 <- sum(TOS[Arm==0])
    out$pu1 <- sum(TOS[Arm==1])
    sA <- support(Arm)
    out$n0 <- sum(Arm == sA[1])
    out$n1 <- sum(Arm == sA[2])
    out$n <- out$n0 + out$n1
    out$Z <- ans$stat/ans$var^0.5
    out$sided <- sided
    class(out) <- "survtest"
    out$call <- .call.
    out
  }

"print.survtest" <- function(x, ...) 
{
    sided <- x$sided
    Z <- x$stat/x$var^0.5
    pval <- (2^(sided == 2)) * (1 - pnorm(abs(Z)))
    cat("call:\n")
    print(x$call)
    z.nm <- as.character(x$call$formula[[3]])
    n0 <- x$n0
    n1 <- x$n1
    pu0 <- x$pu0
    pu1 <- x$pu1
    obs1 <- sum(x$nevent1)
    obs0 <- sum(x$nevent) - obs1
    exp0 <- sum(x$nevent * (1 - x$nrisk1/x$nrisk))
    exp1 <- sum(x$nevent * x$nrisk1/x$nrisk)
    tbl <- data.frame(n = c(n0, n1), pu = c(pu0, pu1), Obs = c(obs0, obs1), Exp = c(exp0, 
        exp1))
    dimnames(tbl) <- list(z.nm %,% c("=0:", "=1:"), names(tbl))
    print(tbl)
    ans <- c(Z, pval)
    names(ans) <- c("z", "p-val")
    print(ans)
    invisible(x)
}

"summary.survtest" <- function(object, ...) 
{
    sided <- object$sided
    Z <- object$stat/object$var^0.5
    pval <- (2^(sided == 2)) * (1 - pnorm(abs(Z)))
    cat("call:\n")
    print(object$call)
    z.nm <- as.character(object$call$formula[[3]])
    n0 <- object$n0
    n1 <- object$n1
    pu0 <- object$pu0
    pu1 <- object$pu1
    obs1 <- sum(object$nevent1)
    obs0 <- sum(object$nevent) - obs1
    exp0 <- sum(object$nevent * (1 - object$nrisk1/object$nrisk))
    exp1 <- sum(object$nevent * object$nrisk1/object$nrisk)
    tbl <- data.frame(n = c(n0, n1), pu = c(pu0, pu1), Obs = c(obs0, obs1), Exp = c(exp0, 
        exp1))
    dimnames(tbl) <- list(z.nm %,% c("=0:", "=1:"), names(tbl))
    object$tbl <- tbl
    ans <- c(Z, pval)
    names(ans) <- c("z", "p-val")
    object$ans
    object
}

"mysurvfit" <-
function (formula = formula(data), data = parent.frame(), subset,
    na.action = na.fail)
{
    m <- .call. <- match.call()
    m[[1]] <- as.name("model.frame")
    m <- eval(m, parent.frame())
    mt <- attr(m, "terms")
    R <- model.extract(m, "response")
    ind.too.small <- (R[,1]<1e-10)
    n.too.small <- sum(ind.too.small)
    TOS <- R[, 1]
    Event <- R[, 2]
    ntimes <- length(unique(TOS[Event!=0]))
    n <- length(Event)
    Arm <- 0 * Event
    Arm.levs <- 0
    nb <- 1
    int.only <- (length(attr(mt, "factors"))==0)
    if(!int.only){
        Arm <- model.matrix(mt, m)[,-1]
        Arm.f <- as.factor(as.character(Arm))
        Arm.levs <- levels(Arm.f)
        nb <- length(Arm.levs)
        Arm <- as.numeric(Arm.f) - 1
    }
    ans <- .C("mysurvfit",
              TOS = as.double(TOS),
              Event = as.integer(Event),
              Arm = as.integer(Arm),
              pn = as.integer(n),
              time = as.double(rep(0, ntimes)),
              nrisk = as.integer(rep(0, ntimes * nb)),
              nevent = as.integer(rep(0, ntimes * nb)),
              pntimes = as.integer(ntimes),
              pnblocks = as.integer(nb),
              PACKAGE = "PwrGSD")

    tbl <- as.data.frame(cbind(ans$time, t(matrix(ans$nrisk, nb, ntimes)),
                               t(matrix(ans$nevent, nb, ntimes))))
    nms <- c("time", "nrisk" %,% Arm.levs, "nevent" %,% Arm.levs)
    names(tbl) <- nms
    out <- list(call=.call., Table=tbl)
    class(out) <- "blkdcp"
    out
}

"print.blkdcp" <-
function(x, ...)
{
    cat("Call:\n",deparse(x$call),"\n\n")
    cat("Table:\n")
    print(x$Table)
    invisible(x$Table)
}

"plot.blkdcp" <-
function(x, event.name="Incidence", colors=NULL, ...)
{
    M <- as.matrix(x$Table)
    nb <- (ncol(M)-1)/2
    if(missing(colors)) colors <- 1:nb
    ti <- M[,1]
    ri.ind <- 1+(1:nb)
    ev.ind <- 1+nb+(1:nb)
    A <- apply(M[,ev.ind]/M[,ri.ind], 2, FUN=cumsum)
    rng.ti <- range(ti)
    rng.A <- range(A)
    plot(rng.ti, rng.A, type="n", xlab="Time on Study", ylab="Cummulative " %,% event.name)
    for(k in 1:nb) lines(ti, A[,k], type="s", col=colors[k])
    invisible(x)
}

"CondPower" <-
function(Z, frac, drift, drift.end, err.I, sided=1)
{
  mu.c.H0 <- Z * frac^0.5
  mu.c.HA <- mu.c.H0 + drift.end - drift
  rt <- drift.end >0
  z.c <- rt * qnorm(1 - err.I/sided) + (1 - rt) * qnorm(err.I/sided)
  Pr.cond.typeIIerr <- rt * pnorm((z.c - mu.c.HA)/(1 - frac)^0.5) +
    (1 - rt) * (1-pnorm((z.c - mu.c.HA)/(1 - frac)^0.5))
  Pr.cond.typeIerr <- rt * (1 - pnorm((z.c - mu.c.H0)/(1 -
        frac)^0.5)) + (1 - rt) * pnorm((z.c - mu.c.H0)/(1 - frac)^0.5)
  cbind(Pr.cond.typeIerr=Pr.cond.typeIerr, Pr.cond.typeIIerr=Pr.cond.typeIIerr)
}

"stpplt" <-
function (x, y, ...) 
{
    d.y <- dim(y)
    n <- d.y[1]
    d <- d.y[2]
    cls <- c("aquamarine", "magenta")
    polygon(c(0, x, 1), c(0, y[, 1], 0), border = NA, col = cls[1])
    for (j in 2:d) polygon(c(x[n:1], x), c(y[n:1, j - 1], y[, 
        j]), border = NA, col = cls[2 - (j%%2)])
    polygon(c(x, 0), c(y[, d], 1), border = NA, col = cls[2 - 
        ((d + 1)%%2)])
    for(j in 1:d)
      lines(x,y[,j], ...)
}

"plot.cpd.PwrGSD" <-
function(x, formula, subset, na.action, ...)
{
  ow <- options("warn")
  options(warn = -1)
  .call. <- match.call(expand=FALSE)
  dots <- .call.$...
  nms.dots <- names(dots)
  given.values <- rows <- columns <- show.given <- col <-
             pch <- bar.bg <- fac <- xlab <- ylab <- subscripts <- axlabels <-
             number <- overlap <- xlim <- ylim <- NULL  
  if("given.values" %in% nms.dots) given.values <- dots[["given.values"]]
  if("rows" %in% nms.dots) rows <- dots[["rows"]]
  if("columns" %in% nms.dots) columns <- dots[["columns"]]
  if("show.given" %in% nms.dots) show.given <- dots[["show.given"]]
  if("col" %in% nms.dots) col <- dots[["col"]]
  if("pch" %in% nms.dots) pch <- dots[["pch"]]
  if("bar.bg" %in% nms.dots) bar.bg <- dots[["bar.bg"]]
  if("fac" %in% nms.dots) fac <- dots[["fac"]]
  if("xlab" %in% nms.dots) xlab <- dots[["xlab"]]
  if("ylab" %in% nms.dots) ylab <- dots[["ylab"]]
  if("subscripts" %in% nms.dots) subscripts <- dots[["subscripts"]]
  if("axlabels" %in% nms.dots) axlabels <- dots[["axlabels"]]
  if("number" %in% nms.dots) number <- dots[["number"]]
  if("overlap" %in% nms.dots) overlap <- dots[["overlap"]]
  if("xlim" %in% nms.dots) xlim <- dots[["xlim"]]
  if("ylim" %in% nms.dots) ylim <- dots[["ylim"]]
  .call.[[1]] <- as.name("costopplot")
  .call.$data <- as.call(expression(as.data.frame.cpd.PwrGSD))
  .call.$data$x <- .call.$x
  .call.$x <- NULL
  .call.$given.values <- given.values
  .call.$rows  <- rows
  .call.$columns <- columns
  .call.$show.given <- show.given
  .call.$col <- col
  .call.$pch <- pch
  .call.$bar.bg <- bar.bg
  .call.$fac <- fac
  .call.$xlab <- xlab
  .call.$ylab <- ylab
  .call.$subscripts <- subscripts
  .call.$axlabels <- axlabels
  .call.$number <- number
  .call.$overlap <- overlap
  .call.$xlim <- xlim 
  .call.$ylim <- ylim
  eval(.call.)
  options(ow)
}

"costopplot" <- 
function (formula, data, given.values, rows, columns, show.given = TRUE, 
    col = par("fg"), pch = par("pch"), bar.bg = c(num = gray(0.8), 
        fac = gray(0.95)), xlab = c(x.name, paste("Given :", 
        a.name)), ylab = c("Type II Error Prob (Aqua) & Power (Magenta)", paste("Given :", b.name)),
    subscripts = FALSE, axlabels = function(f) abbreviate(levels(f)), 
    number = 6, overlap = c(0,0), xlim, ylim, subset, na.action, ...) 
{
    formula.new <- (y ~ x)
    formula.new[[3]] <- formula[[2]]
    formula.new[[2]] <- (I(cbind(dE, dP)) ~ x)[[2]]
    formula <- formula.new
    is.subset <- !missing(subset)
    is.na.action <- !missing(na.action)
    .call. <- match.call()
    .call.$formula <- formula
    deparen <- function(expr) {
        while (is.language(expr) && !is.name(expr) && deparse(expr[[1]]) == 
            "(") expr <- expr[[2]]
        expr
    }
    bad.formula <- function() stop("invalid conditioning formula")
    bad.lengths <- function() stop("incompatible variable lengths")
    formula <- deparen(formula)
    if (!inherits(formula, "formula")) 
        bad.formula()
    y <- deparen(formula[[2]])
    rhs <- deparen(formula[[3]])
    if (deparse(rhs[[1]]) != "|") 
        bad.formula()
    x <- deparen(rhs[[2]])
    rhs <- deparen(rhs[[3]])
    if (is.language(rhs) && !is.name(rhs) && (deparse(rhs[[1]]) == 
        "*" || deparse(rhs[[1]]) == "+")) {
        have.b <- TRUE
        a <- deparen(rhs[[2]])
        b <- deparen(rhs[[3]])
    }
    else {
        have.b <- FALSE
        a <- rhs
    }
    if (missing(data)) 
        data <- parent.frame()
    form <- (yv ~ xv + av)
    form[[2]] <- y
    form[[3]][[2]] <- x
    form[[3]][[3]] <- a
    if (have.b) {
        form <- (yv ~ xv + av * bv)
        form[[2]] <- y
        form[[3]][[2]] <- x
        form[[3]][[3]][[2]] <- a
        form[[3]][[3]][[3]] <- b
    }
    mdl <- as.call(expression(model.frame))
    mdl$formula <- form
    mdl$data <- data
    if (is.subset) 
        mdl$subset <- .call.$subset
    if (is.na.action) 
        mdl$na.action <- .call.$na.action
    mdl <- eval(mdl, parent.frame())
    y.name <- deparse(y)
    y <- model.extract(mdl, "response")
    nlook <- attr(data, "detail")["nlook"]
    d.y <- dim(y)
    ntot.orig <- d.y[1]
    n.conds <- ntot.orig/nlook
    y.new <- matrix(0, n.conds, 2 * nlook)
    for (j in 1:n.conds) {
        y.new[j, 2 * (1:nlook) - 1] <- y[nlook * (j - 1) + (1:nlook), 
            1]
        y.new[j, 2 * (1:nlook)] <- y[nlook * (j - 1) + (1:nlook), 
            2]
        y.new[j, ] <- cumsum(y.new[j, ])
    }
    y <- y.new
    mdl <- mdl[nlook * (0:((ntot.orig - 1)/nlook)) + 1, ]
    x.name <- deparse(x)
    x <- mdl[, x.name]
    nobs <- length(x)
    d.y <- c(nobs, length(y)/nobs)
    d <- d.y[2]
    if (length(y)%%nobs) 
        bad.lengths()
    a.name <- deparse(a)
    a <- as.factor(mdl[, a.name])
    if (length(a) != nobs) 
        bad.lengths()
    if (is.character(a)) 
        a <- as.factor(a)
    a.is.fac <- is.factor(a)
    if (have.b) {
        b.name <- deparse(b)
        b <- as.factor(mdl[, b.name])
        if (length(b) != nobs) 
            bad.lengths()
        if (is.character(b)) 
            b <- as.factor(b)
        b.is.fac <- is.factor(b)
        missingrows <- which(is.na(x) | is.na(y) | is.na(a) | 
            is.na(b))
    }
    else {
        missingrows <- which(is.na(x) | is.na(y) | is.na(a))
        b <- NULL
        b.name <- ""
    }
    number <- as.integer(number)
    if (length(number) == 0 || any(number < 1)) 
        stop("number must be integer >= 1")
    if (any(overlap >= 1)) 
        stop("overlap must be < 1 (and typically >= 0).")
    bad.givens <- function() stop("invalid given.values")
    if (missing(given.values)) {
        a.intervals <- if (a.is.fac) {
            i <- seq(along = a.levels <- levels(a))
            a <- as.numeric(a)
            cbind(i - 0.5, i + 0.5)
        }
        else co.intervals(a, number = number[1], overlap = overlap[1])
        b.intervals <- if (have.b) {
            if (b.is.fac) {
                i <- seq(along = b.levels <- levels(b))
                b <- as.numeric(b)
                cbind(i - 0.5, i + 0.5)
            }
            else {
                if (length(number) == 1) 
                  number <- rep.int(number, 2)
                if (length(overlap) == 1) 
                  overlap <- rep.int(overlap, 2)
                co.intervals(b, number = number[2], overlap = overlap[2])
            }
        }
    }
    else {
        if (!is.list(given.values)) 
            given.values <- list(given.values)
        if (length(given.values) != (if (have.b) 
            2
        else 1)) 
            bad.givens()
        a.intervals <- given.values[[1]]
        if (a.is.fac) {
            a.levels <- levels(a)
            if (is.character(a.intervals)) 
                a.intervals <- match(a.intervals, a.levels)
            a.intervals <- cbind(a.intervals - 0.5, a.intervals + 
                0.5)
            a <- as.numeric(a)
        }
        else if (is.numeric(a)) {
            if (!is.numeric(a.intervals)) 
                bad.givens()
            if (!is.matrix(a.intervals) || ncol(a.intervals) != 
                2) 
                a.intervals <- cbind(a.intervals - 0.5, a.intervals + 
                  0.5)
        }
        if (have.b) {
            b.intervals <- given.values[[2]]
            if (b.is.fac) {
                b.levels <- levels(b)
                if (is.character(b.intervals)) 
                  b.intervals <- match(b.intervals, b.levels)
                b.intervals <- cbind(b.intervals - 0.5, b.intervals + 
                  0.5)
                b <- as.numeric(b)
            }
            else if (is.numeric(b)) {
                if (!is.numeric(b.intervals)) 
                  bad.givens()
                if (!is.matrix(b.intervals) || ncol(b.intervals) != 
                  2) 
                  b.intervals <- cbind(b.intervals - 0.5, b.intervals + 
                    0.5)
            }
        }
    }
    if (any(is.na(a.intervals)) || (have.b && any(is.na(b.intervals)))) 
        bad.givens()
    if (have.b) {
        rows <- nrow(b.intervals)
        columns <- nrow(a.intervals)
        nplots <- rows * columns
        if (length(show.given) < 2) 
            show.given <- rep.int(show.given, 2)
    }
    else {
        nplots <- nrow(a.intervals)
        if (missing(rows)) {
            if (missing(columns)) {
                rows <- ceiling(round(sqrt(nplots)))
                columns <- ceiling(nplots/rows)
            }
            else rows <- ceiling(nplots/columns)
        }
        else if (missing(columns)) 
            columns <- ceiling(nplots/rows)
        if (rows * columns < nplots) 
            stop("rows * columns too small")
    }
    total.columns <- columns
    total.rows <- rows
    f.col <- f.row <- 1
    if (show.given[1]) {
        total.rows <- rows + 1
        f.row <- rows/total.rows
    }
    if (have.b && show.given[2]) {
        total.columns <- columns + 1
        f.col <- columns/total.columns
    }
    mar <- if (have.b) 
        rep.int(0, 4)
    else c(0.5, 0, 0.5, 0)
    oma <- c(5, 6, 5, 4)
    if (have.b) {
        oma[2] <- 5
        if (!b.is.fac) 
            oma[4] <- 5
    }
    if (a.is.fac && show.given[1]) 
        oma[3] <- oma[3] - 1
    opar <- par(mfrow = c(total.rows, total.columns), oma = oma, 
        mar = mar, xaxs = "r", yaxs = "r", new = FALSE)
    on.exit(par(opar))
    plot.new()
    if (missing(xlim)) 
        xlim <- range(as.numeric(x), finite = TRUE)
    if (missing(ylim)) 
        ylim <- range(as.numeric(y), finite = TRUE)
    pch <- rep(pch, length.out = nobs)
    col <- rep(col, length.out = nobs)
    do.panel <- function(index, subscripts = FALSE, id) {
        Paxis <- function(side, x) {
            if (nlevels(x)) {
                lab <- axlabels(x)
                axis(side, labels = lab, at = seq(lab), xpd = NA)
            }
            else axis(side, xpd = NA)
        }
        istart <- (total.rows - rows) + 1
        i <- total.rows - ((index - 1)%/%columns)
        j <- (index - 1)%%columns + 1
        par(mfg = c(i, j, total.rows, total.columns))
        plot.new()
        plot.window(xlim, ylim)
        if (any(is.na(id))) 
            id[is.na(id)] <- FALSE
        if (any(id)) {
            n.id <- sum(id)
            grid(lty = "solid")
            if (subscripts) 
                stpplt(x[id], y[id, ], subscripts = id, ...)
            else stpplt(x[id], y[id, ], ...)
        }
        if ((i == total.rows) && (j%%2 == 0)) 
            Paxis(1, x)
        else if ((i == istart || index + columns > nplots) && 
            (j%%2 == 1)) 
            Paxis(3, x)
        if ((j == 1) && ((total.rows - i)%%2 == 0)) 
            Paxis(2, y)
        else if ((j == columns || index == nplots) && ((total.rows - 
            i)%%2 == 1)) 
            Paxis(4, y)
        box()
    }
    if (have.b) {
        count <- 1
        for (i in 1:rows) {
            for (j in 1:columns) {
                id <- ((a.intervals[j, 1] <= a) & (a <= a.intervals[j, 
                  2]) & (b.intervals[i, 1] <= b) & (b <= b.intervals[i, 
                  2]))
                do.panel(count, subscripts, id)
                count <- count + 1
            }
        }
    }
    else {
        for (i in 1:nplots) {
            id <- ((a.intervals[i, 1] <= a) & (a <= a.intervals[i, 
                2]))
            do.panel(i, subscripts, id)
        }
    }
    mtext(xlab[1], side = 1, at = 0.5 * f.col, outer = TRUE, 
        line = 3.5, xpd = NA)
    mtext(ylab[1], side = 2, at = 0.5 * f.row, outer = TRUE, 
        line = 3.5, xpd = NA)
    if (length(xlab) == 1) 
        xlab <- c(xlab, paste("Given :", a.name))
    if (show.given[1]) {
        par(fig = c(0, f.col, f.row, 1), mar = mar + c(3 + (!a.is.fac), 
            0, 0, 0), new = TRUE)
        plot.new()
        nint <- nrow(a.intervals)
        a.range <- range(a.intervals, finite = TRUE)
        plot.window(a.range + c(0.03, -0.03) * diff(a.range), 
            0.5 + c(0, nint))
        rect(a.intervals[, 1], 1:nint - 0.3, a.intervals[, 2], 
            1:nint + 0.3, col = bar.bg[if (a.is.fac) 
                "fac"
            else "num"])
        if (a.is.fac) {
            text(apply(a.intervals, 1, mean), 1:nint, a.levels)
        }
        else {
            axis(3, xpd = NA)
            axis(1, labels = FALSE)
        }
        box()
        mtext(xlab[2], 3, line = 3 - a.is.fac, at = mean(par("usr")[1:2]), 
            xpd = NA)
    }
    else {
        mtext(xlab[2], 3, line = 3.25, outer = TRUE, at = 0.5 * 
            f.col, xpd = NA)
    }
    if (have.b) {
        if (length(ylab) == 1) 
            ylab <- c(ylab, paste("Given :", b.name))
        if (show.given[2]) {
            par(fig = c(f.col, 1, 0, f.row), mar = mar + c(0, 
                3 + (!b.is.fac), 0, 0), new = TRUE)
            plot.new()
            nint <- nrow(b.intervals)
            b.range <- range(b.intervals, finite = TRUE)
            plot.window(0.5 + c(0, nint), b.range + c(0.03, -0.03) * 
                diff(b.range))
            rect(1:nint - 0.3, b.intervals[, 1], 1:nint + 0.3, 
                b.intervals[, 2], col = bar.bg[if (b.is.fac) 
                  "fac"
                else "num"])
            if (b.is.fac) {
                text(1:nint, apply(b.intervals, 1, mean), b.levels, 
                  srt = 90)
            }
            else {
                axis(4, xpd = NA)
                axis(2, labels = FALSE)
            }
            box()
            mtext(ylab[2], 4, line = 3 - b.is.fac, at = mean(par("usr")[3:4]), 
                xpd = NA)
        }
        else {
            mtext(ylab[2], 4, line = 3.25, at = 0.5 * f.row, 
                outer = TRUE, xpd = NA)
        }
    }
    if (length(missingrows) > 0) {
        cat("\nMissing rows:", missingrows, "\n")
        invisible(missingrows)
    }
}

"agghaz" <- function(t.agg, time, nrisk, nevent)  
{
  n.t <- length(time)
  n.blocks <- ncol(nevent)
  result <-.C("agghaz",
              tagg=as.double(t.agg),
              time=as.double(time),
              nrisk=as.integer(nrisk),
              nevent=as.integer(nevent),
              pndth=as.integer(n.t),
              pnb=as.integer(n.blocks),
              timea=as.double(rep(0,n.t)),
              nriska=as.integer(rep(0,n.blocks*n.t)),
              neventa=as.integer(rep(0,n.blocks*n.t)),
              pnagg=as.integer(0),
              PACKAGE="PwrGSD")
  n.agg <- result$pnagg
  time.a <- result$timea[1:n.agg]
  nrisk.a <- matrix(result$nriska[1:(n.blocks*n.agg)], n.agg, n.blocks)
  nevent.a <- matrix(result$neventa[1:(n.blocks*n.agg)], n.agg, n.blocks)
  out <- as.data.frame(cbind(time.a, nrisk.a, nevent.a))
  names(out) <- c("time", "nrisk" %,% (1:n.blocks), "nevent" %,% (1:n.blocks))
  out
}

"mystack" <- 
function (object, fu.vars, create.idvar = FALSE) 
{
    d.object <- dim(object)
    n <- d.object[1]
    p <- d.object[2]
    for(k in 1:p)
      if(is.character(object[,k]))
        object[,k] <- as.factor(object[,k])

    if (create.idvar) 
        object$id <- 1:n
    out <- NULL
    nms <- names(object)
    n.fu.vars <- length(fu.vars)
    n.fus <- length(grep(fu.vars[1], nms))
    fu.var.pos <- matrix(0, n.fus, n.fu.vars)
    for (k in 1:n.fu.vars) fu.var.pos[, k] <- grep(fu.vars[k], 
        nms)
    base.vars <- object[, -c(fu.var.pos)]
    nbv <- dim(base.vars)[2]
    bvfactors <- NULL
    bvfactorlevs <- list()
    l <- 1
    for (k in 1:nbv) {
        isf <- is.factor(base.vars[, k])
        if (isf) {
            bvfactors <- c(bvfactors, k)
            bvfactorlevs[[l]] <- levels(base.vars[, k])
            base.vars[,k] <- as.numeric(base.vars[,k])
            l <- l + 1
        }
    }
    base.vars <- as.matrix(base.vars)
    nbvf <- l - 1
    fuvars <- object[, fu.var.pos]
    fufactors <- NULL
    fufactorlevs <- list()
    l <- 1
    for (k in 1:n.fu.vars) {
        isf <- is.factor(fuvars[, n.fus * (k - 1) + 1])
        if (isf) {
            fufactors <- c(fufactors, k)
            fufactorlevs[[l]] <- levels(fuvars[, n.fus * (k - 
                1) + 1])
            fuvars <- as.numeric(fuvars[,k])
        }
    }
    fuvars <- as.matrix(fuvars)
    nfuf <- l - 1
    out <- 
    .C("mystack", pn = as.integer(n), pnfus = as.integer(n.fus), 
        pnfuvars = as.integer(n.fu.vars), pnbasevars = as.integer(nbv), 
        basevars = as.double(base.vars), fuvars = as.double(fuvars), 
        out = double(n.fus * n * (nbv + 1 + n.fu.vars)), NAOK = TRUE, 
        PACKAGE = "PwrGSD")$out
    out <- as.data.frame(matrix(out, n.fus * n, nbv + 1 + n.fu.vars))
    names(out) <- c(nms[-fu.var.pos], "fu", fu.vars)
    if (nbvf > 0) 
        for (l in 1:nbvf) out[, bvfactors[l]] <- bvfactorlevs[[l]][out[,bvfactors[l]]]
    if (nfuf > 0) 
        for (l in 1:nfuf) out[, nbv + 1 + fufactors[l]] <- fufactorlevs[[l]][out[,nbv + 1 + fufactors[l]]]
    out
}

.onAttach <- function(libname, pkgname)
{
    require(survival)
}
