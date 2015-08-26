


###### pbsytest.panelmodel: ##################
##
## fixed: Degrees of freedom in the joint test (test="j") of Baltagi/Li (1991). Should be chisquare(2) instead of chisquare(1),
##        see Baltagi/Li (1991), p. 279 and again in Baltagi/Li (1995), p. 136
##
## added: unbalanced versions of tests as in Sosa-Escudero/Bera (2008), Tests for unbalanced error-components models under local misspecification,
##                                                                       The Stata Journal (2008), Vol. 8, Number 1, pp. 68-78
##
##        Implementation follows the formulas for unbalanced panels, which reduce in for balanced data to the formulas for balanced panels.
##        Notation in code largly follows Baltagi's; m in Sosa-Escudero/Bera is total number of observations
## 
## changed: chisq test (two-sided) is standard; for one-sided alternative (normalized version) use argument normal = TRUE
##
## added: Check if we are working with a pooled model in the panelmodel interface, as this is necessary for the test.
##        (The formula interface has such a check already, but check was missing for panelmodel interface. Thus, for the old implementation,
##        if a different model type is passed, the user gets no warning and the statistics are way off.)
##
#
# Reference for the balanced tests="ar"|"re": Bera/Sosa-Escudero/Yoon (2001), Tests for the error component model in the presence of local misspecifcation,
#                                                                             Journal of Econometrics 101 (2001), pp. 1-23
#
#                      for original test="j": Baltagi/Li (1991), A joint test for serial correlation and random individual effects,
#                                                                Statistics & Probability Letters 11 (1991), pp. 277-280
#
# Concise treatment in Baltagi (2005), Econometric Analysis of Panel Data, 3rd edition, pp. 96-97
#                   or Baltagi (2013), Econometric Analysis of Panel Data, 5th edition, pp. 108:

pbsytest_mod_unbalanced.panelmodel <- function(x, test = c("ar","re","j"), normal = FALSE, ...){
  
  test <- match.arg(test)
  if (plm:::describe(x, "model") != "pooling") stop("pbsytest only relevant for pooling models") # added
  if (normal & test %in% c("ar", "j")) stop("normalized tests only for test=\"re\"")
  
  poolres <- resid(x)
  data <- model.frame(x)
  ## extract indices
  index <- attr(data, "index")
  tindex <- index[[2]]
  index <- index[[1]]
  
  ## till here. 
  ## ordering here if needed.
  
  ## this needs ordering of obs. on time, regardless 
  ## whether before that on groups or after
  
  ## and numerosity check
  
  ## order by group, then time
  oo <- order(index,tindex)
  ind <- index[oo]
  tind <- tindex[oo]
  poolres <- poolres[oo]
  pdim <- pdim(x)
  n <- max(pdim$Tint$n) ## det. number of groups
  T_i <- pdim$Tint$Ti
  N_t <- pdim$Tint$nt
  t <- max(T_i) ## det. max. group numerosity
  N_obs <- pdim$nT$N ## det. total number of obs.
  
  ## calc. A and B:
  S1 <- sum( tapply(poolres,ind,sum)^2 )
  S2 <- sum( poolres^2 )
  
  A <- S1/S2 -1
  
  unind <- unique(ind)
  uu <- rep(NA,length(unind))
  uu1 <- rep(NA,length(unind))
  for(i in 1:length(unind)) {
    u.t <- poolres[ind==unind[i]]
    u.t.1 <- u.t[-length(u.t)]
    u.t <- u.t[-1]
    uu[i] <- crossprod(u.t)
    uu1[i] <- crossprod(u.t,u.t.1)
  }
  
  B <- sum(uu1)/sum(uu)
  
  a <- sum(T_i^2)
  
  switch(test,
         ar={
             stat <- (B - (((N_obs - n)/(a - N_obs))* A))^2 * (a - N_obs)*N_obs^2 / ((N_obs - n)*(a - 3*N_obs + 2*n))  # RS*_lambda from Sosa-Escudero/Bera (2008), p. 73 in slightly different notation
             df <- c(df=1)
             names(stat) <- "chisq"
             pstat <- pchisq(stat, df=df, lower.tail=FALSE)
             tname <- "Bera, Sosa-Escudero and Yoon locally robust test"
             myH0_alt <- "AR(1) errors sub random effects"
         },
         
         re={
             stat <- ((N_obs^2) * (2*B-A)^2) / (2*(a -3*N_obs + 2* n)) # RS*_u from Sosa-Escudero/Bera (2008), p. 73 in slightly different notation
             names(stat) <- "chisq"
             df <- c(df=1)
             pstat <- pchisq(stat, df=df, lower.tail=FALSE)
             tname <- "Bera, Sosa-Escudero and Yoon locally robust test"
             myH0_alt <- "random effects sub AR(1) errors"
         },
         
         j={
            stat <- N_obs^2 * ((A^2 - 4*A*B + 4*B^2) / (2*(a - 3*N_obs + 2*n)) + (B^2/(N_obs - n))) # RS_lamba_u in Sosa-Escudero/Bera (2008), p. 74 
            df <- c(df=2)
            names(stat) <- "chisq"
            pstat <- pchisq(stat, df=df, lower.tail=FALSE) # fixed: df=df
            tname <- "Baltagi and Li AR-RE joint test"
            myH0_alt <- "AR(1) errors or random effects"
         }
  )
  
  dname <- paste(deparse(substitute(formula)))
  balanced.type <- ifelse(pdim$balanced, "balanced", "unbalanced")
  tname <- paste(tname, "-", balanced.type, "panel", collapse = " ")
  
  if (normal & test == "re") {
    # for re test: return normalized version with one-sided alternative
    myH0_alt <- paste(myH0_alt, "- one-sided", collapse = " ")
    stat <- sqrt(stat) # convert to normalized version
    names(stat) <- "normal"
    pstat <- pnorm(stat, lower.tail=FALSE) # alternative is one-sided
    parameters = NULL
    
    RVAL <- list(statistic = stat,
                 method = tname,
                 alternative = myH0_alt,
                 p.value = pstat,
                 data.name = dname)
  } else {
    if (normal) stop("normalized test only for test=\"re\"")
    # return chisq statistic
    myH0_alt <- paste(myH0_alt, "- two-sided", collapse = " ")
    
    RVAL <- list(statistic = stat,
                 parameter = df,
                 method = tname,
                 alternative = myH0_alt,
                 p.value = pstat,
                 data.name = dname)
  }
  
  class(RVAL) <- "htest"
  return(RVAL)
} ###### END pbsytest_mod_unbalanced.panelmodel: ##################