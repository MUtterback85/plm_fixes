# Own quick'n'dirty fixes (?)/enhancements to plm version 1.4-0
# no warranty
# License: GPL
#
# Version of this file 0.4-1
#
#
# Please find this file also at https://github.com/helix123/plm_fixes
# Updated documentation for plm (mainly new text book editions) is at https://github.com/helix123/plm/tree/master/man

#
# Instructions:
# load this file after package plm is loaded
# The following functions are then masked:
#   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
#   - summary.plm() and print.summary.plm() (used by summary())
#   - plmtest() [additionally Baltagi/Li (1990) is directly available by baltagi_li_re()]
#   - pbgtest() allows to pass on type="F" to lmtest::bgtest(), thus offering the small sample test (F test)
#   - pbltest(): added panelmodel interface
#   - pbsytest(): fixed degrees of freedom error when test="j" (Baltagi/Li (1991)); added warning if wrong input model
#   - pbltest_lm5(): new function added to compute test statistic LM5 from Baltagi/Li (1995)
#   - pbltest.panelmodel(): panelmodel interface added for convenience
#   - lag.pseries() can handle negative lags (leading values); lead.pseries() is added for convenience


#
## F statistic when a user specified variance-covariance matrix is supplied (for robust inference)
#  [see http://stackoverflow.com/questions/31163353/computing-f-statistic-with-user-supplied-covariance-matrix]
#
#  Note: Formula interface not fixed

## Durbin-Watson test respecting panel structure
#  [see http://stackoverflow.com/questions/31894055/pdwtest-from-plm-with-wrong-p-value-for-pooled-ols-durbin-watson-test-for-autoc]
#
# Note: There are some points noted down below which I am not sure how to handle
#       (marked with:  Statisticians: Can someone please look into this?)
#

## Breusch-Pagan test for random effects for unbalanced panels as in Baltagi/Li (1990)
# [see http://stackoverflow.com/questions/31988449/implementation-of-breusch-pagan-test-for-random-effects-in-plm-with-unbalanced-p]
#
# References: 
# Baltagi/Li (1990), A lagrange multiplier test for the error components model with incomplete panels,
#                    Econometric Reviews, 9:1, 103-107, DOI: 10.1080/07474939008800180
#        
# see how STATA does it and this produces the same statistic: http://www.stata.com/manuals14/xtxtregpostestimation.pdf

# Breusch-Godfrey test for autocorrelation: pbgtest() now passes on (if supplied) type="F" to lmtst:bgtest() for small sample test



################## pdwtest.panelmodel() adapted from pseries.R [Durbin-Watson test] ##################
pdwtest.panelmodel <- function(x,...) {
  ## residual serial correlation test based on the residuals of the demeaned
  ## model and the regular dwtest() in {lmtest}
  ## reference Baltagi (page 98) for FE application, Wooldridge page 288 for
  ## the general idea.
  
  
  model <- plm:::describe(x, "model")
  effect <- plm:::describe(x, "effect")
  theta <- x$ercomp$theta
  
  
  ## check package availability and load if necessary
  lm.ok <- require("lmtest")
  if(!lm.ok) stop("package lmtest is needed but not available")
  ## ARtest is the lmtest, exception made for the method attribute
  dots <- match.call(expand.dots=FALSE)[["..."]]
  if (is.null(dots$order.by)) order.by <- NULL else order.by <- dots$order.by
  if (is.null(dots$alternative)) alternative <- "greater" else alternative <- dots$alternative
  if (is.null(dots$iterations)) iterations <- 15 else iterations <- dots$iterations
  if (is.null(dots$exact)) exact <- NULL else exact <- dots$exact
  if (is.null(dots$tol)) tol <- 1e-10 else tol <- dots$tol
  
  
  #### Uncomment, if panel structure should not be taken into account #####
  #### Then, for pooled OLS, we can use the normal Durbin Watson test (as implemented in lmtest)
  #### But: see below
  
  #     if(model == "pooling") {
  
  #     # pass lm($formula) to lm::dwtest() rathar than the model.matrix as it is in plm_v1.4-0
  #     # => this avoids a doubled intercept which in turn causes the p-value to be way off [I guess due to wrong # degrees of freedom]
  #     ARtest <- lmtest::dwtest(lm(x$formula), order.by = order.by, alternative = alternative,
  #                    iterations = iterations, exact = exact, tol = tol)
  # 
  #     ARtest$method <- "Durbin-Watson test for serial correlation in panel models (pooled OLS - panel structure not taken into account)"
  #     ARtest$alternative <- "serial correlation in idiosyncratic errors"
  #     ARtest$data.name <- paste(deparse(x$formula))
  #  }
  
  
  #### Take panel structure into account, also for pooled OLS model ###
  #### => Results deviates from just applying lmtest's dwtest (regular Durbin Watson test) to the residuals
  #### (gretl also respects the panel structure for Durbin-Watson test for pooled OLS models)
  #
  # Implemented: approach of Bhargava et al. (1982) for (balanced) panels
  # reference: Bhargava/Franzini/Narendranathan, Serial Correlation and the Fixed Effects Model, Review of Economic Studies (1982), XLIX, pp. 533-549
  #
  # Not implemented:
  # Statisticians: Can someone please look into this?
  # There seems to be an modified version of Bhargava et al. (1982) which accounts for unbalanced and unequally spaced data
  # and an additional test Baltagi/Wu_LBI in this reference
  # Baltagi, B. H., and P. X. Wu. 1999. Unequally spaced panel data regressions with AR(1) disturbances. Econometric Theory 15, pp 814-823.
  # STATA has modified.Bhargava et al. (1982) and Baltagi/Wu LBI: http://www.stata.com/manuals14/xtxtregar.pdf
  # 
  
  # For FE and RE, we  need to take the panel structure into account
  # we can do the calculation based on the FE model as the assumptions are also true for the RE model
  # 
  # Statisticians: Can someone please look into this?
  #
  # references: FIXME/TODO [Verbeek 4th edition, papers?]
  #             FE: 

    # if RE model => calculate associated FE model based on formula
    # if passed model is already FE or pooled OLS => go ahead
    mod_to_test <- if (model == "random") fe_mod <- plm(x$formula, data = model.frame(x), model="within") else x
    
    if (pdim(mod_to_test)$balanced != TRUE) warning("Applying Bhargava et al. (1982) Durbin-Watson test for balanced panels to an unbalanced panel.")
  
    # residuals are now class pseries, so diff.pseries is used and the differences are computed within observational units
    # (not across as it would be the case if base::diff() is used and as it is done for lm-objects)
    # NAs are introduced by the differencing as one observation is lost per observational unit
    dw <- sum(plm:::diff.pseries(residuals(mod_to_test))^2, na.rm = T) / sum(residuals(mod_to_test)^2)
    
    # p-value computation seems to be difficult, so no p-value is calculated
    # maybe someone with more statistical knowledge can look into this one
    
    # constuct htest object
    names(dw) <- "DW"
    ARtest <- list(statistic = dw,
                   method = "Bhargava et al. (1982): Durbin-Watson test for serial correlation in balanced panel models (pooled OLS, within (fixed) and random effects) \n
                   For random effects models, the Durbin-Watson statistic is calculated to the corresponding fixed effect model. \n
                   No p-value computed as distubution of statistic under panel assumptions is difficult to calculate.",
                   alternative = NULL, p.value = NULL, data.name = paste(deparse(x$formula)))
    
    class(ARtest) <- "htest"

  return(ARtest)
}



################## pdwtest.formula() adapted from pseries.R ##################
# Formula interface for pdwtest() => not implemented, just prints a message
# 
pdwtest.formula <- function(x, data, ...) {
  ## formula method for pdwtest;
  ## defaults to a RE model
  
  print("Formula interface to pdwtest() not supported. Pass the estimated model to pdwtest().")
}


################## summary.plm() and print.summary.plm() adapted from plm.methods.R ##################

# fixed: use .vcov (if passed) for F test calculation [use package car::linearHypothesis() for calculation]


summary.plm <- function(object, .vcov = NULL, ...){
  object$fstatistic <- plm:::Ftest(object, test = "F")
  model <- plm:::describe(object, "model")
  effect <- plm:::describe(object, "effect")
  object$r.squared <- c(rsq  = r.squared(object),
                        adjrsq = r.squared(object, dfcor = TRUE))
  # construct the table of coefficients
  if (!is.null(.vcov)){
    std.err <- sqrt(diag(.vcov))
    
    # overwrite already calculated F statistic as we have a user supplied .vcov which needs to be respected
    # use car::linearHypothesis() for calculation
    
    car.ok <- require("car")
    if(!car.ok) stop("package car is needed but not available")

    # Need to check if there is an intercept in the model.
    # Intercept should not be passed to car::linearHypothesis(), because if so, the wrong # degrees of freedom is calculated
    return_car_lH <- car::linearHypothesis(object, names(coef(object))[if ("(Intercept)" %in% names(coef(object))) -1 else TRUE], test="F", vcov. = .vcov)

    # extract values from returned object from car::linearHypothesis
    object$fstatistic$statistic <- c("F" = return_car_lH[3][2, ]) # f statistic
    object$fstatistic$p.value <- c("F" = return_car_lH[4][2, ]) # p-value for F statistic
    object$fstatistic$parameter <- c("df1" = return_car_lH[2][2, ], "df2" = return_car_lH[1][2, ])  # Dfs

  }
  else{
    std.err <- sqrt(diag(vcov(object)))
  }
  b <- coefficients(object)
  z <- b / std.err
  p <- 2 * pt(abs(z), df = object$df.residual, lower.tail = FALSE)
  object$coefficients <- cbind("Estimate"   = b,
                               "Std. Error" = std.err,
                               "t-value"    = z,
                               "Pr(>|t|)"   = p)
  class(object) <- c("summary.plm", "plm", "panelmodel")
  object
}



print.summary.plm <- function(x,digits= max(3, getOption("digits") - 2),
                              width=getOption("width"), subset = NULL, ...){
  formula <- formula(x)
  has.instruments <- (length(formula)[2] == 2)
  effect <- plm:::describe(x, "effect")
  model <- plm:::describe(x, "model")
  cat(paste(plm:::effect.plm.list[effect]," ",sep=""))
  cat(paste(plm:::model.plm.list[model]," Model",sep=""))
  
  if (model=="random"){
    ercomp <- plm:::describe(x, "random.method")
    cat(paste(" \n   (",
              plm:::random.method.list[ercomp],
              "'s transformation)\n",
              sep=""))
  }
  else{
    cat("\n")
  }
  if (has.instruments){
    ivar <- plm:::describe(x, "inst.method")
    cat(paste("Instrumental variable estimation\n   (",
              inst.method.list[ivar],
              "'s transformation)\n",
              sep=""))
  }
  cat("\nCall:\n")
  print(x$call)
  cat("\n")
  pdim <- pdim(x)
  print(pdim)
  if (model == "random"){
    cat("\nEffects:\n")
    print(x$ercomp)
  }
  cat("\nResiduals :\n")
  save.digits <- unlist(options(digits = digits))
  on.exit(options(digits = save.digits))
  print(plm:::sumres(x))
  
  cat("\nCoefficients :\n")
  if (is.null(subset)) printCoefmat(coef(x), digits = digits)
  else printCoefmat(coef(x)[subset, , drop = FALSE], digits = digits)
  cat("\n")
  cat(paste("Total Sum of Squares:    ",signif(plm:::tss.plm(x),digits),"\n",sep=""))
  cat(paste("Residual Sum of Squares: ",signif(deviance(x),digits),"\n",sep=""))
  cat(paste("R-Squared      : ", signif(x$r.squared[1], digits),"\n"))
  cat("      Adj. R-Squared : ", signif(x$r.squared[2], digits),"\n")
  fstat <- x$fstatistic
  if (names(fstat$statistic) == "F"){
    cat(paste("F-statistic: ",signif(fstat$statistic),
              " on ",fstat$parameter["df1"]," and ",fstat$parameter["df2"],
              " DF, p-value: ",format.pval(fstat$p.value,digits=digits),"\n",sep=""))
  }
  else{
    cat(paste("Chisq: ",signif(fstat$statistic),
              " on ",fstat$parameter,
              " DF, p-value: ",format.pval(fstat$p.value,digits=digits),"\n",sep=""))
    
  }
  invisible(x)
}


################  Breusch-Pagan test for random effects model for unbalanced data ################  
################  run instead of original plmtest() - but see below for a wrapper version of plmtest()
# Breusch-Pagan test for random effects model for unbalanced data 
# Baltagi/Li (1990), A lagrange multiplier test for the error components model with incomplete panels,
#                    Econometric Reviews, 9:1, 103-107, DOI: 10.1080/07474939008800180
#
# for individual effects
baltagi_li_re <- function (x,
                        effect = c("individual", "time", "twoways"),
                        type = c("honda", "bp", "ghm","kw"),
                        ...) {
  
  ### this is from original plmtest()
  effect <- match.arg(effect)
  type <- match.arg(type)
  if (plm:::describe(x, "model") != "pooling") x <- update(x, model = "pooling")
  pdim <- pdim(x)
  n <- pdim$nT$n
  T <- pdim$nT$T
  balanced <- pdim$balanced
  index <- attr(model.frame(x), "index")
  id <- index[[1]]
  time <- index[[2]]
  res <- resid(x)
  ### END this is from original plmtst() ###
  
  # Implementation for unbalanced panels
  if (effect != "individual") {stop("Not implemented. Only individual effects are implemented.")}  

  T_i <- as.numeric(table(id))
  T_mean <- mean(T_i)
  
  A1 <-  1 - sum(tapply(res,id,sum)^2)/sum(res^2)
  term <- sum(T_i^2) - (n * T_mean)
  bp_stat_baltagi_li_re <- (((n * T_mean)^2) / 2) * ( A1^2 / term)
  
  parameter <- 1
  names(parameter) <- "df"

  res <- list(statistic = c(chisq  = bp_stat_baltagi_li_re),
              p.value = pchisq(bp_stat_baltagi_li_re, df = 1, lower.tail = FALSE),
              parameter = parameter,
              method = "Lagrange Multiplier Test - individual effects - Breusch-Pagan Test for unbalanced Panels as in Baltagi/Li (1990)",
              data.name = plm:::data.name(x),
              alternative = "significant effects")
  class(res) <- "htest"
  return(res)
} # END baltagi_li_re




############## plmtest() ############################################
## modified to handle unbalanced panels as in Baltagi/li (1990) #####

plmtest <- function(x,...){
  UseMethod("plmtest")
}

plmtest.plm <- function(x,
                        effect = c("individual", "time", "twoways"),
                        type = c("honda", "bp", "ghm","kw"),
                        ...){
  
  effect <- match.arg(effect)
  type <- match.arg(type)
  if (plm:::describe(x, "model") != "pooling") x <- update(x, model = "pooling")
  pdim <- pdim(x)
  n <- pdim$nT$n
  T <- pdim$nT$T
  balanced <- pdim$balanced
  index <- attr(model.frame(x), "index")
  id <- index[[1]]
  time <- index[[2]]
  res <- resid(x)
  
  
  if (balanced == F) { # for unbalanced panels, we need Baltagi/Li (1990)
    
    # Implementation for unbalanced panels
    return(baltagi_li_re(x))
    
  } else { ### balanced panel => use original implementation of plm [I copied it in here]
    
    if (effect != "twoways"){
      if (!type %in% c("honda", "bp"))
        stop("type must be one of honda or bp for a one way model")
      if(effect == "individual"){ condvar <- id ; card.cond <- n ; card.other <- T}
      else{condvar <- time ; card.cond <- T ; card.other <- n}
      stat <-  sqrt(card.other*card.cond/(2*(card.other-1)))*
        (crossprod(tapply(res,condvar,mean))*card.other^2/sum(res^2)-1)
      stat <- switch(type,
                     honda = c(normal = stat),
                     bp    = c(chisq  = stat^2))
      parameter <- switch(type,
                          honda = NULL,
                          bp = 1)
      pval <- switch(type,
                     honda = pnorm(abs(stat), lower.tail = FALSE)*2,
                     bp    = pchisq(stat, df = 1, lower.tail = FALSE))
    }
    else{
      stat1 <-  sqrt(n*T/(2*(T-1)))*(crossprod(tapply(res,id,mean))*T^2/sum(res^2)-1)
      stat2 <-  sqrt(n*T/(2*(n-1)))*(crossprod(tapply(res,time,mean))*n^2/sum(res^2)-1)
      stat <- switch(type,
                     ghm   = c(chisq = max(0,stat1)^2+max(0,stat2)^2),
                     bp    = c(chisq = stat1^2+stat2^2),
                     honda = c(normal = (stat1+stat2)/sqrt(2)),
                     kw    = c(normal = sqrt((T-1)/(n+T-2))*stat1+sqrt((n-1)/(n+T-2))*stat2))
      parameter <- 2
      pval <- switch(type,
                     ghm   = pchisq(stat,df=2,lower.tail=FALSE),
                     honda = pnorm(abs(stat),lower.tail=FALSE)*2,
                     bp    = pchisq(stat,df=2,lower.tail=FALSE),
                     kw    = pnorm(abs(stat),lower.tail=FALSE)*2)
    }
    
    method.type <- switch(type,
                          honda  = "Honda",
                          bp     = "Breusch-Pagan",
                          ghm    = "Gourieroux, Holly and Monfort",
                          kw     = "King and Wu")
    method.effect <- switch(effect,
                            id      = "individual effects",
                            time    = "time effects",
                            twoways = "two-ways effects")
    method <- paste("Lagrange Multiplier Test - ",method.effect,
                    " (",method.type,")\n",sep="")
    
    if(type == "honda"){
      res <- list(statistic = stat,
                  p.value   = pval,
                  method    = method,
                  data.name = plm:::data.name(x))
    }
    else{
      names(parameter) <- "df"
      res <- list(statistic = stat,
                  p.value   = pval,
                  method    = method,
                  parameter = parameter,
                  data.name = plm:::data.name(x))
    }
    res$alternative <- "significant effects"
    class(res) <- "htest"
    res
    
    
  }
} # END plmtest()


############### Breusch-Godfrey test ##################################
### fixed pbgtest(), copied over from https://r-forge.r-project.org/scm/viewvc.php/pkg/R/pserial.R?view=markup&root=plm&pathrev=127
# pbgtest() suffered from the same problem as pdwtest() [intercept passed twice to lmtest::bgtest()] - incorporated that from r-forge
#
# additional fix: match arguments, so that type="F" (and order.by=) and is passed on to lmtest::bgtest(), thus enabling the small sample
#                 variant of the test offered by lmtest::bgtest()


pbgtest.panelmodel<-function(x, order = NULL, ...) {
  ## residual serial correlation test based on the residuals of the demeaned
  ## model (see Wooldridge p.288) and the regular bgtest() in {lmtest}
  
  ## structure:
  ## 1: take demeaned data from 'plm' object
  ## 2: est. auxiliary model by OLS on demeaned data
  ## 3: apply bgtest() to auxiliary model and return the result
  
  model <- plm:::describe(x, "model")
  effect <- plm:::describe(x, "effect")
  theta <- x$ercomp$theta
  
  ## retrieve demeaned data
  demX <- model.matrix(x, model = model, effect = effect, theta=theta)
  demy <- pmodel.response(model.frame(x), model = model, effect = effect, theta=theta)
  
  ## ...and group numerosities
  Ti <- pdim(x)$Tint$Ti
  ## set lag order to minimum group numerosity if not specified by user
  ## (check whether this is sensible)
  
  if(is.null(order)) order <- min(Ti)
  ## bg test on the demeaned model:
  
  ## check package availability and load if necessary
  #lm.ok <- require("lmtest")
  #if(!lm.ok) stop("package lmtest is needed but not available")
  
  ## bgtest is the bgtest, exception made for the method attribute
  dots <- match.call(expand.dots=FALSE)[["..."]]      # fixed: added expand.dots=FALSE
  if (!is.null(dots$type)) type <- dots$type else type <- "Chisq"
  if (!is.null(dots$order.by)) order.by <- dots$order.by else order.by <- NULL
  
  auxformula <- demy~demX-1 #if(model == "within") demy~demX-1 else demy~demX
  lm.mod <- lm(auxformula)
  bgtest <- bgtest(lm.mod, order = order, type = type, order.by = order.by)
  bgtest$method <- "Breusch-Godfrey/Wooldridge test for serial correlation in panel models"
  bgtest$alternative <- "serial correlation in idiosyncratic errors"
  bgtest$data.name <- paste(deparse(x$call$formula))
  names(bgtest$statistic) <- if(length(bgtest$parameter)==1) "chisq" else "F"
  return(bgtest)
}




######### pbltest(): added panelmodel interface ########
# Baltagi and Li Serial Dependence Test For Random Effects Models
#
# Add panelmodel interface
# This is a wrapper: using the formula interface is a bit cumbersome as it makes many
# assumptions on how the arguments should be formed. Calling the panelmodel interface is easier.

pbltest.formula <- pbltest

pbltest <- function (x, alternative = c("twosided", "onesided"), index = NULL, ...) {
  UseMethod("pbltest")
}

pbltest.panelmodel <- function(x, ...) {

  # only continue if random effects model
  if (plm:::describe(x, "model") != "random") stop("Test only for random effects models.")
  
  # call pbltest.formula in the right way
  pbltest.formula(formula(x$formula), data=cbind(index(x), x$model), index=names(index(x)), ...)
}


############# Baltagi/Li (1995), LM5
# LM5 (p. 138-139: An LM test for first-order serial correlation in a fixed effects model):
# LM5 is LM3 (from p. 136) only with the residuals of the FE model instead of the OLS model (in matrix B)
# Thus, matrix B calculated in pblsytest is already in the right form for LM5 but need to
# subsitute the residuals for the FE residuals.

# p. 136: "LM1, is exactly the same as the joint test statistic derived by Baltagi and Li (1991) for
# AR(l) residual disturbances and random individual effects." => matrix B calculated in pblsytest(),
# adapt for use here
#
# see equivalently Baltagi (2005), Econometric Analysis of Panel Data, 3rd edition, pp. 97-98 or
#                  Baltagi (2013), Econometric Analysis of Panel Data, 5th edition, pp. 108-109
#
# Baltagi (2005), p. 98 [=Baltagi (2013), p. 109]:
# "Since the [w]ithin transformation wipes out the individual effects whether fixed or random,
#  one can also use [LM5] to test for serial correlation in the random effects models."

pbltest_lm5 <- function(x, ...) {
  
  if (plm:::describe(x, "model") == "pooling") stop("Test only for within effects (=fixed effects) or random effect models.")
  
##### following code adapted from pblsytest() #####
  poolfe <- resid(x)
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
  poolfe <- poolfe[oo]
  ## det. number of groups and df
  n <- length(unique(index))
  k <- ncol(model.matrix(x))
  ## det. max. group numerosity
  t <- max(pdim(x)$Tint$Ti)
  ## det. total number of obs. (robust vs. unbalanced panels)
  nT <- length(ind)
    
  ## calc. B (with FE residuals:
  unind <- unique(ind)
  uu <- rep(NA,length(unind))
  uu1 <- rep(NA,length(unind))
  for(i in 1:length(unind)) {
    u.t <- poolfe[ind==unind[i]]
    u.t.1 <- u.t[-length(u.t)]
    u.t <- u.t[-1]
    uu[i] <- crossprod(u.t)
    uu1[i] <- crossprod(u.t,u.t.1)
  }
    
  B <- sum(uu1)/sum(uu)
##### END code adapted from pblsytest() #####

  # p. 138. onesided. Under H0, this is asymptotically distributed (for large T) as N(0, 1).
  LM5_statsitic <- sqrt((n * t^2) / (t-1)) * B
  names(LM5_statsitic) <- "LM5"
  
  # not needed: p. 138: This LM statistic is asymptotically distributed (for large T) as chisquare(1) [under the null hypothesis]
  # LM5_statistic_squared <- LM5_statsitic^2

  dname <- paste(deparse(substitute(x$formula)))
  RVAL <- list(statistic = LM5_statsitic,
               parameter = NULL,
               method = "LM test for first-order serial correlation in a fixed effects model \n
                         LM5 in Baltagi/Li (1995), p. 138-139",
               alternative = "AR(1) errors (rho > 0) sub fixed effects",
               p.value = pnorm(LM5_statsitic,lower.tail=FALSE),
               data.name = dname)
  class(RVAL) <- "htest"
  return(RVAL)
} ## END: pbltest_lm5


###### pbsystest.panelmodel: ##################
##
## fixed: Degrees of freedom in the joint test (test="j") of Baltagi/Li (1991). Should be chisquare(2) instead of chisquare(1),
##        see Baltagi/Li (1991), p. 279 and again in Baltagi/Li (1995), p. 136
##
## added: Check to verify we are working with a pooled model in the panelmodel interface, as this is necessary for the test.
##        (The formula interface has such check already, but it is missing for panelmodel interface. If a different model
##         type is passed, the user gets no warning and the statistics are way off.)

pbsytest.panelmodel <- function(x, test=c("ar","re","j"), ...){
  
  if (plm:::describe(x, "model") != "pooling") stop("pbsytest only relevant for pooling models") # added
  
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
  ## det. number of groups and df
  n <- length(unique(index))
  k <- ncol(model.matrix(x))
  ## det. max. group numerosity
  t <- max(pdim(x)$Tint$Ti)
  ## det. total number of obs. (robust vs. unbalanced panels)
  nT <- length(ind)
  
  ## calc. A and B:
  S1 <- sum( tapply(poolres,ind,sum)^2 )
  S2 <- sum( poolres^2 )
  
  A <- S1/S2-1
  
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
  
  switch(match.arg(test),
           ar ={LM <- (n * t^2 * (B - (A/t))^2) / ((t-1)*(1-(2/t)))
           df <- c(df=1)
           names(LM) <- "chisq"
           pLM <- pchisq(LM,df=1,lower.tail=FALSE)
           tname <- "Bera, Sosa-Escudero and Yoon locally robust test"
           myH0 <- "AR(1) errors sub random effects"
         },
           re={LM <- (A - 2*B) * sqrt( (n * t) / (2*(t-1)*(1-(2/t))) )
           names(LM) <- "z"
           df <- NULL
           pLM <- pnorm(LM,lower.tail=FALSE)
           tname <- "Bera, Sosa-Escudero and Yoon locally robust test"
           myH0 <- "random effects sub AR(1) errors"
         },              
           j={LM <- (n * t^2) / (2*(t-1)*(t-2)) * (A^2 - 4*A*B + 2*t*B^2) 
           df <- c(df=2)
           names(LM) <- "chisq"
           pLM <- pchisq(LM,df=df,lower.tail=FALSE) # fixed df=df
           tname <- "Baltagi and Li AR-RE joint test"
           myH0 <- "AR(1) errors or random effects"
         }
  )
  
  dname <- paste(deparse(substitute(formula)))
  RVAL <- list(statistic = LM,
               parameter = df,
               method = tname,
               alternative = myH0,
               p.value = pLM,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)
  
} ###### END pbsytest.panelmodel: ##################

############# lag.pseries: able to create negative lags (=leading values) (use k < 0)
## for convenience: method lead.pseries is also added
#
# also as diff against plm_v1.4-0 on github: https://github.com/cran/plm/compare/cd00e7ce878f8ffa0e0bc53ab8692778cb8aaecc...1ef2620b2851ae6d7374d52bfed8141e11eaff76
# 
lag.pseries <- function(x, k = 1, ...) {
  nx <- names(x)
  index <- attr(x, "index")
  id <- index[[1]]
  time <- index[[2]]
  
  # catch the case when an index of pdata.frame shall be lagged
  if (is.factor(x)) if (all(as.character(x) == as.character(id)) | all(as.character(x)==as.character(time))) stop("Lagged vector cannot be index.")
  
  alag <- function(x, ak){
    if (round(ak) != ak) stop("Lagging value 'k' must be whole-numbered (positive, negative or zero)")
    if (ak > 0){
      # delete first ak observations for each unit
      isNAtime <- c(rep(T,ak), diff(as.numeric(time), lag = ak)) != ak
      isNAid <- c(rep(T,ak), diff(as.numeric(id), lag = ak)) != 0
      isNA <- as.logical(isNAtime + isNAid)
      if (is.factor(x)) levs <- levels(x)
      result <- c(rep(NA, ak), x[1:(length(x)-ak)])
      result[isNA] <- NA
      if (is.factor(x)) result <- factor(result, labels = levs)
      structure(result,
                names = nx,
                class = class(x),
                index = index)
    } else if (ak < 0) { # => compute leading values
      
      # delete last ak observations for each unit
      isNAtime <- c(as.numeric(time) - c(tail(as.numeric(time), length(time) + ak), rep(T, -ak))) != ak
      isNAid   <- c(as.numeric(id) - c(tail(as.numeric(id), length(id) + ak) , rep(T, -ak))) != 0
      isNA <- as.logical(isNAtime + isNAid)
      result <- c(x[(1-ak):(length(x))], rep(NA, -ak))
      result[isNA] <- NA
      if (is.factor(x)) levs <- levels(x)
      if (is.factor(x)) result <- factor(result, labels = levs)
      structure(result,
                names = nx,
                class = class(x),
                index = index)
      
    } else return(x) # ak == 0 => nothing to do (no lagging/no leading)
  }
  
  if (length(k) > 1){
    rval <- sapply(k, function(i) alag(x, i))
    colnames(rval) <- k
  }
  else {
    rval <- alag(x, k)
  }
  return(rval)
}

lead <- function (x, k = 1, ...) {
  UseMethod("lead")
}

lead.pseries <- function(x, k = 1, ...) {
  ret <- lag.pseries(x, k = -k)
  if (length(k) > 1) colnames(ret) <- k
  return(ret)
}


