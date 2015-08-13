# Own quick'n'dirty fixes (?) to plm version 1.4-0
# no warranty
#
# this version 0.3
#
# Instructions
# load this file after package plm is loaded
# The following functions are then masked:
#   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
#   - summary.plm() and print.summary.plm() (used by summary())
#   - plmtest() [additionally Baltagi/Li (1990) is directly available by baltagi_li()]



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
# see how STATA does it and this produces the same statistic: http://www.stata.com/manuals13/xtxtregpostestimation.pdf





################## pdwtest.panelmodel() adapted from pseries.R ##################
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
  # Implemented: approach of Bhargava et al. (1982)
  # reference: Bhargava/Franzini/Narendranathan, Serial Correlation and the Fixed Effects Model, Review of Economic Studies (1982), XLIX, pp. 533-549
  #
  # Not implemented:
  # Statisticians: Can someone please look into this?
  # There seems to be an modified version which accounts for unbalanced and unequally spaced data:
  # Baltagi, B. H., and P. X. Wu. 1999. Unequally spaced panel data regressions with AR(1) disturbances. Econometric Theory 15, pp 814-823.
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
    
    # residuals are now class pseries, so diff.pseries is used and the differences are computed within observational units
    # (not across as it would be the case if base::diff() is used and as it is done for lm-objects)
    # NAs are introduced by the differencing as one observation is lost per observational unit
    dw <- sum(plm:::diff.pseries(residuals(mod_to_test))^2, na.rm = T) / sum(residuals(mod_to_test)^2)
    
    # p-value computation seems to be difficult, so no p-value is calculated
    # maybe someone with more statistical knowledge can look into this one
    
    # constuct htest object
    names(dw) <- "DW"
    ARtest <- list(statistic = dw,
                   method = "Durbin-Watson test for serial correlation in panel models (pooled OLS, within (fixed) and random effects) \n
                   For random effects models, the Durbin-Watson statistic is calculated to the coresponding fixed effect model. \n
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
################  run instead of original plmtest() 
# Breusch-Pagan test for random effects model for unbalanced data 
# Baltagi/Li (1990), A lagrange multiplier test for the error components model with incomplete panels,
#                    Econometric Reviews, 9:1, 103-107, DOI: 10.1080/07474939008800180
#
# for individual effects
baltagi_li <- function (x,
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
  
  # Implementation for unblanaced panels
  if (effect != "individual") {stop("Not implemented. Only individual effects are implemented.")}  

  T_i <- as.numeric(table(id))
  T_mean <- mean(T_i)
  
  A1 <-  1 - sum(tapply(res,id,sum)^2)/sum(res^2)
  term <- sum(T_i^2) - (n * T_mean)
  bp_stat_baltagi_li <- (((n * T_mean)^2) / 2) * ( A1^2 / term)
  
  names(bp_stat_baltagi_li) <- "BP_unbalanced"
  parameter <- 1
  names(parameter) <- "df"
  
  
  res <- list(statistic = bp_stat_baltagi_li,
              p.value = pchisq(bp_stat_baltagi_li, df = 1, lower.tail = FALSE),
              parameter = parameter,
              method = "Lagrange Multiplier Test - individual effects - Breusch-Pagan Test for unbalanced Panels as in Baltagi/Li (1990) \n",
              data.name = plm:::data.name(x),
              alternative = "significant effects")
  class(res) <- "htest"
  return(res)
} # END baltagi_li




############## plmtest() ############################################
## modified to handle unblanced panels as in Baltagi/li (1990) ######

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
    
    # Implementation for unblanaced panels
 
    return(baltagi_li(x))
    
  } else { ### balanced panel => use original implementation of plm
    
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
                  data.name = data.name(x))
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
}

