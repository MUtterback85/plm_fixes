# Own quick'n'dirty fixes (?) to plm version 1.4-0
# no warranty
#
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

# Instruction
# load this file after package plm is loaded
# The following functions are then masked:
#   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
#   - summary.plm() and print.summary.plm() (used by summary())
#


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
              random.method.list[ercomp],
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
