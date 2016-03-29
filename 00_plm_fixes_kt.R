# Own quick'n'dirty fixes (?)/enhancements to plm version 1.5-16/rev. 222 as on r-forge (development version)
# and lmtest as on CRAN (2015-12-30)
#
# plm's development version on r-forge (see there for how to install): https://r-forge.r-project.org/R/?group_id=406
# See original package plm on CRAN: https://cran.r-project.org/package=plm
# See also orginal package lmtest: https://cran.r-project.org/package=lmtest
# See there for original authors of these packages.
#
# In this file, some routines are copied over from the original packages and are modified.
# Some routines are new.
#
# Version of this file 0.9
#
# no warranty
#
# License: GPL v2, v3
#          http://www.gnu.org/licenses/gpl-2.0.html
#          http://www.gnu.org/licenses/gpl-3.0.html
#
#
# Find this file also at https://github.com/helix123/plm_fixes
#
# Instructions:
# Load this file after package plm is loaded. Modified functions are then (mostly) loaded into the global environment;
# new functions become available.
#
#   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest()): Durbin-Watson test respecting panel structure,
#                see http://stackoverflow.com/questions/31894055/pdwtest-from-plm-with-wrong-p-value-for-pooled-ols-durbin-watson-test-for-autoc]
#                References:
#                Bhargava, Franzini, Narendranathan, Serial Correlation and the Fixed Effects Model, Review of Economic Studies (1982), XLIX, pp. 533-549.
#                Baltagi, B. H., and P. X. Wu. 1999. Unequally spaced panel data regressions with AR(1) disturbances. Econometric Theory 15, pp 814-823.
#
#                Note: There are some points noted down below which I am not sure how to handle
#                (marked with:  Statisticians: Can someone please look into this?)
#
#   - summary.plm() and print.summary.plm() (used by summary()): F statistic corrected when a user specified variance-covariance matrix
#                                                                is supplied (for robust inference)
#                             see http://stackoverflow.com/questions/31163353/computing-f-statistic-with-user-supplied-covariance-matrix
#
#   - plmtest():  [not needed anymore as of (at least) v1.5-16]
#                * all tests (bp, honda, kw, ghm) of original plmtest implemented for unbalanced panels;
#                * use correct mixed chisquare distribution (=chibarsquare) for type="ghm";
#                * fixed p-values for some tests (now everything according to Baltagi's testdata)
#
#   - pbgtest() [Breusch-Godfrey test]: [not needed anymore as of (at least) v1.5-13]
#                           allows to pass on type="F" to lmtest::bgtest(), thus offering the small sample test (F test)
#
#   - pbsytest(): * Fixed degrees of freedom error when test="j" (test of Baltagi/Li (1991), A joint test for serial correlation and random individual effects).
#                 * Added unbalanced panel version from Sosa-Escudero/Bera (2008)
#                 * Added warning if wrong input model.
#
#                 [code is also in seperate branch on r-forge: https://r-forge.r-project.org/scm/viewvc.php/branches/kt_unbalanced/pbsytest/?root=plm]
#
#   - pbltest_lm5() [added]: new function added to compute test statistic LM5 from Baltagi/Li (1995):
#                                               An LM test for first-order serial correlation in a fixed effects model
#                      [can also be used in a random effects model]
#
#   - lag.pseries() [not needed anymore as of (at least) v1.5-14/rev. 176]
#                   * can handle negative lags (leading values); lead.pseries() is added as a wrapper for convenience
#
#   - pbltest(): added: panelmodel interface for convenience
#
#   - pwtest(): [not needed anymore as of (at least) v1.5-16 (rev. 200), 2016-02-15]
#               pwtest.panelmodel: fixed: respect effect argument for panelmodel interface of test (formula interface was not affected)
#
#   - pbptest(): added: Breusch-Pagan test against heteroskedasticity for panelmodels (wrapper which uses lmtest::bptest())
#
#   - pgqtest(): added: Goldfeld-Quandt test against heteroskedasticity for panelmodels (wrapper which uses lmtest::gqtest())
#                original lmtest::gqtest (CRAN v0.9-34) slightly modified to return alternative hypothesis in returned htest object
#
#   - r.squared(): * Adjusted R-squared corrected for pooling models (before, it did not match lm's adj. R-squared)
#                  NB: For pooling models without intercept, the regular R-squared and the adjusted R-squared still diverge
#                      from lm's (adj.) R-squared, a warning is printed.
#
#   - nobs.plm(): [not needed anymore as of (at least) v1.5-13]
#                  added nobs() function for convenience to extract number of total observations used for estimated of plm model (like nobs() for lm models)
#
#   - phtest.formula(): [not needed anymore as of (at least) v1.5-13]
#                      backported regression-based Hausman test from SVN repository; fixed data handling etc.
#
#   - fitted.plm(): issue warning if model has silently dropped coefficients (likely due to linear dependency) which
#                   will give a additional information on why the method fails
#                   [NB: for plm package development: Change implementation of fitted.plm to accomodate such models with dropped coefficients?
#                                                     Align behaviour of plm() with lm() [lm keeps coefficients and sets them NA]?]
#
#   - mylm() [not exported from package]: added warning if coefficients get dropped during model estimation

#### load package plm first ########
library(plm)
library(lmtest)
if (packageVersion("plm")    != "1.5.16") stop("This fixes/enhancements are against plm version 1.5-16 from r-forge (published 2016-02-15)")
if (packageVersion("lmtest") != "0.9.34") stop("This fixes/enhancements are against lmtest version 0.9-34 from CRAN (published 2015-06-06)")

options(warnPartialMatchAttr   = TRUE,
        warnPartialMatchDollar = TRUE,
        warnPartialMatchArgs   = TRUE)
options(showWarnCalls = TRUE)

################## pdwtest.panelmodel() adapted from pseries.R [Durbin-Watson test] ##################
pdwtest.panelmodel <- function(x, ...) {
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
  # STATA has modified. Bhargava et al. (1982) and Baltagi/Wu LBI: http://www.stata.com/manuals14/xtxtregar.pdf
  # 
  
  # For FE and RE, we  need to take the panel structure into account
  # Verbeek (2004, 2nd edition), p. 358:
  # "Because the fixed effects estimator is also consistent in the random effects model, it is
  # also possible to use this panel data Durbin-Watson test in the latter model."
  
  # 
  # Statisticians: Can someone please look into this?
  #

    if (pdim(x)$balanced != TRUE) warning("Applying Bhargava et al. (1982) Durbin-Watson test for balanced panels to an unbalanced panel.")
  
    # residuals are now class pseries, so diff.pseries is used and the differences are computed within observational units
    # (not across as it would be the case if base::diff() is used and as it is done for lm-objects)
    # NAs are introduced by the differencing as one observation is lost per observational unit
    dw <- sum(plm:::diff.pseries(residuals(x))^2, na.rm = T) / sum(residuals(x)^2)
    
    # p-value computation seems to be difficult, so no p-value is calculated
    # maybe someone with more statistical knowledge can look into this one
    
    # constuct htest object
    names(dw) <- "DW_balance_panel"
    ARtest <- list(statistic = dw,
                   method = "Bhargava et al. (1982): Durbin-Watson test for serial correlation in balanced panel models (pooled OLS, within (fixed) and random effects)\n
                   No p-value computed as distubution of statistic under panel assumptions is difficult to calculate.",
                   alternative = NULL,
                   p.value = NULL 
                   ) # data.name = paste(deparse(x$formula))
    
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


################## summary.plm() adapted from plm.methods.R ##################

# For F test: in the robust case, there seems to be a need to adjust the degrees of freedom to cluster_size -1 for panel models (at least FE models)
# Cameron/Miller: paper and data to test: http://cameron.econ.ucdavis.edu/research/papers.html
#  see Stata example 3 on p 16:  http://www.stata.com/manuals14/xtxtreg.pdf
# Cameron, A.C., J.B. Gelbach and D.L. Miller (2011),
#  Robust inference with multiway clustering, Journal of Business & Economic Statistics 29 (2011), no. 2, 238–249. http://dx.doi.org/10.1198/jbes.2010.07136
# and: http://www.stata.com/manuals14/p_robust.pdf#p_robust
# Eviews (middle of page: wald test with df(1,53): http://www.eviews.com/help/helpintro.html#page/EViews%209%20Help/panel.058.3.html
# http://stats.stackexchange.com/questions/93787/f-test-formula-under-robust-standard-error
#
# However, Gretl (2016a - official) currently does use the normal degrees of freedom
# need at least the development version of 2016b  (built date 2016-03-26)

summary.plm <- function(object, .vcov = NULL, ...){
  object$fstatistic <- plm:::Ftest(object, test = "F", .vcov = .vcov) # Ftest(object, test = "F") # TODO: "activate" robust F test by: Ftest(object, test = "F", .vcov = .vcov)
  model <- plm:::describe(object, "model")
  effect <- plm:::describe(object, "effect")
  object$r.squared <- c(rsq  = r.squared(object),
                        adjrsq = r.squared(object, dfcor = TRUE))
  # construct the table of coefficients
  if (!is.null(.vcov)) {
    if (is.matrix(.vcov))   rvcov <- .vcov
    if (is.function(.vcov)) rvcov <- .vcov(object)
    std.err <- sqrt(diag(rvcov))
  }
  else {
    std.err <- sqrt(diag(vcov(object)))
  }
  b <- coefficients(object)
  z <- b / std.err
  p <- 2 * pt(abs(z), df = object$df.residual, lower.tail = FALSE)
  
  # construct the object of class summary.plm
    object$coefficients <- cbind("Estimate"   = b,
                                 "Std. Error" = std.err,
                                 "t-value"    = z,
                                 "Pr(>|t|)"   = p)
    if (!is.null(.vcov)) {
      # put the robust vcov in summary.plm object (next to "normal" vcov)
      object$rvcov <- rvcov
      attr(object$rvcov, which = "rvcov.name") <- paste0(deparse(substitute(.vcov)))
    }
    class(object) <- c("summary.plm", "plm", "panelmodel")
  object
} ## END summary.plm
assignInNamespace("summary.plm", summary.plm, envir = as.environment("package:plm"))
rm(summary.plm)

# adjusted R-squared corrected (now meets lm's adjusted R-squared for pooling models)
# still: pooling models without intercept: R-squared and adj. R-squared are not correct
r.squared <- function(object, model = NULL,
                      type = c('cor', 'rss', 'ess'), dfcor = FALSE){
  if (is.null(model)) model <- plm:::describe(object, "model")
  effect <- plm:::describe(object, "effect")
  type <- match.arg(type)
  if (type == 'cor'){
    y <- pmodel.response(object, model = model, effect = effect)
    haty <- fitted.plm(object, model = model, effect = effect)
    R2 <- cor(y, haty)^2
  }
  if (type == 'rss'){
    R2 <- 1 - deviance(object, model = model) / tss(object, model = model)
  }
  if (type == 'ess'){
    haty <- fitted(object, model = model)
    mhaty <- mean(haty)
    ess <- sum( (haty - mhaty)^2)
    R2 <- ess / tss(object, model = model)
  }
  if (dfcor) {
    R2 <- 1-(1-R2) * (length(resid(object))-1) / df.residual(object) # [does not account for models without intercept!] was:  R2 * df.residual(object) / length(resid(object) - 1)
    if (!"(Intercept)" %in% attr(object$coefficients, "names") & model == "pooling") warning("R-Squared (regular and adjusted) for pooling models without intercept is not correct. Use lm() to get the correct value.")
  }
  return(R2)
}
# END r.squared

### fitted.plm(): added test for existence of all coefficients (dropped coefficients due to linear dependency)
fitted.plm <- function(object, model = NULL, ...){
  # there are two 'models' used ; the fitted model and the
  # transformation used for the fitted values
  fittedmodel <- plm:::describe(object, "model")
  if (is.null(model)) model <- fittedmodel
  effect <- plm:::describe(object, "effect")
  X <- model.matrix(object, model = model)
  y <- pmodel.response(object, model = model)
  beta <- coef(object)
  if (model == "within" & fittedmodel != "within"){
    Xw <- model.matrix(object, model = "within", effect = effect)
    varwith <- colnames(Xw)
    beta <- beta[varwith]
  }
  
  # Test if all coefficients could be estimated by plm
  # [plm silently drops non-estimatable coefficients [v1.5-13]]
  # With this test, we provide an additional warning message to
  # the user to enhance the error message from failing crossprod later in the code
  # which relies on non-dropped coefficients; see also testfile tests/fitted.plm.R
  # Test could be computationally/space expensive due to creation of model.matrix.
  if (!setequal(names(object$coefficients), colnames(model.matrix(object)))) {
     warning("Coefficients of estimated model do not match variables in its specified model.matrix.
            This is likely due to non-estimatable coefficients (compare object$formula with object$coefficients).")
  }
  
  if (fittedmodel == "within"){
    if (model == "pooling"){
      if (has.intercept(object)) X <- X[,-1]
      index <- attr(model.frame(object), "index")
      if (effect != "time") id <- index[[1]]
      if (effect != "individual") time <- index[[2]]
      fe <- switch(effect,
                   individual = fixef(object, effect = "individual")[as.character(id)],
                   time = fixef(object, effect="time")[as.character(time)],
                   twoways = fixef(object, effect = "individual")[as.character(id)] +
                             fixef(object, effect = "time")[as.character(time)])
      fv <- as.numeric(crossprod(t(X), beta)) + fe
    }
    if (model == "between"){
      alpha <- mean(y) - crossprod(apply(X[, -1], 2, mean), beta)
      beta <- c(alpha, beta)
      fv <- as.numeric(crossprod(t(X), beta))
    }
    if (model == "within"){
      fv <- as.numeric(crossprod(t(X), beta))
    }
  }
  else{
    fv <- as.numeric(crossprod(t(X), beta))
  }
  return(structure(fv, index =  index(object), class = "pseries"))
}
## END fitted.plm()


# mylm: added: warning if coefficents are dropped during estimation
mylm <- function(y, X, W = NULL){
  names.X <- colnames(X)
  if (is.null(W))
      result <- lm(y ~ X - 1)
  else
      result <- twosls(y, X, W)
  if (any(is.na(coef(result)))){
    na.coef <- is.na(coef(result))
    warning("Coefficient(s) '", paste((names.X)[na.coef], collapse = ", "), "' could not be estimated and is (are) dropped.")
    X <- X[, !na.coef, drop = FALSE]
    if (is.null(W)) result <- lm(y ~ X - 1)
    else result <- twosls(y, X, W)
  }
  result$vcov <- vcov(result)
  result$X <- X
  result$y <- y
  result$W <- W
  names(result$coefficients) <- colnames(result$vcov) <-
    rownames(result$vcov) <- colnames(X)
  return(result)
}
# assign modified mylm to plm's namespace, needed because mylm is not an exported function in package plm
assignInNamespace("mylm", mylm, envir = as.environment("package:plm"))
rm(mylm)
# End: mylm()





##### lag.pseries with warning if non-consecutive time series detected ####
lag.pseries <- function(x, k = 1, ...) {
  index <- attr(x, "index")
  id <- index[[1]]
  time <- index[[2]]

  # catch the case when an index of pdata.frame shall be lagged (index variables are always factors)
    # NB: this catches - unintentionally - also the case when a factor variable is the same "on the character level"
    # as one of the corresponding index variables but not the index variable itself
    #
    # -> shall we prevent lagging of index variables at all? -> turned off for now, 2016-03-03
    # if (is.factor(x)) if (all(as.character(x) == as.character(id)) | all(as.character(x)==as.character(time))) stop("Lagged vector cannot be index.")

  alag <- function(x, ak){
    if (round(ak) != ak) stop("Lagging value 'k' must be whole-numbered (positive, negative or zero)")
    if (ak > 0) {

      # NB: this code assumes consecutive time periods and produces wrong results
      #     for lag > 1 and non-consecutive time periods

      # delete first ak observations for each unit
      isNAtime <- c(rep(T, ak), (diff(as.numeric(time), lag = ak) != ak))
      isNAid   <- c(rep(T, ak), (diff(as.numeric(id),   lag = ak) != 0))
      isNA <- (isNAtime | isNAid)

      result <- x                                             # copy x first ...
      result[1:ak] <- NA                                      # ... then make first ak obs NA ...
      result[(ak+1):length(result)] <- x[1:(length(x)-ak)]    # ... shift and ...
      result[isNA] <- NA                                      # ... make more NAs in between: this way, we keep: all factor levels, names, classes

    } else if (ak < 0) { # => compute leading values

      # NB: this code assumes consecutive time periods and produces wrong results
      #     for lag > 1 and non-consecutive time periods

      # delete last |ak| observations for each unit
      num_time <- as.numeric(time)
      num_id   <- as.numeric(id)
      isNAtime <- c(c((num_time[1:(length(num_time)+ak)] - num_time[(-ak+1):length(num_time)]) != ak), rep(T, -ak))
      isNAid   <- c(c((num_id[1:(length(num_id)+ak)]     - num_id[(-ak+1):length(num_id)])     != 0),  rep(T, -ak))
      isNA <- (isNAtime | isNAid)

      result <- x                                            # copy x first ...
      result[(length(result)+ak+1):length(result)] <- NA     # ... then make last |ak| obs NA ...
      result[1:(length(result)+ak)] <- x[(1-ak):(length(x))] # ... shift and ...
      result[isNA] <- NA                                     # ... make more NAs in between: this way, we keep: all factor levels, names, classes

    } else { # ak == 0 => nothing to do, return original pseries (no lagging/no leading)
      result <- x
    }

    return(result)
  } # END function alag

  if (length(k) > 1) {
    rval <- sapply(k, function(i) alag(x, i))
    colnames(rval) <- k
  }
  else {
    rval <- alag(x, k)
  }
  if(!all(is.pconsecutive(x, na.rm.tindex = TRUE))) warning("non-consecutive time periods")
  if(!all(is.pconsecutive(x, na.rm.tindex = TRUE)) & abs(k) > 1) warning("non-consecutive time periods and k > 1")
  return(rval)
}
assignInNamespace("lag.pseries", lag.pseries, envir = as.environment("package:plm"))
rm(lag.pseries)


############## plmtest() ############################################
## not needed anymore as of plm v1.15-16 on r-forge
# [code is also in seperate branch on r-forge: https://r-forge.r-project.org/scm/viewvc.php/branches/kt_unbalanced/plmtest/?root=plm]
# modified to handle unbalanced panels for all test statistics 
#
# For a concise overview with original references see
# Baltagi (2013), Econometric Analysis of Panel Data, 5th edition, pp. 68-76 (balanced), pp. 200-203 (unbalanced).
#
# balanced (original) version of Breusch-Pagan test: T.S. Breusch & A.R. Pagan (1979),
#                                           A Simple Test for Heteroscedasticity and Random Coefficient Variation, Econometrica 47, pp. 1287-1294
#
# unbalanced version: Baltagi/Li (1990), A lagrange multiplier test for the error components model with incomplete panels,
#                    Econometric Reviews, 9:1, pp. 103-107,

# plmtest <- function(x,...){
#   UseMethod("plmtest")
# }
# 
# plmtest.plm <- function(x,
#                         effect = c("individual", "time", "twoways"),
#                         type = c("honda", "bp", "ghm", "kw"),
#                         ...) {
#   
#   effect <- match.arg(effect)
#   type <- match.arg(type)
#   if (plm:::describe(x, "model") != "pooling") x <- update(x, model = "pooling")
#   pdim <- pdim(x)
#   n <- pdim$nT$n
#   T <- pdim$nT$T
#   N_obs <- pdim$nT$N
#   balanced <- pdim$balanced
#   index <- attr(model.frame(x), "index")
#   id <- index[[1]]
#   time <- index[[2]]
#   T_i <- pdim$Tint$Ti
#   N_t <- pdim$Tint$nt
#   res <- resid(x)
#   
#   ### calc of parts of test statistic ##
#   # calc. is done w/o using matrix calculation, see e.g. Baltagi/Li (1990), p. 106
#   A1 <- as.numeric(crossprod(tapply(res,id,sum))/sum(res^2) - 1)   # == A1 <- sum(tapply(res,id,sum)^2)/sum(res^2) - 1
#   A2 <- as.numeric(crossprod(tapply(res,time,sum))/sum(res^2) - 1) # == A2 <- sum(tapply(res,time,sum)^2)/sum(res^2) - 1
#   
#   M11 <- sum(T_i^2)
#   M22 <- sum(N_t^2)
#   
#   LM1 <- N_obs * (1/sqrt(2*(M11 - N_obs))) * A1 # == sqrt( (((N_obs)^2) / 2) * ( A1^2 / (M11 - N_obs)) ) [except sign due to positive sqrt]
#   LM2 <- N_obs * (1/sqrt(2*(M22 - N_obs))) * A2 # == sqrt( (((N_obs)^2) / 2) * ( A2^2 / (M22 - N_obs)) ) [except sign due to positive sqrt]
#   ### END calc of parts of test statistic ##
#   
#   
#   if (effect != "twoways"){
#     # oneway
#     if (!type %in% c("honda", "bp", "kw"))
#       stop("type must be one of \"honda\", \"bp\" or \"kw\" for a one way model") # kw oneway coincides with honda
#     
#     ifelse(effect == "individual", stat <- LM1, stat <- LM2)
#     stat <- switch(type,
#                    honda = c(normal = stat),
#                    bp    = c(chisq  = stat^2),
#                    kw    = c(normal = stat))
#     parameter <- switch(type,
#                           honda = NULL,
#                           bp = c(df = 1),  # df = 1 in the oneway case (Baltagi (2013), p. 70)
#                           kw = NULL)
#     pval <- switch(type,
#                      honda = pnorm(stat, lower.tail = FALSE),  # honda oneway ~ N(0,1), alternative is one-sided (Baltagi (2013), p. 71/202)
#                      bp    = pchisq(stat, df = parameter, lower.tail = FALSE), # is df=1 in the one-way case, alternative is two-sided (Baltagi (2013), p. 70/201)
#                      kw    = pnorm(stat, lower.tail = FALSE)) # kw oneway ~ N(0,1), alternative is one-sided (Baltagi (2013), p. 71/202)
#   }
#   else { # twoways
#     stat <- switch(type,
#                      ghm   = c(chibarsq = max(0,LM1)^2+max(0,LM2)^2),
#                      bp    = c(chisq = LM1^2+LM2^2),
#                      honda = c(normal = (LM1+LM2)/sqrt(2)),
#                      kw    = c(normal = (sqrt(M11-N_obs)/sqrt(M11+M22-2*N_obs))*LM1+(sqrt(M22-N_obs)/sqrt(M11+M22-2*N_obs))*LM2))
#     
#     parameter <- switch(type,
#                           ghm   = c(df0 = 0L, df1=1L, df2=2L, w0=1/4, w1=1/2, w2=1/4),
#                           bp    = c(df = 2),
#                           honda = NULL,
#                           kw    = NULL)
#     
#     pval <- switch(type,
#                      ghm   = (1/4)*pchisq(stat, df=0, lower.tail = F) + (1/2) * pchisq(stat, df=1, lower.tail = F) + (1/4) * pchisq(stat, df=2, lower.tail = F), # mixed chisq (also called chi-bar-square), see Baltagi (2013), pp. 71-72, 74, 88, 202-203, 209
#                      honda = pnorm(stat,lower.tail = FALSE),  # honda two-ways ~ N(0,1), alternative is one-sided (Baltagi (2013), p. 71/202)
#                      bp    = pchisq(stat, df = parameter,lower.tail = FALSE), # is df = 2 in the twoway case, alternative is two-sided (Baltagi (2013), p. 70/201)
#                      kw    = pnorm(stat, lower.tail = FALSE)) # kw twoways ~ N(0,1), alternative is one-sided (Baltagi (2013), p. 71/202)
#   }
#   
#   method.type <- switch(type,
#                           honda  = "Honda",
#                           bp     = "Breusch-Pagan",
#                           ghm    = "Gourieroux, Holly and Monfort",
#                           kw     = "King and Wu")
#   
#   method.effect <- switch(effect,
#                             id      = "individual effects",
#                             time    = "time effects",
#                             twoways = "two-ways effects")
#   
#   balanced.type <- ifelse(balanced, "balanced", "unbalanced")
#   
#   method <- paste("Lagrange Multiplier Test - ", method.effect,
#                   " (",method.type,") for ", balanced.type, " panels", sep="")
#   
#   if(type %in% c("honda", "kw")) {
#     RVAL <- list(statistic = stat,
#                  p.value   = pval,
#                  method    = method,
#                  data.name = plm:::data.name(x))
#   }
#   else {
#     RVAL <- list(statistic = stat,
#                  p.value   = pval,
#                  method    = method,
#                  parameter = parameter,
#                  data.name = plm:::data.name(x))
#   }
#   RVAL$alternative <- "significant effects" # TODO: maybe distinguish be b/w one-sided and two-sided alternatives? (bp: two-sided alt.; all others: one-sided alt.?)
#   class(RVAL) <- "htest"
#   return(RVAL)
# } ## END plmtest.plm()
# 
# 
# plmtest.formula <- function(x, data, ...,
#                             effect = c("individual", "time", "twoways"),
#                             type = c("honda", "bp", "ghm", "kw")) {
#   
#   cl <- match.call(expand.dots = TRUE)
#   cl$model <- "pooling" # plmtest is performed on the pooling model...
#   cl$effect <- NULL     # ... and pooling model has no argument effect...
#   cl$type <- NULL       # ... and no argument type => see below: pass on args effect and type to plmtest.plm()
#   names(cl)[2] <- "formula"
#   m <- match(plm.arg, names(cl), 0)
#   cl <- cl[c(1,m)]
#   cl[[1]] <- as.name("plm")
#   plm.model <- eval(cl, parent.frame())
#   plmtest(plm.model, effect = effect, type = type) # pass on args. effect and type
# }


############### Breusch-Godfrey test ##################################
### Fixes for pbgtest() not needed anymore as of (at least) plm v1.5-13
### fixed pbgtest(), copied over from https://r-forge.r-project.org/scm/viewvc.php/pkg/R/pserial.R?view=markup&root=plm&pathrev=127
# pbgtest() suffered from the same problem as pdwtest() [intercept passed twice to lmtest::bgtest()] - incorporated that from r-forge
#
# additional fix: match arguments, so that type="F" (and order.by=) and is passed on to lmtest::bgtest(), thus enabling the small sample
#                 variant of the test offered by lmtest::bgtest()
# 
# 
# pbgtest.panelmodel<-function(x, order = NULL, ...) {
#   ## residual serial correlation test based on the residuals of the demeaned
#   ## model (see Wooldridge p.288) and the regular bgtest() in {lmtest}
# 
#   ## structure:
#   ## 1: take demeaned data from 'plm' object
#   ## 2: est. auxiliary model by OLS on demeaned data
#   ## 3: apply bgtest() to auxiliary model and return the result
# 
#   model <- plm:::describe(x, "model")
#   effect <- plm:::describe(x, "effect")
#   theta <- x$ercomp$theta
# 
#   ## retrieve demeaned data
#   demX <- model.matrix(x, model = model, effect = effect, theta=theta)
#   demy <- pmodel.response(model.frame(x), model = model, effect = effect, theta=theta)
# 
#   ## ...and group numerosities
#   Ti <- pdim(x)$Tint$Ti
#   ## set lag order to minimum group numerosity if not specified by user
#   ## (check whether this is sensible)
# 
#   if(is.null(order)) order <- min(Ti)
#   ## bg test on the demeaned model:
# 
#   ## check package availability and load if necessary
#   #lm.ok <- require("lmtest")
#   #if(!lm.ok) stop("package lmtest is needed but not available")
# 
#   ## bgtest is the bgtest, exception made for the method attribute
#   dots <- match.call(expand.dots=FALSE)[["..."]]      # fixed: added expand.dots=FALSE
#   if (!is.null(dots$type)) type <- dots$type else type <- "Chisq"
#   if (!is.null(dots$order.by)) order.by <- dots$order.by else order.by <- NULL
# 
#   auxformula <- demy~demX-1 #if(model == "within") demy~demX-1 else demy~demX
#   lm.mod <- lm(auxformula)
#   bgtest <- bgtest(lm.mod, order = order, type = type, order.by = order.by)
#   bgtest$method <- "Breusch-Godfrey/Wooldridge test for serial correlation in panel models [KT]"
#   bgtest$alternative <- "serial correlation in idiosyncratic errors"
#   bgtest$data.name <- paste(deparse(x$call$formula))
#   names(bgtest$statistic) <- if(length(bgtest$parameter)==1) "chisq" else "F"
#   return(bgtest)
# }




######### pbltest(): added panelmodel interface ########
# Baltagi and Li Serial Dependence Test For Random Effects Models
#
# from source code: test statistic (one-sided) is LM_4 [see also original paper: Baltagi/Li (1995), p. 137]
# [or Baltagi (2013), p. 108 mentions it, p. 113: LM_4 ~ N(0,1), alternative is one-sided]
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

### KT: copy formula interface here, otherwise the wrapper won't work
pbltest.formula <- function(x, data, alternative = c("twosided", "onesided"), index=NULL, ...) {
  ## this version (pbltest0) based on a "formula, pdataframe" interface
  
  
  ## reduce X to model matrix value (no NAs)
  X<-model.matrix(x,data=data)
  ## reduce data accordingly
  data <- data[which(row.names(data)%in%row.names(X)),]
  
  
  data <- pdata.frame(data,index=index)
  
  ## need name of individual index
  gindex <- dimnames(attr(data, "index"))[[2]][1]
  
  ## make random effects formula
  rformula <- NULL
  eval(parse(text=paste("rformula <- ~1|",gindex,sep="")))
  
  ## est. MLE model
  mymod <- nlme::lme(x,data=data,random=rformula,method="ML")
  
  nt. <- mymod$dims$N
  n. <- as.numeric(mymod$dims$ngrps[1])
  t. <- nt./n.
  Jt <- matrix(1,ncol=t.,nrow=t.)/t.
  Et <- diag(1,t.)-Jt
  ## make 'bidiagonal' matrix (see BL, p.136)
  G <- matrix(0,ncol=t.,nrow=t.)
  for(i in 2:t.) {
    G[i-1,i] <- 1
    G[i,i-1] <- 1
  }
  
  ## retrieve composite (=lowest level) residuals
  uhat <- residuals(mymod,level=0)
  
  ## sigma2.e and sigma2.1 as in BL
  ## break up residuals by group to get rid of Kronecker prod.
  ## data have to be balanced and sorted by group/time, so this works
  uhat.i <- vector("list",n.)
  for(i in 1:n.) {
    uhat.i[[i]] <- uhat[t.*(i-1)+1:t.]
  }
  s2e <- rep(NA,n.)
  s21 <- rep(NA,n.)
  for(i in 1:n.) {
    u.i <- uhat.i[[i]]
    s2e[i] <- as.numeric(crossprod(u.i,Et) %*% u.i)
    s21[i] <- as.numeric(crossprod(u.i,Jt) %*% u.i)
  }
  sigma2.e <- sum(s2e) / (n.*(t.-1))
  sigma2.1 <- sum(s21) / n.
  
  ## calc. score under the null:
  star1 <- (Jt/sigma2.1 + Et/sigma2.e) %*% G %*% (Jt/sigma2.1 + Et/sigma2.e)
  star2 <- rep(NA,n.)
  ## again, do this group by group to avoid Kronecker prod.
  for(i in 1:n.) {
    star2[i] <- as.numeric(crossprod(uhat.i[[i]],star1) %*% uhat.i[[i]])
  }
  star2 <- sum(star2)
  Drho <- (n.*(t.-1)/t.) * (sigma2.1-sigma2.e)/sigma2.1 + sigma2.e/2 * star2
  ## star2 is (crossprod(uhat, kronecker(In, star1)) %*% uhat)
  
  ## components for the information matrix
  a <- (sigma2.e-sigma2.1)/(t.*sigma2.1)
  j.rr <- n. * (2 * a^2 * (t.-1)^2 + 2*a*(2*t.-3) + (t.-1))
  j.12 <- n.*(t.-1)*sigma2.e / sigma2.1^2
  j.13 <- n.*(t.-1)/t. * sigma2.e * (1/sigma2.1^2 - 1/sigma2.e^2)
  j.22 <- (n. * t.^2) / (2 * sigma2.1^2)
  j.23 <- (n. * t.) / (2 * sigma2.1^2)
  j.33 <- (n./2) * (1/sigma2.1^2 + (t.-1)/sigma2.e^2)
  
  ## build up information matrix
  Jmat <- matrix(nrow=3,ncol=3)
  Jmat[1,] <- c(j.rr,j.12,j.13)
  Jmat[2,] <- c(j.12,j.22,j.23)
  Jmat[3,] <- c(j.13,j.23,j.33)
  
  J11 <- n.^2 * t.^2 * (t.-1) / (det(Jmat) * 4*sigma2.1^2 * sigma2.e^2)
  ## this is the same as J11 <- solve(Jmat)[1,1], see BL page 73
  
  switch(match.arg(alternative),
         onesided = {
           LMr.m <- Drho * sqrt(J11)
           pval <- pnorm(LMr.m,lower.tail=FALSE)
           names(LMr.m) <- "z"
           method1 <- "one-sided"
           method2 <- "H0: rho = 0, HA: rho > 0"
           parameter <- NULL
         },
         twosided = {
           LMr.m <- Drho^2 * J11
           pval <- pchisq(LMr.m,1,lower.tail=FALSE)
           names(LMr.m) <- "chisq"
           parameter <- c(df=1)
           method1 <- "two-sided"
           method2 <- "H0: rho = 0, HA: rho != 0"
         }
  )
  dname <- paste(deparse(substitute(x)))
  method <- paste("Baltagi and Li", method1,"LM test")
  alternative <- paste0("AR(1)/MA(1) errors in RE panel models. ", method2)
  
  res <- list(statistic = LMr.m,
              p.value = pval,
              method = method,
              alternative = alternative,
              parameter = parameter,
              data.name = dname)
  
  class(res) <- "htest"
  res
}
# END pbltest







############# Baltagi/Li (1995), LM5 - An LM test for first-order serial correlation in a fixed effects model (can also be used in RE models)
# LM5 (p. 138-139: An LM test for first-order serial correlation in a fixed effects model):
# LM5 is LM3 (from p. 136) with the residuals of the FE model instead of the OLS model (in matrix B)
# Thus, matrix B calculated in pblsytest is already in the right form for LM5 but need to
# subsitute the residuals for the FE residuals.

# p. 136: "LM1, is exactly the same as the joint test statistic derived by Baltagi and Li (1991) for
# AR(1) residual disturbances and random individual effects." => matrix B calculated in pblsytest(),
# adapted for use here
#
# see (equivalent) Baltagi (2005), Econometric Analysis of Panel Data, 3rd edition, pp. 97-98 or
#                  Baltagi (2013), Econometric Analysis of Panel Data, 5th edition, pp. 108-109
#
# Baltagi (2005), p. 98 [=Baltagi (2013), p. 109]:
# "Since the [w]ithin transformation wipes out the individual effects whether fixed or random,
#  one can also use [LM5] to test for serial correlation in the random effects models."
#
# Baltagi/Li (1995), p. 138 or Baltagi (2013), pp. 109, 113:
#       onesided. Under H0, this is asymptotically distributed (for large T) as N(0, 1).

pbltest_lm5 <- function(x, ...) { # panelmodel interface (plm interface)
  
  if (plm:::describe(x, "model") == "pooling") stop("Test only for within effects (=fixed effects) or random effects models.")
  
##### following code adapted from pbsytest() #####
  resfe <- resid(x)
  data <- model.frame(x)
  
  ## extract indices
  indexvars <- attr(data, "index")
  index <- indexvars[[1]]
  tindex <- indexvars[[2]]
    
  ## till here. 
  ## ordering here if needed.
    
  ## this needs ordering of obs. on time, regardless 
  ## whether before that on groups or after
    
  ## and numerosity check
    
  ## order by group, then time
  oo <- order(index,tindex)
  ind <- index[oo]
  tind <- tindex[oo]
  resfe <- resfe[oo]
  ## det. number of groups and df
  n <- length(unique(index))
  k <- ncol(model.matrix(x))
  ## det. max. group numerosity
  t <- max(pdim(x)$Tint$Ti)
  ## det. total number of obs. (robust vs. unbalanced panels)
  nT <- length(ind)
    
  ## calc. B (with FE residuals):
  unind <- unique(ind)
  uu <- rep(NA,length(unind))
  uu1 <- rep(NA,length(unind))
  for(i in 1:length(unind)) {
    u.t <- resfe[ind==unind[i]]
    u.t.1 <- u.t[-length(u.t)]
    u.t <- u.t[-1]
    uu[i] <- crossprod(u.t)
    uu1[i] <- crossprod(u.t,u.t.1)
  }

  B <- sum(uu1)/sum(uu)
##### END code adapted from pblsytest() #####

  LM5_statsitic <- sqrt((n * t^2) / (t-1)) * B # see calculation in Baltagi/Li (1995), p. 138
  names(LM5_statsitic) <- "LM5"
  
  dname <- paste(deparse(substitute(x$formula)))
  RVAL <- list(statistic = LM5_statsitic,
               parameter = NULL,
               method = "LM test for first-order serial correlation in a fixed effects model \n
                         LM5 in Baltagi/Li (1995), p. 138-139",
               alternative = "AR(1) errors (rho > 0) sub fixed effects",
               p.value = pnorm(LM5_statsitic,lower.tail=FALSE), # Baltagi/Li (1995), p. 138 or Baltagi (2013), pp. 109, 113: onesided.
               data.name = dname)                               #  Under H0, this is asymptotically distributed (for large T) as N(0, 1).
  class(RVAL) <- "htest"
  return(RVAL)
} ## END: pbltest_lm5


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
## changed: chisq test (two-sided) is standard; for one-sided alternative (normalized version) of test="re", use argument normal = TRUE
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

pbsytest.panelmodel <- function(x, test = c("ar","re","j"), normal = FALSE, ...) {
  
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
  S1 <- sum(tapply(poolres,ind,sum)^2)
  S2 <- sum(poolres^2)
  
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
           pstat <- pchisq(stat, df=df, lower.tail=FALSE) # # Degrees of freedom in the joint test (test="j") of Baltagi/Li (1991) should be 2 (chisquare(2) distributed, see Baltagi/Li (1991), p. 279 and again in Baltagi/Li (1995), p. 136
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
} ###### END pbsytest.panelmodel: ##################




############# lag.pseries: able to create negative lags (=leading values) (use k < 0)
# [not needed anymore as of (at least) v1.5-14/rev. 176]
## for convenience: method lead.pseries is also added
#
# also as diff against plm_v1.4-0 on github: https://github.com/cran/plm/compare/cd00e7ce878f8ffa0e0bc53ab8692778cb8aaecc...1ef2620b2851ae6d7374d52bfed8141e11eaff76
# 
# lag.pseries <- function(x, k = 1, ...) {
#   nx <- names(x)
#   index <- attr(x, "index")
#   id <- index[[1]]
#   time <- index[[2]]
#   
#   # catch the case when an index of pdata.frame shall be lagged
#   if (is.factor(x)) if (all(as.character(x) == as.character(id)) | all(as.character(x)==as.character(time))) stop("Lagged vector cannot be index.")
#   
#   alag <- function(x, ak){
#     if (round(ak) != ak) stop("Lagging value 'k' must be whole-numbered (positive, negative or zero)")
#     if (ak > 0){
#       # delete first ak observations for each unit
#       isNAtime <- c(rep(T,ak), diff(as.numeric(time), lag = ak)) != ak
#       isNAid <- c(rep(T,ak), diff(as.numeric(id), lag = ak)) != 0
#       isNA <- as.logical(isNAtime + isNAid)
#       if (is.factor(x)) levs <- levels(x)
#       result <- c(rep(NA, ak), x[1:(length(x)-ak)])
#       result[isNA] <- NA
#       if (is.factor(x)) result <- factor(result, labels = levs)
#       structure(result,
#                 names = nx,
#                 class = class(x),
#                 index = index)
#     } else if (ak < 0) { # => compute leading values
#       
#       # delete last ak observations for each unit
#       isNAtime <- c(as.numeric(time) - c(tail(as.numeric(time), length(time) + ak), rep(T, -ak))) != ak
#       isNAid   <- c(as.numeric(id) - c(tail(as.numeric(id), length(id) + ak) , rep(T, -ak))) != 0
#       isNA <- as.logical(isNAtime + isNAid)
#       result <- c(x[(1-ak):(length(x))], rep(NA, -ak))
#       result[isNA] <- NA
#       if (is.factor(x)) levs <- levels(x)
#       if (is.factor(x)) result <- factor(result, labels = levs)
#       structure(result,
#                 names = nx,
#                 class = class(x),
#                 index = index)
#       
#     } else return(x) # ak == 0 => nothing to do (no lagging/no leading)
#   }
#   
#   if (length(k) > 1){
#     rval <- sapply(k, function(i) alag(x, i))
#     colnames(rval) <- k
#   }
#   else {
#     rval <- alag(x, k)
#   }
#   return(rval)
# }
# 
# lead <- function (x, k = 1, ...) {
#   UseMethod("lead")
# }
# 
# lead.pseries <- function(x, k = 1, ...) {
#   ret <- lag.pseries(x, k = -k)
#   if (length(k) > 1) colnames(ret) <- k
#   return(ret)
# }

#### pwfdtest.panelmodel
##
## added: informative informative message if wrong model is passed

pwfdtest.panelmodel <- function(x, ..., h0=c("fd","fe")) {
  ## first-difference-based serial correlation test for panel models
  ## ref.: Wooldridge (2003), par. 10.6 
  if (plm:::describe(x, "model") != "fd") stop("pwfdtest only relevant for first-differenced models") ## added this line
  
  if(!require(car)) stop("Library 'car' is needed")
  
  ## fetch fd residuals
  FDres <- resid(x)
  ## indices (full length! must reduce by 1st time period)
  ## this is an ad-hoc solution for the fact that the 'fd' model
  ## carries on the full indices while losing the first time period
  index <- attr(model.frame(x), "index")
  time <- as.numeric(index[[2]])
  id <- as.numeric(index[[1]])
  
  ## fetch dimensions and adapt to those of indices
  pdim <- pdim(x)
  n <- pdim$nT$n
  
  
  ## (re)create groupwise-separated index from 1 to nT 
  ## - dropping first time period
  ## - correcting Ti=Ti+1
  Ti <- pdim$Tint$Ti-1
  
  redind <- vector("list",n)
  tfirst <- 0
  for(i in 1:n) {
    redind[[i]] <- (tfirst+2):(tfirst+Ti[i]+1)
    tfirst <- max(redind[[i]])
  }
  ## reduce indices by 1st time period
  redind <- unlist(redind)
  time <- time[redind]
  id <- id[redind]
  
  N <- length(FDres)
  
  FDres.1 <- c(NA,FDres[1:(N-1)])
  
  lagid  <-  id-c(NA,id[1:(N-1)])
  
  FDres.1[lagid!=0] <- NA
  
  
  ## make (panel) dataframe for auxiliary regression
  auxdata <- as.data.frame(cbind(id,time))
  auxdata$FDres <- FDres
  auxdata$FDres.1 <- FDres.1
  ## pooling model FDres vs. lag(FDres),
  ## with intercept (might as well do it w.o.)
  auxmod <- plm(FDres ~ FDres.1, na.omit(auxdata), model = "pooling")
  
  switch(match.arg(h0), 
         fd = {h0des<-"differenced"
         ## theoretical rho under H0: no serial 
         ## corr. in differenced errors is 0
         rho.H0 <- 0},
         fe = {h0des<-"original"
         ## theoretical rho under H0: no serial 
         ## corr. in original errors is -0.5
         rho.H0 <- -0.5})
  
  ## test H0: rho=rho.H0 with HAC t-test (HC0-3 parm may be passed)
  myvcov <- function(x) vcovHC(x, method="arellano", ...)
  
  myH0 <- paste("FDres.1 = ", as.character(rho.H0), sep="")
  lhtest <- linearHypothesis(model=auxmod, myH0, vcov.=myvcov, ...)
  
  ##(insert usual htest features)  
  FDARstat <- lhtest[2,3]
  names(FDARstat) <- dimnames(lhtest)[[2]][3] 
  if (names(FDARstat)=="Chisq") names(FDARstat) <- "chisq"
  ## this is either 'F' or 'Chisq' and is the name of 3rd
  ## column because we are supplying a vcov matrix
  pFDAR<-lhtest[2,4]
  
  dname <- paste(deparse(substitute(x)))
  RVAL <- list(statistic = FDARstat, parameter = NULL,
               method = "Wooldridge's first-difference test for serial correlation in panels",
               alternative = paste("serial correlation in", h0des, "errors"),
               p.value = pFDAR,
               data.name =   dname)
  class(RVAL) <- "htest"
  return(RVAL)
  
}
## END pwfdtest


### pwtest(): ### # [not needed anymore as of (at least) v1.5-16 (rev. 200), 2016-02-15]
###           fixed panelmodel interface which did not respect the effect parameter, i. e.
###           for a supplied panelmode effect="individual" and effect="time" deliver the same result for CRAN version 1.4-0
###           formula interface is not affected
# pwtest <- function(x, ...){
#   UseMethod("pwtest")
# }
# 
# pwtest.formula <- function(x, data, ...) {
#   cl <- match.call(expand.dots = TRUE)
#   if (names(cl)[3] == "") names(cl)[3] <- "data"
#   if (is.null(cl$model)) cl$model <- "pooling"
#   if (cl$model != "pooling") stop("pwtest only relevant for pooling models")
#   names(cl)[2] <- "formula"
#   m <- match(plm:::plm.arg,names(cl),0)
#   cl <- cl[c(1,m)]
#   cl[[1]] <- as.name("plm")
#   plm.model <- eval(cl,parent.frame())
#   effect <- plm:::describe(plm.model, "effect")
#   pwtest.panelmodel(plm.model, effect=effect)
#   
#   ## "RE" test ? la Wooldridge, see 10.4.4
#   ## (basically the scaled and standardized estimator for sigma from REmod)
#   ## does not rely on normality or homoskedasticity; 
#   ## H0: composite errors uncorrelated
#   
#   ## ref. Wooldridge, p.264
#   
#   ######### from here generic testing interface from
#   ######### plm to my code
# }
# 
# pwtest.panelmodel <- function(x, effect=c("individual", "time"), ...){
#   ## tind is actually not needed here
#   if (plm:::describe(x, "model") != "pooling") stop("pwtest only relevant for pooling models")
#   
#   effect <- match.arg(effect) # effect <- plm:::describe(x, "effect") # here we want the effect as in the call of pwtest(), not of the already estimated model
#                               # or do we?
#   data <- model.frame(x)
#   ## extract indices
#   
#   ## if effect="individual" std., else swap
#   index <- attr(data, "index")
#   if (effect == "individual"){
#     index <- index[[1]]
#     tindex <- index[[2]]
#   }
#   else {
#     index <- index[[2]]
#     tindex <- index[[1]]
#   }
#   ## det. number of groups and df
#   n <- length(unique(index))
#   X <- model.matrix(x)
#   
#   k <- ncol(X)
#   ## det. total number of obs. (robust vs. unbalanced panels)
#   nT <- nrow(X)
#   ## det. max. group numerosity
#   t <- max(tapply(X[,1],index,length))
#   
#   ## ref. Wooldridge, p.264
#   
#   ## extract resids
#   u <- resid(x)
#   
#   ## est. random effect variance
#   ## "pre-allocate" an empty list of length n
#   tres <- vector("list", n)
#   
#   ## list of n "empirical omega-blocks"
#   ## with averages of xproducts of t(i) residuals
#   ## for each group 1..n 
#   ## (possibly different sizes if unbal., thus a list
#   ## and thus, unlike Wooldridge (eq.10.37), ve divide 
#   ## every block by *his* t(t-1)/2)
#   #  unind <- unique(ind)
#   unind <- unique(index) # ????
#   
#   for(i in 1:n) {
#     ut <- u[index == unind[i]]
#     tres[[i]] <- ut%o%ut
#   }
#   
#   ## sum over all upper triangles of emp. omega blocks:
#   ## define aux. function
#   uptrisum <- function(x) {
#     uts <- sum(x[upper.tri(x,diag=FALSE)])
#     return(uts)}
#   
#   ## det. # of upper triangle members (n*t(t-1)/2 if balanced)
#   ti <- sapply(tres, function(x) dim(x)[[1]])
#   uptrinum <- sum(ti*(ti-1)/2)  # don't need this!!
#   
#   ## ...apply to list and sum over resulting vector (df corrected)
#   W <- sum(sapply(tres,uptrisum)) # /sqrt(n) simplifies out
#   
#   ## calculate se(Wstat) as in 10.40
#   seW <- sqrt( sum( sapply(tres,uptrisum)^2 ) )
#   
#   ## NB should we apply a df correction here, maybe that of the standard
#   ## RE estimator? (see page 261) 
#   
#   Wstat <- W/seW
#   names(Wstat) <- "z"
#   pW <- 2*pnorm(abs(Wstat),lower.tail=FALSE) # unlike LM, test is two-tailed!
#   
#   ##(insert usual htest features)
#   dname <- paste(deparse(substitute(formula)))
#   RVAL <- list(statistic = Wstat, parameter = NULL,
#                method = paste("Wooldridge's test for unobserved ",
#                               effect,"effects "),
#                alternative = "unobserved effect",
#                p.value = pW,
#                data.name =   dname)
#   class(RVAL) <- "htest"
#   return(RVAL)
#   
# }
# # END pwtest



## Breusch-Pagan test against heteroskedasticity for panelmodels
#
# only panelmodel interface implemented, not formula interface
# Code skeleton adapted from pbgtest() 

pbptest <-function(x, ...) {
  ## residual heteroskedasticity test based on the residuals of the demeaned
  ## model and the regular bptest() in {lmtest}
  
  ## structure:
  ## 1: take demeaned data from 'plm' object
  ## 2: est. auxiliary model by OLS on demeaned data
  ## 3: apply bptest() to auxiliary model and return the result
  
  if (!inherits(x, "plm")) stop("need to supply a panelmodel estimated with plm()")
  model <- plm:::describe(x, "model")
  effect <- plm:::describe(x, "effect")
  theta <- x$ercomp$theta
  
  ## retrieve demeaned data
  demX <- model.matrix(x, model = model, effect = effect, theta = theta)
  demy <- pmodel.response(model.frame(x), model = model, effect = effect, theta = theta)
  
  ## ...and group numerosities
  Ti <- pdim(x)$Tint$Ti
  ## set lag order to minimum group numerosity if not specified by user
  ## (check whether this is sensible)
  
  if (is.null(order)) order <- min(Ti)
  ## bg test on the demeaned model:
  
  ## check package availability and load if necessary
  lm.ok <- require("lmtest")
  if(!lm.ok) stop("package lmtest is needed but not available")
  
  ## bptest is the bptest, exception made for the method attribute
  dots <- match.call(expand.dots=FALSE)[["..."]]      # fixed: added expand.dots=FALSE
  if (!is.null(dots$type)) type <- dots$type else type <- "Chisq"
  if (!is.null(dots$order.by)) order.by <- dots$order.by else order.by <- NULL
  
  auxformula <- demy~demX-1
  lm.mod <- lm(auxformula)
  return(lmtest::bptest(lm.mod, ...)) # call and return bptest from package 'lmtest'
} # END pbptest()

## Goldfeld-Quandt test against heteroskedasticity for panelmodels
#
# only panelmodel interface implemented, not formula interface
# Code skeleton from pbgtest() adapted
# lmtest::gqtest on https://cran.r-project.org/package=lmtest

pgqtest <- function(x, ...) {
  ## residual heteroskedasticity test based on the residuals of the demeaned
  ## model and the regular gqtest() in {lmtest}
  
  ## structure:
  ## 1: take demeaned data from 'plm' object
  ## 2: est. auxiliary model by OLS on demeaned data
  ## 3: apply gqtest() to auxiliary model and return the result
  
  if (!inherits(x, "plm")) stop("need to supply a panelmodel estimated with plm()")
  model <- plm:::describe(x, "model")
  effect <- plm:::describe(x, "effect")
  theta <- x$ercomp$theta
  
  ## retrieve demeaned data
  demX <- model.matrix(x, model = model, effect = effect, theta = theta)
  demy <- pmodel.response(model.frame(x), model = model, effect = effect, theta = theta)
  
  ## ...and group numerosities
  Ti <- pdim(x)$Tint$Ti
  ## set lag order to minimum group numerosity if not specified by user
  ## (check whether this is sensible)
  
  if (is.null(order)) order <- min(Ti)
  ## bg test on the demeaned model:
  
  ## check package availability and load if necessary
  lm.ok <- require("lmtest")
  if(!lm.ok) stop("package lmtest is needed but not available")
  
  ## bptest is the bptest, exception made for the method attribute
  dots <- match.call(expand.dots=FALSE)[["..."]]      # fixed: added expand.dots=FALSE
  if (!is.null(dots$type)) type <- dots$type else type <- "Chisq"
  if (!is.null(dots$order.by)) order.by <- dots$order.by else order.by <- NULL
  
  auxformula <- demy~demX-1
  lm.mod <- lm(auxformula)
  return(gqtest(lm.mod, ...)) # call and return gqtest from package 'lmtest'
} # END pgqtest()

######## Original Goldfeld-Quandt test (gqtest) from lmtest
# CRAN v0.9-34
# added: return alternative
# fixed: when order.by is specified by a formula, make sure the right observations for ordering are used
gqtest <- function(formula, point = 0.5, fraction = 0,
                   alternative = c("greater", "two.sided", "less"), order.by = NULL, data = list())
{
  dname <- paste(deparse(substitute(formula)))
  alternative <- match.arg(alternative)
  
  if(!inherits(formula, "formula")) {
    X <- if(is.matrix(formula$x))
      formula$x
    else model.matrix(terms(formula), model.frame(formula))
    y <- if(is.vector(formula$y))
      formula$y
    else model.response(model.frame(formula))
  } else {
    mf <- model.frame(formula, data = data)
    y <- model.response(mf)
    X <- model.matrix(formula, data = data)
  }  
  
  k <- ncol(X)
  n <- nrow(X)
  if(point > 1) {
    if(fraction < 1) fraction <- floor(fraction * n)
    point1 <- point - ceiling(fraction/2)
    point2 <- point + ceiling(fraction/2 + 0.01)
  } else {
    if(fraction >= 1) fraction <- fraction/n
    point1 <- floor((point-fraction/2) * n)
    point2 <- ceiling((point+fraction/2) * n + 0.01)
  }
  if (point2 > n-k+1 | point1 < k) stop("inadmissable breakpoint/too many central observations omitted")
  
  if(!is.null(order.by))
  {
    if(inherits(order.by, "formula")) {
    
      if(inherits(formula, "formula")) {
        z <- model.matrix(order.by, data = data)
        z <- z[rownames(z) %in% rownames(X), ] # added this line
        z <- as.vector(z[,ncol(z)])
      } else { # is lm-model
        z <- X[ , as.character(terms(order.by)[[2]])]
      }
      

    } else {
      z <- order.by
    }
    X <- as.matrix(X[order(z),])
    y <- y[order(z)]
  }
  
  rss1 <- sum(lm.fit(as.matrix(X[1:point1,]),y[1:point1])$residuals^2)
  rss2 <- sum(lm.fit(as.matrix(X[point2:n,]),y[point2:n])$residuals^2)
  mss <- c(rss1/(point1-k), rss2/(n-point2+1-k))
  
  gq <- mss[2]/mss[1]
  df <- c(n-point2+1-k, point1-k)
  names(df) <- c("df1", "df2")
  
  PVAL <- switch(alternative,
                 "two.sided" = (2*min(pf(gq, df[1], df[2]), pf(gq, df[1], df[2], lower.tail = FALSE))),
                 "less" = pf(gq, df[1], df[2]),
                 "greater" = pf(gq, df[1], df[2], lower.tail = FALSE))
  
  alternative <- switch(alternative,
                        "two.sided" = "variance changes from segment 1 to 2",
                        "less" = "variance decreases from segment 1 to 2",
                        "greater" = "variance increases from segment 1 to 2")
  
  method <- "Goldfeld-Quandt test"
  names(gq) <- "GQ"
  RVAL <- list(statistic = gq,
               parameter = df,
               method = method,
               p.value= PVAL,
               alternative = alternative, # added: alternative = alternative
               data.name=dname)
  
  class(RVAL) <- "htest"
  return(RVAL)
}
# END gqtest




### Imported from r-forge, development version of plm rev. 125
### Regression-based Hausman test which can be made robust
# [not needed anymore as of (at least) v1.5-13]
#
# Wooldridge (2010), pp. 328-334 (for robust test esp. p. 332-333)
# Baltagi (2013), pp. 76-77
# phtest.formula <- function(x, data, model = c("within","random"),
#                            method = c("chisq", "aux"),
#                            index=NULL, vcov=NULL, ...){
#   if(length(model)!=2) stop("two models should be indicated")
#   for (i in 1:2){
#     model.name <- model[i]
#     if(!(model.name %in% names(plm:::model.plm.list))){
#       stop("model must be one of ",oneof(model.plm.list))
#     }
#   }
#   switch(match.arg(method),
#          chisq={
#            cl <- match.call(expand.dots = TRUE)
#            cl$model <- model[1]
#            names(cl)[2] <- "formula"
#            m <- match(plm:::plm.arg,names(cl),0)
#            cl <- cl[c(1,m)]
#            cl[[1]] <- as.name("plm")
#            plm.model.1 <- eval(cl,parent.frame())
#            plm.model.2 <- update(plm.model.1, model = model[2])
#            return(phtest(plm.model.1, plm.model.2))
#          },
#          aux={
#  ## some interface checks here
#                if(model[1]!="within") {
#                    stop("Please supply 'within' as first model type")
#                }
#              
#                ## set pdata
#                if (!inherits(data, "pdata.frame")) data <- plm.data(data, indexes=index) #, ...)
#                
#                row.names(data) <- NULL # reset rownames of original data set (number rownames in clean sequence) to make rownames
#                                        # comparable for later comparision to obs used in estimation of models (get rid of NA values)
#                                        # [needed becausepmodel.response() and model.matrix() do not retain fancy rownames, but rownames]
#                
#                # calculatate FE and RE model
#                fe_mod <- plm(formula=x, data=data, model=model[1])
#                re_mod <- plm(formula=x, data=data, model=model[2])
#                
#                reY <- pmodel.response(re_mod)
#                reX <- model.matrix(re_mod)[ , -1] # intercept not needed
#                feX <- model.matrix(fe_mod)
#                dimnames(feX)[[2]] <- paste(dimnames(feX)[[2]],
#                                            "tilde", sep=".")
#                
#                ## estimated models could have fewer obs (due droping of NAs) compared to the original data
#                ## => match original data and observations used in estimated models
#                ## routine adapted from lmtest::bptest
#                commonrownames <- intersect(intersect(intersect(row.names(data), names(reY)), row.names(reX)), row.names(feX))
#                if (!(all(c(row.names(data) %in% commonrownames, commonrownames %in% row.names(data))))) {
#                  data <- data[commonrownames, ]
#                  reY  <- reY[commonrownames]
#                  reX  <- reX[commonrownames, ]
#                  feX  <- feX[commonrownames, ]
#                }
#                
#                # Tests of correct matching of obs (just for safety ...)
#                 if (!all.equal(length(reY), nrow(data), nrow(reX), nrow(feX)))
#                   stop("number of cases/observations do not match, most likely due to NAs in \"data\"")
#                 if (any(c(is.na(names(reY)), is.na(row.names(data)), is.na(row.names(reX)), is.na(row.names(feX)))))
#                     stop("one (or more) rowname(s) is (are) NA")
#                 if (!all.equal(names(reY), row.names(data), row.names(reX), row.names(feX)))
#                   stop("row.names of cases/observations do not match, most likely due to NAs in \"data\"")
# 
#                ## fetch indices here, check pdata
#                ## construct data set and formula for auxiliary regression
#                data <- data.frame(cbind(data[, 1:2], reY, reX, feX))
#                auxfm <- as.formula(paste("reY~",
#                                          paste(dimnames(reX)[[2]],
#                                                collapse="+"), "+",
#                                          paste(dimnames(feX)[[2]],
#                                                collapse="+"), sep=""))
#                auxmod <- plm(formula=auxfm, data=data, model="pooling")
#                nvars <- dim(feX)[[2]]
#                R <- diag(1, nvars)
#                r <- rep(0, nvars) # here just for clarity of illustration
#                omega0 <- vcov(auxmod)[(nvars+2):(nvars*2+1),
#                                       (nvars+2):(nvars*2+1)]
#                Rbr <- R %*% coef(auxmod)[(nvars+2):(nvars*2+1)] - r
# 
#                h2t <- crossprod(Rbr, solve(omega0, Rbr))
#                ph2t <- pchisq(h2t, df=nvars, lower.tail=FALSE)
# 
#                df <- nvars
#                names(df) <- "df"
#                names(h2t) <- "chisq"
# 
#                if(!is.null(vcov)) {
#                    vcov <- paste(", vcov: ",
#                                   paste(deparse(substitute(vcov))),
#                                   sep="")
#                }
# 
#                haus2 <- list(statistic   = h2t,
#                              p.value     = ph2t,
#                              parameter   = df,
#                              method      = paste("Regression-based Hausman test",
#                                               vcov, sep=""),
#                              alternative = "one model is inconsistent",
#                              data.name   = paste(deparse(substitute(x))))
#                class(haus2) <- "htest"
#                return(haus2)
#          })
# }

# [not needed anymore as of (at least) v1.5-13]
# nobs.plm <- function(x) {
#   if (inherits(x, "plm") | inherits(x, "panelmodel")) return(plm::pdim(x)$nT$N) else stop("Input x needs to be of class 'plm' (or 'panelmodel'), i. e. a panel model estimated by plm()")
# }


## texreg: avoid partial matching
library(texreg)
extract.plm <- function(model, include.rsquared = TRUE, include.adjrs = TRUE, 
    include.nobs = TRUE, ...) {
  s <- summary(model, ...)
  
  coefficient.names <- rownames(coef(s))
  coefficients <- coef(s)[, 1]
  standard.errors <- coef(s)[, 2]
  significance <- coef(s)[, 4]
  
  rs <- s$r.squared[1]
  adj <- s$r.squared[2]
  n <- length(residuals(s))
  
  gof <- numeric()
  gof.names <- character()
  gof.decimal <- logical()
  if (include.rsquared == TRUE) {
    gof <- c(gof, rs)
    gof.names <- c(gof.names, "R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.adjrs == TRUE) {
    gof <- c(gof, adj)
    gof.names <- c(gof.names, "Adj.\ R$^2$")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (include.nobs == TRUE) {
    gof <- c(gof, n)
    gof.names <- c(gof.names, "Num.\ obs.")
    gof.decimal <- c(gof.decimal, FALSE)
  }
  
  tr <- createTexreg(
      coef.names = coefficient.names, 
      coef = coefficients, 
      se = standard.errors, 
      pvalues = significance, 
      gof.names = gof.names, 
      gof = gof, 
      gof.decimal = gof.decimal
  )
  return(tr)
}
## end
assignInNamespace("extract.plm", extract.plm, envir = as.environment("package:texreg"))
rm(extract.plm)
