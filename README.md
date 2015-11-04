# Own quick'n'dirty fixes/enhancements to plm version 1.4-0
see https://cran.r-project.org/package=plm

Some code is new by me, some is adapted or just slightly modified from the indicated source

License: GPL

no warranty



##### F statistic (summary())
 F statistic when a user specified variance-covariance matrix is supplied (for robust inference)
 [see http://stackoverflow.com/questions/31163353/computing-f-statistic-with-user-supplied-covariance-matrix]

 Note: Formula interface not fixed

##### Durbin-Watson test (pdwtest())
 Durbin-Watson test respecting panel structure
  [see http://stackoverflow.com/questions/31894055/pdwtest-from-plm-with-wrong-p-value-for-pooled-ols-durbin-watson-test-for-autoc]

Note: There are some points in the file which I am not sure how to handle
       (marked with:  Statisticians: Can someone please look into this?)

##### Breusch-Pagan test for random effects (plmtest())
Breusch-Pagan test for random effects for unbalanced panels as in Baltagi/Li (1991)
 [see http://stackoverflow.com/questions/31988449/implementation-of-breusch-pagan-test-for-random-effects-in-plm-with-unbalanced-p]

##### Breusch-Godfrey test for autocorrelation (serial correlation) - enabling small sample test
Allow further arguments to be passed on to lmtest::bgtest(), esp. type="F", which enables the small sample test (F test)

##### Baltagi/Li (1995) test: added panelmodel interface for convenience
plm 1.4-0 only offers the formula interface for this test. Using the formula interface is a bit cumbersome as it makes many assumptions on how the arguments should be formed. Calling the panelmodel interface is easier.

##### plmtest(): fixed p-value, unbalanced version of tests

##### pbsytest(..., test="j") (Baltagi/Li (1991)): fixed degrees of freedom, unbalanced version of tests
also added check for correct input

##### Baltagi/Li (1995): implemented statistic and test for LM5 [pbltest_lm5()]
LM test for first-order serial correlation in a fixed effects model

##### Updated documentation
For some updates to plm's documentation (mainly new text book editions) see https://github.com/helix123/plm/tree/master/man

##### lag.pseries() can handle negative lags (leading values)
lead.pseries() is also added for convenience

##### plmtest(): all tests (bp, honda, kw, ghm) of original plmtest implemented for unbalanced panels;
 use correct mixed chisquare distribution (=chibarsquare) for type="ghm"

##### pbsytest(): Fixed degrees of freedom error when test="j" (test of Baltagi/Li (1991), A joint test for serial correlation and random individual effects).
 Added unbalanced panel version from Sosa-Escudero/Bera (2008)
 Added warning if wrong input model.

##### pbltest(): added panelmodel interface is added for convenience

##### pwtest(): pwtest.panelmodel: fixed: respect effect argument for panelmodel interface of test (formula interface was not affected)

##### pbptest(): added Breusch-Pagan test against heteroskedasticity for panelmodels (wrapper which uses lmtest::bptest())

#####  pgqtest(): added Goldfeld-Quandt test against heteroskedasticity for panelmodels (wrapper which uses lmtest::gqtest())

#####  lmtest::gqtest (CRAN v0.9-34) slightly modified to return alternative hypothesis in returned htest object
 
#####  r.squared(): Adjusted R-squared corrected for pooling models, now matches lm's adj. R-squared.
NB: For pooling models _without_ intercept, the regular R-squared and the adjusted R-squared still diverge from lm's (adj.) R-squared, a warning is printed. See also testfile: https://github.com/helix123/plm_fixes/blob/master/testfiles/test_adj_R2.R

##### nobs(): added nobs() function for convenience to extract number of total observations used for estimating plm model (like nobs() for lm models)



## How to use
 
 Load this file after package plm is loaded. The following functions are then masked or new:
   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
   - summary.plm() and print.summary.plm() (used by summary())
   - plmtest()
   - pbgtest.panelmodel() 
   - pbltest()
   - pbsytest()
   - pbltest.panelmodel()
   - pbltest_lm5()
   - lag.pseries()
   - lead.pseries() is added for convenience
   - pbgtest()
   - plmtest()
   - pbptest()
   - pgqtest()
   - r.squared()
   - nobs()

