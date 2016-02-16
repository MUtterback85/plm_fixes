# Own quick'n'dirty fixes/enhancements to plm version 1.5-16
see https://cran.r-project.org/package=plm and https://r-forge.r-project.org/projects/plm/

Some code is new by me, some is adapted or just slightly modified from the indicated source

License: GPL

no warranty

## How to use
 
 Load this file after package plm is loaded. The following functions are then loaded into the global environment and thus used instead of their respective counterparts of package plm:
   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
   - summary.plm() and print.summary.plm() (used by summary())
   - plmtest() - now in dev version of plm (v1.15-16)
   - pbgtest.panelmodel() 
   - pbltest()
   - pbsytest()
   - pbltest.panelmodel()
   - pbltest_lm5()
   - lag.pseries() - now in dev version of plm (v1.15-16)
   - lead.pseries() is added for convenience  - now in dev version of plm (v1.15-16)
   - pbgtest()
   - pbptest()
   - pgqtest()
   - r.squared()
   - nobs()  - now in dev version of plm (v1.15-16)

##### F statistic (summary())
 F statistic when a user specified variance-covariance matrix is supplied (for robust inference)
 [see http://stackoverflow.com/questions/31163353/computing-f-statistic-with-user-supplied-covariance-matrix]


##### Durbin-Watson test (pdwtest())
 Durbin-Watson test respecting panel structure
  [see http://stackoverflow.com/questions/31894055/pdwtest-from-plm-with-wrong-p-value-for-pooled-ols-durbin-watson-test-for-autoc]

Note: There are some points in the file which I am not sure how to handle
       (marked with:  Statisticians: Can someone please look into this?)


##### Baltagi/Li (1995) test: added panelmodel interface for convenience
plm 1.4-0 only offers the formula interface for this test. Using the formula interface is a bit cumbersome as it makes many assumptions on how the arguments should be formed. Calling the panelmodel interface is easier.

##### Baltagi/Li (1995): implemented statistic and test for LM5 [pbltest_lm5()]
LM test for first-order serial correlation in a fixed effects model

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
