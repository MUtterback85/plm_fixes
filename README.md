#Own quick'n'dirty fixes (?)/enhancements to plm version 1.4-0

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

##### pbsytest(..., test="j") (Baltagi/Li (1991)): fixed degrees of freedom
also added check for correct input

##### Baltagi/Li (1995): implemented statistic and test for LM5
LM test for first-order serial correlation in a fixed effects model

##### Updated documentation
For some updates to plm's documentation (mainly new text book editions) see https://github.com/helix123/plm/tree/master/man

##### lag.pseries() can handle negative lags (leading values)
lead.pseries() is also added for convenience

## How to use
 
 Load this file after package plm is loaded. The following functions are then masked:
   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
   - summary.plm() and print.summary.plm() (used by summary())
   - plmtest()
   - pbgtest.panelmodel() 
   - pbltest(): added panelmodel interface
   - pbsytest(): fixed degrees of freedom error when test="j" (Baltagi/Li (1991))
   - pbltest_lm5: new function added to compute test statistic LM5 from Baltagi/Li (1995)
   - lag.pseries() can handle negative lags (leading values); lead.pseries() is added for convenience

