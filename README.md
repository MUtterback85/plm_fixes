#Own quick'n'dirty fixes (?) to plm version 1.4-0

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


 ## Instructions:
 
 Load this file after package plm is loaded. The following functions are then masked:
   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
   - summary.plm() and print.summary.plm() (used by summary())
   - plmtest()

