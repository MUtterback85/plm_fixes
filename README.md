Own quick'n'dirty fixes (?) to plm version 1.4-0
no warranty


# F statistic when a user specified variance-covariance matrix is supplied (for robust inference)
 [see http://stackoverflow.com/questions/31163353/computing-f-statistic-with-user-supplied-covariance-matrix]

 Note: Formula interface not fixed

# Durbin-Watson test respecting panel structure
  [see http://stackoverflow.com/questions/31894055/pdwtest-from-plm-with-wrong-p-value-for-pooled-ols-durbin-watson-test-for-autoc]

Note: There are some points noted down below which I am not sure how to handle
       (marked with:  Statisticians: Can someone please look into this?)


# Instructions
 load this file after package plm is loaded
 The following functions are then masked:
   - pdwtest.panelmodel() and pdwtest.formula (used by pdwtest())
   - summary.plm() and print.summary.plm() (used by summary())

