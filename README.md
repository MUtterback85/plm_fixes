# plm_fixes
little fixes to package plm (panel model estimation) for R

## F statistic when a user specified variance-covariance matrix is supplied (for robust inference)
#  [see http://stackoverflow.com/questions/31163353/computing-f-statistic-with-user-supplied-covariance-matrix]
#
#  Forumula interface not fixed

## Durbin-Watson test respecting panel structure
#  [see http://stackoverflow.com/questions/31894055/pdwtest-from-plm-with-wrong-p-value-for-pooled-ols-durbin-watson-test-for-autoc]


# Instruction:
# load this file after package plm is loaded
# The following functions are then masked:
#   - pdwtest.panelmodel() (used by pdwtest())
#   - summary.plm() (used by summary())
