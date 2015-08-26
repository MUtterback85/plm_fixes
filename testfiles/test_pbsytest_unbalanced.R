# test pbsytest()



################### Bera, Sosa-Escudero and Yoon (2001) and joint test of Baltagi/Li (1991) ###############
# see Baltagi (2005), Econometric Analysis of Panel Data, 3rd edition, pp. 96-97
#  or Baltagi (2013), Econometric Analysis of Panel Data, 5th edition, pp. 108:
# "For the Grunfeld data, we computed LM_mu = 798.162 in Table 4.2 using the xttest0 command
# in Stata. Using TSP, LM_p = 143.523, LM*_mu = 664.948, LM*_p = 10.310 and the joint LM1 statistic
# in (5.36) is 808.471.
# The joint test rejects the null of no first-order serial correlation and no random firm effects.
# The one-directional tests LMp and LM*_p reject the null of no first-order serial correlation,
# while the one-directional tests LM_mu and LM*_mu reject the null of no random firm effects."


# ATTENTION: take pooling model as input!
#            p-value of pbsytest(..., test="j") in v1.4-0 wrong due to wrong degrees of freedom => fixed in plm_fixes_v1.4-0_kt.R
#
require(plm)
data("Grunfeld")
Grunfeldpdata <- pdata.frame(Grunfeld, index = c("firm", "year"), drop.index = FALSE, row.names = TRUE)
pool_grunfeld  <- plm(inv ~ value + capital, data=Grunfeldpdata, model="pooling")

# original implementation (v1.4-0) (plm:::pbsytest.panelmodel() this is meant to call the version 1.4-0 as on CRAN)
plm:::pbsytest.panelmodel(pool_grunfeld, test = "ar") # chisq = 10.31                   => LM*_p in Baltagi's book
plm:::pbsytest.panelmodel(pool_grunfeld, test = "re") # z = 25.787 [25.787^2=664.9694]  => sqrt(LM*_mu) in Baltagi's book
plm:::pbsytest.panelmodel(pool_grunfeld, test = "j")  # chisq = 808.47                  => LM1 statistic in Baltagi's book

### results from Bera et al. (2001), p. 13:
# To replicate, a special version of the Grunfeld data set is needed: only 5 selected firms (total of 100 obs)
# from http://pages.stern.nyu.edu/~wgreene/Text/tables/TableF13-1.txt
# or   http://statmath.wu.ac.at/~zeileis/grunfeld/TableF13-1.txt

Grunfeld_greene_5firms <- read.csv("M:/Stat_tests/datasets/Greene_Grunfeld_TableF13-1.txt", sep="")
pGrunfeld_greene_5firms <- pdata.frame(Grunfeld_greene_5firms, index = c("Firm", "Year"), drop.index = FALSE, row.names = TRUE)
pool_grunfeld_half <- plm(I ~ F + C, data=pGrunfeld_greene_5firms, model="pooling")
re_grunfeld_half <- plm(I ~ F + C, data=pGrunfeld_greene_5firms, model="random")
plm:::pbsytest.panelmodel(pool_grunfeld_half, test = "ar") # chisq = 3.7125 => RS*_p
plm:::pbsytest.panelmodel(pool_grunfeld_half, test = "re") # z = 19.601 => ()^2 => chisq = 384.199 => RS*_u 
plm:::pbsytest.panelmodel(pool_grunfeld_half, test = "j")  # chisq = 457.53 => RS_up
plmtest(pool_grunfeld_half, type = "bp")                   # chisq = 453.82 => RS_u


############### unbalanced version ###################

################ load new implemtation first, for this test:
################ need to provide it under the name pbsytest_mod_unbalanced.panelmodel
################ if already loaded under the real name, use e. g. the follwing for convenience:
#     pbsytest_mod_unbalanced.panelmodel <- pbsytest.panelmodel 

# results for balanced data set grunfeld
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld, test = "ar")
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld, test = "re")
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld, test = "j")

# parts of the grunfeld data (for tests as in Bera et al. (2001), p. 13)
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld_half, test = "ar") # chisq = 3.7125; p = 0.054 => RS*_p
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld_half, test = "re") # chisq = 384.18; p = 0     => RS*_u
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld_half, test = "j")  # chisq = 457.53; p = 0     => RS_up

pbsytest_mod_unbalanced.panelmodel(pool_grunfeld_half, test = "re", normal = T) # normal = 19.601; p = 0 => RSO*_u (tiny diff due to rounding) in Bera et al. (2001), p. 13

# should result in an error
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld_half, test = "j",  normal = T)
pbsytest_mod_unbalanced.panelmodel(pool_grunfeld_half, test = "ar", normal = T)



############ Replicate tests from original paper Sosa-Escudero/Bera (2008) ####################
## data set for test from Sosa-Escudero/Bera (2008), pp. 75-77
## available as STATA .dta at http://www.stata-journal.com/software/sj8-1/sg164_1/ginipanel5.dta
require(haven)
ginipanel5 <- read_dta("http://www.stata-journal.com/software/sj8-1/sg164_1/ginipanel5.dta")
pginipanel5 <- pdata.frame(ginipanel5, index = c("naglo", "ano"), drop.index = FALSE, row.names = TRUE)

# STATA command for RE model: xtreg gini ie ie2 indus adpubedsal desempleo tactiv invipib apertura pyas4 e64 supc tamfam, re
# use pooling model in R:
pool_gini <- plm(gini ~ ie + ie2 + indus + adpubedsal + desempleo + tactiv + invipib + apertura + pyas4 + e64 + supc + tamfam, data=pginipanel5, model="pooling")

# STATA's Output of xttest1, unadjusted (Sosa-Escudero/Bera (2008), p. 77):
#
# Random Effects, Two Sided:
#   LM(Var(u)=0) = 13.50 Pr>chi2(1) = 0.0002
# ALM(Var(u)=0) = 6.03 Pr>chi2(1) = 0.0141               # test="re"
#
# Random Effects, One Sided:
#   LM(Var(u)=0) = 3.67 Pr>N(0,1) = 0.0001
# ALM(Var(u)=0) = 2.46 Pr>N(0,1) = 0.0070                # test="re", normal = T
#
# Serial Correlation:
#   LM(lambda=0) = 9.32 Pr>chi2(1) = 0.0023
# ALM(lambda=0) = 1.86 Pr>chi2(1) = 0.1732               # test="ar"
#
# Joint Test:
#   LM(Var(u)=0,lambda=0) = 15.35 Pr>chi2(2) = 0.0005    # test="j"

pbsytest_mod_unbalanced.panelmodel(pool_gini, test = "ar")
pbsytest_mod_unbalanced.panelmodel(pool_gini, test = "re")
pbsytest_mod_unbalanced.panelmodel(pool_gini, test = "j")
pbsytest_mod_unbalanced.panelmodel(pool_gini, test = "re", normal = T)

# compare old implementation vs. new implementation (unbalanced) should result in a sufficently large difference [TRUE]
stat_re_normal            <- plm:::pbsytest.panelmodel(pool_gini, test = "re") # original version 1.4-0, test="re" already normalized
stat_re_normal_unbalanced <- pbsytest_mod_unbalanced.panelmodel(pool_gini, test = "re", normal = T)
abs(stat_re_normal$statistic - stat_re_normal_unbalanced$statistic) > 0.01


