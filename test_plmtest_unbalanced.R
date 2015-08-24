# Test of new plmtest implementation (handeling unbalanced panels)
#
# compare to Baltagi (2013), p. 74/75 (Table 4.1/4.2)


# Table 4.1
############ [(critical values at %5 level)]
# BP 798.162   6.454  804.615
#   (3.841)   (3.841) (5.991)
#
# HO 28.252   ???2.540    18.181
#   (1.645)   (1.645) (1.645)
#
# KW 28.252   ???2.540    21.832
#   (1.645)   (1.645)   (1.645)
#
# SLM 32.661  ???2.433      - 
#   (1.645)   (1.645)     -
#
# GHM   -       -     798.162 
#       -       -     (4.231)


# Table 4.2 [EViews]
########### [(p-values)]
# BP 798.162   6.454  804.615
#   (0.000)   (0.0111) (0.0000)
#
# HO 28.252   ???2.540    18.181
#   (0.000)   (0.9945) (0.0000)
#
# KW 28.252   ???2.540    21.832
#   (0.000)   (0.9945)   (0.0000)
#
# SLM 32.661  ???2.433      - 
#   (0.000)   (0.9925)    -
#
# GHM   -       -     798.162 
#       -       -     (0.0000)


require(plm)
data("Grunfeld")
Grunfeldpdata <- pdata.frame(Grunfeld, index = c("firm", "year"), drop.index = FALSE, row.names = TRUE)

fe_grunfeld  <- plm(inv ~ value + capital, data=Grunfeldpdata, model="within")
re_grunfeld  <- plm(inv ~ value + capital, data=Grunfeldpdata, model="random")
pool_grunfeld  <- plm(inv ~ value + capital, data=Grunfeldpdata, model="pooling")



Grunfeldpdata_unbalanced <- Grunfeld[1:(nrow(Grunfeld)-1), ] 
Grunfeldpdata_unbalanced <- pdata.frame(Grunfeldpdata_unbalanced, index=c("firm"), drop.index = F)
fe_grunfeld_unbalanced <- plm(inv ~ value + capital, data=Grunfeldpdata_unbalanced, model="within")
re_grunfeld_unbalanced <- plm(inv ~ value + capital, data=Grunfeldpdata_unbalanced, model="random")
pool_grunfeld_unbalanced <- plm(inv ~ value + capital, data=Grunfeldpdata_unbalanced, model="pooling")




########## call original plmtest.plm() to get the test from original v1.4-0 ####
# individual
honda_orig            <- plm:::plmtest.plm(pool_grunfeld, type="honda")
honda_orig_unbalanced <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="honda")
bp_orig               <- plm:::plmtest.plm(pool_grunfeld, type="bp")
bp_orig_unbalanced    <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="bp")

# time
honda_orig            <- plm:::plmtest.plm(pool_grunfeld, type="honda", effect="time")
honda_orig_unbalanced <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="honda", effect="time")
bp_orig               <- plm:::plmtest.plm(pool_grunfeld, type="bp", effect="time")
bp_orig_unbalanced    <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="bp", effect="time")


# twoways
honda_orig_tw            <- plm:::plmtest.plm(pool_grunfeld, type="honda", effect="twoways")
honda_orig_unbalanced_tw <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="honda", effect="twoways")
bp_orig_tw               <- plm:::plmtest.plm(pool_grunfeld, type="bp", effect="twoways")
bp_orig_unbalanced_tw    <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="bp", effect="twoways")
kw_orig_tw               <- plm:::plmtest.plm(pool_grunfeld, type="kw", effect="twoways")
kw_orig_unbalanced_tw    <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="kw", effect="twoways")
ghm_orig_tw              <- plm:::plmtest.plm(pool_grunfeld, type="ghm", effect="twoways")
ghm_orig_unbalanced_tw   <- plm:::plmtest.plm(pool_grunfeld_unbalanced, type="ghm", effect="twoways")

#################### load modified version first (eg. source(plm_fixes_v1.4-0_kt.R)) ############
### generalized Version to handle also unbalanced panels
# individual
honda_mod            <- plmtest(pool_grunfeld, type="honda")
honda_mod_unbalanced <- plmtest(pool_grunfeld_unbalanced, type="honda")
bp_mod               <- plmtest(pool_grunfeld, type="bp")
bp_mod_unbalanced    <- plmtest(pool_grunfeld_unbalanced, type="bp")
kw_mod               <- plmtest(pool_grunfeld, type="kw")

# time
honda_mod_time            <- plmtest(pool_grunfeld, type="honda", effect="time")
honda_mod_time_unbalanced <- plmtest(pool_grunfeld_unbalanced, type="honda", effect="time")
bp_mod_time               <- plmtest(pool_grunfeld, type="bp", effect="time")
bp_mod_time_unbalanced    <- plmtest(pool_grunfeld_unbalanced, type="bp", effect="time")
kw_mod_time               <- plmtest(pool_grunfeld, type="kw", effect="time")
# twoways
honda_mod_tw            <- plmtest(pool_grunfeld, type="honda", effect="twoways")
honda_mod_unbalanced_tw <- plmtest(pool_grunfeld_unbalanced, type="honda", effect="twoways")
bp_mod_tw               <- plmtest(pool_grunfeld, type="bp", effect="twoways")
bp_mod_unbalanced_tw    <- plmtest(pool_grunfeld_unbalanced, type="bp", effect="twoways")
kw_mod_tw               <- plmtest(pool_grunfeld, type="kw", effect="twoways")
kw_mod_unbalanced_tw    <- plmtest(pool_grunfeld_unbalanced, type="kw", effect="twoways")
ghm_mod_tw              <- plmtest(pool_grunfeld, type="ghm", effect="twoways")
ghm_mod_unbalanced_tw   <- plmtest(pool_grunfeld_unbalanced, type="ghm", effect="twoways")


# Tests - balanced - should be TRUE
  # individual
  abs(honda_mod$statistic) - honda_orig$statistic < 0.00000000001
  abs(bp_mod$statistic)    - bp_orig$statistic    < 0.00000000001
  
  # twoways
  abs(kw_mod_tw$statistic)  - kw_orig_tw$statistic    < 0.00000000001
  abs(ghm_mod_tw$statistic) - ghm_orig_tw$statistic   < 0.00000000001


# Tests - unbalanced - should be FALSE
  # individual
  abs(abs(honda_mod_unbalanced$statistic) - honda_orig_unbalanced$statistic) < 0.001
  abs(abs(bp_mod_unbalanced$statistic)    - bp_orig_unbalanced$statistic) < 0.001
  
  # twoways
  abs(abs(kw_mod_unbalanced_tw$statistic)  - kw_orig_unbalanced_tw$statistic)    < 0.001
  abs(abs(ghm_mod_unbalanced_tw$statistic) - ghm_orig_unbalanced_tw$statistic)   < 0.001
  
  
  

  
  
  
  
  
  
  


