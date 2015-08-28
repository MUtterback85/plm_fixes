### Testfile for pwtest()


###  panelmodel interface which did not respect the effect parameter, i. e.
###  for a supplied panelmodel effect="individual" and effect="time" deliver the same result for CRAN version 1.4-0
###  formula interface is not affected



require(plm)
data("Produc", package="plm")
formula <- log(gsp)~log(pcap)+log(pc)+log(emp)+unemp
pwtest(formula, data=Produc)
pwtest(formula, data=Produc, effect="individual")
pwtest(formula, data=Produc, effect="time")

pool_prodc <- plm(formula, data=Produc, model="pooling")
pwtest(pool_prodc)

# for CRAN version 1.4-0, this following two tests results are the same, but should be different (see above use of formula interface)
pwtest(pool_prodc, effect="individual")
pwtest(pool_prodc, effect="time")

