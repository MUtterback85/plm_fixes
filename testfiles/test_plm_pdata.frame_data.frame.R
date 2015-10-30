### (1) Comparision of extraction of data.frame and pdata.frame and (2) class 'pseries' of estimated_model$model
### applys to plm v1.4-0 (CRAN as of 2015-10-28)



### (1) Comparision of extraction of data.frame and pdata.frame ###
# from ?pdata.frame: "The "[" method behaves as for data.frame, except that the extraction is also applied to the index attribute."
data(Grunfeld, package="plm")
class(Grunfeld)
pGrunfeld <- pdata.frame(Grunfeld, index = c("firm", "year"), drop.index = F)
class(pGrunfeld)

nrow(Grunfeld[Grunfeld$inv == 317.60, ])    # 1 row and ...
class(Grunfeld[Grunfeld$inv == 317.60, ])   # ... it is a data.frame

class(pGrunfeld[pGrunfeld$inv == 317.60, ]) # this is classes 'pseries' and 'list'
nrow(pGrunfeld[pGrunfeld$inv == 317.60, ])  # then, operations like this give unbehaved answers



### (2) class 'pseries' of estimated_model$model [cosmetic] ###
mod <- plm(inv ~ value + capital, data=pGrunfeld, model = "pooling")
class(mod$model)
class(mod$model$inv) # 'pseries' appears twice
df <- as.data.frame(mod$model)
class(df)
class(df$inv) # 'pseries' is still here twice


