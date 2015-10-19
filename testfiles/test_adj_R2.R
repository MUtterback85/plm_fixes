## Comparision of lm()'s adjusted R-squared and plm()'s R-squared
## for pooling-Models!


require(plm)
data("Grunfeld")

### lm
a_lm <- lm(inv ~ value + capital, Grunfeld)
r2_lm     <- summary(a_lm)$r.squared
r2_adj_lm <- summary(a_lm)$adj.r.squared

### plm
a_plm <- plm(inv ~ value + capital , data=Grunfeld, model = "pooling")
r2_plm     <- summary(a_plm)$r.squared[1] # R^2
r2_adj_plm <- summary(a_plm)$r.squared[2] # adj. R^2

# current formula in plm::r.squared()
# R2 * df.residual(object) / length(resid(object) - 1)
# should be (???): R2 * df.residual(object) / (length(resid(object)) - 1)
r2_plm * df.residual(a_plm) / length(resid(a_plm)-1)
r2_plm * df.residual(a_plm) / (length(resid(a_plm))-1)

# to meet summary.lm() for pooling models - no checking for intercept
names(r2_plm) <- NULL
1-(1-r2_plm) * (length(resid(a_plm))-1) / df.residual(a_plm)


##### without intercept - regular R-squared diverges
## lm - see summary.lm's source how the checking for the presence of an intercept is done
a2_lm <- lm(inv ~ value + capital -1 , Grunfeld)
r2_lm2     <- summary(a2_lm)$r.squared
r2_adj_lm2 <- summary(a2_lm)$adj.r.squared

## plm
a2_plm <- plm(inv ~ value + capital -1, data=Grunfeld, model = "pooling")
r2_plm2     <- summary(a2_plm)$r.squared[1] # R^2
r2_adj_plm2 <- summary(a2_plm)$r.squared[2] # adj. R^2

