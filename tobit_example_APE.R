# -------------------------------------
# example using R to calculate the APE
# of a tobit model, from wooldridge
# -------------------------------------

library(AER)
library(haven)

setwd("c:/users/tomas/documents/Wooldridge_R_code/data")

# read in mroz data set
mroz <- read_dta("mroz.dta")

# run the tobit model
modl.tobit <- tobit(hours ~ nwifeinc + educ + exper + expersq + age + kidslt6 +
                      kidsge6, data=mroz)

# calcualte the APEs
sigma <- modl.tobit$scale

vars <- c("nwifeinc", "educ", "exper", "expersq", "age", "kidslt6", "kidsge6")
x <- as.matrix(cbind(1, mroz[, vars ]))
n <- nrow(mroz)
(1/n) * sum(pnorm(((x %*% coef(modl.tobit)))/sigma, 0, 1))

# -------------------------------------
# example of a heckit model from
# wooldridge pp. 807
# -------------------------------------

# first OLS
# note that woman with no wage have . change to NA

summary(lm(lwage ~ educ + exper + expersq, data=mroz))

# heckit model

modl.particip <- glm(inlf ~ nwifeinc + educ + exper + expersq + age + kidslt6 + kidsge6,
                     family=binomial(link="probit"), data=mroz)

# calcualte the mills ratios

vars <- c("nwifeinc", "educ", "exper", "expersq", "age", "kidslt6", "kidsge6")
x <- as.matrix(cbind(1, mroz[, vars ]))
c <- x %*% coef(modl.particip)
mills <- dnorm(c)/pnorm(c)

# do an OLS including the mills ration in regression

summary(lm(lwage ~ educ + exper + expersq + mills, data=mroz))
