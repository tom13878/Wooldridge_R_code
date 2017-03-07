# wooldridge (2010) control function example
library(haven)
card <- read_dta("D:/wooldridge_datasets/CARD.DTA")

# test whether the variable educ is endogenous

# 
ols1 <- lm(lwage ~ educ + black:educ + exper + expersq + black + smsa + south +
     reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 +
     smsa66, data=card)

# reduced form 1
ols2 <- lm(educ ~ exper + expersq + black + smsa + south +
     reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 +
     smsa66 + nearc4 + black:nearc4, data = card)

v21 <- ols2$residuals

# make dependent variable for second reduced form model
card$black_educ <- card$black * card$educ

# reduced form 2
ols3 <- lm(black_educ ~ exper + expersq + black + smsa + south +
     reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 +
     smsa66 + nearc4 + black:nearc4, data = card)

v22 <- ols3$residuals

# add residuals to ols1
ols4 <- lm(lwage ~ educ + black:educ + exper + expersq + black + smsa + south +
             reg661 + reg662 + reg663 + reg664 + reg665 + reg666 + reg667 + reg668 +
             smsa66 + v21 + v22, data=card)

# F test to determine whether v21 and v22 are jointly significant
RSS1 <- deviance(ols4) # full model
RSS0 <- deviance(ols1) # reduced model (no residuals)

q <- 2 # number of coefficients set to zero
n <- nrow(card) # number of observations
k <- length(coef(ols4)) # number of regressors in full model

# n - (k + 1) is the residual degrees of freedom

F <- ((RSS0 - RSS1)/q)/((RSS1)/(n-(k+1))) # 30.9

# to calculate the p value use the f quantile
# distribution

pval <- 1 - pf(F, q, n - (k+1)) # 1.73e-06

# The p.value is greater than 0.05 and so we do not
# reject exogeneity of educ and black:educ