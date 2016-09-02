# -------------------------------------
# R code for the replication of
# analysis in:
# 
# Panel data methods for fractional
# response variables with an
# application to test pass rates
# Leslie E. Papke, Jeffrey M. Wooldridge
# Journal of Econometrics (2008)
# 
# for data see the following URL:
# http://www.econ.msu.edu/faculty/papke/
# -------------------------------------

setwd("C:/Users/Tomas/Documents/Wooldridge_R_code/data/Papke_Wooldridge2008statafiles")

library(haven)
library(plm)

PWdata <- read_dta("meap92_01.dta")

# tables 
mean(PWdata$math4[PWdata$year %in% 1995], na.rm=TRUE)
mean(PWdata$math4[PWdata$year %in% 2001], na.rm=TRUE)
mean(PWdata$rexppp[PWdata$year %in% 1995], na.rm=TRUE)
mean(PWdata$rexppp[PWdata$year %in% 2001], na.rm=TRUE)
mean(PWdata$lunch[PWdata$year %in% 1995], na.rm=TRUE)
mean(PWdata$lunch[PWdata$year %in% 2001], na.rm=TRUE)
mean(PWdata$enroll[PWdata$year %in% 1995], na.rm=TRUE)
mean(PWdata$enroll[PWdata$year %in% 2001], na.rm=TRUE)

# -------------------------------------
# estimations
# -------------------------------------

# next try some estimations.
# 1. straightforward linear model
# 2. A fixed effects "within" estimation with plm package
# 3. Time averages (equivalent to fixed effects)
# all models have year dummies for 1996-2001

# only look at 1995-2001
PWdata2 <- PWdata[PWdata$year >= 1995,]

# linear model
lm0 <- lm(math4 ~ lavgrexp + lunch + lenroll + y01 + y00 + y99 + y98 + y97 + y96, data=PWdata2)

# fixed effects model
fixed0 <- plm(math4 ~ lavgrexp + lunch + lenroll + y01 + y00 + y99 + y98 + y97 + y96,
             data=PWdata2, index=c("distid", "year"), model="within")

# cre model - equivalent in this case
# to a fixed effects model
cre0 <- lm(math4 ~ lavgrexp + lunch + lenroll + y96 + y97 + y98 + y99 + y00 + y01 + alunch + alenroll + alavgrexp, data=PWdata2)

# -------------------------------------
# and now the endogeneity bit
# -------------------------------------

# following line indictes exppp is in real (2001) values
quantile(PWdata$exppp[PWdata$year %in% 1994]/PWdata$cpi[PWdata$year %in% 1994], c(0.1, 0.25, 0.5, 0.75, 0.9), na.rm=TRUE

iv_stage1 <- lm(lavgrexp ~ lfound + lexppp94 + lunch + lenroll + alunch + alenroll +
                  y96 + y97 + y98 + y99 + y00 + y01 +
                  lfndy96 + lfndy97 + lfndy98 + lfndy99 + lfndy00 + lfndy01 +
                  le94y96 + le94y97 + le94y98 + le94y99 + le94y00 + le94y01
                  , data=PWdata2)

v <- residuals(iv_stage1)

iv_stage2 <- lm(math4 ~ lavgrexp + lunch + lenroll +
                  y96 + y97 + y98 + y99 + y00 + y01 +
                  lexppp94 +
                  le94y96 + le94y97 + le94y98 + le94y99 + le94y00 + le94y01 +
                  alunch + alenroll + v, data=PWdata2)

test <- PWdata[PWdata$year %in% 1995, ]
quantile(test$exppp/test$cpi, c(0.1, 0.25, 0.5, 0.75, 0.9))
