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

# table 1
# table 2
# table 3
table3 <- PWdata %>%
  filter(year %in% c(1995, 2001)) %>%
  group_by(year) %>%
  summarise(math4_mean=mean(math4), math4_sd=sd(math4),
            rexppp_mean=mean(rexppp), rexppp_sd=sd(rexppp),
            found_mean=mean(found/cpi), found_sd=sd(found/cpi),
            lunch_mean=mean(lunch), lunch_sd=sd(lunch),
            enroll_mean=mean(enroll), enroll_sd=sd(enroll),
            n=n())

# -------------------------------------
# estimations
# -------------------------------------

# next try some estimations.
# 1. straightforward linear model
# 2. A fixed effects "within" estimation with plm package
# 3. Time averages (equivalent to fixed effects)
# all models have year dummies for 1996-2001

# linear model
lm0 <- lm(math4 ~ lavgrexp + lunch + lenroll +
            y01 + y00 + y99 + y98 + y97 + y96,
          data=PWdata[PWdata$year > 1994,])

# fixed effects model
fixed0 <- plm(math4 ~ lavgrexp + lunch + lenroll +
                y01 + y00 + y99 + y98 + y97 + y96,
              data=PWdata[PWdata$year > 1994,], index=c("distid", "year"), model="within")

# cre model - equivalent in this case
# to a fixed effects model
cre0 <- lm(math4 ~ lavgrexp + lunch + lenroll +
             y96 + y97 + y98 + y99 + y00 + y01 +
             alunch + alenroll + alavgrexp, data=PWdata[PWdata$year > 1994,])

# -------------------------------------
# and now the endogeneity bit
# -------------------------------------

iv_stage1 <- lm(lavgrexp ~ lfound + lexppp94 + lunch + lenroll + alunch + alenroll +
                  y96 + y97 + y98 + y99 + y00 + y01 +
                  lfndy96 + lfndy97 + lfndy98 + lfndy99 + lfndy00 + lfndy01 +
                  le94y96 + le94y97 + le94y98 + le94y99 + le94y00 + le94y01
                  , data=PWdata[PWdata$year > 1994,])

v2hat <- residuals(iv_stage1)

iv_stage2 <- lm(math4 ~ lavgrexp + lunch + lenroll +
                  y96 + y97 + y98 + y99 + y00 + y01 +
                  lexppp94 +
                  le94y96 + le94y97 + le94y98 + le94y99 + le94y00 + le94y01 +
                  alunch + alenroll + v2hat, data=PWdata[PWdata$year > 1994,])

# SEs should be fully robust - we can 
# try this using the sandwich or
# multiwayvcov packages

library(sandwich)
cluster <- PWdata[PWdata$year > 1994,]$distid
X <- model.matrix(iv_stage2)
e <- residuals(iv_stage2)
uj <- apply(e*X, 2, function(x) tapply(x, cluster, sum))
uj  <- apply(estfun(iv_stage2),2, function(x) tapply(x, cluster, sum))

# small sample correction
M <- length(unique(cluster))
N <- length(cluster)
K <- iv_stage2$rank
dfc <- (M/(M-1))*((N-1)/(N-K))

# 
# vcovCL_sandwich <- dfc*sandwich(iv_stage2, meat=crossprod(uj)/N)
vcovCL_sandwich <- sandwich(iv_stage2, meat=crossprod(uj)/N)
VcovCL <- vcovCL_sandwich

# results are close but not exactly
# 
library(lmtest)
coeftest(iv_stage2, VcovCL)


# bootstrap the SE for the probit
# and make them robust. 

# first bootstrap the SE
library(multiwayvcov)
library(lmtest)

vcov_distid <- cluster.vcov(iv_stage2, PWdata[PWdata$year > 1994,]$distid) 
iv_stage2_SE <- coeftest(iv_stage2, vcov_distid)

PWdata2 <- PWdata[PWdata$year > 1994,]

# fractional probit QMLE - essentially
# just a probit model
fpr <- glm(math4 ~ lavgrexp + v2hat + lunch + alunch + lenroll + alenroll +
      y96 + y97 + y98 + y99 + y00 + y01 +
      lexppp94 + le94y96 + le94y97 + le94y98 + le94y99 + le94y00 + le94y01,
    data=PWdata2, family=quasibinomial(link=probit))


library(multiwayvcov)
library(frmpd)
X <- model.matrix(fpr)
frm(y=PWdata2$math4,x=X, linkfrac="probit")

cluster.boot(fpr, )

# work out the APE
sigma <- fpr$


form <- formula(fpr)

# make a function for bootstrapping
refit <- function(data, i){
  coef(glm(form, data=data[i,], family=binomial(link=probit)))
}

# boot function
library(boot)
boot <- boot(PWdata[PWdata$year > 1994,], refit, R=500)
SE <- apply(boot$t, 2, sd)

# but we also need the SEs to be robust
# so we need to take account of the
# clustering as well, The idea is to
# draw from the 501 districts. But again
# replace = TRUE
# the idea is then to boostrap the 501
# districts by using 500 bootstrap
# replications

# for each bootstrap sample there should be
# a two step process apparently. According to
# http://ftp.iza.org/dp6887.pdf

PWdata2 <- PWdata[PWdata$year > 1994,]

# start of cluster part
cluster <- PWdata2$distid
clusters <- unique(PWdata2$distid)
reps <- 500
boot_out <- matrix(nrow=reps, ncol=length(fpr$coefficients))

# first run the first stage regression
# in order to get the v2hat
# then do the second stage regression.
form1 <- formula(iv_stage1)# formula for stage 1

form2 <- formula(fpr)

for (i in 1:reps){
  index <- sample(1:length(clusters), length(clusters), replace=TRUE)
  clusters_sampled <- clusters[index]
  cluster_count <- table(clusters_sampled)
  boot_data <- NULL
  for(j in 1:max(cluster_count)) { # check for those clusters sampled 1, 2, 3, 4, 5, ....
    x <- PWdata2[cluster %in% names(cluster_count[cluster_count %in% j]),]
    for(k in 1:j){
      boot_data <- rbind(boot_data, x)
    } 
  }
  iv_stage <- lm(form1, data=boot_data)
  boot_data$v2hat <- residuals(iv_stage)
  boot_out[i,] <- coef(glm(form2, data=boot_data, family=binomial(link=probit)))
}

# stata has an inflation factor thing 
M <- length(unique(cluster))
N <- length(cluster)
K <- fpr$rank
dfc <- (M/(M-1))*((N-1)/(N-K))

SE_boot_cluster <- apply(boot_out, 2, sd)
SE_boot_cluster_stata <- dfc*SE_boot_cluster

# standard errors from the paper for key variables
# larexppp = .759
# lunch .202
# lenroll = .209
# v2 = .811

SE_boot_cluster[c(2, 4, 6, 3)]
SE_boot_cluster_stata[c(2, 4, 6, 3)]

# 0.7392311 0.2273359 0.2176823 0.8099235
# 

# 0.7447456 0.2290317 0.2193062 0.8159653
#

# work ou the t values
b <- coefficients(fpr)
df <- N - length(b)
t <- b/SE_boot_cluster
p.value = 2*pt(abs(t), df=df, lower=FALSE)
