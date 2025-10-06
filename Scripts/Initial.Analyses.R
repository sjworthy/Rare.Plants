# Plotting and Analyzing Rare Plant Data 

library(tidyverse)
library(popbio)

# read in data
# data is for 4 species, including year of survey and stem counts

dat = read.csv("./Formatted.Data/Rare.Plant.Data.csv")

# split the data

ARHI = dat %>%
  filter(Species == "Aralia.hispida")
CYRE = dat %>%
  filter(Species == "Cypripedium.reginae")
COCA = dat %>%
  filter(Species == "Cornus.canadensis")
SPLU = dat %>%
  filter(Species == "Spiranthes.lucida")

# make a plot for each species

ggplot(ARHI, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1)+
  geom_point()+
  theme_classic(base_size = 20)+
  #geom_smooth(color = "pink")+
  ylab("Number of Stems")
ggplot(CYRE, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1)+
  geom_point()+
  #geom_smooth(color = "pink")+
  theme_classic(base_size = 20)+
  ylab("Number of Stems")
ggplot(COCA, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1)+
  geom_point()+
  #geom_smooth(color = "pink")+
  theme_classic(base_size = 20)+
  ylab("Number of Stems")
ggplot(SPLU, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1)+
  geom_point()+
  #geom_smooth(color = "pink")+
  theme_classic(base_size = 20)+
  ylab("Number of Stems")

# Need to remove extremely small counts within each population
# (<=10), threshold serves as the minimum for reliable modeling results
# Count-based PVAs assume that the variation
# in population size observed is strictly due to
# environmental stochasticity. However, very small
# populations fluctuate mostly due to demographic,
# not environmental stochasticity (Bernardo et al. 219)

# An issue for ARHI (6 time points) and SPLU (4 time points)

ARHI.2 = ARHI %>%
  filter(Stem.Count > 10)
SPLU.2 = SPLU %>%
  filter(Stem.Count > 10)

#### Calculating log response ratios ####
## calculate  log(Nt+1/Nt)
nt = length(ARHI.2$Stem.Count)
LRR = log(ARHI.2$Stem.Count[-1]/ARHI.2$Stem.Count[-nt])
## transformation for unequal variances
Year.transform <- sqrt(ARHI.2$Year[-1] - ARHI.2$Year[-length(ARHI.2$Year)])
ARHI.3 = as.data.frame(cbind(LRR,Year.transform))
ARHI.3$Corrected.LRR = ARHI.3$LRR/ARHI.3$Year.transfom
ARHI.3$Stems.year.T = ARHI.2$Stem.Count[-nt]

## calculate  log(Nt+1/Nt)
nt = length(COCA$Stem.Count)
LRR = log(COCA$Stem.Count[-1]/COCA$Stem.Count[-nt])
## transformation for unequal variances
Year.transform <- sqrt(COCA$Year[-1] - COCA$Year[-length(COCA$Year)])
COCA.2 = as.data.frame(cbind(LRR,Year.transform))
COCA.2$Corrected.LRR = COCA.2$LRR/COCA.2$Year.transform
COCA.2$Stems.year.T = COCA$Stem.Count[-nt]

## calculate  log(Nt+1/Nt)
nt = length(CYRE$Stem.Count)
LRR = log(CYRE$Stem.Count[-1]/CYRE$Stem.Count[-nt])
## transformation for unequal variances
Year.transform <- sqrt(CYRE$Year[-1] - CYRE$Year[-length(CYRE$Year)])
CYRE.2 = as.data.frame(cbind(LRR,Year.transform))
CYRE.2$Corrected.LRR = CYRE.2$LRR/CYRE.2$Year.transform
CYRE.2$Stems.year.T = CYRE$Stem.Count[-nt]

## calculate  log(Nt+1/Nt)
nt = length(SPLU.2$Stem.Count)
LRR = log(SPLU.2$Stem.Count[-1]/SPLU.2$Stem.Count[-nt])
## transformation for unequal variances
Year.transform <- sqrt(SPLU.2$Year[-1] - SPLU.2$Year[-length(SPLU.2$Year)])
SPLU.3 = as.data.frame(cbind(LRR,Year.transform))
SPLU.3$Corrected.LRR = SPLU.3$LRR/SPLU.3$Year.transform
SPLU.3$Stems.year.T = SPLU.2$Stem.Count[-nt]


#### Estimate mean and variance of growth rates using linear regression ####
# force the intercept to be 0: enforcing the rule that there can be no change
# in population size if no time has elapsed
# mu is the slope coefficient
# variance is the mean square residual in anova table

# a positive mu indicates an environment in which most realizations tend to 
# grow, whereas a negative mu indicates that declining realizations predominate.
# The more the population growth rate (lambda) varies from year to year as a result
# of environmental stochasiticy, the greater will be the value of sigma-squared
# and the greater the range of possible population sizes in the future.

ARHI.mod = lm(Corrected.LRR ~ 0 + Year.transform, data = ARHI.3)
summary(ARHI.mod)
ARHI.mu = coef(ARHI.mod) # mu = -0.02003551 
ARHI.sig2 = anova(ARHI.mod)[["Mean Sq"]][2] # variance = 0.1217645
confint(ARHI.mod,1) # mu CIs = -0.1590133, 0.1189423
df1 = length(ARHI.3$LRR)-1
(df1*anova(ARHI.mod)[["Mean Sq"]][2])/qchisq(c(0.975,0.025), df = df1)
# var CIs = 0.06526691 0.30285769

COCA.mod = lm(Corrected.LRR ~ 0 + Year.transform, data = COCA.2)
summary(COCA.mod)
COCA.mu = coef(COCA.mod) # mu = -0.0500674 
COCA.sig2 = anova(COCA.mod)[["Mean Sq"]][2] # variance = 0.1230193
confint(COCA.mod,1) # CIs = -0.2014644 0.1013296
df1 = length(COCA.2$LRR)-1
(df1*anova(COCA.mod)[["Mean Sq"]][2])/qchisq(c(0.975,0.025), df = df1)
# var CIs = 0.06173402 0.35463891

CYRE.mod = lm(Corrected.LRR ~ 0 + Year.transform, data = CYRE.2)
summary(CYRE.mod)
CYRE.mu = coef(CYRE.mod) # mu = -0.01633807  
CYRE.sig2 = anova(CYRE.mod)[["Mean Sq"]][2] # variance = 0.1234609
confint(CYRE.mod,1) # CIs = -0.1274633 0.09478714
df1 = length(CYRE.2$LRR)-1
(df1*anova(CYRE.mod)[["Mean Sq"]][2])/qchisq(c(0.975,0.025), df = df1)
# var CIs = 0.07384683 0.24731923

SPLU.mod = lm(Corrected.LRR ~ 0 + Year.transform, data = SPLU.3)
summary(SPLU.mod)
SPLU.mu = coef(SPLU.mod) # mu = 0.06933172  
SPLU.sig2 = anova(SPLU.mod)[["Mean Sq"]][2] # variance = 0.2846666
confint(SPLU.mod,1) # CIs = -0.2212901 0.3599535
df1 = length(SPLU.3$LRR)-1
(df1*anova(SPLU.mod)[["Mean Sq"]][2])/qchisq(c(0.975,0.025), df = df1)
# var CIs = 0.1463791 0.7756956

#### Plotting LRRs ####

ggplot(ARHI.3, aes(x = Year.transform, y = Corrected.LRR))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, color = "blue")+
  xlab(expression(paste("Sqrt time between censuses ", (t[t + 1] - t[i])^{
    1 / 2})))+
  ylab(expression(log(N[t + 1] / N[t]) / (t[t + 1] - t[i])^{
    1 / 2}))+
  theme_classic(base_size = 20)

ggplot(COCA.2, aes(x = Year.transform, y = Corrected.LRR))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, color = "blue")+
  xlab(expression(paste("Sqrt time between censuses ", (t[t + 1] - t[i])^{
    1 / 2})))+
  ylab(expression(log(N[t + 1] / N[t]) / (t[t + 1] - t[i])^{
    1 / 2}))+
  theme_classic(base_size = 20)

ggplot(CYRE.2, aes(x = Year.transform, y = Corrected.LRR))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, color = "blue")+
  xlab(expression(paste("Sqrt time between censuses ", (t[t + 1] - t[i])^{
    1 / 2})))+
  ylab(expression(log(N[t + 1] / N[t]) / (t[t + 1] - t[i])^{
    1 / 2}))+
  theme_classic(base_size = 20)

ggplot(SPLU.3, aes(x = Year.transform, y = Corrected.LRR))+
  geom_point()+
  geom_hline(yintercept = 0)+
  geom_smooth(method = "lm", formula = y ~ 0 + x, se = FALSE, color = "blue")+
  xlab(expression(paste("Sqrt time between censuses ", (t[t + 1] - t[i])^{
    1 / 2})))+
  ylab(expression(log(N[t + 1] / N[t]) / (t[t + 1] - t[i])^{
    1 / 2}))+
  theme_classic(base_size = 20)

#### test for temporal autocorrelation ####
install.packages("lmtest")
library(lmtest)
?dwtest

# Using Durbin-Watson test to measure the strength of autocorrelation in the 
# regression residuals. Tests the assumptions of the diffusion approximation
# that the environmental conditions are uncorrelated from one inter-census interval
# to the next. That is whether a particular interval was good or bad for birth
# or death is independent of wheter preceding or succeeding intervals were
# good or bad

dwtest(ARHI.mod, alternative = "two.sided") # p = 0.027, DW = 0.95
dwtest(CYRE.mod, alternative = "two.sided") # p = 0.019, DW = 2.89
dwtest(COCA.mod, alternative = "two.sided") # p = 0.133, DW = 1.18
dwtest(SPLU.mod, alternative = "two.sided") # p = 0.281, DW = 1.43

# test for autocorrelation using the residuals the first-order autocorrelation
# of the residuals

# Get residuals
res <- resid(ARHI.mod)
# Create lagged residuals (drop first observation for lag)
res_t   <- res[-1]         # residuals t = 2, ..., n
res_t_1 <- res[-length(res)]  # residuals t-1 = 1, ..., n-1
# Compute Pearson correlation
cor.test(res_t, res_t_1) # r = 0.42, p = 0.13

# Get residuals
res <- resid(CYRE.mod)
# Create lagged residuals (drop first observation for lag)
res_t   <- res[-1]         # residuals t = 2, ..., n
res_t_1 <- res[-length(res)]  # residuals t-1 = 1, ..., n-1
# Compute Pearson correlation
cor.test(res_t, res_t_1) # r = -0.48, p = 0.02

# Get residuals
res <- resid(COCA.mod)
# Create lagged residuals (drop first observation for lag)
res_t   <- res[-1]         # residuals t = 2, ..., n
res_t_1 <- res[-length(res)]  # residuals t-1 = 1, ..., n-1
# Compute Pearson correlation
cor.test(res_t, res_t_1) # r = 0.11, p = 0.74

# Get residuals
res <- resid(SPLU.mod)
# Create lagged residuals (drop first observation for lag)
res_t   <- res[-1]         # residuals t = 2, ..., n
res_t_1 <- res[-length(res)]  # residuals t-1 = 1, ..., n-1
# Compute Pearson correlation
cor.test(res_t, res_t_1) # r = 0.23, p = 0.47

#### test for outliers ####
# rstudent gives the studentized residual for each (x,y) pair in the regression
# data points with studentized residuals greater than 2 are suspected outliers
# studentized residuals greater than 3 is often used for small datasets

# dffits is a statistic that measures the influence each data point has on the
# regression parameter estimates. For a linear regression with no intercept,
# a value of Dffits greater than 2*sqrt(1/q) where q is the number of data points,
# or transitions, suggests high influence.

influence.measures(ARHI.mod)
influence.measures(CYRE.mod)
influence.measures(COCA.mod)
influence.measures(SPLU.mod)

# test for outlier in y
rstudent(ARHI.mod) # one > 2, 13, none greater than 3
rstudent(COCA.mod) # one > 2, 12, yes it's greater than 3
rstudent(CYRE.mod) # three > 2, 2,9,10, none greater than 3
rstudent(SPLU.mod) # one > 2, 8, none greater than 3

# test for impact on model prediction
dffits(ARHI.mod)[dffits(ARHI.mod) > 2*sqrt(1/15)] # none
dffits(CYRE.mod)[dffits(CYRE.mod) > 2*sqrt(1/23)] # 2
dffits(COCA.mod)[dffits(COCA.mod) > 2*sqrt(1/12)] # 12
dffits(SPLU.mod)[dffits(SPLU.mod) > 2*sqrt(1/13)] # none

#### REDO models removing outliers ####
# remove transition 12 for COCA
COCA.3 = COCA.2[c(1:11),]

# rerun models
COCA.mod = lm(Corrected.LRR ~ 0 + Year.transform, data = COCA.3)
summary(COCA.mod)
COCA.mu.2 = coef(COCA.mod) # mu = -0.08939685
COCA.sig2.2 = anova(COCA.mod)[["Mean Sq"]][2] # variance = 0.03477893
confint(COCA.mod,1) # CIs = -0.1725025 -0.006291231
df1 = length(COCA.2$LRR)-1
(df1*anova(COCA.mod)[["Mean Sq"]][2])/qchisq(c(0.975,0.025), df = df1)
# var CIs = 0.01745289 0.10026035

# remove transition 2 for CYRE
CYRE.3 = CYRE.2[c(1,3:23),]

CYRE.mod = lm(Corrected.LRR ~ 0 + Year.transform, data = CYRE.3)
summary(CYRE.mod)
CYRE.mu.2 = coef(CYRE.mod) # mu = -0.04095235  
CYRE.sig2.2 = anova(CYRE.mod)[["Mean Sq"]][2] # variance = 0.1039081
confint(CYRE.mod,1) # CIs = -0.1456448 0.06374014
df1 = length(CYRE.2$LRR)-1
(df1*anova(CYRE.mod)[["Mean Sq"]][2])/qchisq(c(0.975,0.025), df = df1)
# var CIs = 0.06215157 0.20815080

#### extinction time cumulative distribution function ####

ARHI.cdf = extCDF(ARHI.mu,ARHI.sig2,Nc = 33, Ne = 10, tmax = 50)
ARHI.cdf = as.data.frame(ARHI.cdf)
ARHI.cdf$Years = 1:50

ggplot(ARHI.cdf, aes(y = ARHI.cdf, x = Years))+
  geom_line()+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

n <- seq(20, 100, 2)
exts <- numeric(length(n))
for (i in 1:length(n) ){
  ex <- extCDF(ARHI.mu, ARHI.sig2, Nc=n[i], Ne=20)
  exts[i] <- ex[50]
}

ARHI.countCDF = countCDFxt(mu = ARHI.mu, sig2 = ARHI.sig2, nt = 15,
                           Nc = 33, Ne = 10, tq = 29,
                           tmax = 50, Nboot = 500, plot = TRUE)

ggplot(ARHI.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

plot(n, exts, type='l', las=1,
     xlab="Current population size",
     ylab="Probability of quasi-extinction by year 50")

CYRE.cdf = extCDF(CYRE.mu,CYRE.sig2,Nc = 53, Ne = 10, tmax = 50)
CYRE.cdf = as.data.frame(CYRE.cdf)
CYRE.cdf$Years = 1:50

ggplot(CYRE.cdf, aes(y = CYRE.cdf, x = Years))+
  geom_line()+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

n <- seq(20, 100, 2)
exts <- numeric(length(n))
for (i in 1:length(n) ){
  ex <- extCDF(CYRE.mu, CYRE.sig2, Nc=n[i], Ne=20)
  exts[i] <- ex[50]
}

plot(n, exts, type='l', las=1,
     xlab="Current population size",
     ylab="Probability of quasi-extinction by year 50")

CYRE.countCDF = countCDFxt(mu = CYRE.mu, sig2 = CYRE.sig2, nt = 23,
                           Nc = 53, Ne = 10, tq = 43,
                           tmax = 50, Nboot = 500, plot = TRUE)

ggplot(CYRE.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

COCA.cdf = extCDF(COCA.mu,COCA.sig2,Nc = 478, Ne = 10, tmax = 50)
COCA.cdf = as.data.frame(COCA.cdf)
COCA.cdf$Years = 1:50

ggplot(COCA.cdf, aes(y = COCA.cdf, x = Years))+
  geom_line()+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

n <- seq(20, 100, 2)
exts <- numeric(length(n))
for (i in 1:length(n) ){
  ex <- extCDF(COCA.mu, COCA.sig2, Nc=n[i], Ne=20)
  exts[i] <- ex[50]
}

plot(n, exts, type='l', las=1,
     xlab="Current population size",
     ylab="Probability of quasi-extinction by year 50")

COCA.countCDF = countCDFxt(mu = COCA.mu, sig2 = COCA.sig2, nt = 12,
                           Nc = 478, Ne = 10, tq = 26,
                           tmax = 50, Nboot = 500, plot = TRUE)

ggplot(COCA.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

SPLU.cdf = extCDF(SPLU.mu,SPLU.sig2,Nc = 94, Ne = 10, tmax = 50)
SPLU.cdf = as.data.frame(SPLU.cdf)
SPLU.cdf$Years = 1:50

ggplot(SPLU.cdf, aes(y = SPLU.cdf, x = Years))+
  geom_line()+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

n <- seq(20, 100, 2)
exts <- numeric(length(n))
for (i in 1:length(n) ){
  ex <- extCDF(SPLU.mu, SPLU.sig2, Nc=n[i], Ne=20)
  exts[i] <- ex[50]
}

plot(n, exts, type='l', las=1,
     xlab="Current population size",
     ylab="Probability of quasi-extinction by year 50")

SPLU.countCDF = countCDFxt(mu = SPLU.mu, sig2 = SPLU.sig2, nt = 13,
                           Nc = 94, Ne = 10, tq = 29,
                           tmax = 50, Nboot = 500, plot = TRUE)

ggplot(SPLU.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

CYRE.cdf.2 = extCDF(CYRE.mu.2,CYRE.sig2.2,Nc = 53, Ne = 10, tmax = 50)
CYRE.cdf.2 = as.data.frame(CYRE.cdf.2)
CYRE.cdf.2$Years = 1:50

ggplot(CYRE.cdf.2, aes(y = CYRE.cdf.2, x = Years))+
  geom_line()+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

CYRE.countCDF = countCDFxt(mu = CYRE.mu.2, sig2 = CYRE.sig2.2, nt = 22,
                           Nc = 53, Ne = 10, tq = 43,
                           tmax = 50, Nboot = 500, plot = TRUE)

ggplot(CYRE.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")


COCA.cdf.2 = extCDF(COCA.mu.2,COCA.sig2.2,Nc = 188, Ne = 10, tmax = 50)

COCA.countCDF = countCDFxt(mu = COCA.mu.2, sig2 = COCA.sig2.2, nt = 15,
                           Nc = 188, Ne = 10, tq = 11,
                           tmax = 50, Nboot = 500, plot = TRUE)

ggplot(COCA.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")


#### Plots to investigate density dependence ####
ggplot(ARHI.3, aes(x = Stems.year.T, y = LRR))+
  geom_point()+
  xlab("Number of Stems in Year T")+
  ylab("Log Response Ratio (LRR)")+
  theme_classic(base_size = 20)

ggplot(COCA.2, aes(x = Stems.year.T, y = LRR))+
  geom_point()+
  xlab("Number of Stems in Year T")+
  ylab("Log Response Ratio (LRR)")+
  theme_classic(base_size = 20)

ggplot(CYRE.2, aes(x = Stems.year.T, y = LRR))+
  geom_point()+
  xlab("Number of Stems in Year T")+
  ylab("Log Response Ratio (LRR)")+
  theme_classic(base_size = 20)

ggplot(SPLU.3, aes(x = Stems.year.T, y = LRR))+
  geom_point()+
  xlab("Number of Stems in Year T")+
  ylab("Log Response Ratio (LRR)")+
  theme_classic(base_size = 20)

#### Testing for density dependence ####
# some papers do a straight regression
# others log transform initial population size

cor.test(ARHI.3$LRR,ARHI.3$Stems.year.T) # r = -0.12, p = 0.67
ARHI.DD.mod = lm(LRR ~ log(Stems.year.T), data = ARHI.3)
summary(ARHI.DD.mod) # not significant
check_normality(ARHI.DD.mod) # not normal
check_heteroskedasticity(ARHI.DD.mod)

cor.test(COCA.2$LRR,COCA.2$Stems.year.T) # r = -0.27, p = 0.40
COCA.DD.mod = lm(LRR ~ log(Stems.year.T), data = COCA.2)
summary(COCA.DD.mod) # not significant
check_normality(COCA.DD.mod) # not normal
check_heteroskedasticity(COCA.DD.mod)

cor.test(COCA.3$LRR,COCA.3$Stems.year.T) # r = 0.06, p = 0.87
COCA.DD.mod = lm(LRR ~ log(Stems.year.T), data = COCA.3)
summary(COCA.DD.mod) # not significant
check_normality(COCA.DD.mod) # not normal
check_heteroskedasticity(COCA.DD.mod)

cor.test(CYRE.2$LRR,CYRE.2$Stems.year.T) # r = -0.48, p = 0.02
CYRE.DD.mod = lm(LRR ~ log(Stems.year.T), data = CYRE.2)
summary(CYRE.DD.mod) # not significant
check_normality(CYRE.DD.mod)
check_heteroskedasticity(CYRE.DD.mod)

cor.test(CYRE.3$LRR,CYRE.3$Stems.year.T) # r = -0.45, p = 0.04
CYRE.DD.mod = lm(LRR ~ log(Stems.year.T), data = CYRE.3)
summary(CYRE.DD.mod) # not significant
check_normality(CYRE.DD.mod)
check_heteroskedasticity(CYRE.DD.mod)

cor.test(SPLU.3$LRR,SPLU.3$Stems.year.T) # r = -0.40, p = 0.17
SPLU.DD.mod = lm(LRR ~ log(Stems.year.T), data = SPLU.3)
summary(SPLU.DD.mod) # not significant
check_normality(SPLU.DD.mod)
check_heteroskedasticity(SPLU.DD.mod)


