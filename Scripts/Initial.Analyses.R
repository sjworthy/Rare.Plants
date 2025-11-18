# Plotting and Analyzing Rare Plant Data 

library(tidyverse)
library(popbio)
library(wesanderson) # for pop colors
library(cowplot)
library(lmtest)
library(measurements)
library(lmerTest)
library(datawizard)
library(performance)
library(MASS)

wes_palette("AsteroidCity1")
wes_palette("Rushmore1")[3]
"#0A9F9D"
"#E54E21"
"#6C8645"
"#35274A"

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

ARHI.plot = ggplot(ARHI, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1, col = "#0A9F9D")+
  geom_point()+
  theme_classic(base_size = 20)+
  ylab("Number of Stems")
CYRE.plot = ggplot(CYRE, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1, col = "#6C8645")+
  geom_point()+
  theme_classic(base_size = 20)+
  ylab("Number of Stems")
COCA.plot = ggplot(COCA, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1, col = "#E54E21")+
  geom_point()+
  theme_classic(base_size = 20)+
  ylab("Number of Stems")
SPLU.plot = ggplot(SPLU, aes(x = Year, y = Stem.Count))+
  geom_line(linewidth = 1, col = "#35274A")+
  geom_point()+
  theme_classic(base_size = 20)+
  ylab("Number of Stems")

stem.count.plots = plot_grid(ARHI.plot,COCA.plot,CYRE.plot,SPLU.plot, labels = c("A.","B.","C.","D."))

#ggsave(stem.count.plots, file = "./Plots/stem.counts.plot.pdf", height = 10, width = 10)

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
ARHI.3$Corrected.LRR = ARHI.3$LRR/ARHI.3$Year.transform
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
# Using Durbin-Watson test to measure the strength of autocorrelation in the 
# regression residuals. Tests the assumptions of the diffusion approximation
# that the environmental conditions are uncorrelated from one inter-census interval
# to the next. That is whether a particular interval was good or bad for birth
# or death is independent of whether preceding or succeeding intervals were
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

# testing for autocorrelation again
dwtest(COCA.mod, alternative = "two.sided") # p = 0.2832, DW = 2.57


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

# testing for autocorrelation again
dwtest(CYRE.mod, alternative = "two.sided") # p = 0.03, DW = 2.86


#### Finite population growth rate ####
# not sure if this is correct
exp(ARHI.mu) # 0.9801639
exp(CYRE.mu) # 0.9837947 
exp(CYRE.mu.2) # 0.9598749
exp(COCA.mu) # 0.9511653 
exp(COCA.mu.2) # 0.9144826 
exp(SPLU.mu) # 1.071792 

#### Continuous rate of increase and average finite rate of increase ####

# continuous rate of increase (rbar), this is per unit time
# if r = 0.02 this mean 2% increase per year
ARHI.r = ARHI.mu + ARHI.sig2/2 # 0.04084674 
ARHI.r.low = ARHI.r + qnorm(0.025)*sqrt(ARHI.sig2 * ((1/29) + (ARHI.sig2 / (2 * (15 - 1)))))
ARHI.r.high = ARHI.r - qnorm(0.025)*sqrt(ARHI.sig2 * ((1/29) + (ARHI.sig2 / (2 * (15 - 1)))))
CYRE.r = CYRE.mu + CYRE.sig2/2 # 0.04539237 
CYRE.r.low = CYRE.r + qnorm(0.025)*sqrt(CYRE.sig2 * ((1/43) + (CYRE.sig2 / (2 * (23 - 1)))))
CYRE.r.high = CYRE.r - qnorm(0.025)*sqrt(CYRE.sig2 * ((1/43) + (CYRE.sig2 / (2 * (23 - 1)))))
CYRE.r.2 = CYRE.mu.2 + CYRE.sig2.2/2 # 0.01100171 
CYRE.r.low.2 = CYRE.r.2 + qnorm(0.025)*sqrt(CYRE.sig2.2 * ((1/43) + (CYRE.sig2.2 / (2 * (22 - 1)))))
CYRE.r.high.2 = CYRE.r.2 - qnorm(0.025)*sqrt(CYRE.sig2.2 * ((1/43) + (CYRE.sig2.2 / (2 * (22 - 1)))))
COCA.r = COCA.mu + COCA.sig2/2 # 0.01144227 
COCA.r.low = COCA.r + qnorm(0.025)*sqrt(COCA.sig2 * ((1/26) + (COCA.sig2 / (2 * (12 - 1)))))
COCA.r.high = COCA.r - qnorm(0.025)*sqrt(COCA.sig2 * ((1/26) + (COCA.sig2 / (2 * (12 - 1)))))
COCA.r.2 = COCA.mu.2 + COCA.sig2.2/2 # -0.07200738
COCA.r.low.2 = COCA.r.2 + qnorm(0.025)*sqrt(COCA.sig2.2 * ((1/25) + (COCA.sig2.2 / (2 * (11 - 1)))))
COCA.r.high.2 = COCA.r.2 - qnorm(0.025)*sqrt(COCA.sig2.2 * ((1/25) + (COCA.sig2.2 / (2 * (11 - 1)))))
SPLU.r = SPLU.mu + SPLU.sig2/2 # 0.211665 
SPLU.r.low = SPLU.r + qnorm(0.025)*sqrt(SPLU.sig2 * ((1/16) + (SPLU.sig2 / (2 * (13 - 1)))))
SPLU.r.high = SPLU.r - qnorm(0.025)*sqrt(SPLU.sig2 * ((1/16) + (SPLU.sig2 / (2 * (13 - 1)))))


# Average finite rate of increase (lambdabar), discrete time
# average population growth rate

# lambda 1.02 = 2% increase per time step
exp(ARHI.r) # 1.041692
exP(ARHI.r.low)
exp(ARHI.r.high)
exp(CYRE.r) # 1.046438
exp(CYRE.r.low)
exp(CYRE.r.high)
exp(CYRE.r.2) # 1.011062 
exp(CYRE.r.low.2)
exp(CYRE.r.high.2)
exp(COCA.r) # 1.011508 
exp(COCA.r.low)
exp(COCA.r.high)
exp(COCA.r.2) # 0.930524 
exp(COCA.r.low.2)
exp(COCA.r.high.2)
exp(SPLU.r) # 1.235734 
exp(SPLU.r.low)
exp(SPLU.r.high)

#### probability of reaching extinction ####
# mu (mean) is positive for SPLU so need to calculate this probability
# probability is 1 for populations where mean is negative

SPLU.new.sig = (13-1)*SPLU.sig2/13
SPLU.prob = (10/94)^(2*SPLU.mu/SPLU.new.sig)

# interpretation of CDF when mean is positive:
# median time to extinction from the CDF is 105 years. This does NOT mean that half 
# of all realizations will have reached the extinction threshold by 105 years, but
# instead that half of all realizations that will eventually hit the threshold will
# have done so by 105 years. 
# can still calculate the total probability that the population has gone extinct,
# accounting for ALL possible realizations, if we calculate both the probability that
# the extinction threshold is eventually reaches and the conditional extinction time CDF.

# probability that extinction will occur multiplied by the conditional probability that
# extinction will have occurred by 100 years given that it will occur eventually

# formula is probability of extinction multiplied by CDF value at 20 years
SPLU.prob*1.851535e-01 # 0.05675601 

#### extinction time cumulative distribution function ####

ARHI.cdf = extCDF(ARHI.mu,ARHI.sig2,Nc = 33, Ne = 10, tmax = 50)
ARHI.cdf = as.data.frame(ARHI.cdf)
ARHI.cdf$Years = 1:50

x_at_half <- ARHI.cdf %>%
  filter(abs(ARHI.cdf - 0.5) == min(abs(ARHI.cdf - 0.5))) %>%
  pull(Years)

ggplot(ARHI.cdf, aes(y = ARHI.cdf, x = Years))+
  geom_line()+
  geom_segment(aes(x = 0, y = 0.5, xend = x_at_half, yend = 0.5), size = 0.8) +
  geom_segment(aes(x = x_at_half, y = 0, xend = x_at_half, yend = 0.5), size = 0.8) +
  ylim(0.0,1.0)+
  theme_classic(base_size = 20)+
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
                           tmax = 50, Nboot = 10000, plot = TRUE)

ggplot(ARHI.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

min(ARHI.countCDF$Gbest[ARHI.countCDF$Gbest != 0
                        ])

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
  theme_classic(base_size = 20)+
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
                           Nc = 94, Ne = 10, tq = 16,
                           tmax = 50, Nboot = 10000, plot = TRUE)

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

x_at_half <- CYRE.cdf.2 %>%
  filter(abs(CYRE.cdf.2 - 0.5) == min(abs(CYRE.cdf.2 - 0.5))) %>%
  pull(Years)

ggplot(CYRE.cdf.2, aes(y = CYRE.cdf.2, x = Years))+
  geom_line()+
  geom_segment(aes(x = 0, y = 0.5, xend = x_at_half, yend = 0.5), size = 0.8) +
  geom_segment(aes(x = x_at_half, y = 0, xend = x_at_half, yend = 0.5), size = 0.8) +
  ylim(0.0,1.0)+
  theme_classic(base_size = 20)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

CYRE.countCDF = countCDFxt(mu = CYRE.mu.2, sig2 = CYRE.sig2.2, nt = 22,
                           Nc = 53, Ne = 10, tq = 43,
                           tmax = 50, Nboot = 10000, plot = TRUE)

ggplot(CYRE.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")


COCA.cdf.2 = extCDF(COCA.mu.2,COCA.sig2.2,Nc = 188, Ne = 10, tmax = 50)
COCA.cdf.2 = as.data.frame(COCA.cdf.2)
COCA.cdf.2$Years = 1:50

x_at_half <- COCA.cdf.2 %>%
  filter(abs(COCA.cdf.2 - 0.5) == min(abs(COCA.cdf.2 - 0.5))) %>%
  pull(Years)

ggplot(COCA.cdf.2, aes(y = COCA.cdf.2, x = Years))+
  geom_line()+
  geom_segment(aes(x = 0, y = 0.5, xend = x_at_half, yend = 0.5), size = 0.8) +
  geom_segment(aes(x = x_at_half, y = 0, xend = x_at_half, yend = 0.5), size = 0.8) +
  ylim(0.0,1.0)+
  theme_classic(base_size = 20)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

COCA.countCDF = countCDFxt(mu = COCA.mu.2, sig2 = COCA.sig2.2, nt = 11,
                           Nc = 188, Ne = 10, tq = 25,
                           tmax = 50, Nboot = 10000, plot = TRUE)

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

#### 50-year simulation to calculate the probability of extinction ####
# follows code from Bernardo et al. 2018

# LRR mean and sd for distribution
ARHI.LRR.mean = mean(ARHI.3$Corrected.LRR)
ARHI.LRR.sd = sd(ARHI.3$Corrected.LRR)

# output files
ARHI.sim = matrix(NA,nrow = 50, ncol = 10000)
ARHI.ext = rep(NA,10000)

# loop; 10000 replicate 50-year projections
for (t in 1:10000){

  # distribution of values to pick from
  ARHI.LRR.dist=rnorm(50, ARHI.LRR.mean, ARHI.LRR.sd)
  
  # starting population size 
  ARHI.Nt=33

  # loop; 50-year projection
  for (n in 1:50) 
  {
    Nt1 <- ARHI.Nt*exp(ARHI.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(Nt1<1) {
      Nt1=0 # to set a population as zero if it drops below 1 individual
    } else {
      Nt1=Nt1
    }
    ARHI.sim[n,t]=Nt1
    } #end of 50 year projection
  
  #capture extinctions
  if (min(ARHI.sim[,t])<=10){ #quasi-extinction threshold  = 10 Ind
    ARHI.ext[t]=1
  } else {
    ARHI.ext[t]=0
  }
} #end of 1,0000 replicates

# calculate the probability of extinction
ARHI.probext=sum(ARHI.ext)/10000 

#extract results; each column is one replicate 20-year projection
write.csv(ARHI.sim,"./Formatted.Data/ARHI.sim.graphing.csv")

# LRR mean and sd for distribution
COCA.LRR.mean = mean(COCA.2$Corrected.LRR)
COCA.LRR.sd = sd(COCA.2$Corrected.LRR)

# output files
COCA.sim = matrix(NA,nrow = 50, ncol = 10000)
COCA.ext = rep(NA,10000)

# loop; 10000 replicate 50-year projections
for (t in 1:10000){
  
  # distribution of values to pick from
  COCA.LRR.dist=rnorm(50, COCA.LRR.mean, COCA.LRR.sd)
  
  # starting population size 
  COCA.Nt=478
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    Nt1 <- COCA.Nt*exp(COCA.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(Nt1<1) {
      Nt1=0 # to set a population as zero if it drops below 1 individual
    } else {
      Nt1=Nt1
    }
    COCA.sim[n,t]=Nt1
  } #end of 50 year projection
  
  #capture extinctions
  if (min(COCA.sim[,t])<=10){ #quasi-extinction threshold  = 10 Ind
    COCA.ext[t]=1
  } else {
    COCA.ext[t]=0
  }
} #end of 1,0000 replicates

# calculate the probability of extinction
COCA.probext=sum(COCA.ext)/10000 

#extract results; each column is one replicate 20-year projection
write.csv(COCA.sim,"./Formatted.Data/COCA.sim.graphing.csv")

# LRR mean and sd for distribution
CYRE.LRR.mean = mean(CYRE.2$Corrected.LRR)
CYRE.LRR.sd = sd(CYRE.2$Corrected.LRR)

# output files
CYRE.sim = matrix(NA,nrow = 50, ncol = 10000)
CYRE.ext = rep(NA,10000)

# loop; 10000 replicate 50-year projections
for (t in 1:10000){
  
  # distribution of values to pick from
  CYRE.LRR.dist=rnorm(50, CYRE.LRR.mean, CYRE.LRR.sd)
  
  # starting population size 
  CYRE.Nt=53
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    Nt1 <- CYRE.Nt*exp(CYRE.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(Nt1<1) {
      Nt1=0 # to set a population as zero if it drops below 1 individual
    } else {
      Nt1=Nt1
    }
    CYRE.sim[n,t]=Nt1
  } #end of 50 year projection
  
  #capture extinctions
  if (min(CYRE.sim[,t])<=10){ #quasi-extinction threshold  = 10 Ind
    CYRE.ext[t]=1
  } else {
    CYRE.ext[t]=0
  }
} #end of 1,0000 replicates

# calculate the probability of extinction
CYRE.probext=sum(CYRE.ext)/10000 

#extract results; each column is one replicate 20-year projection
write.csv(CYRE.sim,"./Formatted.Data/CYRE.sim.graphing.csv")

# generate SPLU for Nt = 94 and Nt = 2 (the real last recorded value)

# LRR mean and sd for distribution
SPLU.LRR.mean = mean(SPLU.3$Corrected.LRR)
SPLU.LRR.sd = sd(SPLU.3$Corrected.LRR)

# output files
SPLU.sim = matrix(NA,nrow = 50, ncol = 10000)
SPLU.ext = rep(NA,10000)

# loop; 10000 replicate 50-year projections
for (t in 1:10000){
  
  # distribution of values to pick from
  SPLU.LRR.dist=rnorm(50, SPLU.LRR.mean, SPLU.LRR.sd)
  
  # starting population size 
  SPLU.Nt=94
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    Nt1 <- SPLU.Nt*exp(SPLU.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(Nt1<1) {
      Nt1=0 # to set a population as zero if it drops below 1 individual
    } else {
      Nt1=Nt1
    }
    SPLU.sim[n,t]=Nt1
  } #end of 50 year projection
  
  #capture extinctions
  if (min(SPLU.sim[,t])<=10){ #quasi-extinction threshold  = 10 Ind
    SPLU.ext[t]=1
  } else {
    SPLU.ext[t]=0
  }
} #end of 1,0000 replicates

# calculate the probability of extinction
SPLU.probext=sum(SPLU.ext)/10000 

#extract results; each column is one replicate 20-year projection
write.csv(SPLU.sim,"./Formatted.Data/SPLU.sim.graphing.csv")

#### Climate Data Prep ####
# annual climate year is the year of the survey, so counts in 1994 have 1994 annual climate
# monthly climate year is the year before the survey so counts in 1994 have climate from July 1993-June 1994

# merge LRR values for all pops into one file
ARHI.3$Species = "ARHI"
ARHI.3$Year = ARHI.2[c(1:15),2]
ARHI.3$Climate.Year = ARHI.3$Year-1
COCA.2$Species = "COCA"
COCA.2$Year = COCA[c(1:12),2]
COCA.2$Climate.Year = COCA.2$Year-1
CYRE.2$Species = "CYRE"
CYRE.2$Year = CYRE[c(1:23),2]
CYRE.2$Climate.Year = CYRE.2$Year-1
SPLU.3$Species = "SPLU"
SPLU.3$Year = SPLU[c(1:13),2]
SPLU.3$Climate.Year = SPLU.3$Year-1

LRR.all = rbind(ARHI.3,COCA.2,CYRE.2,SPLU.3)

# merge all count values for all pops
counts.all = rbind(ARHI.2,COCA,CYRE,SPLU.2)
counts.all$Climate.Year = counts.all$Year-1

# read in climate data
annual.clim = read.csv("./Formatted.Data/annual.climate.csv")
monthly.clim = read.csv("./Formatted.Data/monthly.climate.csv")

# changing units
annual.clim$ppt.mm = conv_unit(annual.clim$ppt..inches.,"inch","mm")
annual.clim$temp.C = conv_unit(annual.clim$tmean..degrees.F.,"F","C")
colnames(annual.clim)[1]="Year"

monthly.clim$ppt.mm = conv_unit(monthly.clim$ppt..inches.,"inch","mm")
monthly.clim$temp.C = conv_unit(monthly.clim$tmean..degrees.F.,"F","C")

# calculate yearly values from monthly values. Year is july 1 to June 30
monthly.clim = monthly.clim %>% 
  mutate(Date = ym(Date),
         Year = year(Date),
         Month = month(Date))

# Assign July-June as climate year
monthly.clim = monthly.clim %>% 
  mutate(Climate.Year = ifelse(Month >= 7, Year, Year - 1))

# summarize
climate.year = monthly.clim %>% 
  group_by(Species,Climate.Year) %>% 
  summarise(mean.ppt = mean(ppt.mm),
            mean.temp = mean(temp.C),
            sd.temp = sd(temp.C),
            cv.ppt = (sd(ppt.mm)/mean(ppt.mm))*100)

# merge climate data with LRR data
LRR.all.annual.clim = left_join(LRR.all,annual.clim, by = c("Species","Year"))
LRR.all.climate.year = left_join(LRR.all,climate.year, by = c("Species","Climate.Year"))

# merge climate data with count data
counts.all.annual.clim = left_join(counts.all,annual.clim, by = c("Species","Year"))
counts.all.climate.year = left_join(counts.all,climate.year, by = c("Species","Climate.Year"))

# scale ppt and temp
LRR.all.annual.clim = standardise(LRR.all.annual.clim, select = c("ppt.mm","temp.C"), append = TRUE)
LRR.all.climate.year = standardise(LRR.all.climate.year, select = c("mean.ppt","mean.temp"), append = TRUE)
counts.all.annual.clim = standardise(counts.all.annual.clim, select = c("ppt.mm","temp.C"), append = TRUE)
counts.all.climate.year = standardise(counts.all.climate.year, select = c("mean.ppt","mean.temp"), append = TRUE)

#### Climate Models ####

LRR.annual.mod = lmer(Corrected.LRR ~ ppt.mm_z + temp.C_z + (1|Species), data = LRR.all.annual.clim)
summary(LRR.annual.mod)
r2(LRR.annual.mod)
check_model(LRR.annual.mod)

LRR.annual.mod.2 = lmer(Corrected.LRR ~ ppt.mm_z + temp.C_z + I(temp.C_z^2) +(1|Species), data = LRR.all.annual.clim)
summary(LRR.annual.mod.2)
AIC(LRR.annual.mod,LRR.annual.mod.2)

LRR.climate.year.mod = lmer(Corrected.LRR ~ mean.ppt_z + mean.temp_z + (1|Species), data = LRR.all.climate.year)
summary(LRR.climate.year.mod)
r2(LRR.climate.year.mod)
check_model(LRR.climate.year.mod)

LRR.climate.year.mod.2 = lmer(Corrected.LRR ~ mean.ppt_z + mean.temp_z + I(mean.temp_z^2) +(1|Species), data = LRR.all.climate.year)
summary(LRR.annual.mod.2)
AIC(LRR.climate.year.mod,LRR.annual.mod.2)

counts.annual.mod = glmer.nb(Stem.Count ~ ppt.mm_z + temp.C_z + (1|Species),
                          data = counts.all.annual.clim)
summary(counts.annual.mod) # ppt is significant, negative
r2(counts.annual.mod)
check_model(counts.annual.mod)

counts.climate.year.mod = glmer.nb(Stem.Count ~ mean.ppt_z + mean.temp_z + (1|Species), data = counts.all.climate.year)
summary(counts.climate.year.mod)
r2(counts.climate.year.mod)
check_model(counts.climate.year.mod)

# species specific models
split.LRR.annual = split(LRR.all.annual.clim, LRR.all.annual.clim$Species)
split.LRR.climate.year = split(LRR.all.climate.year, LRR.all.climate.year$Species)
split.counts.annual = split(counts.all.annual.clim, counts.all.annual.clim$Species)
split.counts.climate.year = split(counts.all.climate.year, counts.all.climate.year$Species)

ARHI.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$ARHI)
summary(ARHI.LRR.annual.mod)
check_normality(ARHI.LRR.annual.mod)
check_heteroscedasticity(ARHI.LRR.annual.mod)
check_outliers(ARHI.LRR.annual.mod)
check_collinearity(ARHI.LRR.annual.mod)
check_autocorrelation(ARHI.LRR.annual.mod)

ARHI.LRR.climate.year.mod = lm(Corrected.LRR ~ mean.ppt_z + mean.temp_z, data = split.LRR.climate.year$ARHI)
summary(ARHI.LRR.climate.year.mod)
check_normality(ARHI.LRR.climate.year.mod)
check_heteroscedasticity(ARHI.LRR.climate.year.mod)
check_outliers(ARHI.LRR.climate.year.mod)
check_collinearity(ARHI.LRR.climate.year.mod)
check_autocorrelation(ARHI.LRR.climate.year.mod)# autocorrelated residuals

COCA.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$COCA)
summary(COCA.LRR.annual.mod)
check_normality(COCA.LRR.annual.mod) # not normal
check_heteroscedasticity(COCA.LRR.annual.mod)
check_outliers(COCA.LRR.annual.mod)
check_collinearity(COCA.LRR.annual.mod)
check_autocorrelation(COCA.LRR.annual.mod)

COCA.LRR.climate.year.mod = lm(Corrected.LRR ~ mean.ppt_z + mean.temp_z, data = split.LRR.climate.year$COCA)
summary(COCA.LRR.climate.year.mod)
check_normality(COCA.LRR.climate.year.mod)
check_heteroscedasticity(COCA.LRR.climate.year.mod)
check_outliers(COCA.LRR.climate.year.mod)
check_collinearity(COCA.LRR.climate.year.mod)
check_autocorrelation(COCA.LRR.climate.year.mod)

CYRE.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$CYRE)
summary(CYRE.LRR.annual.mod)

CYRE.LRR.climate.year.mod = lm(Corrected.LRR ~ mean.ppt_z + mean.temp_z, data = split.LRR.climate.year$CYRE)
summary(CYRE.LRR.climate.year.mod)
SPLU.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$SPLU)
summary(SPLU.LRR.annual.mod) # temperature is significant
SPLU.LRR.climate.year.mod = lm(Corrected.LRR ~ mean.ppt_z + mean.temp_z, data = split.LRR.climate.year$SPLU)
summary(SPLU.LRR.climate.year.mod)

ARHI.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$ARHI)
summary(ARHI.counts.annual.mod)
check_overdispersion()
check_zeroinflation()
check_outliers()
check_collinearity()
check_distribution()



ARHI.counts.climate.year.mod = glm.nb(Stem.Count ~ mean.ppt_z + mean.temp_z, data = split.counts.climate.year$ARHI)
summary(ARHI.counts.climate.year.mod)
COCA.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$COCA)
summary(COCA.counts.annual.mod)
COCA.counts.climate.year.mod = glm.nb(Stem.Count ~ mean.ppt_z + mean.temp_z, data = split.counts.climate.year$COCA)
summary(COCA.counts.climate.year.mod)
CYRE.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$CYRE)
summary(CYRE.counts.annual.mod)
CYRE.counts.climate.year.mod = glm.nb(Stem.Count ~ mean.ppt_z + mean.temp_z, data = split.counts.climate.year$CYRE)
summary(CYRE.counts.climate.year.mod)
SPLU.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$SPLU)
summary(SPLU.counts.annual.mod)
SPLU.counts.climate.year.mod = glm.nb(Stem.Count ~ mean.ppt_z + mean.temp_z, data = split.counts.climate.year$SPLU)
summary(SPLU.counts.climate.year.mod)


#### Simulate restoration practices ####
### Relationship between climate and counts
# Interesting to evaluate sensitivity of cdf to 
# Initial population size, extinction threshold,
# amount of environmental variability, number of time steps

