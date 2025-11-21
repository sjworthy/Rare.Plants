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
library(corrplot)

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
sort(rstudent(ARHI.mod)) # one > 2, 13, none greater than 3
sort(rstudent(COCA.mod)) # one > 2, 12, yes it's greater than 3
sort(rstudent(CYRE.mod)) # three > 2, 2,9,10, none greater than 3
sort(rstudent(SPLU.mod)) # one > 2, 8, none greater than 3

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
exp(ARHI.r.low)
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

# formula is probability of extinction multiplied by CDF value at 50 years
SPLU.prob*2.757946e-01 # 0.08454067  

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
  ex <- extCDF(ARHI.mu, ARHI.sig2, Nc=n[i], Ne=10)
  exts[i] <- ex[50]
}

ARHI.countCDF = countCDFxt(mu = ARHI.mu, sig2 = ARHI.sig2, nt = 15,
                           Nc = 33, Ne = 10, tq = 29,
                           tmax = 50, Nboot = 10000, plot = TRUE)
ARHI.countCDF$Year = 1:50
write.csv(ARHI.countCDF, file = "./Formatted.Data/ARHI.diff.graphing.csv")

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
  ex <- extCDF(CYRE.mu, CYRE.sig2, Nc=n[i], Ne=10)
  exts[i] <- ex[50]
}

plot(n, exts, type='l', las=1,
     xlab="Current population size",
     ylab="Probability of quasi-extinction by year 50")

CYRE.countCDF = countCDFxt(mu = CYRE.mu, sig2 = CYRE.sig2, nt = 23,
                           Nc = 53, Ne = 10, tq = 43,
                           tmax = 50, Nboot = 500, plot = TRUE)
CYRE.countCDF$Year = 1:50
write.csv(CYRE.countCDF, file = "./Formatted.Data/CYRE.diff.graphing.csv")

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
  ex <- extCDF(COCA.mu, COCA.sig2, Nc=n[i], Ne=10)
  exts[i] <- ex[50]
}

plot(n, exts, type='l', las=1,
     xlab="Current population size",
     ylab="Probability of quasi-extinction by year 50")

COCA.countCDF = countCDFxt(mu = COCA.mu, sig2 = COCA.sig2, nt = 12,
                           Nc = 478, Ne = 10, tq = 26,
                           tmax = 50, Nboot = 500, plot = TRUE)
COCA.countCDF$Year = 1:50
write.csv(COCA.countCDF, file = "./Formatted.Data/COCA.diff.graphing.csv")


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
  ex <- extCDF(SPLU.mu, SPLU.sig2, Nc=n[i], Ne=10)
  exts[i] <- ex[50]
}

plot(n, exts, type='l', las=1,
     xlab="Current population size",
     ylab="Probability of quasi-extinction by year 50")

SPLU.countCDF = countCDFxt(mu = SPLU.mu, sig2 = SPLU.sig2, nt = 13,
                           Nc = 94, Ne = 10, tq = 16,
                           tmax = 50, Nboot = 10000, plot = TRUE)
SPLU.countCDF$Year = 1:50
write.csv(SPLU.countCDF, file = "./Formatted.Data/SPLU.diff.graphing.csv")

ggplot(SPLU.countCDF, aes(x = 1:50, y = Gbest))+
  geom_line()+
  geom_ribbon(aes(ymin = Glo, ymax = Gup), fill = "red", alpha = 0.2)+
  ylim(0.0,1.0)+
  theme_classic(base_size = 22)+
  xlab("Years into the Future")+
  ylab("Cumulative Probability of Quasi-Extinction")

# using real last census size of 2 in 2023
SPLU.countCDF = countCDFxt(mu = SPLU.mu, sig2 = SPLU.sig2, nt = 13,
                           Nc = 2, Ne = 10, tq = 29,
                           tmax = 50, Nboot = 10000, plot = TRUE)
SPLU.countCDF$Year = 1:50
write.csv(SPLU.countCDF, file = "./Formatted.Data/SPLU.diff.graphing.n2023.csv")

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

CYRE.countCDF$Year = 1:50
write.csv(CYRE.countCDF, file = "./Formatted.Data/CYRE.outlier.diff.graphing.csv")

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

COCA.countCDF$Year = 1:50
write.csv(COCA.countCDF, file = "./Formatted.Data/COCA.outlier.diff.graphing.csv")

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
set.seed(13)
for (t in 1:10000){

  # distribution of values to pick from
  ARHI.LRR.dist=rnorm(50, ARHI.LRR.mean, ARHI.LRR.sd)
  
  # starting population size 
  ARHI.Nt=33

  # loop; 50-year projection
  for (n in 1:50) 
  {
    ARHI.Nt <- ARHI.Nt*exp(ARHI.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(ARHI.Nt<1) {
      ARHI.Nt=0 # to set a population as zero if it drops below 1 individual
    } else {
      ARHI.Nt=ARHI.Nt
    }
    ARHI.sim[n,t]=ARHI.Nt
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
binom.test(7186,10000)

#extract results; each column is one replicate 20-year projection
write.csv(ARHI.sim,"./Formatted.Data/ARHI.sim.graphing.csv")

# output files
ARHI.sim.cdf = matrix(NA,nrow = 50, ncol = 10000)
ARHI.ext.cdf = rep(0,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  ARHI.LRR.dist=rnorm(50, ARHI.LRR.mean, ARHI.LRR.sd)
  
  # starting population size 
  ARHI.Nt=33
  extinct = FALSE
  
  # loop; 50-year projection
  for (n in 1:50) {
    
    # if already extinct in a previous year, keep at zero
    if (extinct) {
      ARHI.Nt <- 0
    } else {
      # project abundance
      ARHI.Nt <- ARHI.Nt * exp(ARHI.LRR.dist[n])
      
      # enforce integer quasi-extinction threshold
      if (ARHI.Nt < 1) ARHI.Nt <- 0
    }
    
    # store value
    ARHI.sim.cdf[n, t] <- ARHI.Nt
    
    # check for quasi-extinction (threshold = 10)
    if (ARHI.Nt <= 10 & !extinct) {
      extinct <- TRUE
      ARHI.ext.cdf[t] <- 1
    }
  }
}

yearly_ext_prob <- rowMeans(ARHI.sim.cdf <= 10)
yearly_ext_prob

n.reps <- ncol(ARHI.sim.cdf)

yearly_counts <- rowSums(ARHI.sim.cdf <= 10)

yearly_CI <- t(sapply(yearly_counts, function(k)
  binom.test(k, n.reps)$conf.int
))

colnames(yearly_CI) <- c("lower", "upper")
yearly_CI

ARHI.sim.cdf = as.data.frame(yearly_CI)
ARHI.sim.cdf$prob = yearly_ext_prob

write.csv(ARHI.sim.cdf, file = "./Formatted.Data/ARHI.sim.cdf.csv")

# LRR mean and sd for distribution
COCA.LRR.mean = mean(COCA.2$Corrected.LRR)
COCA.LRR.sd = sd(COCA.2$Corrected.LRR)

# output files
COCA.sim = matrix(NA,nrow = 50, ncol = 10000)
COCA.ext = rep(NA,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  COCA.LRR.dist=rnorm(50, COCA.LRR.mean, COCA.LRR.sd)
  
  # starting population size 
  COCA.Nt=478
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    COCA.Nt <- COCA.Nt*exp(COCA.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(COCA.Nt<1) {
      COCA.Nt=0 # to set a population as zero if it drops below 1 individual
    } else {
      COCA.Nt=COCA.Nt
    }
    COCA.sim[n,t]=COCA.Nt
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
binom.test(1193,10000)

#extract results; each column is one replicate 20-year projection
write.csv(COCA.sim,"./Formatted.Data/COCA.sim.graphing.csv") 

# output files
COCA.sim.cdf = matrix(NA,nrow = 50, ncol = 10000)
COCA.ext.cdf = rep(0,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  COCA.LRR.dist=rnorm(50, COCA.LRR.mean, COCA.LRR.sd)
  
  # starting population size 
  COCA.Nt=478
  extinct = FALSE
  
  # loop; 50-year projection
  for (n in 1:50) {
    
    # if already extinct in a previous year, keep at zero
    if (extinct) {
      COCA.Nt <- 0
    } else {
      # project abundance
      COCA.Nt <- COCA.Nt * exp(COCA.LRR.dist[n])
      
      # enforce integer quasi-extinction threshold
      if (COCA.Nt < 1) COCA.Nt <- 0
    }
    
    # store value
    COCA.sim.cdf[n, t] <- COCA.Nt
    
    # check for quasi-extinction (threshold = 10)
    if (COCA.Nt <= 10 & !extinct) {
      extinct <- TRUE
      COCA.ext.cdf[t] <- 1
    }
  }
}

yearly_ext_prob <- rowMeans(COCA.sim.cdf <= 10)
yearly_ext_prob

n.reps <- ncol(COCA.sim.cdf)

yearly_counts <- rowSums(COCA.sim.cdf <= 10)

yearly_CI <- t(sapply(yearly_counts, function(k)
  binom.test(k, n.reps)$conf.int
))

colnames(yearly_CI) <- c("lower", "upper")
yearly_CI

COCA.sim.cdf = as.data.frame(yearly_CI)
COCA.sim.cdf$prob = yearly_ext_prob

write.csv(COCA.sim.cdf, file = "./Formatted.Data/COCA.sim.cdf.csv")

# LRR mean and sd for distribution
CYRE.LRR.mean = mean(CYRE.2$Corrected.LRR)
CYRE.LRR.sd = sd(CYRE.2$Corrected.LRR)

# output files
CYRE.sim = matrix(NA,nrow = 50, ncol = 10000)
CYRE.ext = rep(NA,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  CYRE.LRR.dist=rnorm(50, CYRE.LRR.mean, CYRE.LRR.sd)
  
  # starting population size 
  CYRE.Nt=53
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    CYRE.Nt <- CYRE.Nt*exp(CYRE.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(CYRE.Nt<1) {
      CYRE.Nt=0 # to set a population as zero if it drops below 1 individual
    } else {
      CYRE.Nt=CYRE.Nt
    }
    CYRE.sim[n,t]=CYRE.Nt
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
binom.test(4641,10000)

#extract results; each column is one replicate 20-year projection
write.csv(CYRE.sim,"./Formatted.Data/CYRE.sim.graphing.csv")

# output files
CYRE.sim.cdf = matrix(NA,nrow = 50, ncol = 10000)
CYRE.ext.cdf = rep(0,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  CYRE.LRR.dist=rnorm(50, CYRE.LRR.mean, CYRE.LRR.sd)
  
  # starting population size 
  CYRE.Nt=53
  extinct = FALSE
  
  # loop; 50-year projection
  for (n in 1:50) {
    
    # if already extinct in a previous year, keep at zero
    if (extinct) {
      CYRE.Nt <- 0
    } else {
      # project abundance
      CYRE.Nt <- CYRE.Nt * exp(CYRE.LRR.dist[n])
      
      # enforce integer quasi-extinction threshold
      if (CYRE.Nt < 1) CYRE.Nt <- 0
    }
    
    # store value
    CYRE.sim.cdf[n, t] <- CYRE.Nt
    
    # check for quasi-extinction (threshold = 10)
    if (CYRE.Nt <= 10 & !extinct) {
      extinct <- TRUE
      CYRE.ext.cdf[t] <- 1
    }
  }
}

yearly_ext_prob <- rowMeans(CYRE.sim.cdf <= 10)
yearly_ext_prob

n.reps <- ncol(CYRE.sim.cdf)

yearly_counts <- rowSums(CYRE.sim.cdf <= 10)

yearly_CI <- t(sapply(yearly_counts, function(k)
  binom.test(k, n.reps)$conf.int
))

colnames(yearly_CI) <- c("lower", "upper")
yearly_CI

CYRE.sim.cdf = as.data.frame(yearly_CI)
CYRE.sim.cdf$prob = yearly_ext_prob

write.csv(CYRE.sim.cdf, file = "./Formatted.Data/CYRE.sim.cdf.csv")

# generate SPLU for Nt = 94 and Nt = 2 (the real last recorded value)

# LRR mean and sd for distribution
SPLU.LRR.mean = mean(SPLU.3$Corrected.LRR)
SPLU.LRR.sd = sd(SPLU.3$Corrected.LRR)

# output files
SPLU.sim = matrix(NA,nrow = 50, ncol = 10000)
SPLU.ext = rep(NA,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  SPLU.LRR.dist=rnorm(50, SPLU.LRR.mean, SPLU.LRR.sd)
  
  # starting population size 
  SPLU.Nt=94
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    SPLU.Nt <- SPLU.Nt*exp(SPLU.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(SPLU.Nt<1) {
      SPLU.Nt=0 # to set a population as zero if it drops below 1 individual
    } else {
      SPLU.Nt=SPLU.Nt
    }
    SPLU.sim[n,t]=SPLU.Nt
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
binom.test(2049,10000)

#extract results; each column is one replicate 20-year projection
write.csv(SPLU.sim,"./Formatted.Data/SPLU.sim.graphing.csv")

# output files
SPLU.sim.cdf = matrix(NA,nrow = 50, ncol = 10000)
SPLU.ext.cdf = rep(0,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  SPLU.LRR.dist=rnorm(50, SPLU.LRR.mean, SPLU.LRR.sd)
  
  # starting population size 
  SPLU.Nt=94
  extinct = FALSE
  
  # loop; 50-year projection
  for (n in 1:50) {
    
    # if already extinct in a previous year, keep at zero
    if (extinct) {
      SPLU.Nt <- 0
    } else {
      # project abundance
      SPLU.Nt <- SPLU.Nt * exp(SPLU.LRR.dist[n])
      
      # enforce integer quasi-extinction threshold
      if (SPLU.Nt < 1) SPLU.Nt <- 0
    }
    
    # store value
    SPLU.sim.cdf[n, t] <- SPLU.Nt
    
    # check for quasi-extinction (threshold = 10)
    if (SPLU.Nt <= 10 & !extinct) {
      extinct <- TRUE
      SPLU.ext.cdf[t] <- 1
    }
  }
}

yearly_ext_prob <- rowMeans(SPLU.sim.cdf <= 10)
yearly_ext_prob

n.reps <- ncol(SPLU.sim.cdf)

yearly_counts <- rowSums(SPLU.sim.cdf <= 10)

yearly_CI <- t(sapply(yearly_counts, function(k)
  binom.test(k, n.reps)$conf.int
))

colnames(yearly_CI) <- c("lower", "upper")
yearly_CI

SPLU.sim.cdf = as.data.frame(yearly_CI)
SPLU.sim.cdf$prob = yearly_ext_prob

write.csv(SPLU.sim.cdf, file = "./Formatted.Data/SPLU.sim.cdf.csv")

# LRR mean and sd for distribution
COCA.LRR.mean.2 = mean(COCA.3$Corrected.LRR)
COCA.LRR.sd.2 = sd(COCA.3$Corrected.LRR)

# output files
COCA.sim.2 = matrix(NA,nrow = 50, ncol = 10000)
COCA.ext.2 = rep(NA,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  COCA.LRR.dist.2=rnorm(50, COCA.LRR.mean.2, COCA.LRR.sd.2)
  
  # starting population size 
  COCA.Nt=188
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    COCA.Nt <- COCA.Nt*exp(COCA.LRR.dist.2[n]) #computes N(t+1) = Nt*exp(LRR);
    if(COCA.Nt<1) {
      COCA.Nt=0 # to set a population as zero if it drops below 1 individual
    } else {
      COCA.Nt=COCA.Nt
    }
    COCA.sim.2[n,t]=COCA.Nt
  } #end of 50 year projection
  
  #capture extinctions
  if (min(COCA.sim.2[,t])<=10){ #quasi-extinction threshold  = 10 Ind
    COCA.ext.2[t]=1
  } else {
    COCA.ext.2[t]=0
  }
} #end of 1,0000 replicates

# calculate the probability of extinction
COCA.probext.2=sum(COCA.ext.2)/10000 
binom.test(8768,10000)

#extract results; each column is one replicate 20-year projection
write.csv(COCA.sim.2,"./Formatted.Data/COCA.outlier.sim.graphing.csv")

# output files
COCA.sim.2.cdf = matrix(NA,nrow = 50, ncol = 10000)
COCA.ext.2.cdf = rep(0,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  COCA.LRR.dist.2=rnorm(50, COCA.LRR.mean.2, COCA.LRR.sd.2)
  
  # starting population size 
  COCA.Nt=188
  extinct = FALSE
  
  # loop; 50-year projection
  for (n in 1:50) {
    
    # if already extinct in a previous year, keep at zero
    if (extinct) {
      COCA.Nt <- 0
    } else {
      # project abundance
      COCA.Nt <- COCA.Nt * exp(COCA.LRR.dist.2[n])
      
      # enforce integer quasi-extinction threshold
      if (COCA.Nt < 1) COCA.Nt <- 0
    }
    
    # store value
    COCA.sim.2.cdf[n, t] <- COCA.Nt
    
    # check for quasi-extinction (threshold = 10)
    if (COCA.Nt <= 10 & !extinct) {
      extinct <- TRUE
      COCA.ext.2.cdf[t] <- 1
    }
  }
}

yearly_ext_prob <- rowMeans(COCA.sim.2.cdf <= 10)
yearly_ext_prob

n.reps <- ncol(COCA.sim.2.cdf)

yearly_counts <- rowSums(COCA.sim.2.cdf <= 10)

yearly_CI <- t(sapply(yearly_counts, function(k)
  binom.test(k, n.reps)$conf.int
))

colnames(yearly_CI) <- c("lower", "upper")
yearly_CI

COCA.2.sim.cdf = as.data.frame(yearly_CI)
COCA.2.sim.cdf$prob = yearly_ext_prob

write.csv(COCA.2.sim.cdf, file = "./Formatted.Data/COCA.2.sim.cdf.csv")

# LRR mean and sd for distribution
CYRE.LRR.mean.2 = mean(CYRE.3$Corrected.LRR)
CYRE.LRR.sd.2 = sd(CYRE.3$Corrected.LRR)

# output files
CYRE.sim.2 = matrix(NA,nrow = 50, ncol = 10000)
CYRE.ext.2 = rep(NA,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  CYRE.LRR.dist.2=rnorm(50, CYRE.LRR.mean.2, CYRE.LRR.sd.2)
  
  # starting population size 
  CYRE.Nt=53
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    CYRE.Nt <- CYRE.Nt*exp(CYRE.LRR.dist.2[n]) #computes N(t+1) = Nt*exp(LRR);
    if(CYRE.Nt<1) {
      CYRE.Nt=0 # to set a population as zero if it drops below 1 individual
    } else {
      CYRE.Nt=CYRE.Nt
    }
    CYRE.sim.2[n,t]=CYRE.Nt
  } #end of 50 year projection
  
  #capture extinctions
  if (min(CYRE.sim.2[,t])<=10){ #quasi-extinction threshold  = 10 Ind
    CYRE.ext.2[t]=1
  } else {
    CYRE.ext.2[t]=0
  }
} #end of 1,0000 replicates

# calculate the probability of extinction
CYRE.probext.2=sum(CYRE.ext.2)/10000 
binom.test(6819,10000)

#extract results; each column is one replicate 20-year projection
write.csv(CYRE.sim.2,"./Formatted.Data/CYRE.outlier.sim.graphing.csv")

# output files
CYRE.sim.2.cdf = matrix(NA,nrow = 50, ncol = 10000)
CYRE.ext.2.cdf = rep(0,10000)

# loop; 10000 replicate 50-year projections
set.seed(13)
for (t in 1:10000){
  
  # distribution of values to pick from
  CYRE.LRR.dist.2=rnorm(50, CYRE.LRR.mean.2, CYRE.LRR.sd.2)
  
  # starting population size 
  CYRE.Nt=53
  extinct = FALSE
  
  # loop; 50-year projection
  for (n in 1:50) {
    
    # if already extinct in a previous year, keep at zero
    if (extinct) {
      CYRE.Nt <- 0
    } else {
      # project abundance
      CYRE.Nt <- CYRE.Nt * exp(CYRE.LRR.dist.2[n])
      
      # enforce integer quasi-extinction threshold
      if (CYRE.Nt < 1) CYRE.Nt <- 0
    }
    
    # store value
    CYRE.sim.2.cdf[n, t] <- CYRE.Nt
    
    # check for quasi-extinction (threshold = 10)
    if (CYRE.Nt <= 10 & !extinct) {
      extinct <- TRUE
      CYRE.ext.2.cdf[t] <- 1
    }
  }
}

yearly_ext_prob <- rowMeans(CYRE.sim.2.cdf <= 10)
yearly_ext_prob

n.reps <- ncol(CYRE.sim.2.cdf)

yearly_counts <- rowSums(CYRE.sim.2.cdf <= 10)

yearly_CI <- t(sapply(yearly_counts, function(k)
  binom.test(k, n.reps)$conf.int
))

colnames(yearly_CI) <- c("lower", "upper")
yearly_CI

CYRE.2.sim.cdf = as.data.frame(yearly_CI)
CYRE.2.sim.cdf$prob = yearly_ext_prob

write.csv(CYRE.2.sim.cdf, file = "./Formatted.Data/CYRE.2.sim.cdf.csv")

# SPLU with 2023 pop size

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
  SPLU.Nt=2
  
  # loop; 50-year projection
  for (n in 1:50) 
  {
    SPLU.Nt <- SPLU.Nt*exp(SPLU.LRR.dist[n]) #computes N(t+1) = Nt*exp(LRR);
    if(SPLU.Nt<1) {
      SPLU.Nt=0 # to set a population as zero if it drops below 1 individual
    } else {
      SPLU.Nt=SPLU.Nt
    }
    SPLU.sim[n,t]=SPLU.Nt
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
write.csv(SPLU.sim,"./Formatted.Data/SPLU.sim.graphing.n2023.csv")

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
SPLU.3$Year = SPLU.2[c(1:13),2]
SPLU.3$Climate.Year = SPLU.3$Year-1

LRR.all = rbind(ARHI.3,COCA.2,CYRE.2,SPLU.3)

# read in climate data
monthly.clim = read.csv("./Formatted.Data/monthly.climate.csv")

# changing units
monthly.clim$ppt.mm = conv_unit(monthly.clim$ppt..inches.,"inch","mm")
monthly.clim$temp.C = conv_unit(monthly.clim$tmean..degrees.F.,"F","C")

# calculate yearly values from monthly values. Year is September 1 to August 31
monthly.clim = monthly.clim %>% 
  mutate(Date = ym(Date),
         Year = year(Date),
         Month = month(Date))

# Assign September-August as climate year
monthly.clim = monthly.clim %>% 
  mutate(Climate.Year = ifelse(Month >= 9, Year, Year - 1))

# summarize
climate.year = monthly.clim %>% 
  group_by(Species,Climate.Year) %>% 
  summarise(total.ppt = sum(ppt.mm),
            mean.temp = mean(temp.C),
            sd.temp = sd(temp.C),
            cv.ppt = (sd(ppt.mm)/mean(ppt.mm))*100)

season.climate = monthly.clim %>% 
  group_by(Species,Climate.Year) %>% 
  filter(Month %in% c(12,1,2)) %>% 
  summarise(winter.temp = mean(temp.C))

early.spring.temp = monthly.clim %>% 
  group_by(Species,Climate.Year) %>% 
  filter(Month %in% c(3,4)) %>% 
  summarise(early.spring.temp = mean(temp.C))

late.spring.temp = monthly.clim %>% 
  group_by(Species,Climate.Year) %>% 
  filter(Month %in% c(5,6)) %>% 
  summarise(late.spring.temp = mean(temp.C))

spring.ppt = monthly.clim %>% 
  group_by(Species,Climate.Year) %>% 
  filter(Month %in% c(3,4,5,6)) %>% 
  summarise(spring.ppt = sum(ppt.mm))

non.spring.ppt =  monthly.clim %>% 
  group_by(Species,Climate.Year) %>% 
  filter(Month %in% c(7,8,9,10,11,12,1,2,7)) %>% 
  summarise(non.spring.ppt = sum(ppt.mm))

dfs <- list(season.climate, early.spring.temp, late.spring.temp, spring.ppt, non.spring.ppt)

all.season.climate <- reduce(dfs, full_join, by = c("Species", "Climate.Year"))

# merge climate data with LRR data
LRR.all.climate.year = left_join(LRR.all,climate.year, by = c("Species","Climate.Year"))

# scale ppt and temp
LRR.all.climate.year = standardise(LRR.all.climate.year, select = c("total.ppt","mean.temp"), append = TRUE)
all.season.climate = standardise(all.season.climate, select = c("winter.temp","early.spring.temp",
                                                                "late.spring.temp","spring.ppt",
                                                                "non.spring.ppt"), append = TRUE)

# correlation of seasonal variables
cor.season = cor(all.season.climate[,c(8:12)], use = "pairwise")
corrplot(cor.season,method = "number")
# winter temp and early spring temp only variables semi-correlated (r = 0.47)

LRR.all.season = left_join(LRR.all,all.season.climate, by = c("Species","Climate.Year"))

#### Climate Year Models ####

LRR.climate.year.mod = lmer(Corrected.LRR ~ total.ppt_z + mean.temp_z + (1|Species), data = LRR.all.climate.year)
summary(LRR.climate.year.mod)
r2(LRR.climate.year.mod)
check_model(LRR.climate.year.mod)

# species specific models
split.LRR.climate.year = split(LRR.all.climate.year, LRR.all.climate.year$Species)
ARHI.climate.year = split.LRR.climate.year$ARHI[,c(1:11)]
COCA.climate.year = split.LRR.climate.year$COCA[,c(1:11)]
CYRE.climate.year = split.LRR.climate.year$CYRE[,c(1:11)]
SPLU.climate.year = split.LRR.climate.year$SPLU[,c(1:11)]

ARHI.climate.year = standardise(ARHI.climate.year, select = c("total.ppt","mean.temp"), append = TRUE)
COCA.climate.year = standardise(COCA.climate.year, select = c("total.ppt","mean.temp"), append = TRUE)
CYRE.climate.year = standardise(CYRE.climate.year, select = c("total.ppt","mean.temp"), append = TRUE)
SPLU.climate.year = standardise(SPLU.climate.year, select = c("total.ppt","mean.temp"), append = TRUE)

ARHI.LRR.climate.year.mod = lm(Corrected.LRR ~ total.ppt_z + mean.temp_z, data = ARHI.climate.year)
summary(ARHI.LRR.climate.year.mod)
check_normality(ARHI.LRR.climate.year.mod)
check_heteroscedasticity(ARHI.LRR.climate.year.mod)
check_collinearity(ARHI.LRR.climate.year.mod)
check_autocorrelation(ARHI.LRR.climate.year.mod)

COCA.LRR.climate.year.mod = lm(Corrected.LRR ~ total.ppt_z + mean.temp_z, data = COCA.climate.year)
summary(COCA.LRR.climate.year.mod)
check_normality(COCA.LRR.climate.year.mod)
check_heteroscedasticity(COCA.LRR.climate.year.mod)
check_collinearity(COCA.LRR.climate.year.mod)
check_autocorrelation(COCA.LRR.climate.year.mod)

CYRE.LRR.climate.year.mod = lm(Corrected.LRR ~ total.ppt_z + mean.temp_z, data = CYRE.climate.year)
summary(CYRE.LRR.climate.year.mod)
check_normality(CYRE.LRR.climate.year.mod)
check_heteroscedasticity(CYRE.LRR.climate.year.mod)
check_collinearity(CYRE.LRR.climate.year.mod)
check_autocorrelation(CYRE.LRR.climate.year.mod) # autocorrelated residuals

SPLU.LRR.climate.year.mod = lm(Corrected.LRR ~ total.ppt_z + mean.temp_z, data = SPLU.climate.year)
summary(SPLU.LRR.climate.year.mod) # mean temp is significant 0.012 (negative)
check_normality(SPLU.LRR.climate.year.mod)
check_heteroscedasticity(SPLU.LRR.climate.year.mod)
check_collinearity(SPLU.LRR.climate.year.mod)
check_autocorrelation(SPLU.LRR.climate.year.mod)

#### Seasonal climate Models ####
LRR.seasonal.mod = lmer(Corrected.LRR ~ early.spring.temp_z + late.spring.temp_z +
                          spring.ppt_z + non.spring.ppt_z + (1|Species), data = LRR.all.season)
summary(LRR.seasonal.mod)
r2(LRR.seasonal.mod)
check_model(LRR.seasonal.mod)

# species specific models
split.LRR.season = split(LRR.all.season, LRR.all.season$Species)
ARHI.LRR.season = split.LRR.season$ARHI[,c(1:12)]
COCA.LRR.season = split.LRR.season$COCA[,c(1:12)]
CYRE.LRR.season = split.LRR.season$CYRE[,c(1:12)]
SPLU.LRR.season = split.LRR.season$SPLU[,c(1:12)]

ARHI.LRR.season = standardise(ARHI.LRR.season, select = c("winter.temp","early.spring.temp",
                                                          "late.spring.temp","spring.ppt",
                                                          "non.spring.ppt"), append = TRUE)
COCA.LRR.season = standardise(COCA.LRR.season, select = c("winter.temp","early.spring.temp",
                                                          "late.spring.temp","spring.ppt",
                                                          "non.spring.ppt"), append = TRUE)
CYRE.LRR.season = standardise(CYRE.LRR.season, select = c("winter.temp","early.spring.temp",
                                                          "late.spring.temp","spring.ppt",
                                                          "non.spring.ppt"), append = TRUE)
SPLU.LRR.season = standardise(SPLU.LRR.season, select = c("winter.temp","early.spring.temp",
                                                          "late.spring.temp","spring.ppt",
                                                          "non.spring.ppt"), append = TRUE)


summary(lm(Corrected.LRR ~ early.spring.temp_z, data = ARHI.LRR.season))
summary(lm(Corrected.LRR ~ early.spring.temp_z, data = COCA.LRR.season))
summary(lm(Corrected.LRR ~ early.spring.temp_z, data = CYRE.LRR.season))
summary(lm(Corrected.LRR ~ early.spring.temp_z, data = SPLU.LRR.season))

summary(lm(Corrected.LRR ~ late.spring.temp_z, data = ARHI.LRR.season))
summary(lm(Corrected.LRR ~ late.spring.temp_z, data = COCA.LRR.season))
summary(lm(Corrected.LRR ~ late.spring.temp_z, data = CYRE.LRR.season))
summary(lm(Corrected.LRR ~ late.spring.temp_z, data = SPLU.LRR.season))

summary(lm(Corrected.LRR ~ spring.ppt_z, data = ARHI.LRR.season))
summary(lm(Corrected.LRR ~ spring.ppt_z, data = COCA.LRR.season))
summary(lm(Corrected.LRR ~ spring.ppt_z, data = CYRE.LRR.season))
summary(lm(Corrected.LRR ~ spring.ppt_z, data = SPLU.LRR.season))

summary(lm(Corrected.LRR ~ non.spring.ppt_z, data = ARHI.LRR.season))
summary(lm(Corrected.LRR ~ non.spring.ppt_z, data = COCA.LRR.season))
summary(lm(Corrected.LRR ~ non.spring.ppt_z, data = CYRE.LRR.season))
summary(lm(Corrected.LRR ~ non.spring.ppt_z, data = SPLU.LRR.season))

summary(lm(Corrected.LRR ~ winter.temp_z, data = ARHI.LRR.season))
summary(lm(Corrected.LRR ~ winter.temp_z, data = COCA.LRR.season))
summary(lm(Corrected.LRR ~ winter.temp_z, data = CYRE.LRR.season))
summary(lm(Corrected.LRR ~ winter.temp_z, data = SPLU.LRR.season))
# significant, negative 0.0013

#### Figure of CDFS and sims ####

# Read in CDF data
ARHI.cdf = read.csv("./Formatted.Data/ARHI.diff.graphing.csv", row.names = 1)
ARHI.cdf$Species = "ARHI"
COCA.cdf = read.csv("./Formatted.Data/COCA.diff.graphing.csv", row.names = 1)
COCA.cdf$Species = "COCA"
CYRE.cdf = read.csv("./Formatted.Data/CYRE.diff.graphing.csv", row.names = 1)
CYRE.cdf$Species = "CYRE"
SPLU.cdf = read.csv("./Formatted.Data/SPLU.diff.graphing.csv", row.names = 1)
SPLU.cdf$Species = "SPLU"
COCA.out.cdf = read.csv("./Formatted.Data/COCA.outlier.diff.graphing.csv", row.names = 1)
COCA.out.cdf$Species = "COCA.out"
CYRE.out.cdf = read.csv("./Formatted.Data/CYRE.outlier.diff.graphing.csv", row.names = 1)
CYRE.out.cdf$Year = 1:50
CYRE.out.cdf$Species = "CYRE.out"
SPLU.n.cdf = read.csv("./Formatted.Data/SPLU.diff.graphing.n2023.csv", row.names = 1)
SPLU.n.cdf$Species = "SPLU.n"

# merge
dfs <- list(ARHI.cdf, COCA.cdf, CYRE.cdf, SPLU.cdf, COCA.out.cdf, CYRE.out.cdf)

all.cdfs <- reduce(dfs, full_join)

# plot
cdf.outlier.plot = ggplot(all.cdfs, aes(x = Year, y = Gbest, color = Species))+
  geom_line(aes(linetype = Species))+
  scale_color_manual(values = c("#0A9F9D","#E54E21","#E54E21","#6C8645","#6C8645","#35274A"),
                     labels = c("A. hispida","C. canadensis", "C. canadensis",
                                "C. reginae", "C. reginae", "S. lucida"))+
  scale_linetype_manual(values = c("solid", "solid", "dashed",
                                   "solid", "dashed", "solid"),
                        labels = c("A. hispida","C. canadensis", "C. canadensis",
                                   "C. reginae", "C. reginae", "S. lucida")) +
  theme_classic(base_size = 20)+
  labs(x = "Years into the Future", y = "Cumulative Probability of Quasi-Extinction",
       title = "Diffusion Approximation")+
 # theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(0,1)
cdf.outlier.plot

#ggsave(cdf.outlier.plot, file = "./Plots/cdf.outlier.plot.legend.pdf", height = 8, width = 8)
#ggsave(cdf.outlier.plot, file = "./Plots/cdf.outlier.plot.pdf", height = 8, width = 8)


# read in simulations
ARHI.sim = read.csv("./Formatted.Data/ARHI.sim.cdf.csv",row.names = 1)
ARHI.sim$Species = "ARHI"
ARHI.sim$Year = 1:50
COCA.sim = read.csv("./Formatted.Data/COCA.sim.cdf.csv",row.names = 1)
COCA.sim$Species = "COCA"
COCA.sim$Year = 1:50
CYRE.sim = read.csv("./Formatted.Data/CYRE.sim.cdf.csv",row.names = 1)
CYRE.sim$Species = "CYRE"
CYRE.sim$Year = 1:50
SPLU.sim = read.csv("./Formatted.Data/SPLU.sim.cdf.csv",row.names = 1)
SPLU.sim$Species = "SPLU"
SPLU.sim$Year = 1:50
COCA.out.sim = read.csv("./Formatted.Data/COCA.2.sim.cdf.csv",row.names = 1)
COCA.out.sim$Species = "COCA.out"
COCA.out.sim$Year = 1:50
CYRE.out.sim = read.csv("./Formatted.Data/CYRE.2.sim.cdf.csv",row.names = 1)
CYRE.out.sim$Species = "CYRE.out"
CYRE.out.sim$Year = 1:50

# merge
dfs.sim <- list(ARHI.sim, COCA.sim, CYRE.sim, COCA.out.sim, CYRE.out.sim, SPLU.sim)

all.sims <- reduce(dfs.sim, full_join)

# plot
ggplot(all.sims, aes(x = Year, y = prob, color = Species))+
  geom_line()

# plot
sim.outlier.plot = ggplot(all.sims, aes(x = Year, y = prob, color = Species))+
  geom_line(aes(linetype = Species))+
  scale_color_manual(values = c("#0A9F9D","#E54E21","#E54E21","#6C8645","#6C8645","#35274A"),
                     labels = c("A. hispida","C. canadensis", "C. canadensis",
                                "C. reginae", "C. reginae", "S. lucida"))+
  scale_linetype_manual(values = c("solid", "solid", "dashed",
                                   "solid", "dashed", "solid"),
                        labels = c("A. hispida","C. canadensis", "C. canadensis",
                                   "C. reginae", "C. reginae", "S. lucida")) +
  theme_classic(base_size = 20)+
  labs(x = "Years into the Future", y = "Cumulative Probability of Quasi-Extinction", 
       title = "Stochastic Projection")+
  theme(legend.position = "none")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylim(0,1)
sim.outlier.plot

#ggsave(sim.outlier.plot, file = "./Plots/sim.outlier.plot.legend.pdf", height = 8, width = 8)
#ggsave(sim.outlier.plot, file = "./Plots/sim.outlier.plot.pdf", height = 8, width = 8)



