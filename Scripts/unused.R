# merge all count values for all pops
counts.all = rbind(ARHI.2,COCA,CYRE,SPLU.2)
counts.all$Climate.Year = counts.all$Year-1

# read in climate data
annual.clim = read.csv("./Formatted.Data/annual.climate.csv")

annual.clim$ppt.mm = conv_unit(annual.clim$ppt..inches.,"inch","mm")
annual.clim$temp.C = conv_unit(annual.clim$tmean..degrees.F.,"F","C")
colnames(annual.clim)[1]="Year"

# merge climate data with LRR data
LRR.all.annual.clim = left_join(LRR.all,annual.clim, by = c("Species","Year"))

# merge climate data with count data
counts.all.annual.clim = left_join(counts.all,annual.clim, by = c("Species","Year"))
counts.all.climate.year = left_join(counts.all,climate.year, by = c("Species","Climate.Year"))

# scale ppt and temp
LRR.all.annual.clim = standardise(LRR.all.annual.clim, select = c("ppt.mm","temp.C"), append = TRUE)
counts.all.annual.clim = standardise(counts.all.annual.clim, select = c("ppt.mm","temp.C"), append = TRUE)
counts.all.climate.year = standardise(counts.all.climate.year, select = c("total.ppt","mean.temp"), append = TRUE)

counts.all.season = left_join(counts.all,all.season.climate, by = c("Species","Climate.Year"))


LRR.annual.mod = lmer(Corrected.LRR ~ ppt.mm_z + temp.C_z + (1|Species), data = LRR.all.annual.clim)
summary(LRR.annual.mod)
r2(LRR.annual.mod)
check_model(LRR.annual.mod)

LRR.annual.mod.2 = lmer(Corrected.LRR ~ ppt.mm_z + temp.C_z + I(temp.C_z^2) +(1|Species), data = LRR.all.annual.clim)
summary(LRR.annual.mod.2)
AIC(LRR.annual.mod,LRR.annual.mod.2)

counts.annual.mod = glmer.nb(Stem.Count ~ ppt.mm_z + temp.C_z + (1|Species),
                             data = counts.all.annual.clim)
summary(counts.annual.mod) # ppt is significant, negative
r2(counts.annual.mod)
check_model(counts.annual.mod)

counts.climate.year.mod = glmer.nb(Stem.Count ~ total.ppt_z + mean.temp_z + (1|Species), data = counts.all.climate.year)
summary(counts.climate.year.mod)
r2(counts.climate.year.mod)
check_model(counts.climate.year.mod)

split.LRR.annual = split(LRR.all.annual.clim, LRR.all.annual.clim$Species)
split.counts.annual = split(counts.all.annual.clim, counts.all.annual.clim$Species)
split.counts.climate.year = split(counts.all.climate.year, counts.all.climate.year$Species)

ARHI.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$ARHI)
summary(ARHI.LRR.annual.mod)
check_normality(ARHI.LRR.annual.mod)
check_heteroscedasticity(ARHI.LRR.annual.mod)
check_outliers(ARHI.LRR.annual.mod)
check_collinearity(ARHI.LRR.annual.mod)
check_autocorrelation(ARHI.LRR.annual.mod)

COCA.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$COCA)
summary(COCA.LRR.annual.mod)
check_normality(COCA.LRR.annual.mod) # not normal
check_heteroscedasticity(COCA.LRR.annual.mod)
check_outliers(COCA.LRR.annual.mod)
check_collinearity(COCA.LRR.annual.mod)
check_autocorrelation(COCA.LRR.annual.mod)

CYRE.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$CYRE)
summary(CYRE.LRR.annual.mod)

SPLU.LRR.annual.mod = lm(Corrected.LRR ~ ppt.mm_z + temp.C_z, data = split.LRR.annual$SPLU)
summary(SPLU.LRR.annual.mod) # temperature is significant

ARHI.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$ARHI)
summary(ARHI.counts.annual.mod)
check_overdispersion()
check_zeroinflation()
check_outliers()
check_collinearity()
check_distribution()

ARHI.counts.climate.year.mod = glm.nb(Stem.Count ~ total.ppt_z + mean.temp_z, data = split.counts.climate.year$ARHI)
summary(ARHI.counts.climate.year.mod)
COCA.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$COCA)
summary(COCA.counts.annual.mod)
COCA.counts.climate.year.mod = glm.nb(Stem.Count ~ total.ppt_z + mean.temp_z, data = split.counts.climate.year$COCA)
summary(COCA.counts.climate.year.mod)
CYRE.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$CYRE)
summary(CYRE.counts.annual.mod)
CYRE.counts.climate.year.mod = glm.nb(Stem.Count ~ total.ppt_z + mean.temp_z, data = split.counts.climate.year$CYRE)
summary(CYRE.counts.climate.year.mod)
SPLU.counts.annual.mod = glm.nb(Stem.Count ~ ppt.mm_z + temp.C_z, data = split.counts.annual$SPLU)
summary(SPLU.counts.annual.mod)
SPLU.counts.climate.year.mod = glm.nb(Stem.Count ~ total.ppt_z + mean.temp_z, data = split.counts.climate.year$SPLU)
summary(SPLU.counts.climate.year.mod)

counts.seasonal.mod = glmer.nb(Stem.Count ~ early.spring.temp_z + late.spring.temp_z +
                                 spring.ppt_z + non.spring.ppt_z + (1|Species), data = counts.all.season)
summary(counts.seasonal.mod) # non spring ppt is significant, negative
r2(counts.seasonal.mod)
check_model(counts.seasonal.mod)

summary(glmer.nb(Stem.Count ~ early.spring.temp_z + (1|Species), data = counts.all.season))
summary(glmer.nb(Stem.Count ~ late.spring.temp_z + (1|Species), data = counts.all.season))
summary(glmer.nb(Stem.Count ~ spring.ppt_z + (1|Species), data = counts.all.season))
summary(glmer.nb(Stem.Count ~ non.spring.ppt_z + (1|Species), data = counts.all.season)) # significant
summary(glmer.nb(Stem.Count ~ winter.temp_z + (1|Species), data = counts.all.season))

split.counts.season = split(counts.all.season, counts.all.season$Species)

summary(glm.nb(Stem.Count ~ early.spring.temp_z, data = split.counts.season$ARHI))
summary(glm.nb(Stem.Count ~ early.spring.temp_z, data = split.counts.season$COCA))
summary(glm.nb(Stem.Count ~ early.spring.temp_z, data = split.counts.season$CYRE))
summary(glm.nb(Stem.Count ~ early.spring.temp_z, data = split.counts.season$SPLU))

summary(glm.nb(Stem.Count ~ late.spring.temp_z, data = split.counts.season$ARHI))
summary(glm.nb(Stem.Count ~ late.spring.temp_z, data = split.counts.season$COCA))
summary(glm.nb(Stem.Count ~ late.spring.temp_z, data = split.counts.season$CYRE))
summary(glm.nb(Stem.Count ~ late.spring.temp_z, data = split.counts.season$SPLU))

summary(glm.nb(Stem.Count ~ spring.ppt_z, data = split.counts.season$ARHI))
summary(glm.nb(Stem.Count ~ spring.ppt_z, data = split.counts.season$COCA))
summary(glm.nb(Stem.Count ~ spring.ppt_z, data = split.counts.season$CYRE))
summary(glm.nb(Stem.Count ~ spring.ppt_z, data = split.counts.season$SPLU))

summary(glm.nb(Stem.Count ~ non.spring.ppt_z, data = split.counts.season$ARHI))
summary(glm.nb(Stem.Count ~ non.spring.ppt_z, data = split.counts.season$COCA))
summary(glm.nb(Stem.Count ~ non.spring.ppt_z, data = split.counts.season$CYRE))
# significant, negative
summary(glm.nb(Stem.Count ~ non.spring.ppt_z, data = split.counts.season$SPLU))

summary(glm.nb(Stem.Count ~ winter.temp_z, data = split.counts.season$ARHI))
summary(glm.nb(Stem.Count ~ winter.temp_z, data = split.counts.season$COCA))
summary(glm.nb(Stem.Count ~ winter.temp_z, data = split.counts.season$CYRE))
summary(glm.nb(Stem.Count ~ winter.temp_z, data = split.counts.season$SPLU))


