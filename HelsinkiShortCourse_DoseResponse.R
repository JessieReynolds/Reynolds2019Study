## Topic 4: Dose Response Analysis

library(devtools)
install_github("doseResponse/drcData")
install_github("doseResponse/drc")
library(drc)
library(drcData)

## Data for ryegrass are already included in package
head(ryegrass)

# Start by making a simple plot
plot(rootl ~ conc, data=ryegrass)

# Sample plot yt with log-transformed x-axis
plot(rootl ~conc, data=ryegrass, log="x")
# but note warning message that 6 values are omitted from log plot because they are less than 0

# Fitting a 4-parameter log-logistic model to data
ryegrass.LL.4 <- drm(rootl ~conc, data=ryegrass, fct=LL.4())

# Polt the model fit
plot(ryegrass.LL.4, broken =TRUE)

# Same plot, but including all data points
plot(ryegrass.LL.4, broken=TRUE, type = "all")

# Model control

# Residual plot
plot(resid(ryegrass.LL.4) ~ fitted(ryegrass.LL.4))
# note different variation of residuals 

# qqplot
qqnorm(resid(ryegrass.LL.4))
qqline(resid(ryegrass.LL.4))
## qqplot looks fine

# Summary of model
summary(ryegrass.LL.4)

# Confidence intervals
confint(ryegrass.LL.4)

# Plot with confidence bands
plot(ryegrass.LL.4, broken=TRUE, type ="confidence")

# Fitting a model where the lower limit is zero, known beforehand
ryegrass.LL.3 <- drm (rootl ~ conc, data = ryegrass, fct=LL.3())
summary(ryegrass.LL.3)

# Summary of new model
summary(ryegrass.LL.3)

# Estimate effective concentrations
ED(ryegrass.LL.4, c(10,20,50))

# Effective concentrations with confidence intervals
ED(ryegrass.LL.4, c(10,20,50), interval = "delta")

# Dealing with slight model misspecification
library(lmtest)
library(sandwich)

coeftest(ryegrass.LL.4, vcov. = sandwich)
# sandwich allows for standard error to vary; doesn't assume homogenity of variance

## Another data example - binary data
echovirus

# Fitting a dose-response model
echovirus.LL.2 <- drm(infected/total ~ dose, data=echovirus, fct=LL.2(), 
                      weights = total, type = "binomial")

# Plotting the model
plot(echovirus.LL.2, broken=TRUE, xlim=c(0,100000), ylim=c(0,1), bp=1)

# Fitting a 2-parameter log-normal model
echovirus.LN.2 <- drm(infected/total ~ dose, data = echovirus , fct=LN.2(),
                     weights=total, type="binomial")

# Fitting a 2-parameter weibull model
echovirus.W2.2 <- drm(infected/total ~ dose, data=echovirus, fct=W2.2(),
                     weights=total, type="binomial")

# Add fitted curves to existing plot
plot(echovirus.LN.2, broken=TRUE, xlim=c(0,100000, ylim=c(0,1), bp=1, add=TRUE, col="red"))
## NOT WORKING = something must be wrong, check PDF

# Estimate ED10
ED(echovirus.LL.2, 10)
ED(echovirus.LN.2, 10)
ED(echovirus.W2.2, 10)

# Finding the model-averaged ED10
maED(echovirus.LL.2, list(W2.2(), LN.2()), c(10), interval = 'kang')

# Install package bmd from github
install_github("doseResponse/bmd")
library(bmd)
library(bmdMA)


# Estimating BMD and BMDL for BMR =0.1 and added risk definition
bmd(echovirus.LL.2, bmr=0.1, def="additional", backgType = "modelBased")

# Estimate the mode-averaged BMD and BMDL with BMR=0.1 and added risk
bmdMA(list(echovirus.LL.2, echovirus.LN.2, echovirus.W2.2), bmr=0.1, 
      def="additional", backgType = "modelBased", modelWeights = "AIC",
      type = "Kang")

## Another example with binomial data

head(earthworms)

# Fitting a 2-parameter log-logistic model
earthworms.LL.2. <- drm(number/total ~ dose, data=earthworms, fct=LL.2(), 
                        weights=total, type="binomial")
plot(earthworms.LL.2., broken=TRUE)

# Upper limit should be 0.5, because expect worms to be distributed even in both
# halves of the box under control conditions

# Fit a model with the assumption that the upper limit is 0.5
earthworms.LL.3.fixed <- drm(number/total ~ dose, data=earthworms, fct=LL.3(fixed=c(NA,0.5,NA)),
                             weights = total, type ="binomial")
# order of variables: b, d, e

plot(earthworms.LL.3.fixed, broken = TRUE, add=TRUE, col="red")
# saying "add" adds this to the previous plot

## New example - time to event

head(chickweed)
tail(chickweed)
# after 281.5 hours, all plants that did not germinate are specified as Inf 160

# Fitting a model to these data
chickweed.LL.3 <- drm(count~start+end, data=chickweed, fct=LL.3(), type="event")
# lower limit is known = 0, so 3 parameter model is appropriate
# interval specied as X+Y

# Plot the fitted model
plot(chickweed.LL.3, log="")

# Summary of the model fit
summary(chickweed.LL.3)
# e means the time it takes for half of the seeds (.10) to germinate = 196 hr

# The time it takes for 5% of ALL seds to germinate
ED(chickweed.LL.3, 0.05, type="absolute")

# Last example - analysis of 2 dose-response curves - count data

head(decontaminants)
# concentrations, count, group of contaminants

# Tabulate data
with(decontaminants, table(conc,group))
# common control haphazardly assigned to oxalic, so need to code for that

# Fit a model that takes the common control into account
# Fitting a Poison dose response model with a common upper, estimated upper limit (common control)

decontaminants.LL.3 <- drm(count ~ conc, group, data=decontaminants, fct=LL.3(), 
                    type="Poisson", pmodels = list(~ group -1, ~ 1, ~ group-1))

plot(decontaminants.LL.3, broken=TRUE)

# Comparing slopes
compParm(decontaminants.LL.3, "b", "-") # difference in slopes

# Comparing ED50, e parameter
compParm(decontaminants.LL.3, "e", "/") # the relative potency

# Comparing ED10
EDcomp(decontaminants.LL.3,c(10,10))
# hpc looks to be 4 times that of oxalic, but note p value

# Dealing with overdispersion
library(lmtest)
library(sandwich)

coeftest(decontaminants.LL.3, vcov. = sandwich)
# doesn't change estimates, just standard errors; got bigger

# Fitting the model with a negative binomial distribution
decontaminants.LL.3.negbin <- drm(count ~ conc, group, data=decontaminants,
                                  fct=LL.3(), type="negbin2",
                                  pmodels = list(~ group -1, ~1, ~ group -1))

# Summary
summary(decontaminants.LL.3.negbin)

# Comparing ED value with sandwich standard errors
EDcomp(decontaminants.LL.3,c(10,10), vcov. = sandwich)

