#-------------------------------------------------#
# Association Between Enactment of a “Stand Your 
# Ground” Self-defense Law and Unlawful Homicides 
# in Florida
# JAMA Internal Medicine
# DOI:10.1001/jamainternmed.2017.3433
# 
# David K. Humphreys, Antonio Gasparrini, 
# Douglas J. Wiebe (2017b)
# Email: david.humphreys@spi.ox.ac.uk
#-------------------------------------------------#

# This file provides the R code used for the analysis of example dataset used 
# in the tutorial paper.

#install.packages("lmtest") ; install.packages("Epi")
#install.packages("tsModel"); install.packages("vcd")
library(foreign) ; library(tsModel) ; library("lmtest") ; library("Epi")
library("splines") ; library("vcd"); 

##################
# 1. Read the data
##################

# 1.1 Read data from csv file
jhom <- read.csv("Data/All_Florida_Justifiable Homicide.csv")
head(jhom)

##################
# 2. Prerequisites
##################

# 2.1 Create a justifiable homicide rate per 100,000 population
jhom$rate <- (jhom$j_hom/jhom$stdpop)*100000

# 2.2 Create an unlawful homicide count
jhom$unlaw.h <- jhom$cdc_hom - jhom$j_hom 

# 2.3 Create an unlawful homicide rate
jhom$unlaw.h.rte <- (jhom$unlaw.h / jhom$stdpop) * 100000

##################
# 3. Descriptives
##################

# Table
# 3.1 Mean monthy counts: justifiable
summary(jhom$j_hom[jhom$Effective==0])
summary(jhom$j_hom[jhom$Effective==1])

# 3.2 Mean monthy counts: unlawful
summary(jhom$unlaw.h[jhom$Effective==0])
summary(jhom$unlaw.h[jhom$Effective==1])

# 3.3 Mean monthly rate per 100,000 pop: justifiable
summary(jhom$rate[jhom$Effective==0])
summary(jhom$rate[jhom$Effective==1])

# 3.4 Mean monthly rate per 100,000 pop: unlawful
summary(jhom$unlaw.h.rte[jhom$Effective==0])
summary(jhom$unlaw.h.rte[jhom$Effective==1])

#########################################################
# 4. Running a seasonal model
#########################################################

# 4.1 Justifiable homicide

jhom.model3 <- glm(j_hom ~ offset(log(stdpop)) + Effective + time + 
                     harmonic(month,2,12), family=quasipoisson, jhom)

summary(jhom.model3)
summary(jhom.model3)$dispersion
round(ci.lin(jhom.model3,Exp=T),3)

# 4.1.1 Breusch-Godfrey test for autocorrelation
bgtest(jhom.model3)
bgtest(jhom.model3, order=12) # p=0.001

# 4.1.2 Accounting for autocorrelation using a sandwich estimator

# Stand covariance matrix (returned by 'vcov' or as a component of 'summary')
vcov(jhom.model3 )[1:3,1:3]
summary(jhom.model3 )$cov.scaled[1:3,1:3]

# Standard estimates
summary(jhom.model3 )
coeftest(jhom.model3 )

# Covariance matrix, accounting for heteroskedascity and autocorrelation
library(sandwich)
vcov(jhom.model3 )[1:3,1:3]
vcovHAC(jhom.model3 )[1:3,1:3]

# 'Robust' estimates using the function vcovHAC in 'coeftest'
coeftest(jhom.model3 )
coeftest(jhom.model3 ,vcov=vcovHAC)

# 95% CIs must be calculated by hand: step change
coef <- coef(jhom.model3 )["Effective"]
se <- sqrt(vcovHAC(jhom.model3 )["Effective","Effective"])
c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se))

# 95% CIs must be calculated by hand: trend change 
coef2 <- coef(jhom.model3 )["time"]
se2 <- sqrt(vcovHAC(jhom.model3 )["time","time"])
c(RR=exp(coef2),ll=exp(coef2-qnorm(0.975)*se2),ul=exp(coef2+qnorm(0.975)*se2))



# 4.2 Unlawful homicide
jlawhom.m3 <- glm(unlaw.h ~ offset(log(stdpop)) + Effective + time + 
                    harmonic(month,2,12), family=quasipoisson, jhom)
summary(jlawhom.m3)
summary(jlawhom.m3)$dispersion
round(ci.lin(jlawhom.m3,Exp=T),3)

# 4.2.1 Breusch-Godfrey test for autocorrelation
bgtest(jlawhom.m3)
bgtest(jlawhom.m3, order=12) # p= 0.04

# 4.2.2 Accounting for autocorrelation using a sandwich estimator

# Stand covariance matrix (returned by 'vcov' or as a component of 'summary')
vcov(jlawhom.m3)[1:3,1:3]
summary(jlawhom.m3)$cov.scaled[1:3,1:3]

# Standard estimates
summary(jlawhom.m3)
coeftest(jlawhom.m3)

# Covariance matrix, accounting for heteroskedascity and autocorrelation
vcov(jlawhom.m3)[1:3,1:3]
vcovHAC(jlawhom.m3)[1:3,1:3]

# 'Robust' estimates using the function vcovHAC in 'coeftest'
coeftest(jlawhom.m3)
coeftest(jlawhom.m3,vcov=vcovHAC)

# 95% CIs must be calculated by hand: step change
coef <- coef(jlawhom.m3)["Effective"]
se <- sqrt(vcovHAC(jlawhom.m3)["Effective","Effective"])
c(RR=exp(coef),ll=exp(coef-qnorm(0.975)*se),ul=exp(coef+qnorm(0.975)*se))

# 95% CIs must be calculated by hand: trend change
coef2 <- coef(jlawhom.m3 )["time"]
se2 <- sqrt(vcovHAC(jlawhom.m3)["time","time"])
c(RR=exp(coef2),ll=exp(coef2-qnorm(0.975)*se2),ul=exp(coef2+qnorm(0.975)*se2))

#########################################################
# 5. Lawful Homicide as a proportion of Unlawful Homicide
# (In results section)
#########################################################

jhom$propJtohom <- (jhom$j_hom/jhom$cdc_hom) * 100
summary(jhom$propJtohom)
summary(jhom$propJtohom[jhom$Effective==0]) #3.4%
summary(jhom$propJtohom[jhom$Effective==1]) #8.7%


#---------------------------------END-----------------------------------#