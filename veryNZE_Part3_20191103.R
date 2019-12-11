##################################################################
# Titel:      On the waning of forms - a corpus-based analysis of 
#             decline and loss in adjective amplification
# Part:       3
# R version:  3.5.1 (2018-07-02) -- "Feather Spray"
# Autor:      Martin Schweinberger
# Date:       2019-10-21
# Contact:    martin.schweinberger.hh@gmail.com
# Disclaimer: If you have questions,suggestions or you found errors
#             or in case you would to provide feedback, questions
#             write an email to martin.schweinberger.hh@gmail.com.
# Citation:   If you use this script or results thereof, please cite it as:
#             Schweinberger, Martin. 2019. On the waning of forms - a 
#             corpus-based analysis of decline and loss in the New 
#             Zealand English amplifier system.
#             Unpublished R script, The University of Queensland.
###############################################################
#                   START
###############################################################
# remove all lists from the current workspace
rm(list=ls(all=T))
# set wd
setwd("D:\\Uni\\Projekte\\02-Intensification\\AmpNZE\\LossAndDeclineAdjAmp")
# load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(stringr)
# set options
options(stringsAsFactors = F)
options(scipen = 999)
options(max.print=10000)
# define image directory
imageDirectory<-"images"
###############################################################
# load data
verywsc <- read.delim("datatables/ampwsc_regdat.txt", sep = "\t", header = T, skipNul = T)
# create a very dep. variable
verywsc <- verywsc %>%
  dplyr::select(-ID, - FileSpeaker, -Speaker, -File)
# simplify adjective
table(verywsc$Adjective)[order(table(verywsc$Adjective), decreasing = T)]

freqadj <- names(table(verywsc$Adjective))[which(table(verywsc$Adjective) > 20)]
verywsc$Adjective <- ifelse(verywsc$Adjective %in% freqadj, verywsc$Adjective, "other")
# factorize variables 
fctr <- c("Adjective", "Ethnicity", "Gender", "Education", "L1", "Age", 
          "Function", "Priming", "SemanticCategory", "Emotionality", "very")
verywsc[fctr] <- lapply(verywsc[fctr], factor)
# inspect data
nrow(verywsc); str(verywsc); head(verywsc); summary (verywsc)

###########################################################################
#                   SCALING
# frequency
summary(verywsc$Frequency)

verywsc$Frequency <- as.vector(scale(verywsc$Frequency))
summary(verywsc$Frequency)

plot(verywsc$Frequency, verywsc$very)
abline(lm(verywsc$very ~verywsc$Frequency))

# gradability
verywsc$Gradability <- as.vector(scale(verywsc$Gradability-1))
summary(verywsc$Gradability)

plot(verywsc$Gradability, verywsc$very)
abline(lm(verywsc$very ~ verywsc$Gradability))
###########################################################################
#                       STATISTICAL MODELLING
###########################################################################
#                             BORUTA
# perform variable selection
###########################################################################
#                  BORUTA
# load library
library(Boruta)
# create dada for boruta
borutadata <- verywsc %>%
  na.omit()
# save and load data
write.table(borutadata, "datatables/borutadata.txt", sep = "\t", 
            col.names = T, row.names = F, quote = F)
# run 1
boruta.verywsc <- Boruta(very ~.,data=borutadata)
print(boruta.verywsc)

getConfirmedFormula(boruta.verywsc)

png("images/BorutaVeryWsc.png",  width = 1800, height = 300)
plot(boruta.verywsc, cex = .75)
dev.off()
plot(boruta.verywsc)

png("images/BorutaVeryWsc_History1.png",  width = 680, height = 480)
plotImpHistory(boruta.verywsc)
dev.off()
plotImpHistory(boruta.verywsc)

# remove superfluous variables
borutadata$Gender <- NULL
borutadata$Function <- NULL
borutadata$Education <- NULL
# run2
boruta.verywsc <- Boruta(very~.,data=borutadata)
print(boruta.verywsc)

getConfirmedFormula(boruta.verywsc)

png("images/BorutaVeryWsc2.png",  width = 1200, height = 300)
plot(boruta.verywsc, cex = .75)
dev.off()
plot(boruta.verywsc)

png("images/BorutaNSP_History2.png",  width = 680, height = 480)
plotImpHistory(boruta.verywsc)
dev.off()
plotImpHistory(boruta.verywsc)

png("images/BorutaVeryWsc_final.png",  width = 1500, height = 750)
par(mar = c(22, 8, 4, 2) + 0.1)
plot(boruta.verywsc, cex.axis=3, las=2, xlab="", ylab = "", cex = 3, 
     col = c(rep("grey50", 9), rep("grey90",3)))
abline(v = 3.5, lty = "dashed")
mtext("Predictors", 1, line = 21, at = 9, cex = 3)
mtext("Control", 1, line = 21, at = 2, cex = 3)
mtext("Importance", 2, line = 5, at = 25, cex = 3, las = 0)
dev.off()
par(mar = c(5, 4, 4, 2) + 0.1)

###########################################################################
#           MIXED EFFECTS BIONOMIAL LOGISTIC REGRESSION
# load library
library(rms)
# set options
options(contrasts  =c("contr.treatment", "contr.poly"))
verywsc.dist <- datadist(verywsc)
options(datadist = "verywsc.dist")
# generate initial minimal regression model 
m0.glm = glm(very ~ 1, family = binomial, data = verywsc) # baseline model glm
m0.lrm = lrm(very ~ 1, data = verywsc, x = T, y = T) # baseline model lrm
# inspect results
summary(m0.glm)

m0.lrm

###########################################################################
# load libraries
library(lme4)
library(car)
# create model with a random intercept for token
m0.lmer <- lmer(very ~ (1|Adjective), data = verywsc, family = binomial)
# Baayen (2008:278-284) uses the call above but the this call is now longer
# up-to-date becwsce the "family" parameter is deprecated
# we switch to glmer (suggested by R) instead but we will also
# create a lmer object of the final minimal adequate model as some functions
# will not (yet) work on glmer
m0.glmer = glmer(very ~ (1|Adjective), data = verywsc, family = binomial)

# results of the lmer object
print(m0.lmer, corr = F)

# check if including the random effect is permitted by comparing the aic from the glm to aic from the glmer model
aic.glmer <- AIC(logLik(m0.glmer))
aic.glm <- AIC(logLik(m0.glm))
aic.glmer; aic.glm

# the aic of the glmer object is smaller which shows that including the random
# intercepts is justified

# test random effects
null.id = -2 * logLik(m0.glm) + 2 * logLik(m0.glmer)
pchisq(as.numeric(null.id), df=1, lower.tail=F) # sig m0.glmer better than m0.glm

# inspect results
summary(m0.glm)

summary(m0.glmer)

###########################################################################
# model fitting
# fit the model to find the "best" model, i.e. the minimal adequate model
# we will use a step-wise step up procedure
# we need to add "control = glmerControl(optimizer = "bobyqa")" 
# becwsce otherwise R fails to converge
#	manual modelfitting
m0.glmer <- glmer(very ~ 1 + (1|Adjective), family = binomial, data = verywsc, 
                  control=glmerControl(optimizer="bobyqa"))
# add Age
ifelse(min(ftable(verywsc$Age, verywsc$very)) == 0, "not possible", "possible")
m1.glm <- update(m0.glm, .~.+Age)
m1.glmer <- update(m0.glmer, .~.+Age)
anova(m0.glmer, m1.glmer, test = "Chi")                            # SIG! (p = 0.00000000000000022 ***)
Anova(m1.glmer, type = "III", test = "Chi")                        # SIG! (p = 0.00000000000000022 ***)

# add Adjective
ifelse(min(ftable(verywsc$Adjective, verywsc$very)) == 0, "not possible", "possible")
# not possible

# add Frequency
m2.glm <- update(m1.glm, .~.+Frequency)
max(vif(m2.glm)[,1])                                               # VIFs ok
m2.glmer <- update(m1.glmer, .~.+Frequency)
anova(m2.glmer, m1.glmer, test = "Chi")                            # not sig (p=0.618)

# add Emotionality
ifelse(min(ftable(verywsc$Emotionality, verywsc$very)) == 0, "not possible", "possible")
m3.glm <- update(m1.glm, .~.+Emotionality)
max(vif(m3.glm)[,1])                                               # VIFs ok
m3.glmer <- update(m1.glmer, .~.+Emotionality)
anova(m3.glmer, m1.glmer, test = "Chi")                            # SIG! (p = 0.00004315 ***)
Anova(m3.glmer, type = "III", test = "Chi")                        # SIG! (p = 0.000146 ***)

# add Ethnicity
ifelse(min(ftable(verywsc$Ethnicity, verywsc$very)) == 0, "not possible", "possible")
m4.glm <- update(m3.glm, .~.+Ethnicity)
max(vif(m4.glm)[,1])                                               # VIFs ok
m4.glmer <- update(m3.glmer, .~.+Ethnicity)
anova(m4.glmer, m3.glmer, test = "Chi")                            # not sig (p=0.9503)

# add SemanticCategory
ifelse(min(ftable(verywsc$SemanticCategory, verywsc$very)) == 0, "not possible", "possible")
m5.glm <- update(m3.glm, .~.+SemanticCategory)
max(vif(m5.glm)[,1])                                               # VIFs ok
m5.glmer <- update(m3.glmer, .~.+Gradability)
anova(m5.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06859 .)
Anova(m5.glmer, type = "III", test = "Chi")                        # mar sig (p = 0.0683430 .)

# add Gradability
m6.glm <- update(m3.glm, .~.+Gradability)
max(vif(m6.glm)[,1])                                               # VIFs ok
m6.glmer <- update(m3.glmer, .~.+Gradability)
anova(m6.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06859 .)
Anova(m6.glmer, type = "III", test = "Chi")                        # mar sig (p = 0.0683430 .)

# add Priming
ifelse(min(ftable(verywsc$Priming, verywsc$very)) == 0, "not possible", "possible")
m7.glm <- update(m3.glm, .~.+Priming)
max(vif(m7.glm)[,1])                                               # VIFs ok
m7.glmer <- update(m3.glmer, .~.+Priming)
anova(m7.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.8786)

# add L1
ifelse(min(ftable(verywsc$L1, verywsc$very)) == 0, "not possible", "possible")
m8.glm <- update(m3.glm, .~.+L1)
max(vif(m8.glm)[,1])                                               # VIFs ok
m8.glmer <- update(m3.glmer, .~.+L1)
anova(m8.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.09292 .)
Anova(m8.glmer, type = "III", test = "Chi")                        # mar sig (p = 0.0869448 .)

###########################################################################
# find all 2-way interactions
library(utils)
colnames(verywsc)

vars <- c("Age", "Frequency", "Emotionality", "Ethnicity", 
          "SemanticCategory", "Gradability", "Priming", "L1")
intac <- t(combn(vars, 2))
intac

# add Age * Frequency
m9.glm <- update(m3.glm, .~.+Age * Frequency)
max(vif(m9.glm)[,1])                                               # high VIFs
m9.glmer <- update(m3.glmer, .~.+Age * Frequency)
anova(m9.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.1792)

# add Age * Emotionality
m10.glm <- update(m3.glm, .~.+Age * Emotionality)
max(vif(m10.glm)[,1])                                              # VIFs unacceptable

# add Age * Ethnicity
m11.glm <- update(m3.glm, .~.+Age * Ethnicity)
max(vif(m11.glm)[,1])                                               
# Error in vif.default(m11.glm) : there are aliased coefficients in the model

# add Age * SemanticCategory
m12.glm <- update(m3.glm, .~.+Age * SemanticCategory)
max(vif(m12.glm)[,1])                                               
# Error in vif.default(m12.glm) : there are aliased coefficients in the model

# add Age * Gradability
m13.glm <- update(m3.glm, .~.+Age * Gradability)
max(vif(m13.glm)[,1])                                               # VIFs ok
m13.glmer <- update(m3.glmer, .~.+Age * Gradability)
anova(m13.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.3312)

# add Age * Priming
m14.glm <- update(m3.glm, .~.+Age * Priming)
max(vif(m14.glm)[,1])                                               # high VIFs
m14.glmer <- update(m3.glmer, .~.+Age * Priming)
anova(m14.glmer, m3.glmer, test = "Chi")                            # SIG! (p = 0.01175 *)
Anova(m14.glmer, type = "III", test = "Chi")                        # SIG! (p = 0.01175 *)
bic <- anova(m14.glmer, m3.glmer, test = "Chi")
ifelse(bic$BIC[2]-bic$BIC[1] < 0, "ok", "BIC inflation!")           # BIC inflation!
# return to model m3.glmer

# add Age * L1
m15.glm <- update(m3.glm, .~.+Age * L1)
max(vif(m15.glm)[,1])                                               # high VIFs
m15.glmer <- update(m3.glmer, .~.+Age * L1)
anova(m15.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2736)

# add Frequency * Emotionality
m16.glm <- update(m3.glm, .~.+Frequency * Emotionality)
max(vif(m16.glm)[,1])                                               # VIFs unacceptable

# add Frequency * Ethnicity
m17.glm <- update(m3.glm, .~.+Frequency * Ethnicity)
max(vif(m17.glm)[,1])                                               # VIFs ok
m17.glmer <- update(m3.glmer, .~.+Frequency * Ethnicity)
anova(m17.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.7228)

# add Frequency * SemanticCategory
m18.glm <- update(m3.glm, .~.+Frequency * SemanticCategory)
max(vif(m18.glm)[,1])                                               # VIFs unacceptable

# add Frequency * Gradability
m19.glm <- update(m3.glm, .~.+Frequency * Gradability)
max(vif(m19.glm)[,1])                                               # VIFs ok
m19.glmer <- update(m3.glmer, .~.+Frequency * Gradability)
anova(m19.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.1545)

# add Frequency * Priming
m20.glm <- update(m3.glm, .~.+Frequency * Priming)
max(vif(m20.glm)[,1])                                               # VIFs ok
m20.glmer <- update(m3.glmer, .~.+Frequency * Priming)
anova(m20.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.1897)

# add Frequency * L1 
m21.glm <- update(m3.glm, .~.+Frequency * L1)
max(vif(m21.glm)[,1])                                               # VIFs ok
m21.glmer <- update(m3.glmer, .~.+Frequency * L1)
anova(m21.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.3614)

# add Emotionality * Ethnicity  
m22.glm <- update(m3.glm, .~.+Emotionality * Ethnicity)
max(vif(m22.glm)[,1])                                               # VIFs unacceptable

# add Emotionality * SemanticCategory 
m23.glm <- update(m3.glm, .~.+Emotionality * SemanticCategory)
max(vif(m23.glm)[,1])                                               
# Error in vif.default(m23.glm) : there are aliased coefficients in the model

# add Emotionality * Gradability 
m24.glm <- update(m3.glm, .~.+Emotionality * Gradability)
max(vif(m24.glm)[,1])                                               # VIFs ok
m24.glmer <- update(m3.glmer, .~.+Emotionality * Gradability)
anova(m24.glmer, m3.glmer, test = "Chi")                            # SIG! (p = 0.004597 **)
Anova(m24.glmer, type = "III", test = "Chi")                        # SIG! (p = 0.0089823 **)
bic <- anova(m24.glmer, m3.glmer, test = "Chi")
ifelse(bic$BIC[2]-bic$BIC[1] < 0, "ok", "BIC inflation!")           # BIC inflation!
# return to model m3.glmer

# add Emotionality * Priming 
m25.glm <- update(m3.glm, .~.+Emotionality * Priming)
max(vif(m25.glm)[,1])                                               # high VIFs
m25.glmer <- update(m3.glmer, .~.+Emotionality * Priming)
anova(m25.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.9929)

# add Emotionality * L1 
m26.glm <- update(m3.glm, .~.+Emotionality * L1)
max(vif(m26.glm)[,1])                                               # high VIFs 
m26.glmer <- update(m3.glmer, .~.+Emotionality * L1)
anova(m26.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06568 .)
Anova(m26.glmer, type = "III", test = "Chi")                        # not sig (p = 0.1294207)

# add Ethnicity * SemanticCategory 
m27.glm <- update(m3.glm, .~.+Ethnicity * SemanticCategory)
max(vif(m27.glm)[,1])                                               # high VIFs 
# Error in vif.default(m27.glm) : there are aliased coefficients in the model

# add Ethnicity * Gradability 
m28.glm <- update(m3.glm, .~.+Ethnicity * Gradability)
max(vif(m28.glm)[,1])                                               # VIFs ok
m28.glmer <- update(m3.glmer, .~.+Ethnicity * Gradability)
anova(m28.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2035)

# add Ethnicity * Priming  
m29.glm <- update(m3.glm, .~.+Ethnicity * Priming)
max(vif(m29.glm)[,1])                                               # VIFs ok
m29.glmer <- update(m3.glmer, .~.+Ethnicity * Priming)
anova(m29.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.9995)

# add Ethnicity * L1   
m30.glm <- update(m3.glm, .~.+Ethnicity * L1)
max(vif(m30.glm)[,1])                                               # VIFs ok
m30.glmer <- update(m3.glmer, .~.+Ethnicity * L1)
anova(m30.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2482)

# add SemanticCategory * Gradability   
m31.glm <- update(m3.glm, .~.+SemanticCategory * Gradability)
max(vif(m31.glm)[,1])                                               # VIFs unacceptable

# add SemanticCategory * Priming   
m32.glm <- update(m3.glm, .~.+SemanticCategory * Priming)
max(vif(m32.glm)[,1])                                               # VIFs unacceptable

# add SemanticCategory * L1   
m33.glm <- update(m3.glm, .~.+SemanticCategory * L1)
max(vif(m33.glm)[,1])                                               # VIFs unacceptable
# Error in vif.default(m33.glm) : there are aliased coefficients in the model

# add Gradability * Priming  
m34.glm <- update(m3.glm, .~.+Gradability * Priming)
max(vif(m34.glm)[,1])                                               # VIFs ok
m34.glmer <- update(m3.glmer, .~.+Gradability * Priming)
anova(m34.glmer, m3.glmer, test = "Chi")                            # SIG! (p = 0.04327 *)
Anova(m34.glmer, type = "III", test = "Chi")                        # SIG! (p = 0.032951 *)
bic <- anova(m34.glmer, m3.glmer, test = "Chi")
ifelse(bic$BIC[2]-bic$BIC[1] < 0, "ok", "BIC inflation!")           # BIC inflation!
# return to model m3.glmer

# add Gradability * L1  
m35.glm <- update(m3.glm, .~.+Gradability * L1)
max(vif(m35.glm)[,1])                                               # VIFs ok
m35.glmer <- update(m3.glmer, .~.+Gradability * L1)
anova(m35.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06003 .)
Anova(m35.glmer, type = "III", test = "Chi")                        # not sig (p = 0.278465)

# add Priming * L1  
m36.glm <- update(m3.glm, .~.+Priming * L1)
max(vif(m36.glm)[,1])                                               # VIFs ok
m36.glmer <- update(m3.glmer, .~.+Priming * L1)
anova(m36.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2725)

#########################################
# load function for regression table summary
source("D:\\R/meblr.summary.R")
# set up summary table
meblrm_ampwsc <- meblrm.summary(m0.glm, m3.glm, m0.glmer, m3.glmer, verywsc$very) #
meblrm_ampwsc

# save results to disc
write.table(meblrm_ampwsc, "datatables/meblrm_ampwsc.txt", sep="\t")

# load function
library(car)
meblrm_ampwsc_Anova <- Anova(m3.glmer, type = "III", test = "Chi")
meblrm_ampwsc_Anova

# save results to disc
write.table(meblrm_ampwsc_Anova, "datatables/meblrm_ampwsc_Anova.txt", sep="\t")

effectage <- anova(m1.glmer, m0.glmer, test = "Chi")

effectemotionality <- anova(m3.glmer, m1.glmer, test = "Chi")

# use customized model comparison function
# create compariveryns
m1.m0 <- anova(m1.glmer, m0.glmer, test = "Chi")                            # SIG! (p = 0.00000000000000022 ***)
m2.m1 <- anova(m2.glmer, m1.glmer, test = "Chi")                            # not sig (p=0.618)
m3.m1 <- anova(m3.glmer, m1.glmer, test = "Chi")                            # SIG! (p = 0.00004315 ***)
m4.m3 <- anova(m4.glmer, m3.glmer, test = "Chi")                            # not sig (p=0.9503)
m5.m3 <- anova(m5.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06859 .)
m6.m3 <- anova(m6.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06859 .)
m7.m3 <- anova(m7.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.8786)
m8.m3 <- anova(m8.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.09292 .)
m9.m3 <- anova(m9.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.1792)
m13.m3 <- anova(m13.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.3312)
m14.m3 <- anova(m14.glmer, m3.glmer, test = "Chi")                            # SIG! (p = 0.01175 *)
m15.m3 <- anova(m15.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2736)
m17.m3 <- anova(m17.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.7228)
m19.m3 <- anova(m19.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.1545)
m20.m3 <- anova(m20.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.1897)
m21.m3 <- anova(m21.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.3614)
m24.m3 <- anova(m24.glmer, m3.glmer, test = "Chi")                            # SIG! (p = 0.004597 **)
m25.m3 <- anova(m25.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.9929)
m26.m3 <- anova(m26.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06568 .)
m28.m3 <- anova(m28.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2035)
m29.m3 <- anova(m29.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.9995)
m30.m3 <- anova(m30.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2482)
m34.m3 <- anova(m34.glmer, m3.glmer, test = "Chi")                            # SIG! (p = 0.04327 *)
m35.m3 <- anova(m35.glmer, m3.glmer, test = "Chi")                            # mar sig (p = 0.06003 .)
m36.m3 <- anova(m36.glmer, m3.glmer, test = "Chi")                            # not sig (p = 0.2725)

# create a list of the model compariveryns
mdlcmp <- list(m1.m0, m2.m1, m3.m1, m4.m3, m5.m3, m6.m3, m7.m3, m8.m3, 
               m9.m3, m13.m3, m14.m3, m15.m3, m17.m3, m19.m3, m20.m3, 
               m21.m3, m24.m3, m25.m3, m26.m3, m28.m3, m29.m3, m30.m3, 
               m34.m3, m35.m3, m36.m3)
# load function
source("D:\\R/ModelFittingSummarySWSU.R") # for Mixed Effects Model fitting (step-wise step-up): Binary Logistic Mixed Effects Models
# apply function
mdl.cmp.glmersc.swsu.dm <- mdl.fttng.swsu(mdlcmp)
# inspect output
mdl.cmp.glmersc.swsu.dm

write.table(mdl.cmp.glmersc.swsu.dm, "datatables/mdl_cmp_glmersc_swsu_verywsc.txt", sep="\t")
###########################################################
# Post-hoc analysis
library (multcomp)
summary(glht(m3.glmer, mcp(Age="Tukey")))

summary(glht(m3.glmer, mcp(Emotionality="Tukey")))

################################################################
#                 IMPORTANT OBJECTS
################################################################
# inspect very important objects
head(verywsc)

# glmer
effectage

effectemotionality

meblrm_ampwsc

meblrm_ampwsc_Anova

###############################################################
# predict probs of nativelike for effects
verywsc$PredictedFrequency <- predict(m3.glmer, verywsc, type="response")
summary(verywsc$PredictedFrequency)

# create response variable
verywsc$PredictedResponse <- ifelse(verywsc$PredictedFrequency > .5, "very", "other")
summary(verywsc$PredictedResponse)

# load library
library(caret)
# create confusion matrix
confusionMatrix(as.factor(verywsc$PredictedResponse), as.factor(verywsc$very))

# prepare plot data (pd) 
pd <- verywsc

# effect emotionality
p5 <- ggplot(pd, aes(x = Emotionality, y = PredictedFrequency)) +
  stat_summary(fun.y = mean, geom = "point") +          
  stat_summary(fun.data = mean_cl_boot, geom = "errorbar", width = 0.2) + 
  coord_cartesian(ylim = c(0, 1)) +
  theme_set(theme_bw(base_size = 20)) +
  #  theme(legend.position="none") +
  labs(x = "Emotionality", y = "Predicted probability\nof very") +
  ggsave(file = paste(imageDirectory,"PredictedEmotionality.png",sep="/"), 
       height = 5,  width = 7,  dpi = 320)
p5

# convert Age column
Agelbs <- names(table(pd$Age))
pd$Age <- ifelse(pd$Age == "16-19", 6,
                 ifelse(pd$Age == "20-29", 5,
                        ifelse(pd$Age == "30-39", 4, 
                               ifelse(pd$Age == "40-49", 3, 
                                      ifelse(pd$Age == "50-59", 2, 
                                             ifelse(pd$Age == "60+", 1, pd$Age))))))
pd$Age<- as.numeric(pd$Age)
# effect age
p6 <- ggplot(pd, aes(x = Age, y = PredictedFrequency)) +
  geom_smooth(aes(y = PredictedFrequency, x = Age), 
              colour="black", size=1, se = T, method = "loess") +
  coord_cartesian(ylim = c(0, 1)) +
  theme_set(theme_light(base_size = 15)) +
  #  theme(legend.position="none") +
  labs(x = "Frequency of adj. type", y = "Predicted probability\nof very") +
  scale_x_continuous(name = "Age",
                     breaks = c(1, 2, 3, 4, 5, 6),
                     labels=rev(Agelbs))
ggsave(file = paste(imageDirectory,"PredictedAge.png",sep="/"), 
       height = 5,  width = 5,  dpi = 320)
p6

library(effects)
png("images/effectsfinalmodel.png",  width = 960, height = 480) 
plot(allEffects(m3.glmer), type="response", ylim=c(0,1), grid=TRUE, 
     lines = list(col="black",
                  lty = 1,
                  confint=list(style="bars",
                               col = "grey80")), 
     ylab = "Prob (very)")
dev.off()

randomtb <- ranef(m3.glmer)
rndmadj <- as.vector(unlist(randomtb$`Adjective`))
adj <- as.vector(unlist(rownames(randomtb$`Adjective`)))
rndmadjtb <- data.frame(adj, rndmadj)
colnames(rndmadjtb) <- c("Adjective", "Intercept")
rndmadjtb <- rndmadjtb[order(rndmadjtb$Intercept, decreasing = T),]
rndmadjtb

p7 <- ggplot(rndmadjtb, aes(Adjective, Intercept)) +
  geom_point(aes(reorder(Adjective, -Intercept, fun = Intercept), y=Intercept)) +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  theme_set(theme_bw(base_size = 15)) +
  theme(legend.position="none", axis.text.x = element_text(size=15, angle=90)) +
  labs(x = "Adjective type", y = "Adjustment to Intercept")
ggsave(file = paste(imageDirectory,"RanAdjective.png",sep="/"), 
       height = 5,  width = 10,  dpi = 320)
# activate (remove #) to show
p7

###############################################################
write.table(verywsc, "datatables/verywsc.txt", row.names= F, sep = "\t")
###############################################################
#              END PART 3
###############################################################


