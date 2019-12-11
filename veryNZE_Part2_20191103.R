##################################################################
# Titel:      On the waning of forms - a corpus-based analysis of 
#             decline and loss in adjective amplification
# Part:       2
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
# load packages
library(dplyr)
library(ggplot2)
library(stringr)
# set options
options(stringsAsFactors = F)
options(scipen = 999)
options(max.print=10000)
# define image directory
imageDirectory<-"images"
###############################################################
# load data
ampwsc <- read.delim("datatables/wsc_fullclean.txt", sep = "\t", 
                     header = T, skipNul = T)
# inspect data
str(ampwsc)

###############################################################
#                           EXAMPLES
exdata <- ampwsc %>%
  dplyr::select(PreContextLong, Token, PostContext, PreContext, Function, File, Speaker) %>%
  mutate(Example = paste(PreContextLong, Token, PostContext, sep = " ")) %>%
  mutate(ExampleClean = str_replace_all(Example, "/[A-Z]{1,3}", "")) %>%
  dplyr::select(ExampleClean, Function, PreContext, File, Speaker) %>%
  dplyr::rename(Example = ExampleClean, File = File, Amplifier = PreContext, Speaker = Speaker)
# extract predicative contexts
examples_predicative <- exdata[exdata$Function == "Predicative",]
examples_predicative <- examples_predicative[order(sample(nrow(examples_predicative), nrow(examples_predicative))),]
# variants
very_predicative <- examples_predicative[examples_predicative$Amplifier == "very",]
pretty_predicative <- examples_predicative[examples_predicative$Amplifier == "pretty",]
so_predicative <- examples_predicative[examples_predicative$Amplifier == "so",]
really_predicative <- examples_predicative[examples_predicative$Amplifier == "really",]
# extract attributive contexts
examples_attributive <- exdata[exdata$Function == "Attributive",]
examples_attributive <- examples_attributive[order(sample(nrow(examples_attributive),nrow(examples_attributive))),]
# variants
very_attributive <- examples_attributive[examples_attributive$Amplifier == "very",]
pretty_attributive <- examples_attributive[examples_attributive$Amplifier == "pretty",]
so_attributive <- examples_attributive[examples_attributive$Amplifier == "so",]
real_attributive <- examples_attributive[examples_attributive$Amplifier == "real",]
head(very_attributive); head(very_predicative)

###############################################################
#                  CLEAN DATA
ampwsc <- ampwsc %>%
  dplyr::filter(Type == "Conversation") %>%
  dplyr::select(- Text, -CleanText, -TextPOS, -OriginalString, 
                -ExtractMinutes, -Words.1, -PreContext, -PostContext,
                -Minutes, -PreContextLong, -Id, - Topic, -Words, -Type) %>%
  na.omit()
# recode age
ampwsc$Age <- ifelse(ampwsc$Age == "16-19", "16-19",
                     ifelse(ampwsc$Age == "20-24", "20-29",
                            ifelse(ampwsc$Age == "25-29", "20-29",
                                   ifelse(ampwsc$Age == "30-34", "30-39",
                                          ifelse(ampwsc$Age == "35-39", "30-39",
                                                 ifelse(ampwsc$Age == "40-44", "40-49",
                                                        ifelse(ampwsc$Age == "45-49", "40-49",
                                                               ifelse(ampwsc$Age == "50-54", "50-59",
                                                                      ifelse(ampwsc$Age == "55-59", "50-59", "60+"))))))))) 
# create a numeric age variable
ampwsc$AgeNumeric <- ifelse(ampwsc$Age == "16-19", 17.5,
                     ifelse(ampwsc$Age == "20-24", 22,
                            ifelse(ampwsc$Age == "25-29", 27,
                                   ifelse(ampwsc$Age == "30-34", 32,
                                          ifelse(ampwsc$Age == "35-39", 37,
                                                 ifelse(ampwsc$Age == "40-44", 42,
                                                        ifelse(ampwsc$Age == "45-49", 47,
                                                               ifelse(ampwsc$Age == "50-54", 52,
                                                                      ifelse(ampwsc$Age == "55-59", 57,
                                                                             ifelse(ampwsc$Age == "60-54", 62,
                                                                                    ifelse(ampwsc$Age == "65-69", 67,
                                                                                           ifelse(ampwsc$Age == "70-74", 72,
                                                                                                  ifelse(ampwsc$Age == "75-79", 77,
                                                                                                         ifelse(ampwsc$Age == "80-84", 82,
                                                                                                                ifelse(ampwsc$Age == "85-89", 87,
                                                                                                                       ifelse(ampwsc$Age == "90+", 90, NA))))))))))))))))
summary(ampwsc$AgeNumeric)

# recode education
University <- c("1.5 YRS UNIVERSITY", "1ST YEAR BA", "1ST YEAR BA/LLB", 
                "2 YRS UNI", "2ND YEAR UNIVERSITY",  "3RD YEAR UNIVERSITY",  
                "ACA", "ALMOST LLB/BA", "ANZIQS", "ATLL", "B AGRI", 
                "B CAPS DEGREE", "B E", "B ED", "B ED DIP BUS STUDIES", 
                "B MUS (HONS)", "B SOCSCI", "B TECH HONS", "BA", "BA (HONS)", 
                "BA (LATER BA HONS)", "BA BSC", "BA EDUC", 
                "BA EDUCATION & PHILOSOPHY", "BA HONS", "BA/DIP TOWN PLANNING", 
                "BA/POSTGRAD DIP", "BACHELOR DEGREE", "BBSC", "BCA", "BCA / LLB", 
                "BCA HONS", "BCA/BA", "BCOM", "BCOM (HONS)", "BMUS", 
                "BMUS (HONS)", "BP RM", "BSC", "BSC BA", "BSC CHEMISTRY", 
                "BSC HONS", "BSW HONS",  "BURSARY + UNI PAPERS", 
                "BURSARY, 4TH YR UNIV", "FINAL YR ENGLISH STUDENT", 
                "HIGHER DIP TEACH, BA", "HLC", "IST YEAR UNIVERSITY", 
                "LLB BCA", "M AG SCI", "M MUS", "M SOC SCI", "MA", "MA APPLIED", 
                "MA DIP ED", "MA HONS", "MA HONS DIP CLINICAL PSYCH", "MA MBA", 
                "MASTERS", "MATRICULATION", "MBA", "MSC", 
                "NCBC (BUSINESS COMPUTING)", "NZCE + UNI", "NZCE, B E", 
                "NZCS", "PART OF BA", "PARTIAL BA", "PARTIAL SC", "PHD", 
                "POST GRAD DEGREE", "POST GRAD DIP", "POSTGRAD DIP", 
                "PROFICIENCY", "PROFICIENCY, FTCL", "PUBLIC SERVICE", 
                "REGISTERED NURSE", "SCHOLARSHIP",  "SCRNZN", "SFC", 
                "STAGE 2 ACCOUNTING", "STAGE 2 UNI", "STAGE 2 UNIVERSITY", 
                "STAGE 3 UNIVERSITY", "STAGE ONE UNIVERSITY", "TERTIARY DEGREE", 
                "UE, ASTU UNI PAPERS",  "UNDER GRADUATE", 
                "UNI POSTGRAD DIP", "UNIVERSITY DEGREE")
NoUniversity <- c("1.5 YRS SECONDARY", "1ST QUAL IN PASTRY", 
                  "2 YRS SECONDARY SCHOOL", "3 YRS SECONDARY", 
                  "4 YRS PRIMARY SCHOOL", "4TH FORM", "5 YEARS SECONDARY", 
                  "5TH FORM", "6TH FORM", "6TH FORM ENGLISH", 
                  "6TH FORM MAORI", "7TH FORM CERT", "A BURSARY", "A BURSARY/UE",
                  "FORM 2", "FORM 3", "FORM 5", "HSC",  "HSC/BURSARY", "HSC/UE",
                  "NONE", "SC", "SC / TRADE CERT", "SC, CONTINUING EDUC", 
                  "SC, TRADE CERT", "SCHOOL", "BURSARY ENGLISH","BURSARY",
                  "BURSARY/SCHOLARSHIP", "UE / NURSING DIPLOMA","UE/BURSARY",
                  "NIL", "NURSING REGISTRATION", "NZ BUILDING CERT", 
                  "NZ LIB CERT", "CHAMBER OF COMM CERT", "COMMERCIAL EXAM JP", 
                  "DIP DESIGN", "DIP HOTEL & CATERING", "DIP MED LAB TECH", 
                  "DIP MGMT DEVT", "DIP OCCUPATIOANL THERAPY", "DIP PHARMACY", 
                  "DIP PROF DRAMA", "DIP TCHG", "DIP VISUAL COMM DES", "ECE DIP", 
                  "ADV DIP TCHG", "ADV TEACHER DIP", "ADV TRADE CERT", "ENGLISH", 
                  "LIB DIP", "B BURSARY",  "MANAGEMENT CERT", "HIGHER DIP TCHG", 
                  "TOHU MATAURANGA MAORI", "TRADE CERT", "TRADE CERT CATERING", 
                  "TTC", "TTC/TOHU MAAORI", "UE",  "UE EQUIV", "UE LTCL", 
                  "UE, COMM TEACH CERT", "UE, NZ LIBRARY CERT")
ampwsc$Education <- ifelse(ampwsc$Education %in% University, "University",
                           ifelse(ampwsc$Education %in% NoUniversity, "NoUniversity",
                                  NA))
# recode ethnicity
ampwsc$Ethnicity <- ifelse(ampwsc$Ethnicity == "PAKEHA", "PAKEHA",
                           ifelse(ampwsc$Ethnicity == "MAORI", "MAORI",
                                  ifelse(ampwsc$Ethnicity == "MAORI/PAKEHA", "MAORI/PAKEHA",
                                         ifelse(ampwsc$Ethnicity == "PAKEHA/MAORI", "MAORI/PAKEHA", 
                                                ifelse(is.na(ampwsc$Ethnicity) == T, NA, "OTHER")))))
# recode L1
ampwsc$L1 <- ifelse(ampwsc$L1 == "ENGLISH", "ENGLISH", 
                    ifelse(is.na(ampwsc$L1) == T, NA, "OTHER"))
# inspect data
summary(ampwsc)

###############################################################
#                           LD EXAMPLE TABLE
# inspect data
str(ampwsc)
ldtb <- ampwsc %>%
  dplyr::select(Age, Adjective, Variant) %>%
  dplyr::group_by(Age, Variant) %>%
  dplyr::summarize(AdjectiveTypes = length(table(Adjective)), AllTypes = sum(table(Adjective))) %>%
  dplyr::filter(Variant == "very") %>%
  dplyr::mutate(LD = AdjectiveTypes/AllTypes)
head(ldtb)

###############################################################
# prepare data for plotting
# create data frame with relevant variables
pd <- data.frame(ampwsc$Age, ampwsc$Age, ampwsc$Function, ampwsc$Amplified, ampwsc$Variant)
# clean col names
colnames(pd) <- gsub("ampwsc.", "", colnames(pd))
# multiply Amplified * 100 to get percent for Variant
pd$Amplified <- ifelse(pd$Amplified == 1, 100, 0)
table(pd$Variant)[order(table(pd$Variant), decreasing = T)]

famps <- names(table(pd$Variant))[which(table(pd$Variant) >= 100)]
# reclassify Adjectives - infreq. Adjectives are collapsed Amplifiedo category other
pd$Variant <- ifelse(pd$Variant  %in% famps, pd$Variant , "other")
# factorize variables
fctrs <- c("Function", "Variant")
pd[fctrs] <- lapply(pd[fctrs], as.factor)
# save levels of Age
agelevels <- rev(names(table(pd$Age)))
pd$Age <- ifelse(pd$Age == "16-19", 6,
                 ifelse(pd$Age == "20-29", 5,
                        ifelse(pd$Age == "30-39", 4, 
                               ifelse(pd$Age == "40-49", 3, 
                                      ifelse(pd$Age == "50-59", 2, 
                                             ifelse(pd$Age == "60+", 1, pd$Age))))))
pd$Age <- as.numeric(pd$Age)
# convert  Amplified into numeric variables
clnm <- c("Amplified")
pd[clnm] <- lapply(pd[clnm], as.numeric)
# create variables 
pd$other <- ifelse(pd$Variant == "other", 100, 0)
pd$pretty <- ifelse(pd$Variant == "pretty", 100, 0)
pd$really <- ifelse(pd$Variant == "really", 100, 0)
pd$so <- ifelse(pd$Variant == "so", 100, 0) 
pd$very <- ifelse(pd$Variant == "very", 100, 0)
pd$zero <- ifelse(pd$Variant == "0", 100, 0)
###############################################################
# p1
p1d <- pd
# start plot: Amplified
p1 <- ggplot(p1d, aes(x = jitter(Age), y = Amplified)) +
  geom_smooth(aes(y = Amplified), col = "gray30", lty = "solid", se = F, size = .5) +
  facet_grid(vars(Function)) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "Age", y = "Percent of Amplification") +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_x_discrete(name ="Age",
                   limits = c(1:length(agelevels)),
                   labels= agelevels) +
  ggsave(file = paste(imageDirectory,"Amplified_Function.png",sep="/"), 
         width = 15, height = 10, units = c("cm"),  dpi = 320)
p1

###############################################################
# p2
p2d <- pd
# remove non-amplified instances
p2d <- p2d[p2d$Amplified != 0,]
# start plot: all
p2 <- ggplot(p2d, aes(x = jitter(Age), y = very)) +
  facet_grid(vars(Function)) +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), se = F, size = .5) +
  geom_smooth(aes(y = pretty, color = "pretty", linetype = "pretty"), se = F, size = .5) +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), se = F, size = .5) +
  geom_smooth(aes(y = so, color = "so", linetype = "so"), se = F, size = .5) +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), se = F, size = .5) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dashed", "dotted", "longdash", "twodash", "solid"),
                        name="Variant",
                        breaks = c("other", "pretty", "really", "so", "very"), 
                        labels = c("other", "pretty", "really", "so", "very")) +
  scale_colour_manual(values=c("grey20",  "grey20", "grey20",  "grey20", "grey20"),
                      name="Variant", 
                      breaks=c("other", "pretty", "really", "so", "very"), 
                      labels = c("other", "pretty", "really", "so", "very")) +
  theme_bw(base_size = 15) +
  theme(legend.position="top", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE) +
  scale_x_discrete(name ="Age",
                   limits = 1:6,
                   labels= agelevels) +
  ggsave(file = paste(imageDirectory,"Variant_Function.png",sep="/"), 
         width = 15, height = 10, units = c("cm"),  dpi = 320)
p2

###############################################################
# recode variants (very, really, other)
pd$Variant <- ifelse(pd$Variant == "very", "very",
                     ifelse(pd$Variant == "really", "really",
                            ifelse(pd$Variant == "0", "0", "other")))
pd$pretty <- NULL
pd$so <- NULL
pd$other <- ifelse(pd$Amplified == 100 & pd$very == 0 & pd$really == 0, 100, 0)
###############################################################
# p3
p3d <- pd
# remove non-amplified instances
p3d <- p3d %>%
  dplyr::filter(Amplified != 0)
# start plot: all
p3 <- ggplot(p3d, aes(x = jitter(Age), y = very)) +
  facet_wrap(vars(Function)) +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), se = F, size = .5, method = "loess") +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), se = F, size = .5, method = "loess") +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), se = F, size = .5, method = "loess") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dotted", "dashed", "solid"),
                        name="Variant",
                        breaks = c("other",  "really",  "very"), 
                        labels = c("other",  "really",  "very")) +
  scale_colour_manual(values=c("grey20", "grey20", "grey20"),
                      name="Variant", 
                      breaks=c("other", "really", "very"), 
                      labels = c("other", "really",  "very")) +
  theme_bw(base_size = 15) +
  theme(legend.position="top", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent of Amplification") +
  guides(size = FALSE)+
  guides(alpha = FALSE) +
  scale_x_discrete(name ="Age",
                   limits = 1:6,
                   labels= agelevels) +
  ggsave(file = paste(imageDirectory,"Variant_Function.png",sep="/"), 
         width = 20, height = 10, units = c("cm"),  dpi = 320)
p3

###############################################################
#               PLOTTING LEXICAL DIVESRITY
# load library
table(ampwsc$Variant)[order(table(ampwsc$Variant), decreasing = T)]

famps <- names(table(ampwsc$Variant))[which(table(ampwsc$Variant) > 300)]
p4d <- ampwsc %>%
  dplyr::select(Age, Function, Adjective, Variant) %>%
  dplyr::filter(Variant != "0") %>%
  dplyr::mutate(VariantS = ifelse(Variant  %in% famps, Variant , "other")) %>%
  dplyr::group_by(Function, Age, VariantS) %>%
  dplyr::summarize(AdjectiveTypes = length(table(Adjective)), AllTypes = sum(table(Adjective))) %>%
  dplyr::mutate(LD = AdjectiveTypes/AllTypes) %>%
  dplyr::select(Age, Function, VariantS, LD) %>%
  tidyr::spread(VariantS, LD)
# create a numeric variable from Age
p4d$Age <- ifelse(p4d$Age == "16-19", 6,
                  ifelse(p4d$Age == "20-29", 5,
                         ifelse(p4d$Age == "30-39", 4, 
                                ifelse(p4d$Age == "40-49", 3, 
                                       ifelse(p4d$Age == "50-59", 2, 
                                              ifelse(p4d$Age == "60+", 1, p4d$Age))))))
p4d$Age <- as.numeric(p4d$Age)
p4d$Function <- as.factor(p4d$Function)
p4d

# start plot: Amplified
p4 <- ggplot(p4d, aes(x = Age, y = very)) +
  facet_wrap(vars(Function)) +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), se = F, size = .5, method = "loess") +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), se = F, size = .5, method = "loess") +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), se = F, size = .5, method = "loess") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dotted", "dashed", "solid"),
                        name="Variant",
                        breaks = c("other",  "really", "very"), 
                        labels = c("other",  "really",  "very")) +
  scale_colour_manual(values=c("grey20",    "grey20", "grey20"),
                      name="Variant", 
                      breaks=c("other",  "really",  "very"), 
                      labels = c("other",  "really",  "very")) +
  theme_set(theme_bw(base_size = 15)) +
  theme(legend.position="top", 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Age", y = "Lexical Diversity") +
  guides(size = FALSE)+
  guides(alpha = FALSE) +
  scale_x_discrete(name ="Age",
                   limits = 1:6,
                   labels= agelevels) +
  ggsave(file = paste(imageDirectory,"LD_VariantTime.png", sep="/"), 
         width = 15,  height = 10, units = c("cm"),  dpi = 320)
p4

##############################################################
#    COVARYING COLLEXEME ANALYSIS: CHANGE IN COLL. STRENGTH
ccd1 <- ampwsc %>%
  dplyr::filter(Variant != "0")
table(ccd1$Adjective)[order(table(ccd1$Adjective), decreasing = T)]

table(ccd1$Variant)[order(table(ccd1$Variant), decreasing = T)]

ntfrqadj <- names(table(ccd1$Adjective))[which(table(ccd1$Adjective) < 100)]
nfvrnt <- names(table(ccd1$Variant))[which(table(ccd1$Variant) < 300)]
ccd <- ampwsc %>%
  dplyr::select(Adjective, Function, Variant, Age) %>%
  dplyr::filter(Variant != "0") %>%
  dplyr::mutate(VariantS = ifelse(Variant %in% nfvrnt, "other", Variant)) %>%
  dplyr::mutate(AdjectiveS = ifelse(Adjective %in% ntfrqadj, "other", Adjective)) %>%
  dplyr::select(Age, Function, AdjectiveS, VariantS) %>%
  dplyr::group_by(Age, Function, AdjectiveS, VariantS) %>%
  dplyr::summarize(CollocationFrequency = n())
ccd <- ccd %>%
  dplyr::group_by(Age, Function) %>%
  dplyr::mutate(AdjectiveToken = sum(CollocationFrequency))
ccd <- ccd %>%
  dplyr::group_by(Age, Function, VariantS) %>%
  dplyr::mutate(VariantToken = sum(CollocationFrequency))
ccd <- ccd %>%
  dplyr::group_by(Age, Function, AdjectiveS) %>%
  dplyr::mutate(AmplifierTypeByAdjectiveType = sum(CollocationFrequency))
ccd <- ccd %>%
  dplyr::rename(Adjective = AdjectiveS, Variant = VariantS)
ccd <- ccd %>%
  dplyr::mutate(VariantAdjective = CollocationFrequency, 
                VariantOtherAdjective = (VariantToken-CollocationFrequency),
                OtherVariantAdjective = (AmplifierTypeByAdjectiveType-CollocationFrequency),
                OtherVariantOtherAdjective = AdjectiveToken-(VariantAdjective+VariantOtherAdjective+OtherVariantAdjective)) %>%
  dplyr::select(Age, Function, Variant, Adjective, VariantAdjective,
                VariantOtherAdjective, OtherVariantAdjective, 
                OtherVariantOtherAdjective) %>%
  dplyr::mutate(AmplifierFrequency = (VariantAdjective+VariantOtherAdjective),
                AdjectiveFrequency = (VariantAdjective+OtherVariantAdjective),
                TotalFrequency = (VariantAdjective+VariantOtherAdjective+
                                    OtherVariantAdjective+OtherVariantOtherAdjective))
ccd <- ccd %>%
  dplyr::mutate(Expected = AdjectiveFrequency*AmplifierFrequency/TotalFrequency,
                Probability = Expected/AmplifierFrequency,
                OddsRatio = as.vector(unlist(fisher.test(matrix(c(VariantAdjective, VariantOtherAdjective, OtherVariantAdjective, OtherVariantOtherAdjective), nrow = 2, byrow = T), simulate.p.value=TRUE)))[3],
                p = pbinom(VariantAdjective, AmplifierFrequency, Probability),
                Effect = round(abs(log(as.vector(unlist(p, 10)), 10)), 2))
# save data to disc
write.table(ccd, "datatables/ccd.txt", sep = "\t", row.names = F)
# continue data processing
ccd3 <- ccd %>%
  dplyr::rename(Observed = VariantAdjective) %>%
  dplyr::select(Age, Function, Variant, Adjective, 
                Observed, Expected, TotalFrequency, p, Effect) %>%
  dplyr::mutate(Type = ifelse(as.numeric(Observed) > as.numeric(Expected), "attr", "repel"),
                Corrected05 = 0.05/TotalFrequency,
                Corrected01 = 0.01/TotalFrequency,
                Corrected001 = 0.001/TotalFrequency,
                Bonferroni = ifelse(p <= Corrected001, "p<.001",
                                    ifelse(p <= Corrected01, "p<.01",
                                           ifelse(p <= Corrected05, "p<.05", "n.s."))))  %>%
  dplyr::select(Age, Function, Variant, Adjective, 
                Observed, Expected, Effect, p, Bonferroni)
# reshape data and replace missing values and include missing combinations
p5d <- ccd3 %>%
  dplyr::select(Age, Function, Variant, Adjective, Effect) %>%
  tidyr::spread(Variant, Effect) %>%
  #    dplyr::filter(Function == "Predicative") %>%
  tidyr::replace_na(list(so = 0, other = 0)) %>%
  tidyr::gather(Variant, Effect, other:very) %>%
  #  tidyr::spread(Age, Effect) %>%
  #  tidyr::replace_na(list(`16-19` = 0, other = 0)) %>%
  #  tidyr::gather(Age, Effect, `1675-1750`:`1851-1930`) %>%
  tidyr::spread(Variant, Effect)
# create a numeric variable from Age
p5d$Age <- ifelse(p5d$Age == "16-19", 6,
                  ifelse(p5d$Age == "20-29", 5,
                         ifelse(p5d$Age == "30-39", 4, 
                                ifelse(p5d$Age == "40-49", 3, 
                                       ifelse(p5d$Age == "50-59", 2, 
                                              ifelse(p5d$Age == "60+", 1, p5d$Age))))))
p5d$Age <- as.numeric(p5d$Age)
# factorize function and adjective
p5d$Function <- as.factor(p5d$Function)
p5d$Adjective <- as.factor(p5d$Adjective)
# inpsect data
p5d

# start plot
p5 <- ggplot(p5d, aes(x = Age, y = very)) +
  facet_grid(Adjective ~ Function, scales = "free") +
  geom_smooth(aes(y = other, color = "other", linetype = "other"), se = F, size = .5) +
  geom_smooth(aes(y = really, color = "really", linetype = "really"), se = F, size = .5) +
  geom_smooth(aes(y = very, color = "very", linetype = "very"), se = F, size = .5) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dotted", "dashed", "solid"),
                        name="Variant",
                        breaks = c("other",  "really",  "very"), 
                        labels = c("other",  "really",  "very")) +
  scale_colour_manual(values=c("grey20",  "grey20",  "grey20"),
                      name="Variant", 
                      breaks=c("other",  "really",  "very"), 
                      labels = c("other",  "really", "very")) +
  theme_set(theme_bw(base_size = 15)) +
  theme(legend.position="top", 
        plot.margin = unit(c(1,1,1,1), "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(3, "lines")) +
  #  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "Age", y = "Collocation Strength (LOG(p), 10)") +
  guides(size = FALSE)+
  guides(alpha = FALSE) +
  scale_x_discrete(name ="Age",
                   limits = 1:6,
                   labels= agelevels) +
  ggsave(file = paste(imageDirectory,"CVCLA_Function.png",sep="/"), 
         width = 20, height = 15, units = c("cm"),  dpi = 320)
p5

##############################################################
#    COVARYING COLLEXEME ANALYSIS: SIG. COLLOCATIONS
ccad1 <- ampwsc %>%
  dplyr::filter(Variant != "0")
ccad <- ampwsc %>%
  dplyr::select(Adjective, Function, Variant, Age) %>%
  dplyr::filter(Variant != "0") %>%
  dplyr::select(Age, Function, Adjective, Variant) %>%
  dplyr::group_by(Age, Function, Adjective, Variant) %>%
  dplyr::summarize(CollocationFrequency = n())
ccad <- ccad %>%
  dplyr::group_by(Age, Function) %>%
  dplyr::mutate(AdjectiveToken = sum(CollocationFrequency))
ccad <- ccad %>%
  dplyr::group_by(Age, Function, Variant) %>%
  dplyr::mutate(VariantToken = sum(CollocationFrequency))
ccad <- ccad %>%
  dplyr::group_by(Age, Function, Adjective) %>%
  dplyr::mutate(AmplifierTypeByAdjectiveType = sum(CollocationFrequency))
ccad <- ccad %>%
  dplyr::mutate(VariantAdjective = CollocationFrequency, 
                VariantOtherAdjective = (VariantToken-CollocationFrequency),
                OtherVariantAdjective = (AmplifierTypeByAdjectiveType-CollocationFrequency),
                OtherVariantOtherAdjective = AdjectiveToken-(VariantAdjective+VariantOtherAdjective+OtherVariantAdjective)) %>%
  dplyr::select(Age, Function, Variant, Adjective, VariantAdjective,
                VariantOtherAdjective, OtherVariantAdjective, 
                OtherVariantOtherAdjective) %>%
  dplyr::mutate(AmplifierFrequency = (VariantAdjective+VariantOtherAdjective),
                AdjectiveFrequency = (VariantAdjective+OtherVariantAdjective),
                TotalFrequency = (VariantAdjective+VariantOtherAdjective+
                                    OtherVariantAdjective+OtherVariantOtherAdjective))
ccad <- ccad %>%
  dplyr::mutate(Expected = AdjectiveFrequency*AmplifierFrequency/TotalFrequency,
                Probability = Expected/AmplifierFrequency,
                OddsRatio = as.vector(unlist(fisher.test(matrix(c(VariantAdjective, VariantOtherAdjective, OtherVariantAdjective, OtherVariantOtherAdjective), nrow = 2, byrow = T), simulate.p.value=TRUE)))[3],
                p = pbinom(VariantAdjective, AmplifierFrequency, Probability),
                Effect = round(abs(log(as.vector(unlist(p, 10)), 10)), 2))
# continue data processing
ccad3 <- ccad %>%
  dplyr::rename(Observed = VariantAdjective) %>%
  dplyr::select(Age, Function, Variant, Adjective, 
                Observed, Expected, TotalFrequency, p, Effect) %>%
  dplyr::mutate(Type = ifelse(as.numeric(Observed) > as.numeric(Expected), "attr", "repel"),
                Corrected05 = 0.05/TotalFrequency,
                Corrected01 = 0.01/TotalFrequency,
                Corrected001 = 0.001/TotalFrequency,
                Bonferroni = ifelse(p <= Corrected001, "p<.001",
                                    ifelse(p <= Corrected01, "p<.01",
                                           ifelse(p <= Corrected05, "p<.05", "n.s."))))  %>%
  dplyr::select(Age, Function, Variant, Adjective, 
                Observed, Expected, Effect, p, Bonferroni)

# extract significant collocations
ccad3$Bonferroni <- as.factor(ccad3$Bonferroni)
SignificantCollocationsWSC <- ccad3 %>%
  dplyr::filter(p < .05)
SignificantCollocationsWSC

write.table(SignificantCollocationsWSC, "datatables/SignificantCollocationsWSC.txt", sep = "\t", 
            row.names = F, quote = F, eol="\n")
###############################################################
#                  CHANGES IN Adjective FREQ
table(ampwsc$Adjective)[order(table(ampwsc$Adjective), decreasing = T)]

ntfrqadj <- names(table(ampwsc$Adjective))[which(table(ampwsc$Adjective) < 320)]

AdjectiveChangeData <- ampwsc %>%
  dplyr::select(Adjective, Function, Age, Amplified) %>%
  dplyr::mutate(AdjectiveS = ifelse(Adjective %in% ntfrqadj, "other", Adjective)) %>%
  dplyr::group_by(Age, Function, AdjectiveS) %>%
  dplyr::mutate(AdjectiveTypeFrequency = n()) %>%
  dplyr::select(Age, Function, Amplified, AdjectiveS, AdjectiveTypeFrequency)
AdjectiveChangeData <-  AdjectiveChangeData %>%
  dplyr::group_by(Age, Function) %>%
  dplyr::mutate(AdjectiveFrequency = n())
AdjectiveChangeData <- AdjectiveChangeData %>%
  dplyr::mutate(RelativeFrequency = AdjectiveTypeFrequency/AdjectiveFrequency*100) %>%
  dplyr::select(Age, Function, AdjectiveS, RelativeFrequency) %>%
  unique()
AdjectiveChangeData <- AdjectiveChangeData %>%
  tidyr::spread(AdjectiveS, RelativeFrequency)
AdjectiveChangeData$Age <- ifelse(AdjectiveChangeData$Age == "16-19", 6,
                                  ifelse(AdjectiveChangeData$Age == "20-29", 5,
                                         ifelse(AdjectiveChangeData$Age == "30-39", 4, 
                                                ifelse(AdjectiveChangeData$Age == "40-49", 3, 
                                                       ifelse(AdjectiveChangeData$Age == "50-59", 2, 
                                                              ifelse(AdjectiveChangeData$Age == "60+", 1, AdjectiveChangeData$Age))))))
AdjectiveChangeData$Age <- as.numeric(AdjectiveChangeData$Age)
AdjectiveChangeData$Function <- as.factor(AdjectiveChangeData$Function)
AdjectiveChangeData  

# save data to disc
write.table(AdjectiveChangeData, "datatables/AdjectiveChangeData.txt", sep = "\t", row.names = F)

# start plot: Adjective
p7 <- ggplot(AdjectiveChangeData, aes(x = Age, y = other)) +
  facet_grid(vars(Function)) +
  geom_smooth(aes(y = big, color = "big", lty = "big"), size=.5, se = F) +  
  geom_smooth(aes(y = first, color = "first", lty = "first"), size=.5, se = F) +
  geom_smooth(aes(y = good, color = "good", lty = "good"), size=.5, se = F) +
  geom_smooth(aes(y = last, color = "last", lty = "last"), size=.5, se = F) +
  geom_smooth(aes(y = little, color = "little", lty = "little"), size=.5, se = F) +
  geom_smooth(aes(y = new, color = "new", lty = "new"), size=.5, se = F) +
  geom_smooth(aes(y = other, color = "other", lty = "other"), size=.5, se = F) +
  geom_smooth(aes(y = right, color = "right", lty = "right"), size=.5, se = F) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dotted", "dotted", "dashed", "dashed",  "longdash", "longdash", "solid", "solid"),
                        name="Adjective",
                        breaks = c("big", "first", "good", "last", "little", "new", "other", "right"), 
                        labels = c("big", "first", "good", "last", "little", "new", "other", "right")) +
  scale_colour_manual(values=c("grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80"),
                      name="Adjective", 
                      breaks=c("big", "first", "good", "last", "little", "new", "other", "right"), 
                      labels = c("big", "first", "good", "last", "little", "new", "other", "right")) +
  theme_bw(base_size = 15) +
  theme(legend.position="top", 
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(3, "lines")) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(x = "Age", y = "Percent \n (of all adjectives)") +
  scale_x_discrete(name ="Age",
                   limits = 1:6,
                   labels= agelevels) +
  ggsave(file = paste(imageDirectory,"AdjectiveFrequency_Time_Function.png", sep="/"), 
         width = 15,  height = 12, units = c("cm"),  dpi = 320)
p7

###############################################################
#                  CHANGES IN Adjective FREQ
table(ampwsc$Adjective)[order(table(ampwsc$Adjective), decreasing = T)]

ntfrqadj <- names(table(ampwsc$Adjective))[which(table(ampwsc$Adjective) < 320)]

AmplifiedAdjectivesData <- ampwsc %>%
  dplyr::select(Adjective, Function, Age, Amplified) %>%
  dplyr::mutate(AdjectiveS = ifelse(Adjective %in% ntfrqadj, "other", Adjective)) %>%
  dplyr::group_by(Age, Function, AdjectiveS) %>%
  dplyr::select(Age, Function, Amplified, AdjectiveS) %>%
  dplyr::mutate(AdjectiveTypeFrequency = n(), AmplifiedN = sum(Amplified)) %>%
  dplyr::mutate(RateAmplification = AmplifiedN/AdjectiveTypeFrequency*100) %>%
  dplyr::select(Age, Function, AdjectiveS, RateAmplification) %>%
  unique()
AmplifiedAdjectivesData <- AmplifiedAdjectivesData %>%
  tidyr::spread(AdjectiveS, RateAmplification)
AmplifiedAdjectivesData$Age <- ifelse(AmplifiedAdjectivesData$Age == "16-19", 6,
                                      ifelse(AmplifiedAdjectivesData$Age == "20-29", 5,
                                             ifelse(AmplifiedAdjectivesData$Age == "30-39", 4, 
                                                    ifelse(AmplifiedAdjectivesData$Age == "40-49", 3, 
                                                           ifelse(AmplifiedAdjectivesData$Age == "50-59", 2, 
                                                                  ifelse(AmplifiedAdjectivesData$Age == "60+", 1, AmplifiedAdjectivesData$Age))))))
AmplifiedAdjectivesData$Age <- as.numeric(AmplifiedAdjectivesData$Age)
AmplifiedAdjectivesData$Function <- as.factor(AmplifiedAdjectivesData$Function)
AmplifiedAdjectivesData  

# save data to disc
write.table(AmplifiedAdjectivesData, "datatables/AmplifiedAdjectivesData.txt", sep = "\t", row.names = F)

# start plot: Adjective
p8 <- ggplot(AmplifiedAdjectivesData, aes(x = Age, y = other)) +
  facet_grid(vars(Function)) +
  geom_smooth(aes(y = big, color = "big", lty = "big"), size=.5, se = F) +  
  geom_smooth(aes(y = first, color = "first", lty = "first"), size=.5, se = F) +
  geom_smooth(aes(y = good, color = "good", lty = "good"), size=.5, se = F) +
  geom_smooth(aes(y = last, color = "last", lty = "last"), size=.5, se = F) +
  geom_smooth(aes(y = little, color = "little", lty = "little"), size=.5, se = F) +
  geom_smooth(aes(y = new, color = "new", lty = "new"), size=.5, se = F) +
  geom_smooth(aes(y = other, color = "other", lty = "other"), size=.5, se = F) +
  geom_smooth(aes(y = right, color = "right", lty = "right"), size=.5, se = F) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  scale_linetype_manual(values=c("dotted", "dotted", "dashed", "dashed",  "longdash", "longdash", "solid", "solid"),
                        name="Adjective",
                        breaks = c("big", "first", "good", "last", "little", "new", "other", "right"), 
                        labels = c("big", "first", "good", "last", "little", "new", "other", "right")) +
  scale_colour_manual(values=c("grey20", "grey80", "grey20", "grey80", "grey20", "grey80", "grey20", "grey80"),
                      name="Adjective", 
                      breaks=c("big", "first", "good", "last", "little", "new", "other", "right"), 
                      labels = c("big", "first", "good", "last", "little", "new", "other", "right")) +
  theme_bw(base_size = 15) +
  theme(legend.position="top", 
        plot.margin = unit(c(1,1,1,1), "cm"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.spacing.x = unit(3, "lines")) +
  coord_cartesian(ylim = c(0, 50)) +
  labs(x = "Age", y = "Percent (Amplification)") +
  scale_x_discrete(name ="Age",
                   limits = 1:6,
                   labels= agelevels) +
  ggsave(file = paste(imageDirectory,"AmplifiedAdjectives_Time_Function.png", sep="/"), 
         width = 15,  height = 12, units = c("cm"),  dpi = 320)
p8

###############################################################
#                  ADJECTIVE FREQUENCY
table(ampwsc$Adjective)[order(table(ampwsc$Adjective), decreasing = T)]

ntfrqadj <- names(table(ampwsc$Adjective))[which(table(ampwsc$Adjective) < 320)]

AdjectiveFrequencyData <- ampwsc %>%
  dplyr::select(Adjective, Function) %>%
  dplyr::mutate(AdjectiveS = ifelse(Adjective %in% ntfrqadj, "other", Adjective)) %>%
  dplyr::group_by(Function, AdjectiveS) %>%
  dplyr::select(Function, AdjectiveS) %>%
  dplyr::mutate(AdjectiveFrequency = n()) %>%
  unique() %>%
  arrange(-AdjectiveFrequency)
AdjectiveFrequencyData  

###########################################################################
#                  CREATE REGRESSION DATA SET
# clean data
ampwsc <- ampwsc %>%
  dplyr::select(-Token)  %>%
  dplyr::mutate(very = ifelse(Variant == "very", "very", "other"))  %>%
  dplyr::select(-Amplified) %>%
  mutate(Date=str_replace_all(Date, ".*/", "")) %>%
  mutate(Date=str_replace_all(Date, " ", "")) %>%
  mutate(Date=ifelse(Date == "", NA, Date)) %>%
  mutate(FileSpeaker=paste(File, Speaker, sep = "")) %>%
  dplyr::select(-Occupation, -AgeNumeric, -Date) %>%
  na.omit() 
# inspect data
nrow(ampwsc); str(ampwsc)

###############################################################
###             TABULARIZATION
###############################################################
# tb 1
Varianttbwsc <- table(ampwsc$Variant)
Varianttbwsc <- Varianttbwsc[order(table(ampwsc$Variant), decreasing = T)]
Variantnames <- as.vector(names(Varianttbwsc))
Variantn <- as.vector(Varianttbwsc)
Variantprcnt <- round(Variantn/sum(Variantn)*100, 2)
Variantprcnt2 <-  c(0, round(Variantn[2:length(Variantn)]/sum(Variantn[2:length(Variantn)])*100, 2))
tb1 <- data.frame(Variantnames, Variantn, Variantprcnt, Variantprcnt2)
colnames(tb1) <- c("Variant", "TokenFrequency", "PercentageSlots", "PercentAgeAmplifiedensifiers")
tb1 <- rbind(tb1, c("Total", sum(as.vector(tb1$TokenFrequency)), "", ""))
rownames(tb1) <- NULL
# inspect data
head(tb1)

# save data to disc
write.table(tb1, "datatables/Table1.txt", sep = "\t", row.names = F)
# WARNING: remove all non amplified rows
ampwsc <- ampwsc %>%
  dplyr::filter(Variant != "0") %>%
  dplyr::select(-Variant)

# create table
Table2 <- ampwsc %>%
  dplyr::select(Age, Gender, FileSpeaker, very)  %>%
  dplyr::group_by(Age, Gender) %>%
  dplyr::summarize(Speakers = length(names(table(FileSpeaker))),
                   AdjectiveSlots = n(),
                   very = table(very)[2]) %>%
  dplyr::mutate(Percent = round(very/AdjectiveSlots*100, 2))
Table2 <- rbind(as.data.frame(Table2), c("", "", sum(Table2$Speakers), 
                  sum(Table2$AdjectiveSlots), 
                  sum(Table2$very), round(mean(Table2$Percent), 2)))
Table2

# save data to disc
write.table(Table2, "datatables/Table2.txt", sep = "\t",  row.names = F, quote = F, eol="\n")
###########################################################################
# save data to disc
write.table(ampwsc, "datatables/ampwsc_regdat.txt", sep = "\t", 
            row.names = F, quote = F, eol="\n")
###########################################################################
#                   END PART 2
###########################################################################
