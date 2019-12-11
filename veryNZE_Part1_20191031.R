##################################################################
# Titel:      On the waning of forms - a corpus-based analysis of 
#             decline and loss in adjective amplification
# Part:       1
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
library(stringr)
# set options
options(stringsAsFactors = F)
options(scipen = 999)
# define image directory
imageDirectory<-"images"
# specify path to corpra
wscpath <- "D:\\Uni\\Korpora\\Original\\WSC"
biowscpath <- "D:\\Uni\\Korpora\\Original\\WSC\\DOC/PTPNTS.txt"
bioguidewsc <- "D:\\Uni\\Korpora\\Original\\WSC\\DOC/GUIDE.txt"
biolinkwsc <- "D:\\Uni\\Korpora\\Original\\WSC\\DOC/LINK.txt"
bioextractswsc <- "D:\\Uni\\Korpora\\Original\\WSC\\DOC/EXTRACTS.txt"
# define corpus files
corpusfiles = list.files(path = wscpath, pattern = ".TXT", all.files = T,
                         full.names = T, recursive = F, ignore.case = T, 
                         include.dirs = F)
###############################################################
# load and start processing corpus
wsc <- sapply(corpusfiles, function(x) {
  x <- scan(x, what = "char", sep = "", quote = "", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  x <- paste(x, sep = " ", collapse = " ")
  x <- strsplit(gsub("(<WSC#)", "~~~\\1", x), "~~~" )
  x <- as.vector(unlist(x))
  x <- x[2:length(x)]
} ) 
# extract number of text elements per file
idx <- sapply(wsc, function(x) length(x))
# extract file names
File <- names(wsc)
File <- rep(File, idx)
File <- gsub("D:\\Uni\\Korpora\\Original\\WSC/", "", File, fixed = T)
File <- gsub(".TXT", "", File)
# extract speakers
Speaker <- as.vector(unlist(sapply(wsc, function(x) {
  x <- gsub(">.*", "", x)
  x <- gsub(".*#", "", x)
  x <- gsub(".*:", "", x)
  })))
# extract raw text
Text <- as.vector(unlist(sapply(wsc, function(x) {
  x <- gsub("<[A-Z]{1,}#[A-Z]{1,}[0-9]{1,}:[0-9]{1,}:[A-Z]{1,}[0-9]{0,}>", "", x)
  x <- str_trim(x, side = "both")
  })))
# clean text
CleanText <- Text %>%
  str_replace_all("<O> {0,1}[[:alnum:] ]{0,}</O>", " ") %>%
  str_replace_all("<&> {0,1}[[:alnum:] ]{0,}[:punct:]{0,}[[:alnum:] ]{0,}</&>", " ") %>%
  str_replace_all("</{0,1}[:alnum:]{0,}>", " ") %>%
  str_replace_all("<[:punct:]{1,2}[:alnum:]{0,2}>", " ") %>%
  str_replace_all("<indig=[:alnum:]{1,}> {0,1}[[:alnum:] ]{0,}</indig=[:alnum:]{1,}>", " ") %>%
  str_replace_all(fixed("\""), " ") %>%
  str_replace_all("<foreign=[:alnum:]{1,}> {0,1}[[:alnum:] ]{0,}</foreign=[:alnum:]{1,}>", " ") %>%
  str_replace_all("<.*>", " ") %>%
  tolower() %>%
  str_replace_all(" {2,}", " ") %>%
  str_trim(side = "both")
# combine elements into a data frame 
wscdf <- data.frame(File, Speaker, Text, CleanText)
# testing
t1 <- wscdf$CleanText[grep("\\[", wscdf$CleanText)] #test
#t1

# remove problematic elements
#wscdf <- wscdf[-grep("\\[", wscdf$CleanText),]
# remove empty speech units
wscdf <- wscdf[wscdf$CleanText != "",]
nrow(wscdf)

###############################################################
# save raw data to disc
write.table(wscdf, "datatables/wscdf_raw.txt", sep = "\t", row.names = F, col.names = T)
wscdf <- read.delim("datatables/wscdf_raw.txt", sep = "\t", header = T, skipNul = T)
###############################################################
# split data into smaller chunks
pos01 <- wscdf$CleanText[1:10000]
pos02 <- wscdf$CleanText[10001:20000]
pos03 <- wscdf$CleanText[20001:22000]
pos04 <- wscdf$CleanText[22001:23000]
pos05 <- wscdf$CleanText[23001:24000]
pos06 <- wscdf$CleanText[24001:25000]
pos07 <- wscdf$CleanText[25001:30000]
pos08 <- wscdf$CleanText[30001:40000]
pos09 <- wscdf$CleanText[40001:50000]
pos10 <- wscdf$CleanText[50001:60000]
pos11 <- wscdf$CleanText[60001:70000]
pos12 <- wscdf$CleanText[70001:80000]
pos13 <- wscdf$CleanText[80001:90000]
pos14 <- wscdf$CleanText[90001:100000]
pos15 <- wscdf$CleanText[100001:nrow(wscdf)]
# reload libraries
source("D:\\R/POStagObject.R") # for pos-tagging objects in R
library(NLP)
library(openNLP)
library(openNLPmodels.en)
# pos tagging data
wscpos01 <- POStag(object = pos01)
wscpos01 <- as.vector(unlist(wscpos01))
writeLines(wscpos01, con = "datatables/wscpos01.txt", sep = "\n", useBytes = FALSE)
# chunk 2
wscpos02 <- POStag(object = pos02)
wscpos02 <- as.vector(unlist(wscpos02))
writeLines(wscpos02, con = "datatables/wscpos02.txt", sep = "\n", useBytes = FALSE)
# chunk 03
wscpos03 <- POStag(object = pos03)
wscpos03 <- as.vector(unlist(wscpos03))
writeLines(wscpos03, con = "datatables/wscpos03.txt", sep = "\n", useBytes = FALSE)
# chunk 04
wscpos04 <- POStag(object = pos04)
wscpos04 <- as.vector(unlist(wscpos04))
writeLines(wscpos04, con = "datatables/wscpos04.txt", sep = "\n", useBytes = FALSE)
# chunk 05
wscpos05 <- POStag(object = pos05)
wscpos05 <- as.vector(unlist(wscpos05))
writeLines(wscpos05, con = "datatables/wscpos05.txt", sep = "\n", useBytes = FALSE)
# chunk 06
wscpos06 <- POStag(object = pos06)
wscpos06 <- as.vector(unlist(wscpos06))
writeLines(wscpos06, con = "datatables/wscpos06.txt", sep = "\n", useBytes = FALSE)
# chunk 07
wscpos07 <- POStag(object = pos07)
wscpos07 <- as.vector(unlist(wscpos07))
writeLines(wscpos07, con = "datatables/wscpos07.txt", sep = "\n", useBytes = FALSE)
# chunk 08
wscpos08 <- POStag(object = pos08)
wscpos08 <- as.vector(unlist(wscpos08))
writeLines(wscpos08, con = "datatables/wscpos08.txt", sep = "\n", useBytes = FALSE)
# chunk 09
wscpos09 <- POStag(object = pos09)
wscpos09 <- as.vector(unlist(wscpos09))
writeLines(wscpos09, con = "datatables/wscpos09.txt", sep = "\n", useBytes = FALSE)
# chunk 10
wscpos10 <- POStag(object = pos10)
wscpos10 <- as.vector(unlist(wscpos10))
writeLines(wscpos10, con = "datatables/wscpos10.txt", sep = "\n", useBytes = FALSE)
# chunk 11
wscpos11 <- POStag(object = pos11)
wscpos11 <- as.vector(unlist(wscpos11))
writeLines(wscpos11, con = "datatables/wscpos11.txt", sep = "\n", useBytes = FALSE)
# chunk 12
wscpos12 <- POStag(object = pos12)
wscpos12 <- as.vector(unlist(wscpos12))
writeLines(wscpos12, con = "datatables/wscpos12.txt", sep = "\n", useBytes = FALSE)
# chunk 13
wscpos13 <- POStag(object = pos13)
wscpos13 <- as.vector(unlist(wscpos13))
writeLines(wscpos13, con = "datatables/wscpos13.txt", sep = "\n", useBytes = FALSE)
# chunk 14
wscpos14 <- POStag(object = pos14)
wscpos14 <- as.vector(unlist(wscpos14))
writeLines(wscpos14, con = "datatables/wscpos14.txt", sep = "\n", useBytes = FALSE)
# chunk 15
wscpos15 <- POStag(object = pos15)
wscpos15 <- as.vector(unlist(wscpos15))
writeLines(wscpos15, con = "datatables/wscpos15.txt", sep = "\n", useBytes = FALSE)
# list pos tagged elements
postag.files = c("datatables/wscpos01.txt", "datatables/wscpos02.txt", 
                 "datatables/wscpos03.txt", "datatables/wscpos04.txt", 
                 "datatables/wscpos05.txt", "datatables/wscpos06.txt",  
                 "datatables/wscpos07.txt", "datatables/wscpos08.txt", 
                 "datatables/wscpos09.txt",  "datatables/wscpos10.txt",
                 "datatables/wscpos11.txt", "datatables/wscpos12.txt",
                 "datatables/wscpos13.txt", "datatables/wscpos14.txt",
                 "datatables/wscpos15.txt")
# load pos tagged elements
wscpos <- sapply(postag.files, function(x) {
  x <- scan(x, what = "char", sep = "\n", quote = "", quiet = T, skipNul = T)
  x <- gsub(" {2,}", " ", x)
  x <- str_trim(x, side = "both")
  x <- str_replace_all(x, fixed("\n"), " ")
})
# unlist pos tagged elements
wscdf$TextPOS <- unlist(wscpos)
###############################################################
# save raw data to disc
write.table(wscdf, "datatables/wscdf_postagged.txt", sep = "\t", row.names = F, col.names = T)
wscdf <- read.delim("datatables/wscdf_postagged.txt", sep = "\t", header = T, skipNul = T)
###############################################################
# extract number of adjs per line
pstggd <- wscdf$TextPOS
lpstggd <- strsplit(pstggd, " ")
nlpstggd <- as.vector(unlist(sapply(lpstggd, function(x){
  x <- x[grep("[A-Z]{0,1}[a-z]{1,}\\/JJ[A-Z]{0,1}", x)]
  x <- length(x) } )))
rp <- nlpstggd
rp <- ifelse(rp == 0, 1, rp)
# detach dplyr package (clash with plyr)
detach("package:dplyr", unload=TRUE)
# load function for concordancing
source("D:\\R/ConcR_2.3_loadedfiles.R")
# set parameters for concordancing
pattern <- "[A-Z]{0,1}[a-z]{1,}\\/JJ[A-Z]{0,1}"
context <- 50
# extract all adjectives (concordance)
concjjwsc <- ConcR(wscdf$TextPOS, pattern, context, all.pre = FALSE)
# repeat rows in data frame as often as there are adjectives in it 
# (if 0 adj, repeat once)
wscadjdf <- wscdf[rep(seq(nrow(wscdf)), rp),]
# combine data sets
wscadjdf <- data.frame(1:nrow(wscadjdf), wscadjdf, concjjwsc)
# remove rows without Tokens
wscadjdf <- wscadjdf[is.na(wscadjdf$Token) == F,]
# add clean column names
colnames(wscadjdf)[1] <- "ID"
# clean adjectives
wscadjdf$Adjective <- str_replace_all(wscadjdf$Token, "/.*", "")
# add Variant column
wscadjdf$Variant <- gsub(".* ", "", str_trim(wscadjdf$PreContext, side = "both")) 
# inspect data
#nrow(wscadjdf); head(wscadjdf)

###############################################################
# save raw data to disc
write.table(wscadjdf, "datatables/wscadjdf.txt", sep = "\t", row.names = F, col.names = T)
wscadjdf <- read.delim("datatables/wscadjdf.txt", sep = "\t", header = T, skipNul = T)
###############################################################
# load files
biowsc <- read.delim(biowscpath, sep = "\t", header = F, skipNul = T)
guidewsc <- read.delim(bioguidewsc, sep = "\t", header = F, skipNul = T)
linkwsc <- read.delim(biolinkwsc, sep = "\t", header = F, skipNul = T)
extractswsc <- read.delim(bioextractswsc, sep = "\t", header = F, skipNul = T)
# add column names
biowsc <- biowsc %>%
  dplyr::rename(Id = V1, Ethnicity = V2, Gender = V3, Age = V4, 
                Occupation = V5, Education = V6, L1 = V7) %>%
  mutate(Gender=replace(Gender, Gender=="F", "Woman")) %>%
  mutate(Gender=replace(Gender, Gender=="M", "Man"))
# add column names to guide
linkwsc <- linkwsc %>%
  dplyr::rename(File = V1, Id = V2, Speaker = V3, Words = V4)
# add column names to extracts
extractswsc <- extractswsc %>%
  dplyr::rename(File = V1, Words = V2, Date = V3, Topic = V4, 
                ExtractMinutes = V5, Minutes = V6)
# join biowsc and linkwsc
wscspeakerinformation <- join(biowsc, linkwsc, by = "Id")
wscspeakerinformation <- join(wscspeakerinformation, extractswsc, by = "File")
# add infromation about type
wscspeakerinformation$Type <- wscspeakerinformation$File %>%
  str_replace_all("[0-9]", "") %>%
  str_replace_all("MSN", "BroadcastNews") %>%
  str_replace_all("MST", "BroadcastMonologue") %>%
  str_replace_all("MSW", "BroadcastWeather") %>%
  str_replace_all("MUC", "SportsCommentary") %>%
  str_replace_all("MUJ", "JudgesSummation") %>%
  str_replace_all("MUL", "Lecture") %>%
  str_replace_all("MUS", "TeacherMonologue") %>%
  str_replace_all("DPC", "Conversation") %>%
  str_replace_all("DPF", "TelephoneConversation") %>%
  str_replace_all("DPH", "OralHistoryInterview") %>%
  str_replace_all("DPP", "SocialDialectInterview") %>%
  str_replace_all("DGB", "RadioTalkback") %>%
  str_replace_all("DGI", "BroadcastInterview") %>%
  str_replace_all("DGU", "ParliamentaryDebate") %>%
  str_replace_all("DGZ", "TransactionsAndMeetings")
# inspect speaker information
head(wscspeakerinformation); nrow(wscspeakerinformation)

###############################################################
# save raw data to disc
write.table(wscspeakerinformation, "datatables/wscspeakerinformation.txt", sep = "\t", row.names = F, col.names = T)
wscspeakerinformation <- read.delim("datatables/wscspeakerinformation.txt", sep = "\t", header = T, skipNul = T)
wscadjdf <- read.delim("datatables/wscadjdf.txt", sep = "\t", header = T, skipNul = T)
###############################################################
# join wscadjdf and wscspeakerinformation
wsc <- join(wscadjdf, wscspeakerinformation, by = c("File", "Speaker"))
################################################################
# code freq of adj type by date category
frqadjtb <- table(wsc$Age, wsc$Adjective)
relfreqadjtb <- round(prop.table(frqadjtb, margin = 1)*100, 5)
relfreqadjdf <- as.data.frame(relfreqadjtb)
colnames(relfreqadjdf)[1:2] <- c("Age", "Adjective")
# add freq by date to data
wscadjdf <- join(wsc, relfreqadjdf, by=c("Age", "Adjective"))
# relabel Freq
colnames(wscadjdf)[which(colnames(wscadjdf) == "Freq")] <- "Frequency"
# inspect data
head(wscadjdf); nrow(wscadjdf)

###############################################################
# save raw data to disc
write.table(wscadjdf, "datatables/wscadjdf_spk.txt", sep = "\t", row.names = F, col.names = T)
wscadjdf <- read.delim("datatables/wscadjdf_spk.txt", sep = "\t", header = T, skipNul = T)
###############################################################
# define amplifiers
amplifiers <- c("absolutely", "actually", "aggressively", 
                "amazingly", "appallingly", "awful", "awfully", 
                "badly", "bloody", "certainly", "clearly",
                "complete", "dead", "completely", "considerably", 
                "crazy", "decidedly", "definitely",  "distinctly", 
                "dreadfully", "enormously", "entirely", "especially", 
                "exactly", "exceedingly", "exceptionally", 
                "excruciatingly", "extraordinarily", "extremely",
                "fiercely", "firmly", "frightfully", "fucking", 
                "fully", "genuinely", "greatly",
                "grossly", "heavily", "highly", "hopelessly", 
                "horrendously", "hugely",
                "immediately", "immensely", "incredibly", 
                "infinitely", "intensely", "irrevocably",
                "mad", "mega", "mighty", "most", "much", 
                "obviously", "openly", "overwhelmingly", "particularly", 
                "perfectly", "plenty", "positively", "precisely", 
                "pretty", "profoundly", "purely", 
                #"quite", 
                "real", "really", "remarkably", "seriously", 
                "shocking",   "significant", "significantly", "so", 
                "specially", "specifically", "strikingly",
                "strongly", "substantially", "super", "surely", 
                "terribly", "terrifically", 
                #"too",
                "total", "totally", "traditionally", "true", 
                "truly", "ultra", "utterly", "very",
                "viciously", 
                #"well", 
                "wholly", "wicked", "wildly")
# determine function
Function1 <- wscadjdf$PostContext %>%
  str_trim(side = "both") %>%
  tolower() %>%
  str_replace_all(" {2,}", " ") %>%
  str_replace_all(" .*", "") %>%
  str_replace_all(" .*", "") %>%
  str_replace_all(".*/n.*", "Attributive") %>%
  str_replace_all(fixed("/."), "PUNCT") %>%
  str_replace_all(".*PUNCT.*", "Attributive")
Function2 <-  wscadjdf$PostContext %>%
  str_trim(side = "both") %>%
  tolower() %>%
  str_replace(" ", "") %>%
  str_replace_all(" .*", "") %>%
  str_replace_all(".*/n.*", "Attributive") %>%
  str_replace_all(fixed("/."), "PUNCT") %>%
  str_replace_all(".*PUNCT.*", "Attributive")
functiontb <- cbind(Function1, Function2)
wscadjdf$Function <-  as.vector(unlist(apply(functiontb, 1, function(x){
  x <-   ifelse(x[1] == "" | x[1] == "Attributive" | x[2] == "Attributive", "Attributive", "Predicative")
})))
# shorten post Context
wscadjdf$PostContext <- substr(wscadjdf$PostContext, 1, ifelse((nchar(wscadjdf$PostContext)+25) <25, max(nchar(wscadjdf$PostContext)), 25))
# pre Context
wscadjdf$PreContext <- str_trim(wscadjdf$PreContext, side = "both")
wscadjdf$PreContextLong <- wscadjdf$PreContext
wscadjdf$PreContextLong <- substr(wscadjdf$PreContextLong, ifelse(nchar(wscadjdf$PreContextLong)-25 <=0, 1, 
                                                              nchar(wscadjdf$PreContextLong)-25), nchar(wscadjdf$PreContextLong))
wscadjdf$PreContext <- gsub(".* ", "", wscadjdf$PreContext)
# amplifier variant
wscadjdf$PreContext <- gsub("\\/.*", "", wscadjdf$PreContext)
wscadjdf$Variant <- ifelse(wscadjdf$PreContext %in% amplifiers, wscadjdf$PreContext, "0")
# amplified y/n
wscadjdf$Amplified <- ifelse(wscadjdf$Variant == "0", 0, 1) 
# adjective
wscadjdf$Adjective <- tolower(wscadjdf$Adjective)
# inspect data
nrow(wscadjdf); head(wscadjdf); table(wscadjdf$Variant)

# code priming
prim1 <- c(rep(0, 1), wscadjdf$Variant[1:length(wscadjdf$Variant)-1])
prim2 <- c(rep(0, 2), wscadjdf$Variant[1:(length(wscadjdf$Variant)-2)])
prim3 <- c(rep(0, 3), wscadjdf$Variant[1:(length(wscadjdf$Variant)-3)])
primtb <- cbind(wscadjdf$Variant, prim1, prim2, prim3)

wscadjdf$Priming <- as.vector(unlist(apply(primtb, 1, function(x){
  x <- ifelse(x[1]== "0" , "NoPrime",
              ifelse(x[1] == x[2] | x[1] == x[3] | x[1] == x[4], "Prime", "NoPrime"))
})))
# define forms that require removal
sups <- c(".*most.*", ".*more.*") 
negs <- c(".*not.*", ".*never.*", ".*n't.*")
downtoners <- c(".*sort/.*", ".*kind/.*", ".* bit/.*", ".*somewhat.*", ".*fairly.*", 
                ".*rather.*", ".*reasonably.*", ".*slightly.*", ".*comparatively.*", ".*semi.*", 
                ".*relatively.*", ".*little.*", ".*somehow.*", ".*almost.*", ".*partly.*", 
                ".*hardly.*", ".* less.*", ".*barely.*", ".* just/.*")
specialforms <- c(".* too.*", ".*quite.*")
PostContextdowntoners <- c(".*enough.*")
nonpropadj <- c("only", "much", "many", "cheaper", "cheaperr", "bests", "larger",
                "bst", "better", "bigger")
# check length of dataset
str(wscadjdf); head(wscadjdf); nrow(wscadjdf)

# find items to be removed
supsidx <- unique(grep(paste(sups,collapse="|"), wscadjdf$PreContextLong, value=F))
negsidx <- unique(grep(paste(negs,collapse="|"), wscadjdf$PreContextLong, value=F))
downtonersidx <- unique(grep(paste(downtoners,collapse="|"), wscadjdf$PreContextLong, value=F))
specialformsidx <- unique(grep(paste(specialforms,collapse="|"), wscadjdf$PreContextLong, value=F))
PostContextdowntonersidx <- unique(grep(paste(PostContextdowntoners,collapse="|"), wscadjdf$PostContext, value=F))
nonpropadjidx <- unique(grep(paste(nonpropadj,collapse="|"), wscadjdf$Adjective, value=F))
# combine indices
idxs <- unique(c(supsidx, negsidx, downtonersidx, specialformsidx, PostContextdowntonersidx, nonpropadjidx))
# remove forms that require removal
wscadjdf <- wscadjdf[-idxs,]
# remove empty values
wscadjdf <- wscadjdf[!wscadjdf$Variant == "", ]
###############################################################
# save raw data to disc
write.table(wscadjdf, "datatables/wscadjdf_wo_neg.txt", sep = "\t", row.names = F)
wscadjdf <- read.delim("datatables/wscadjdf_wo_neg.txt", sep = "\t", header = T, skipNul = T)
###############################################################
# remove items that were not intensified by a minimum of 2 intensifier variants
pintadjtb <- table(wscadjdf$Adjective, wscadjdf$Variant)
#pintadjtb <- pintadjtb[2:nrow(pintadjtb),]
pintadjtb <- pintadjtb[,2:ncol(pintadjtb)]
pintadjtb2 <- apply(pintadjtb, 1, function(x){
  x <- ifelse(x > 1, 1, x)})
pintadjtb3 <- colSums(pintadjtb2)
pintadjschildes <- names(pintadjtb3)[which(pintadjtb3 >= 2)]
wscadjdf <- wscadjdf[wscadjdf$Adjective %in% pintadjschildes, ]
nrow(wscadjdf)

# inspect adjectives
names(table(wscadjdf$Adjective))

# create vector with false adjectives
rmvadj <- c("cos", "di", "er", "g", "helen", "i", "l", "lower", "more", 
            "most", "nah", "ohakune", "overall", "pete", "super", "t", 
            "tara", "tony", "ultra", "um", "un", "w", "wh", "yep")
wscadjdf$remove <- ifelse(wscadjdf$Adjective %in% rmvadj, "remove", wscadjdf$Adjective)
wscadjdf <- wscadjdf[wscadjdf$remove != "remove",]
wscadjdf$remove <- NULL
# inspecta data
nrow(wscadjdf); length(table(wscadjdf$Adjective)); head(wscadjdf)

###############################################################
# save raw data to disc
write.table(wscadjdf, "datatables/wscadjdf_clean.txt", sep = "\t", row.names = F)
wscadjdf <- read.delim("datatables/wscadjdf_clean.txt", sep = "\t", header = T, skipNul = T)
###############################################################
# code gradability
# load Gradability data (derived from COCA)
gradability <- read.delim("D:\\Uni\\Projekte\\09-GradabilityOfAdjectives/Gradability.txt", sep = "\t", header = T, quote = "", skipNul = T)
wscadjdf$Gradability <- ifelse(wscadjdf$Adjective %in% gradability$Adjective, gradability$Beta, 1)
# inspect data
nrow(wscadjdf); head(wscadjdf)

###############################################################
# add semantic types (tagliamonte 2008, based on dixon 1977)
# dimension = semdim (e.g. big, large, little, small, long, short, wide, narrow, thick)
# difficulty = semdif (e.g. difficult, simple)
# physical property = (e.g. hard, soft, heavy, light, rough, smooth, hot, sweet)
# color = semcol (e.g. black, white, red)
# human propensity: semhup (e.g. jealous, happy, kind, clever, generous, gay, rude)
# age = semage (e.g. new, young, old) 
# value (e.g. good, bad, proper, perfect, excellent, delicious, poor), 
# speed Speed (fast, quick, slow)
# position (e.g. right, left, near, far)
# other
# load data
code1 <- read.delim("datatables/semcodecg1.txt", sep = "\t", header = T, skipNul = T)
code2 <- read.delim("datatables/semcodedjm1.txt", sep = "\t", header = T, skipNul = T)
code3 <- read.delim("datatables/semcodedl1.txt", sep = "\t", header = T, skipNul = T)
code4 <- read.delim("datatables/semcodedm1.txt", sep = "\t", header = T, skipNul = T)
code5 <- read.delim("datatables/semcodedma1.txt", sep = "\t", header = T, skipNul = T)
code6 <- read.delim("datatables/semcodedv1.txt", sep = "\t", header = T, skipNul = T)
code7 <- read.delim("datatables/semcodedjw1.txt", sep = "\t", header = T, skipNul = T)
# order data sets
code1 <- code1[order(code1$Id),]
code2 <- code2[order(code2$Id),]
code3 <- code3[order(code3$Id),]
code4 <- code4[order(code4$Id),]
code5 <- code5[order(code5$Id),]
code6 <- code6[order(code6$Id),]
code7 <- code6[order(code7$Id),]
# repair adjectives in code1
code3$Adjective <- code1$Adjective
# combine tables
semcode <- rbind(code1, code2, code3, code4, code5, code6, code7)
# convert coding into numeric values
semcode[,3:12] <- t(apply(semcode[,3:12], 1, function(x) {
#  x <- ifelse(x == "" | is.na(x) == T, 0, 1)}))
  x <- ifelse(x == "" | x == "?"| is.na(x) == T, 0, 1)}))
# convert into data frame
semcode <- as.data.frame(semcode)
# add column names
colnames(semcode)[3:12] <- c("Dimension", "Difficulty", "PhysicalProperty", "Color",
                             "HumanPropensity", "Age", "Value", "Speed", "Position", 
                             "Other")
# load library
library(dplyr)
AdjectiveSemantics <- semcode %>%
  dplyr::group_by(Adjective) %>%
  na.omit() %>%
  dplyr::summarize(Dimension = sum(Dimension), Difficulty = sum(Difficulty),
                PhysicalProperty = sum(PhysicalProperty), Color = sum(Color),
                HumanPropensity = sum(HumanPropensity), Age = sum(Age),
                Value = sum(Value), Speed = sum(Speed),
                Position = sum(Position), Other = sum(Other)) %>%
  dplyr::mutate(OverallScore = rowSums(.[,2:11])) %>%
  dplyr::mutate(Maximum = do.call(pmax, (.[,2:11]))) %>%
  dplyr::mutate(Certainty = Maximum/OverallScore*100)
AdjectiveSemantics <- AdjectiveSemantics %>%
  dplyr::mutate(Category = colnames(AdjectiveSemantics[2:11])[apply(AdjectiveSemantics[2:11],1,which.max)])
# inspect interrater reliability
summary(AdjectiveSemantics$Certainty)

# create vectors
Age <- as.vector(unlist(AdjectiveSemantics %>%
                   dplyr::filter(Category == "Age") %>% select(Adjective)))
Color <- as.vector(unlist(AdjectiveSemantics %>%
                          dplyr::filter(Category == "Color") %>% select(Adjective)))
Difficulty <- as.vector(unlist(AdjectiveSemantics %>%
                          dplyr::filter(Category == "Difficulty") %>% select(Adjective)))
Dimension <- as.vector(unlist(AdjectiveSemantics %>%
                          dplyr::filter(Category == "Dimension") %>% select(Adjective)))
HumanPropensity <- as.vector(unlist(AdjectiveSemantics %>%
                          dplyr::filter(Category == "HumanPropensity") %>% select(Adjective)))
PhysicalProperty <- as.vector(unlist(AdjectiveSemantics %>%
                          dplyr::filter(Category == "PhysicalProperty") %>% select(Adjective)))
Position <- as.vector(unlist(AdjectiveSemantics %>%
                                       dplyr::filter(Category == "Position") %>% select(Adjective)))
Speed <- as.vector(unlist(AdjectiveSemantics %>%
                                       dplyr::filter(Category == "Speed") %>% select(Adjective)))
Value <- as.vector(unlist(AdjectiveSemantics %>%
                            dplyr::filter(Category == "Value") %>% select(Adjective)))
# add semantic category to data
wscadjdf$SemanticCategory <- ifelse(wscadjdf$Adjective %in% Age, "Age",
                                    ifelse(wscadjdf$Adjective %in% Color, "Color",
                                    ifelse(wscadjdf$Adjective %in% Difficulty, "Difficulty",
                                    ifelse(wscadjdf$Adjective %in% Dimension, "Dimension",
                                    ifelse(wscadjdf$Adjective %in% HumanPropensity, "HumanPropensity",
                                    ifelse(wscadjdf$Adjective %in% PhysicalProperty, "PhysicalProperty",
                                    ifelse(wscadjdf$Adjective %in% Position, "Position",
                                    ifelse(wscadjdf$Adjective %in% Speed, "Speed",
                                    ifelse(wscadjdf$Adjective %in% Value, "Value", "Other")))))))))
# table sem class of tokens
table(wscadjdf$SemanticCategory)

# inspect data
head(wscadjdf)

###############################################################
# code emotion
library(syuzhet)
class_emo <- get_nrc_sentiment(wscadjdf$Adjective)

# process sentiment
wscadjdf$Emotionality <- as.vector(unlist(apply(class_emo, 1, function(x){
  x <- ifelse(x[9] == 1, "NegativeEmotional",
              ifelse(x[10] == 1, "PositiveEmotional", "NonEmotional")) } )))
# revert order of factor Emotionality
wscadjdf$Emotionality <- factor(wscadjdf$Emotionality, levels = c("NonEmotional", "NegativeEmotional", "PositiveEmotional"))
# inspect data
head(wscadjdf); str(wscadjdf); nrow(wscadjdf)

###############################################################
# save raw data to disc
write.table(wscadjdf, "datatables/wsc_fullclean.txt", sep = "\t", row.names = F, col.names = T)
###############################################################
#                        END PART 1
###############################################################
