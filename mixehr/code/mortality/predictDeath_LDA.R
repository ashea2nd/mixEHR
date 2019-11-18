library(topicmodels)
library(tm)
library(dplyr)
library(tidytext)

library(doParallel)
library(glmnet)
library(cvAUC)

rm(list=ls(all=TRUE))

source('/code/main/eval_fun.R')

topicNum <- as.numeric(eval(commandArgs(TRUE)[1]))

outdir <- "/results/mortality"
if(!dir.exists(outdir)) dir.create(outdir)

resrda <- sprintf("/results/predictDeath_LDA_test_K%s.RData", topicNum)

#### train on full data ####
dataFileTrain <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"

df <- read.table(dataFileTrain)

colnames(df) <- c("patId","typeId","pheId","stateId","freq")

df$termId <- sprintf("%s_%s_%s", df$typeId, df$pheId, df$stateId)

dtm <- cast_dtm(df, patId, termId, freq)

ldafit <- LDA(dtm, k = topicNum, control = list(seed = 1234))

#### output results
save(ldafit, file=resrda)

dataTrain <- posterior(ldafit, dtm)$topics

patrda <- "/data/Mimic/version_1_3/processed/PATIENTS.RData"
load(patrda) # age, adm

trainData <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"
ehrinp <- read.table(trainData)
subjId <- unique(ehrinp$V1)
trainDeadLabel <- adm_single$DEATHTIME[match(subjId, adm_single$SUBJECT_ID)]!=""

# cl <- makeCluster(detectCores())
cl <- makeCluster(8)
registerDoParallel(cl)

fitFull <- cv.glmnet(dataTrain, trainDeadLabel, family="binomial", 
                     alpha=1, type.measure="auc",
                     grouped=FALSE, standardize=FALSE, parallel=TRUE)

stopCluster(cl)

#### testing performance ####
dataFileTest <- "/data/Mimic/version_1_3/processed/mimic_testData_lastAdm.txt"
admFileTest <- "/data/Mimic/version_1_3/processed/mimic_testData_lastAdm_deathLabel.txt"

df_test <- read.table(dataFileTest)
admTest <- read.delim(admFileTest)

colnames(df_test) <- c("patId","typeId","pheId","stateId","freq")

df_test$termId <- sprintf("%s_%s_%s", df_test$typeId, df_test$pheId, df_test$stateId)

# df_test <- subset(df_test, termId %in% df$termId)

dtm_test <- cast_dtm(df_test, patId, termId, freq)

dataTest <- posterior(ldafit, dtm_test)$topics

pred_test <- predict(fitFull, dataTest, type="response", s="lambda.min")

perf_roc_test <- auroc(pred_test, admTest$dead)

print(perf_roc_test)

res_test <- cbind(pred_test, admTest$dead)

#### prospective predictions ####
dataFileTest_firstAdm <- "/data/Mimic/version_1_3/processed/mimic_testData_firstAdm.txt"
admFileTest_firstAdm <- "/data/Mimic/version_1_3/processed/mimic_testData_firstAdm_deathLabel.txt"

df_test_firstAdm <- read.table(dataFileTest_firstAdm)
admTest_firstAdm <- read.delim(admFileTest_firstAdm)

colnames(df_test_firstAdm) <- c("patId","typeId","pheId","stateId","freq")

df_test_firstAdm$termId <- sprintf("%s_%s_%s", df_test_firstAdm$typeId, df_test_firstAdm$pheId, df_test_firstAdm$stateId)

# df_test_firstAdm <- subset(df_test_firstAdm, termId %in% df$termId)

dtm_test <- cast_dtm(df_test_firstAdm, patId, termId, freq)

dataTest_firstAdm <- posterior(ldafit, dtm_test)$topics

pred_test_firstAdm <- predict(fitFull, dataTest_firstAdm, type="response", s="lambda.min")

adm_filein <- "/data/Mimic/version_1_3/ADMISSIONS.csv.gz"

admAll <- read.csv(adm_filein, stringsAsFactors = FALSE)

deathIdx <- which(admTest$dead)

admTest$daysBeforeDeath <- NA

admTest$daysBeforeDeath[deathIdx] <-
  as.Date(admAll$DEATHTIME[match(admTest$HADM_ID[deathIdx], admAll$HADM_ID)]) -
  as.Date(admAll$ADMITTIME[match(admTest_firstAdm$HADM_ID[deathIdx], admAll$HADM_ID)])

idx1 <- which(admTest$daysBeforeDeath >= 30 & admTest$daysBeforeDeath < 365/2 & admTest$dead)
idx0 <- which(!admTest$dead)
idx <- c(idx0, idx1)

perf_roc_prospect <- auroc(pred_test_firstAdm[idx], admTest$dead[idx])

print(perf_roc_prospect)

res_prospect <- cbind(pred_test_firstAdm[idx], admTest$dead[idx])

#### output results
save(fitFull, res_test, res_prospect, file=resrda)






