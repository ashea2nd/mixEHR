library(randomForest)
library(doParallel)
library(cvAUC)
library(glmnet)

rm(list=ls(all=TRUE))

K <- as.numeric(eval(commandArgs(TRUE)[1]))

#### cross validation ####
inpdir <- "/results/mortality"

#### train on full data ####
inpdir <- "/data/Mimic/version_1_3/processed"
resdir <- "/results"

finalMetaphe_filein <- list.files(resdir, sprintf("mimic_trainData_mimic_trainData_JCVB0_nmar_K%s", K), full.names = TRUE)

niter <- sub("iter","",unlist(strsplit(basename(finalMetaphe_filein), "_"))[8])
niter <- as.numeric(niter)

metapheTrain <- read.csv(finalMetaphe_filein, header=F)

patrda <- "/data/Mimic/version_1_3/processed/PATIENTS.RData"
load(patrda) # age, adm

trainData <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"
ehrinp <- read.table(trainData)

subjId <- unique(ehrinp$V1)
trainDeadLabel <- adm_single$DEATHTIME[match(subjId, adm_single$SUBJECT_ID)]!=""

cl <- makeCluster(detectCores())
registerDoParallel(cl)

fitFull <- cv.glmnet(as.matrix(metapheTrain), trainDeadLabel, family="binomial", 
                     alpha=1, type.measure="auc",
                     grouped=FALSE, standardize=FALSE, parallel=TRUE)

stopCluster(cl)

#### prospective predictions ####
dataFileTest_firstAdm <- sprintf("%s/mimic_testData_lastAdm_mimic_trainData_JCVB0_nmar_K%s_iter%s_metaphe.csv", 
                                 inpdir, K, niter)
admFileTest_firstAdm <- "/data/Mimic/version_1_3/processed/mimic_testData_firstAdm_deathLabel.txt"

dataTest_firstAdm <- read.csv(dataFileTest_firstAdm, header=F)
admTest_firstAdm <- read.delim(admFileTest_firstAdm)

pred_test_firstAdm <- predict(fitFull, as.matrix(dataTest_firstAdm), type="response", s="lambda.min")

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
outdir <- "/results/mortality"
if(!dir.exists(outdir)) dir.create(outdir)

resrda <- sprintf("%s/predictDeath_mixehr_K%s.RData", outdir, K)

save(fitFull, res_prospect, file=resrda)



























