library(cvTools)
library(plyr)
library(tidytext)
library(parallel)

rm(list=ls(all=TRUE))

metadata <- "/data/Mimic/version_1_3/processed/mimic_meta.txt"
ehrmeta <- read.table(metadata, header=F)
colnames(ehrmeta) <- c("typeId","pheId","stateCnt")

trainData <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"
ehrinp <- read.table(trainData, header=F)
colnames(ehrinp) <- c("patId","typeId","pheId","stateId","freq")
ehrinp <- arrange(ehrinp, patId, typeId, pheId, stateId)

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId) # ehrFeatInfo

nfolds <- 5

set.seed(1234)  # set seed for reproducibility

outdir <- "/results/impmimic"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

outfile_meta <- sprintf("%s/metainfo.txt", outdir)
write.table(ehrmeta, outfile_meta, col.names=F, row.names=F, sep=' ', quote=F)  

#### create folds ####
patId <- unique(ehrinp$patId)

set.seed(1234)  # set seed for reproducibility
cv_folds <- cvFolds(length(patId), nfolds)

cvIdx_list <- split(cv_folds$subsets, cv_folds$which)

# convert from patient idx to patient id
cvIdx_list <- lapply(cvIdx_list, function(cvidx) patId[cvidx])

tmp <- lapply(1:nfolds, function(k) {
  
  #### training set ####
  ehrinp_train <- subset(ehrinp, patId %in% unlist(cvIdx_list[-k]))
  ehrinp_train <- arrange(ehrinp_train, patId, typeId, pheId, stateId)
  
  outfile_trainInp <- sprintf("%s/train%s.txt", outdir, k)
  write.table(ehrinp_train, outfile_trainInp, col.names=F, row.names=F, sep=' ', quote=F)
  
  #### test set ####
  outfile_testInp <- sprintf("%s/test%s.txt", outdir, k)
  
  ehrinp_test <- subset(ehrinp, patId %in% cvIdx_list[[k]])
  ehrinp_test <- arrange(ehrinp_test, patId, typeId, pheId, stateId)
  
  write.table(ehrinp_test, outfile_testInp, col.names=F, row.names=F, sep=' ', quote=F)
})










































