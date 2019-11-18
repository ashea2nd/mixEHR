library(plyr)
library(parallel)

rm(list=ls(all=TRUE))

# purpose of this script: convert type-specific data to a single type data
inpdir <- "/results/impmimic"

outdir <- "/results/impmimic/datatype_flattened"
if(!dir.exists(outdir)) dir.create(outdir)

inpfile_meta <- sprintf("%s/metainfo.txt", inpdir)

metainfo <- read.table(inpfile_meta, header=F)
colnames(metainfo) <- c("typeId", "pheId", "stateCnt")

# make metainfo have only one disctint typeID to "flatten" the data type
metainfo$typeId0 <- 1
metainfo$pheId0 <- 1:nrow(metainfo)

nfolds <- 5 # use already created folds for equal comparison

options(mc.cores=detectCores())

tmp <- mclapply(1:nfolds, function(k) {
  
  print(k)
  
  # use already created folds for equal comparison
  inpfile_train <- sprintf("%s/train%s.txt", inpdir, k)
  inpfile_test <- sprintf("%s/test%s.txt", inpdir, k)
  
  train_fold <- read.table(inpfile_train, header=F)
  test_fold <- read.table(inpfile_test, header=F)
  
  colnames(train_fold) <- colnames(test_fold) <- c("patId","typeId","pheId","stateId","freq")
  
  train_matchidx <- match(sprintf("%s_%s", train_fold$typeId, train_fold$pheId), 
                          sprintf("%s_%s", metainfo$typeId, metainfo$pheId))
  
  train_fold$typeId <- metainfo$typeId0[train_matchidx]
  train_fold$pheId <- metainfo$pheId0[train_matchidx]
  
  test_matchidx <- match(sprintf("%s_%s", test_fold$typeId, test_fold$pheId), 
                         sprintf("%s_%s", metainfo$typeId, metainfo$pheId))
  
  test_fold$typeId <- metainfo$typeId0[test_matchidx]
  test_fold$pheId <- metainfo$pheId0[test_matchidx]
  
  # for lab test use only the lab observation not results
  train_fold$stateId <- 1
  
  selcols <- c("patId", "typeId", "pheId", "stateId", "freq")
  
  train_fold <- aggregate(freq~patId+typeId+stateId+pheId, train_fold, sum)[,selcols]
  test_fold <- aggregate(freq~patId+typeId+stateId+pheId, test_fold, sum)[,selcols]
  
  train_fold <- plyr::arrange(train_fold, patId, typeId, pheId, stateId, freq)[,selcols]
  test_fold <- plyr::arrange(test_fold, patId, typeId, pheId, stateId, freq)[,selcols]
  
  outfile_train <- sprintf("%s/train%s.txt", outdir, k)
  outfile_test <- sprintf("%s/test%s.txt", outdir, k)
  
  write.table(train_fold, outfile_train, col.names=F, row.names=F, sep=' ', quote=F)
  write.table(test_fold, outfile_test, col.names=F, row.names=F, sep=' ', quote=F)
})


metainfo$stateCnt <- 1
metainfo$typeId <- metainfo$typeId0
metainfo$pheId <- metainfo$pheId0
metainfo$typeId0 <- metainfo$pheId0 <- NULL

outfile_meta <- sprintf("%s/metainfo.txt", outdir)
write.table(metainfo, outfile_meta, col.names=F, row.names=F, sep=' ', quote=F)  













