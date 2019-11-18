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

# from select_impute_target_code.R
inpfile <- "/results/impmimic/impute_target_pheId.txt"
stat_info_sel <- read.table(inpfile)

colnames(stat_info_sel) <- c("typeId", "pheId")

target_code_matchidx <- match(sprintf("%s_%s", stat_info_sel$typeId, stat_info_sel$pheId), 
                        sprintf("%s_%s", metainfo$typeId, metainfo$pheId))

stat_info_sel0 <- metainfo[target_code_matchidx, c("typeId0","pheId0")]

outfile <- "/results/impmimic/datatype_flattened/impute_target_pheId.txt"

write.table(stat_info_sel0, file=outfile, quote=F, row.names=F, col.names=F, sep=' ')

