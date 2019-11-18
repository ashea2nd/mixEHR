library(plyr)

rm(list=ls(all=TRUE))

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId) # ehrFeatInfo

trainData <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"
ehrinp <- read.table(trainData, header=F)
colnames(ehrinp) <- c("patId","typeId","pheId","stateId","freq")
ehrinp <- arrange(ehrinp, patId, typeId, pheId, stateId)

ehrinp$bin <- ehrinp$freq > 0

stat <- aggregate(bin~typeId+pheId, ehrinp, sum)

stat <- arrange(stat, bin, decreasing=TRUE)

# stat_notes <- subset(stat, typeId==1)
# stat_notes_info <- merge(stat_notes, ehrFeatInfo$notes, sort = FALSE)

npat <- length(unique(ehrinp$patId))

min_pat <- 100 # 0.01*npat
max_pat <- 0.5*npat

min_pat2 <- 10
max_pat2 <- 500

topD <- 100 # top most observed phe/lab

stat_info_sel <- do.call(rbind, lapply(names(ehrFeatInfo), function(datatype) {
  
  print(datatype)
  
  stat <- subset(stat, typeId==which(names(ehrFeatInfo)==datatype))
  
  if(datatype=="notes") {
    stat <- subset(stat, bin > min_pat2 & bin < max_pat2)
  } else {
    stat <- subset(stat, bin > min_pat & bin < max_pat)  
  }
  
  stat_info <- merge(stat, ehrFeatInfo[[datatype]], sort = FALSE)
  
  if(datatype=="icd_cm") {
    stat_info[,c("typeId","pheId")]
  } else {
    head(stat_info[,c("typeId","pheId")], topD)  
  }
}))

print(table(stat_info_sel$typeId))

print(sum(table(stat_info_sel$typeId)))

outfile <- "/results/impmimic/impute_target_pheId.txt"

write.table(stat_info_sel[,c("typeId","pheId")], file=outfile, quote=F, row.names=F, col.names=F, sep=' ')

stat_info_sel_info <- lapply(split(stat_info_sel, stat_info_sel$typeId), function(x) {
  cbind(x, pheName=as.character(ehrFeatInfo[[as.numeric(x$typeId[1])]]$pheName[
    match(x$pheId, ehrFeatInfo[[as.numeric(x$typeId[1])]]$pheId)]))
})
stat_info_sel_info <- do.call(rbind, stat_info_sel_info[unique(stat_info_sel$typeId)]) # order typeId back

# print(sort(as.character(subset(stat_info_sel_info, typeId==2)$pheName))) # notes

print(lapply(split(stat_info_sel_info,stat_info_sel_info$typeId), function(x) sort(as.character(x$pheName))))











