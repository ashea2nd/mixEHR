library(ggplot2)
library(reshape2)
library(plyr)

rm(list=ls(all=TRUE))

source('/code/main/eval_fun.R')

inpdir <- "/data/results/impmimic"

knn <- 100

# perf_target_info <- lapply(c(10, 25, 50, 100, 200, 300, 400, 500), function(knn) {

perf_target_info <- lapply(100, function(knn) {
  
  print(knn)

  # target code files
  target_phe_pred_fileins <- list.files(inpdir, "target_phe_pred.csv.gz", recursive = TRUE, full.names = TRUE)
  target_phe_true_fileins <- list.files(inpdir, "target_phe_true.csv.gz", recursive = TRUE, full.names = TRUE)
  
  target_phe_pred_fileins <- grep("datatype_flattened", target_phe_pred_fileins, invert=TRUE, value=TRUE)
  target_phe_true_fileins <- grep("datatype_flattened", target_phe_true_fileins, invert=TRUE, value=TRUE)
  
  # lab obs files
  target_lab_obs_pred_fileins <- list.files(inpdir, "target_lab_obs_pred.csv.gz", recursive = TRUE, full.names = TRUE)
  target_lab_obs_true_fileins <- list.files(inpdir, "target_lab_obs_true.csv.gz", recursive = TRUE, full.names = TRUE)
  
  target_lab_obs_pred_fileins <- grep("datatype_flattened", target_lab_obs_pred_fileins, invert=TRUE, value=TRUE)
  target_lab_obs_true_fileins <- grep("datatype_flattened", target_lab_obs_true_fileins, invert=TRUE, value=TRUE)
  
  
  # filter files by knn setting
  target_phe_pred_fileins <- grep(sprintf("_knn%s/", knn), target_phe_pred_fileins, value=TRUE)
  target_phe_true_fileins <- grep(sprintf("_knn%s/", knn), target_phe_true_fileins, value=TRUE)
  target_lab_obs_pred_fileins <- grep(sprintf("_knn%s/", knn), target_lab_obs_pred_fileins, value=TRUE)
  target_lab_obs_true_fileins <- grep(sprintf("_knn%s/", knn), target_lab_obs_true_fileins, value=TRUE)
  
  # read in prediction and true labels files
  target_phe_pred <- do.call(rbind, lapply(target_phe_pred_fileins, read.csv, header=F))
  target_phe_true <- do.call(rbind, lapply(target_phe_true_fileins, read.csv, header=F))
  target_lab_obs_pred <- do.call(rbind, lapply(target_lab_obs_pred_fileins, read.csv, header=F))
  target_lab_obs_true <- do.call(rbind, lapply(target_lab_obs_true_fileins, read.csv, header=F))
  
  # combine phe and lab pred and true labels
  target_obs_pred <- cbind(target_phe_pred, target_lab_obs_pred)
  target_obs_true <- cbind(target_phe_true, target_lab_obs_true)
  
  target_pheId_fileins <- list.files(inpdir, "target_pheid.csv.gz", recursive = TRUE, full.names = TRUE)
  target_labId_fileins <- list.files(inpdir, "target_labid.csv.gz", recursive = TRUE, full.names = TRUE)
  target_patId_fileins <- list.files(inpdir, "target_patid.csv.gz", recursive = TRUE, full.names = TRUE)

  target_pheId_fileins <- grep("datatype_flattened", target_pheId_fileins, invert=TRUE, value=TRUE)
  target_labId_fileins <- grep("datatype_flattened", target_labId_fileins, invert=TRUE, value=TRUE)
  target_patId_fileins <- grep("datatype_flattened", target_patId_fileins, invert=TRUE, value=TRUE)
  
  
  # filter files by knn setting
  target_pheId_fileins <- grep(sprintf("_knn%s/", knn), target_pheId_fileins, value=TRUE)
  target_labId_fileins <- grep(sprintf("_knn%s/", knn), target_labId_fileins, value=TRUE)
  target_patId_fileins <- grep(sprintf("_knn%s/", knn), target_patId_fileins, value=TRUE)
  
  target_pheId <- unique(do.call(rbind, lapply(target_pheId_fileins, read.csv, header=F)))
  target_labId <- unique(do.call(rbind, lapply(target_labId_fileins, read.csv, header=F)))
  target_patId <- unique(as.character(as.matrix(do.call(rbind, lapply(target_patId_fileins, read.csv, header=F)))))
  
  #### DEBUG BEGIN ####
  # print(lapply(target_pheId_fileins, function(x) dim(read.csv(x, header=F))))
  # target_labId <- lapply(target_labId_fileins, read.csv, header=F)
  #### DEBUG END ####
  
  target_obs_id <- rbind(target_pheId, target_labId) # combine phe and lab IDs
  
  featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
  load(featId) # ehrFeatInfo
  
  colnames(target_obs_id) <- c("typeId", "pheId")
  
  target_obs_info <- lapply(split(target_obs_id, target_obs_id$typeId), function(x) {
    
    cbind(x, pheName=as.character(ehrFeatInfo[[as.numeric(x$typeId[1])]]$pheName[
      match(x$pheId, ehrFeatInfo[[as.numeric(x$typeId[1])]]$pheId)]))
  })
  
  target_obs_info <- do.call(rbind, target_obs_info[unique(target_obs_id$typeId)]) # order typeId back
  
  dimnames(target_obs_pred) <- dimnames(target_obs_true) <- 
    list(target_patId, sprintf("%s:%s", target_obs_info$typeId, target_obs_info$pheId))
  
  aurocPerf <- sapply(colnames(target_obs_pred), function(j) {
    auroc(target_obs_pred[,j], target_obs_true[,j]>0)
  })
  
  auprcPerf <- sapply(colnames(target_obs_pred), function(j) {
    auprc(target_obs_pred[,j], target_obs_true[,j]>0)
  })
  
  k <- knn
  
  accuPerf <- sapply(colnames(target_obs_pred), function(j) {
    sum((target_obs_pred[,j]>1/k)==target_obs_true[,j])/length(target_obs_true[,j])
  })
  
  perf_target_info <- target_obs_info
  
  perf_target_info$acc <- accuPerf
  
  perf_target_info$auroc <- aurocPerf
  
  perf_target_info$auprc <- auprcPerf
  
  perf_target_info$datatype <- names(ehrFeatInfo)[perf_target_info$typeId]
  
  datatype_ordered_by_auroc <- as.character(plyr::arrange(aggregate(auroc~datatype, perf_target_info, mean), auroc, decreasing=TRUE)$datatype)
  datatype_ordered_by_auprc <- as.character(plyr::arrange(aggregate(auprc~datatype, perf_target_info, mean), auprc, decreasing=TRUE)$datatype)
  perf_target_info$datatype <- factor(perf_target_info$datatype, levels=datatype_ordered_by_auroc)
  
  perf_target_info$knn <- as.factor(knn)
  
  perf_target_info
})

perf_target_info <- do.call(rbind, perf_target_info)

df <- melt(perf_target_info[,c("datatype", "auroc", "auprc", "knn")], c("datatype","knn"))

colnames(df)[c(3,4)] <- c("perf_type", "perf")

datatype_ordered_by_auprc <- as.character(plyr::arrange(aggregate(perf~datatype, subset(df,perf_type=="auprc"), median), 
                                                        perf, decreasing=F)$datatype)
df$datatype <- factor(df$datatype, levels=datatype_ordered_by_auprc)

# overall performance
ggbox_overall <- ggplot(df, aes(x=knn, y=perf)) + 
  geom_boxplot() + facet_wrap(~perf_type, scales='free') + theme_bw() + 
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.text.y=element_text(colour="black",size=14),
        strip.text.x = element_text(size = 16)) +
  xlab("") + ylab("5-Fold CV performance")

# ggout <- sprintf("/results/impmimic_eval_box_overall.eps", outdir)
# ggsave(ggout, ggbox_overall, width=3.2, height=2.4)

ggbox_datatypes <- ggplot(df, aes(x=datatype, y=perf, fill=datatype)) + 
  geom_boxplot(outlier.shape=NA) + 
  facet_grid(perf_type~., scales='free') + theme_bw() + 
  theme(axis.text.x=element_text(colour="black",size=14,angle=30),
        axis.text.y=element_text(colour="black",size=14),
        strip.text.x = element_text(size = 16)) +
  xlab("Data type") + ylab("5-Fold CV performance")


ggout <- "/results/FigS17_impmimic_eval_box_datatypes.eps"
ggsave(ggout, ggbox_datatypes, width=4.8, height=4.8)


#### icd-9 grouping ####
icd_grouping_filein <- "/data/Mimic/version_1_3/processed/metadata/icd9groupings_processed.txt"
icd_grouping <- read.delim(icd_grouping_filein, header=F)

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId) # ehrFeatInfo

colnames(icd_grouping) <- c("level1", "level2", "level3")

perf_target_info_icd <- subset(perf_target_info, datatype=="icd_cm")

perf_target_info_icd <- merge(perf_target_info_icd, ehrFeatInfo$icd_cm)


groupIdx <- sapply(as.character(perf_target_info_icd$ICD9_CODE), function(icd) {
  
  if(substr(icd,1,1)=="E") {
    match(substr(icd, 1, 4), as.character(icd_grouping$level3)) 
  } else {
    match(substr(icd, 1, 3), as.character(icd_grouping$level3))
  }
})

perf_target_info_icd <- cbind(perf_target_info_icd, icd_grouping[groupIdx,])

level1_ordered_by_auroc <- as.character(plyr::arrange(aggregate(auroc~level1, perf_target_info_icd, mean), auroc, decreasing=F)$level1)
level1_ordered_by_auprc <- as.character(plyr::arrange(aggregate(auprc~level1, perf_target_info_icd, mean), auprc, decreasing=F)$level1)

perf_target_info_icd$level1 <- factor(perf_target_info_icd$level1, levels=level1_ordered_by_auprc)

df <- melt(perf_target_info_icd[,c("level1", "auroc", "auprc", "knn")], c("level1","knn"))

colnames(df)[c(3,4)] <- c("perf_type", "perf")

df$knn <- as.factor(df$knn)

longname <- "Endocrine, nutritional and metabolic diseases, and immunity disorders (240-279)"
levels(df$level1)[which(levels(df$level1)==longname)] <- "Endocrine, nutritional and metabolic diseases,\nimmunity disorders (240-279)"

longname <- "Diseases of the blood and blood-forming organs (280-289)"
levels(df$level1)[which(levels(df$level1)==longname)] <- "Diseases of the blood\nand blood-forming organs (280-289)"

longname <- "Supplementary classification of factors influencing health status and contact with health services (V)"
levels(df$level1)[which(levels(df$level1)==longname)] <- "Supplementary classification of factors (V)"

ggbox <- ggplot(df, aes(x=level1, y=perf)) + geom_boxplot() + 
  coord_flip() + facet_grid(.~perf_type, scales='free') + theme_bw() + 
  theme(axis.text.x=element_text(colour="black"),
        axis.text.y=element_text(colour="black",size=16),
        strip.text.x = element_text(size = 16)) +
  xlab("ICD9 group") + ylab("5-Fold CV performance")

ggout <- "/results/FigS18_impmimic_eval_box_icd9group.eps"

ggsave(ggout, ggbox, width=14, height=8)

# save(perf_target_info, file="/results/impmimic/impmimic_eval/impmimic_eval.RData")





















