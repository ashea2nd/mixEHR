library(ggplot2)
library(reshape2)

rm(list=ls(all=TRUE))

source('/code/main/eval_fun.R')

inpdir <- "/data/results/impmimic"

knn <- 100 # 25

# perf_target_info <- lapply(c(10, 25, 50, 100, 200, 300, 400, 500), function(knn) {
perf_target_info <- lapply(100, function(knn) {
  
  print(knn)
  
  # lab obs files
  target_lab_res_pred_fileins <- list.files(inpdir, "target_lab_res_pred.csv.gz", recursive = TRUE, full.names = TRUE)
  target_lab_res_true_fileins <- list.files(inpdir, "target_lab_res_true.csv.gz", recursive = TRUE, full.names = TRUE)
  
  target_lab_res_pred_fileins <- grep("datatype_flattened", target_lab_res_pred_fileins, invert=TRUE, value=TRUE)
  target_lab_res_true_fileins <- grep("datatype_flattened", target_lab_res_true_fileins, invert=TRUE, value=TRUE)
  
  # filter files by knn setting
  target_lab_res_pred_fileins <- grep(sprintf("_knn%s/", knn), target_lab_res_pred_fileins, value=TRUE)
  target_lab_res_true_fileins <- grep(sprintf("_knn%s/", knn), target_lab_res_true_fileins, value=TRUE)
  
  target_lab_obs_true_fileins <- list.files(inpdir, "target_lab_obs_true.csv.gz", recursive = TRUE, full.names = TRUE)
  target_lab_obs_true_fileins <- grep("datatype_flattened", target_lab_obs_true_fileins, invert=TRUE, value=TRUE)
  target_lab_obs_true_fileins <- grep(sprintf("_knn%s/", knn), target_lab_obs_true_fileins, value=TRUE)
  target_lab_obs_true <- do.call(rbind, lapply(target_lab_obs_true_fileins, read.csv, header=F))
  
  # read in prediction and true labels files
  target_lab_res_pred <- do.call(rbind, lapply(target_lab_res_pred_fileins, read.csv, header=F))
  target_lab_res_true <- do.call(rbind, lapply(target_lab_res_true_fileins, read.csv, header=F))
  
  target_labId_fileins <- list.files(inpdir, "target_labid.csv.gz", recursive = TRUE, full.names = TRUE)
  target_patId_fileins <- list.files(inpdir, "target_patid.csv.gz", recursive = TRUE, full.names = TRUE)
  
  target_labId_fileins <- grep("datatype_flattened", target_labId_fileins, invert=TRUE, value=TRUE)
  target_patId_fileins <- grep("datatype_flattened", target_patId_fileins, invert=TRUE, value=TRUE)
  
  # filter files by knn setting
  target_labId_fileins <- grep(sprintf("_knn%s/", knn), target_labId_fileins, value=TRUE)
  target_patId_fileins <- grep(sprintf("_knn%s/", knn), target_patId_fileins, value=TRUE)
  
  target_labId <- unique(do.call(rbind, lapply(target_labId_fileins, read.csv, header=F)))
  target_patId <- as.character(as.matrix(do.call(rbind, lapply(target_patId_fileins, read.csv, header=F))))
  
  featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
  load(featId) # ehrFeatInfo
  
  colnames(target_labId) <- c("typeId", "labId")
  
  labName <- as.character(ehrFeatInfo$lab$pheName[match(target_labId$labId, ehrFeatInfo$lab$pheId)])
  
  target_labInfo <- cbind(target_labId, labName=labName)
  
  dimnames(target_lab_res_pred) <- dimnames(target_lab_res_true) <- dimnames(target_lab_obs_true) <- 
    list(target_patId, target_labInfo$labName)
  
  accPerf <- sapply(colnames(target_lab_res_pred), function(l) {
    obsIdx <- which(target_lab_obs_true[,l]==1)
    sum(target_lab_res_pred[obsIdx,l] == target_lab_res_true[obsIdx,l])/length(obsIdx)
  })
  
  target_labInfo$perf <- accPerf
  
  target_labInfo$knn <- as.factor(knn)
  
  target_labInfo
})

perf_target_info <- do.call(rbind, perf_target_info)

df <- melt(perf_target_info[,c("perf","knn")], "knn")

colnames(df)[c(2,3)] <- c("accuracy", "perf")

ggbox <- ggplot(df, aes(x=knn, y=perf)) + geom_boxplot() + 
  theme_bw() + geom_point() + geom_jitter() +
  theme(axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        axis.text.y=element_text(colour="black")) +
  xlab("") + ylab("Prediction accuracy") + 
  ggtitle(sprintf("Lab results prediction\n(Median %s)", percent(median(df$perf))))


ggout <- "/results/FigS16_impmimic_eval_labres.eps"
ggsave(ggout, ggbox, width=2.4, height=2.4)





























