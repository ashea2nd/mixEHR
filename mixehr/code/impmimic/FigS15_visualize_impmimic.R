library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(dendsort)

rm(list=ls(all=TRUE))

source('/code/main/eval_fun.R')

inpdir <- "/data/results/impmimic"

knn <- 100

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
target_patId <- as.character(as.matrix(do.call(rbind, lapply(target_patId_fileins, read.csv, header=F))))

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
  list(target_patId, sprintf("%s:%s:%s", target_obs_info$typeId, target_obs_info$pheId, target_obs_info$pheName))

# from impmimic_eval.R
rda <- "/data/results/impmimic/impmimic_eval/impmimic_eval.RData"
load(rda) # perf_target_info

perf_target_info <- plyr::arrange(perf_target_info, auprc, decreasing=TRUE)

topD <- 1e3

topD_pheId <- head(sprintf("%s:%s:%s", perf_target_info$typeId, perf_target_info$pheId, perf_target_info$pheName), topD)

# colidx <- grep(sprintf("^%s:",typeId), colnames(target_obs_true))
colidx <- match(topD_pheId, colnames(target_obs_true))


set.seed(1234)
randN <- 5e3
patidx <- sample(1:nrow(target_obs_true), randN)

pred_x <- target_obs_pred[patidx, colidx] > 1/knn
true_x <- target_obs_true[patidx, colidx]

colnames(pred_x) <- colnames(true_x) <- substr(colnames(pred_x), 1, 20)

dend_rows <- dendsort(hclust(dist(pred_x)))
dend_cols <- dendsort(hclust(dist(t(pred_x))))

hm_pred <- Heatmap(pred_x, name='Predicted', 
                   col = colorRamp2(c(0, 1), c("white", "red")),
                   row_dend_reorder = FALSE,
                   column_dend_reorder = FALSE,
                   cluster_columns = dend_cols,
                   show_row_names= FALSE,
                   show_column_names= FALSE,
                   cluster_rows = dend_rows)

hm_true <- Heatmap(true_x, name='True labels', 
                   # row_order=unlist(row_order(hm_pred)),
                   # column_order=column_order(hm_pred),
                   col = colorRamp2(c(0, 1), c("white", "red")),
                   row_dend_reorder = FALSE,
                   column_dend_reorder = FALSE,
                   show_column_names= FALSE,
                   show_row_names= FALSE,
                   cluster_columns = dend_cols,
                   cluster_rows = dend_rows)

hmpout <- "/results/FigS15_obs_vs_imp.png"

png(hmpout, width=3300, height=1300)

pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
draw(hm_true, newpage=FALSE)
upViewport()

pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
draw(hm_pred, newpage=FALSE)
upViewport()

dev.off()

# sum((target_obs_pred[, colidx]>1/knn) == target_obs_true[, colidx])/prod(dim(target_obs_true[, colidx]))
print(sum((target_obs_pred>1/knn) == target_obs_true)/prod(dim(target_obs_true)))
























