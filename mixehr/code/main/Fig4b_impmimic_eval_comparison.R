library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(openxlsx)

rm(list=ls(all=TRUE))

rda_flat <- "/data/precomputed_results/impmimic/impmimic_eval/datatype_flattened/impmimic_eval.RData"
load(rda_flat) # perf_target_info

perf_target_info_flat <- perf_target_info
perf_target_info_flat$modeltype <- "flattened"

rda <- "/data/precomputed_results/impmimic/impmimic_eval/impmimic_eval.RData"
load(rda) # perf_target_info
perf_target_info_mtmd <- perf_target_info
perf_target_info_mtmd$modeltype <- "multimodal"
fileout <- sub("impmimic_eval\\.RData","TableS3_impmimic_eval\\.xlsx",rda)

wb <- createWorkbook()
addWorksheet(wb, "multimodal")
writeData(wb, "multimodal", perf_target_info_mtmd)
addWorksheet(wb, "flattened")
writeData(wb, "flattened", perf_target_info_flat)
saveWorkbook(wb, fileout, overwrite=TRUE)

wcox_acc <- wilcox.test(perf_target_info_mtmd$acc, perf_target_info_flat$acc)
wcox_auroc <- wilcox.test(perf_target_info_mtmd$auroc, perf_target_info_flat$auroc)
wcox_auprc <- wilcox.test(perf_target_info_mtmd$auprc, perf_target_info_flat$auprc)

perf_target_info <- rbind(perf_target_info_mtmd, perf_target_info_flat)

# df <- melt(perf_target_info[,c("modeltype", "acc", "auroc", "auprc", "knn")], c("modeltype","knn"))
df <- melt(perf_target_info[,c("modeltype", "acc", "auroc", "auprc", "knn")], c("modeltype","knn"))

colnames(df)[c(3,4)] <- c("perf_type", "perf")

modeltype_ordered_by_auprc <- as.character(plyr::arrange(aggregate(perf~modeltype, subset(df,perf_type=="auprc"), median), 
                                                        perf, decreasing=F)$modeltype)
df$modeltype <- factor(df$modeltype, levels=modeltype_ordered_by_auprc)

evalDf <- rbind(data.frame(pval=sprintf("Wilcox test:\np < %s", scientific(wcox_acc$p.value)), 
                           perf_type="acc", modeltype="multimodal", perf=1),
                data.frame(pval=sprintf("Wilcox test:\np < %s", scientific(wcox_auroc$p.value)), 
                           perf_type="auroc", modeltype="multimodal", perf=0.9),
                data.frame(pval=sprintf("Wilcox test:\np < %s", scientific(wcox_auprc$p.value)), 
                           perf_type="auprc", modeltype="multimodal", perf=0.25))

evalDf_average <- aggregate(perf~modeltype+perf_type, df, median)

# overall performance
ggbox_overall <- ggplot(df, aes(x=modeltype, y=perf)) + 
  geom_boxplot() + facet_wrap(~perf_type, scales='free') + theme_bw() + 
  theme(axis.text.x=element_text(colour="black",size=14),
        axis.text.y=element_text(colour="black",size=14),
        strip.text.x = element_text(size = 16)) +
  xlab("") + ylab("5-Fold CV performance") +
  geom_text(aes(x=modeltype, y=perf, label=pval), data=evalDf, hjust=1.2) +
  geom_text(aes(x=modeltype, y=perf, label=percent(perf)), data=evalDf_average, vjust=1.3)

ggout <- "/results/Fig4b_multimodal_vs_flattned_overall.eps"

ggsave(ggout, ggbox_overall, width=8, height=6)

