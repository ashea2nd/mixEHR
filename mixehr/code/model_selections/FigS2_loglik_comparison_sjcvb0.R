library(ggplot2)
library(plyr)
library(reshape)

#### JCVB0 ####
# from eval_summary.R
inpdir <- "/data/results/model_selections/mixmimic_sjcvb0"
perfDf_jcvb0 <- read.table(sprintf("%s/eval_summary_fullbatch.txt", inpdir), header=T)
perfDf_jcvb0 <- melt(as.matrix(perfDf_jcvb0))
perfDf_jcvb0$method <- "JCVB0"

perfDf_jcvb0$B <- 4e4
perfDf_jcvb0 <- perfDf_jcvb0[,c("X1","X2","B","value","method")]
colnames(perfDf_jcvb0)[c(1,2)] <- c("Iter", "K")

#### SJCVB0 ####
inpdir <- "/data/results/model_selections/mixmimic_sjcvb0"
fileins <- list.files(inpdir, "mimic_trainData_logPredLik_SJCVB0",full.names = TRUE, recursive = TRUE)
fileins <- grep("K150_B", fileins, value=TRUE)
perf <- do.call(cbind, lapply(fileins, scan))
dimnames(perf) <- list(1:nrow(perf), sub("mixmimic_sjcvb0_K150_", "", basename(dirname(fileins))))
perfDf_sjcvb0 <- melt(as.matrix(perf))
perfDf_sjcvb0$method <- "SJCVB0"

perfDf_sjcvb0$K <- "K150"
perfDf_sjcvb0 <- perfDf_sjcvb0[,c("X1","K","X2","value","method")]
colnames(perfDf_sjcvb0)[c(1,3)] <- c("Iter", "B")

perfDf_sjcvb0$B <- sub("B","",perfDf_sjcvb0$B)
perfDf_sjcvb0$B <- as.numeric(perfDf_sjcvb0$B)
perfDf_sjcvb0$B <- as.factor(perfDf_sjcvb0$B)

perfDf <- rbind(perfDf_jcvb0, perfDf_sjcvb0)

perfDf <- subset(perfDf, value != 0)

perfDf$B <- factor(perfDf$B, levels=as.character(sort(unique(as.numeric(perfDf$B)))))

perfDf <- subset(perfDf, !B %in% c("100","500","1000"))

gg <- ggplot(subset(perfDf, K=="K150"), 
             aes(x=Iter, y=value, color=B, linetype=method)) +  # shape=method 
  # facet_grid(~method) +
  # geom_point() + 
  geom_line() + theme_bw() +
  geom_hline(yintercept = max(subset(perfDf, K=="K150" & method=="JCVB0")$value),
             color='blue', linetype="dashed") +
  scale_color_brewer(type = "qual", palette="Set1") +
  ylab("Predictive log Likelihood") + xlab("Iterations")  

ggout <- "/results/FigS2_loglik_comparison_sjcvb0.eps"

ggsave(ggout, gg, width=6.4, height=4.8)















