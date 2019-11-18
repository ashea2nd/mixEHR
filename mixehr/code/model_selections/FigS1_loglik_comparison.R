library(ggplot2)
library(plyr)
library(reshape)

rm(list=ls(all=TRUE))

inpdir <- "/data/results/model_selections/mixmimic_cv"

fileins <- list.files(inpdir, "logPredLik.txt", full.names = TRUE)

fileins <- grep("SJCVB0", fileins, invert=T, value=T)

niter <- 200

fileins <- grep(sprintf("iter%s_",niter), fileins, value=TRUE)

configs <- do.call(rbind, strsplit(basename(fileins), "_"))


colnames(configs) <- c("data","fold","method","model","K","iter","type")

nfolds <- 5

configs <- as.data.frame(configs, stringsAsFactors = F)

# evaluate only the ones with all folds completed at niter
hasAllFolds <- names(which(sapply(split(configs, configs[,"K"]), nrow) == nfolds))

if(!all(unique(configs[,"K"]) %in% hasAllFolds)) {
  incomplK <- unique(configs[,"K"])[!unique(configs[,"K"]) %in% hasAllFolds]
  print(subset(configs, K %in% incomplK))
  stopifnot(all(unique(configs[,"K"]) %in% hasAllFolds))
}

configs <- subset(configs, K %in% hasAllFolds)

j <- 1

loglikDf <- do.call(rbind, lapply(1:length(fileins), function(j) {
  
  f <- fileins[j]
  
  loglik <- read.table(f, header=F, col.names = "loglik")[-1,,drop=F]
  
  loglik$K <- as.numeric(sub("K","",configs[j,"K"]))
  
  loglik$fold <- configs[j,"fold"]
  
  loglik$iter <- 1:nrow(loglik)
  
  loglik
}))

loglikDf_avg <- aggregate(loglik~K+iter, loglikDf, mean)

loglikDf_avg$K <- as.factor(loglikDf_avg$K)

gg <- ggplot(loglikDf_avg, aes(x=iter, y=loglik, color=K)) + 
  geom_point() + # theme_bw() +
  scale_color_brewer(type='qual',palette='Set1')

# print(gg)

ggout <- sprintf("/results/FigS1_loglik_comparison.eps")

ggsave(ggout, gg, width=8, height=6)



