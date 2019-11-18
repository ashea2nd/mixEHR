library(ggplot2)
library(plyr)
library(ROCR)
library(gridExtra)

rm(list=ls(all=TRUE))

inpdir <- "/data/precomputed_results/mortality"

source('/code/main/eval_fun.R')

rdalist <- list.files(inpdir, "predictDeath", full.names = TRUE)

names(rdalist) <- sub("predictDeath_","",sub(".RData","",basename(rdalist)))
names(rdalist) <- sub("test_","",names(rdalist))

rdalist <- c(grep("LDA", rdalist, invert=T, value=T), grep("LDA", rdalist, value=T))

x <- names(rdalist)

rocdf <- do.call(rbind, lapply(names(rdalist), function(x) {
  load(rdalist[x])
  roceval(res_prospect[,1], res_prospect[,2], x)    
}))

gg_roc <- evalLinePlot(rocdf, "ROC")

prcdf <- do.call(rbind, lapply(names(rdalist), function(x) {
  load(rdalist[x])
  prceval(res_prospect[,1], res_prospect[,2], x)
}))

gg_prc <- evalLinePlot(prcdf, "PRC")

ggout <- "/results/Fig5a_prospectPerf.eps"

cairo_ps(ggout, width=8, height=6)
grid.arrange(gg_roc, gg_prc)
dev.off()




























