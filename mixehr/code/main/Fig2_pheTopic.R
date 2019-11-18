library(Rcpp)
library(gridExtra)
library(plyr)
library(ggplot2)
library(reshape)
library(ComplexHeatmap)

rm(list=ls(all=TRUE))

inpdir <- "../../data/precomputed_results/mixmimic" # here

featId <- "../../data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId)

K <- 75
niter <- 5497

filein_phi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_phi_normalized.csv", inpdir, K, niter)

inpdir <- dirname(filein_phi)

phi <- read.csv(filein_phi, header=F)

colnames(phi)[c(1,2)] <- c("typeId","pheId")

colnames(phi)[-c(1,2)] <- sprintf("M%s", 1:(ncol(phi)-2))

philist <- split(phi, phi$typeId)

phi <- do.call(rbind, lapply(philist, function(x) {
  
  featId <- x[,c(1,2)]
  
  ehrFeatInfo_t <- ehrFeatInfo[[as.numeric(x[1,"typeId"])]]
  
  x <- as.matrix(x[,-c(1,2)])
  
  rownames(x) <- sprintf("%s:%s",
                         names(ehrFeatInfo)[ehrFeatInfo_t$typeId],
                         ehrFeatInfo_t$pheName[match(featId[,"pheId"], ehrFeatInfo_t$pheId)])
  x
}))

# remove uninformative topics and features
featsel <- which(apply(phi, 1, max) > 0.01)
topicsel <- which(apply(phi, 2, max) > 0.01)

phisel <- phi[featsel,topicsel]

pheOrder <- hclust(dist(phisel))$order
topicOrder <- hclust(dist(t(phisel)))$order

phiTopics <- colnames(phisel)[topicOrder]

rdaout <- "/results/Fig2_phePhiTopicOrder.RData"
save(phiTopics, phisel, file=rdaout)


