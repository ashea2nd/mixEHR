library(Rcpp)
library(gridExtra)
library(plyr)
library(ggplot2)
library(reshape)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)

rm(list=ls(all=TRUE))

## this script will plot only the select topics shown in the Figure 2 of the main text ##

inpdir <- "../../data/precomputed_results/mixmimic" # here

# from pheTopic.R
# phiTopicOrderRda <- sprintf("%s/phePhiTopicOrder.RData", inpdir)
phiTopicOrderRda <- "/results/Fig2_phePhiTopicOrder.RData"
load(phiTopicOrderRda) # phiTopics phisel

featId <- "../../data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId)


#### nmar rates under topics ####
K <- 75
niter <- 5497

filein_phi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_phi_normalized.csv", inpdir, K, niter)
filein_psi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_psi.csv", inpdir, K, niter)

inpdir <- dirname(filein_phi)

psi <- read.csv(filein_psi, header=F)

colnames(psi)[c(1:2)] <- c("typeId","pheId")

colnames(psi)[-c(1:2)] <- sprintf("M%s", 1:(ncol(psi)-2))

featId <- psi[,c(1,2)]

ehrFeatInfo_t <- ehrFeatInfo[[as.numeric(psi[1,"typeId"])]]

psi <- as.matrix(psi[,-c(1:2)])

rownames(psi) <- sprintf("%s:%s",
                         names(ehrFeatInfo)[ehrFeatInfo_t$typeId],
                         ehrFeatInfo_t$pheName[match(featId[,"pheId"], ehrFeatInfo_t$pheId)])

labsel <- which(apply(psi, 1, max) > 0.01)

#### select labels ####
icds <- c("250", "296", "331", "042", "205", "415", "571")
icds <- sprintf("icd_cm:%s", icds)

icd <- icds[7]

phiTopics <- sapply(icds, function(icd) {
  idx <- grep(icd, rownames(phisel), value=T)
  print(idx)
  colnames(phisel)[which.max(colMeans(phisel[idx,,drop=F]))]
})

labTopicSel <- phiTopics <- unique(phiTopics)

phisel <- phisel[, phiTopics]
psisel <- psi[labsel, labTopicSel]

topicsel <- colnames(psisel)

topicOrder <- topicsel[order(as.numeric(sub("M","",topicsel)))]

phisel <- phisel[, topicOrder]
psisel <- psisel[, topicOrder]

topW <- 5
topL <- 5

topphe <- unique(unlist(lapply(1:ncol(phisel), function(k) {
  # c(head(rownames(phisel)[order(phisel[,k], decreasing = T)], topW),
  #   head(rownames(psisel)[order(psisel[,k], decreasing = T)], topL))
  head(rownames(phisel)[order(phisel[,k], decreasing = T)], topW)
})))

# phisel <- rbind(phisel, psisel)

topsel <- phisel[topphe,]  

typeColors <- brewer.pal(length(ehrFeatInfo), "Dark2")
names(typeColors) <- names(ehrFeatInfo)

df <- melt(topsel)
colnames(df) <- c("pheId0", "topic", "freq")
df$pheId0 <- as.character(as.matrix(df$pheId0))
df$typeId <- unlist(sapply(strsplit(df$pheId0,":"), function(x) x[[1]]))
pheIds <- sapply(strsplit(df$pheId0,":"), function(x) x[[2]])
df$pheId <- sapply(strsplit(pheIds,"_"), function(x) tail(x,1))

groupAnnot <- rowAnnotation(data.frame(category=names(typeColors[df$typeId])),
                            col=list(category=typeColors[df$typeId]))

featName <- sapply(strsplit(rownames(topsel),":"), function(x) {
  
  datatype <- x[[1]]
  
  if(datatype %in% c("icd_cpt","icd_cm","presc")) {
    
    tmp <- unlist(strsplit(x[[2]], "_"))
    x1 <- tolower(tmp[2])
    substr(x1,1,1) <- toupper(substr(x1,1,1))
    x1
    
  } else if(datatype == "drg") {
    
    tmp <- unlist(strsplit(x[[2]], "_"))
    
    x1 <- tmp[3]
    x1 <- tolower(x1)
    substr(x1,1,1) <- toupper(substr(x1,1,1))
    x1
  } else if(datatype=="notes") {
    x[[2]]
  }
})


featId <- sapply(strsplit(rownames(topsel),":"), function(x) {
  
  datatype <- x[1]
  
  if(datatype %in% c("icd_cpt","icd_cm")) {
    
    tmp <- unlist(strsplit(x[2], "_"))[1]
    
    sprintf("(%s)", tmp)
    
  } else if(datatype == "drg") {
    
    tmp <- paste(gsub(" ","",unlist(strsplit(unlist(x)[2], "_"))[1:2]),collapse = ",")
    
    # tmp2 <- paste(tail(unlist(strsplit(unlist(x)[2], "_")),2),collapse = ",")
    
    sprintf("(%s)", tmp)
    
  } else if (datatype=="presc") {
    
    tmp <- paste(tail(unlist(strsplit(x[2], "_")),2),collapse=",")
    
    sprintf("(%s)", tmp)
    
  } else if(datatype=="notes") {
    
    ""
  }
})

rownames(topsel) <- sprintf("%s %s", featName, featId)

hmp <- Heatmap(topsel, 
               # column_order = topicOrder, row_order=featOrder, 
               cluster_columns=FALSE, 
               cluster_rows=FALSE,
               name="Prob",
               row_names_side = "left",
               row_dend_side = "right",
               # show_row_dend=FALSE,
               # show_column_dend=FALSE,
               rect_gp = gpar(col = "grey", lwd =0.2),
               row_names_gp = gpar(fontsize=10),
               col=c("white", "red"),
               show_row_names = T,
               row_title="",
               column_title="Disease topics") + groupAnnot

myw <- 8
myh <- 6

hmpout <- "/results/Fig2a_TopicHmp_sel.eps"
cairo_ps(hmpout, width=myw, height=myh)
draw(hmp, heatmap_legend_side = "right", padding=unit(c(5,65,5,20), "mm")) # bottom, left, top, right
dev.off()




#### lab abnormal results under topics ####
filein_eta <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_eta_normalized.csv", inpdir, K, niter)
filein_psi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_psi.csv", inpdir, K, niter)

eta <- read.csv(filein_eta, header=F)
psi <- read.csv(filein_psi, header=F)[,-c(1,2)]

colnames(eta)[c(1:3)] <- c("typeId","labId","stateId")
colnames(eta)[-c(1:3)] <- sprintf("M%s", 1:(ncol(eta)-3))

labId <- subset(eta, stateId==1)[,c(1,2)]

ehrFeatInfo_t <- ehrFeatInfo[[as.numeric(eta[1,"typeId"])]]

eta$labId <- ehrFeatInfo_t$pheName[match(eta$labId, ehrFeatInfo_t$pheId)]

rownames(psi) <- subset(eta, stateId==1)$labId

etaNrm <- subset(eta, stateId==0)[,-c(1:3)]
etaAbn <- subset(eta, stateId==1)[,-c(1:3)]

etaAbnPsi <- etaAbn * psi

etaAbnPsi <- scale(etaAbnPsi, center=F, scale=colSums(etaAbnPsi))

rownames(etaAbnPsi) <- rownames(etaNrm) <- rownames(etaAbn) <- rownames(psi)

topL <- 5

# toplab_idx <- unique(as.numeric(sapply(1:K, function(k) head(order(etaAbnPsi[,k], decreasing = TRUE), topL))))

toplab <- unique(as.numeric(sapply(1:ncol(etaAbnPsi), function(k) 
  head(order(etaAbnPsi[,k], decreasing = TRUE), topL))))

etaAbnPsi <- as.matrix(etaAbnPsi[toplab,topicOrder])

require(circlize) # colorRamp2

topres <- names(which(apply(etaAbnPsi, 1, max) > 0.05))
# toptop <- names(which(apply(etaAbnPsi, 2, max) > 0.01))

idx <- which(rownames(etaAbnPsi) %in% topres)

groupAnnot <- rowAnnotation(data.frame(category=rep("lab",nrow(etaAbnPsi))),
                            col=list(category=rep(typeColors[4],nrow(etaAbnPsi))))


labName <- sapply(strsplit(rownames(etaAbnPsi[idx,]),"_"), function(x) x[2])
labName <- tolower(labName)
substr(labName,1,1) <- toupper(substr(labName,1,1))
labId <- sapply(strsplit(rownames(etaAbnPsi[idx,]),"_"), function(x) x[1])

rownames(etaAbnPsi)[idx] <- sprintf("%s (%s)", labName, labId)

hmp <- Heatmap(etaAbnPsi[idx,],
               cluster_columns=F,
               cluster_rows=F,
               row_names_side = "left",
               row_dend_side = "right",
               rect_gp = gpar(col = "grey", lwd =0.2),
               # col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               col = c("white", "red"),
               name="abnm lab",
               column_title="Disease topics") + groupAnnot

myw <- 7
myh <- 6

hmpout <- "/results/FigS4_labEtaTopicHmp_sel.eps"
cairo_ps(hmpout, width=myw, height=myh)
# draw(hmp, heatmap_legend_side = "left", padding=unit(c(2,10,2,100), "mm")) # bottom, left, top, right
draw(hmp, heatmap_legend_side = "right", padding=unit(c(5,50,5,10), "mm")) # bottom, left, top, right
dev.off()

system("rm -f /results/Fig2_phePhiTopicOrder.RData")

