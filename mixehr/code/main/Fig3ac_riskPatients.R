library(ggplot2)
library(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

rm(list=ls(all=TRUE))

K <- 75
niter <- 5497

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId)

inpdir <- "/data/precomputed_results/mixmimic"

### load phi ###
filein_phi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_phi_normalized.csv", inpdir, K, niter)
filein_eta <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_eta_normalized.csv", inpdir, K, niter)
filein_psi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_psi.csv", inpdir, K, niter)
filein_metaphe <- sprintf("%s/mimic_trainData_mimic_trainData_JCVB0_nmar_K%s_iter%s_metaphe.csv", inpdir, K, niter)
filein_ehrinp <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"

#### load phi ####
phi <- read.csv(filein_phi, header=F)
psi <- read.csv(filein_psi, header=F)

# phi <- rbind(phi, psi)

colnames(phi)[c(1,2)] <- c("typeId","pheId")

colnames(phi)[-c(1,2)] <- sprintf("M%s", 1:(ncol(phi)-2))

phi <- plyr::arrange(phi, typeId, pheId)

pheId <- phi[,c(1,2)]

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

#### get observed data ####
# check patient ID
ehrinp <- read.table(filein_ehrinp, header=F)

colnames(ehrinp) <- c("patId","typeId","pheId","stateId","freq")

#### select labels ####
icds <- c("205", "415", "571")
icds <- sprintf("icd_cm:%s", icds)

icd <- icds[2]

metaphe_sel <- sapply(icds, function(icd) 
  which.max(colMeans(phi[grep(icd, rownames(phi), value=T),,drop=F])))

metaphe_sel <- unique(metaphe_sel)

x <- icds[1]

ehrFeatInfo_sel <- lapply(icds, function(x) {
  idx <- grep(sprintf("%s", x), sprintf("icd_cm:%s", ehrFeatInfo[["icd_cm"]]$ICD9_CODE))
  ehrFeatInfo[["icd_cm"]][idx,]
})

names(ehrFeatInfo_sel) <- icds

topD <- 50

metaphe <- read.csv(filein_metaphe, header=F)

rownames(metaphe) <- unique(ehrinp$patId)

# k=metaphe_sel[1]

metaphe_riskpats <- do.call(rbind, lapply(metaphe_sel, function(k) 
  
  head(metaphe[order(metaphe[,k], decreasing = TRUE),], topD)))

colnames(metaphe_riskpats) <- sub("V","M",colnames(metaphe_riskpats))

ehrinp_riskpats <- subset(ehrinp, patId %in% rownames(metaphe_riskpats))

ehrinp_riskpats_list <- split(ehrinp_riskpats, ehrinp_riskpats$patId)

ehrinp_riskpats_list <- ehrinp_riskpats_list[rownames(metaphe_riskpats)]

riskpats_class <- do.call(rbind, lapply(ehrinp_riskpats_list, function(x) {
  sapply(ehrFeatInfo_sel, function(y) {
    any(x$typeId %in% y$typeId & x$pheId %in% y$pheId)
  })
}))

traitGroupColor <- brewer.pal(length(icds), 'Set1')

names(traitGroupColor) <- icds

# traitGroupAnnot <- rowAnnotation(
#   data.frame(TraitGroup=names(traitGroupColor)), 
#   col = list(TraitGroup=traitGroupColor))

col1 <- rep("grey", nrow(riskpats_class))
col1[riskpats_class[,1]==1] <- traitGroupColor[1]
names(col1) <- as.character(riskpats_class[,1])
traitGroupAnnot1 <- rowAnnotation(
  data.frame(Leukemia=riskpats_class[,1]), 
  col=list(Leukemia=col1))

col2 <- rep("grey", nrow(riskpats_class))
col2[riskpats_class[,2]==1] <- traitGroupColor[2]
names(col2) <- as.character(riskpats_class[,2])
traitGroupAnnot2 <- rowAnnotation(
  data.frame(PulmEmbol=riskpats_class[,2]),
  col=list(PulmEmbol=col2))

col3 <- rep("grey", nrow(riskpats_class))
col3[riskpats_class[,3]==1] <- traitGroupColor[3]
names(col3) <- as.character(riskpats_class[,3])
traitGroupAnnot3 <- rowAnnotation(
  data.frame(Cirrhosis=riskpats_class[,3]),
  col=list(Cirrhosis=col3))

colnames(metaphe_riskpats)[!colnames(metaphe_riskpats) %in% sprintf("M%s", metaphe_sel)] <- ""

metaphe_riskpats <- as.matrix(metaphe_riskpats)

hm <- Heatmap(metaphe_riskpats, 
              name="Patient-topics",
              row_names_side = "left",
              cluster_rows = F,
              cluster_columns = F,
              rect_gp = gpar(col = "grey", lwd =0.1),
              col=c("white", "red"),
              show_row_names = F,
              row_title="Top Risk Patients",
              column_title="Disease topics") + 
  traitGroupAnnot1 + traitGroupAnnot2 + traitGroupAnnot3

ggout <- "/results/Fig3a_riskpats_metaphe.eps"

cairo_ps(ggout, width=8, height=6)
draw(hm)
dev.off()


#### top phe ####
topW <- 10

tophi_idx <- as.numeric(sapply(metaphe_sel, function(k) head(order(phi[,k], decreasing = TRUE), topW)))

topphi <- phi[tophi_idx,]

topPheId <- pheId[tophi_idx,]
topPheIdKey <- sprintf("%s_%s", topPheId[,1], topPheId[,2])

riskpats_imp <- tcrossprod(topphi, metaphe_riskpats)

# riskpats_imp_thres <- quantile(as.numeric(riskpats_imp), 0.8)

riskpats_obs <- matrix(0, nrow(riskpats_imp), ncol(riskpats_imp),
                       dimnames=dimnames(riskpats_imp))

ehrinp_sel <- subset(ehrinp, patId %in% rownames(metaphe_riskpats))

j <- 1

for(j in 1:ncol(riskpats_obs)) {
  
  ehrinp_riskpats <- subset(ehrinp_sel, patId == rownames(metaphe_riskpats)[j])
  
  tmp1 <- sprintf("%s_%s", ehrinp_riskpats[,"typeId"], ehrinp_riskpats[,"pheId"])
  
  riskpats_obs[,j] <- ehrinp_riskpats$freq[match(topPheIdKey, tmp1)]
}

riskpats_obs[is.na(riskpats_obs)] <- 0

riskpats_imp <- riskpats_imp[,rownames(riskpats_class)]
riskpats_obs <- riskpats_obs[,rownames(riskpats_class)]

traitGroupAnnot1 <- HeatmapAnnotation(
  data.frame(Leukemia=riskpats_class[,1]), 
  col=list(Leukemia=col1))

traitGroupAnnot2 <- HeatmapAnnotation(
  data.frame(PulmEmbol=riskpats_class[,2]),
  col=list(PulmEmbol=col2))

traitGroupAnnot3 <- HeatmapAnnotation(
  data.frame(Cirrhosis=riskpats_class[,3]),
  col=list(Cirrhosis=col3))

riskpats_labels <- rep("Unrecorded", nrow(riskpats_class))
riskpats_labels[riskpats_class[,1]==1] <- "Leukemia"
riskpats_labels[riskpats_class[,2]==1] <- "PulmEmbol"
riskpats_labels[riskpats_class[,3]==1] <- "Cirrhosis"

riskpats_colors <- rep("grey", nrow(riskpats_class))
riskpats_colors[riskpats_class[,1]==1] <- traitGroupColor[icds[1]]
riskpats_colors[riskpats_class[,2]==1] <- traitGroupColor[icds[2]]
riskpats_colors[riskpats_class[,3]==1] <- traitGroupColor[icds[3]]

names(riskpats_colors) <- riskpats_labels

traitGroupAnnot <- HeatmapAnnotation(
  data.frame(ICD9=riskpats_labels), 
  col=list(ICD9=riskpats_colors))

riskpats_obs[riskpats_obs>0] <- 1

typeColors <- brewer.pal(length(ehrFeatInfo), "Dark2")
names(typeColors) <- names(ehrFeatInfo)
groupAnnot <- rowAnnotation(data.frame(category=names(typeColors[topPheId$typeId])),
                            col=list(category=typeColors[topPheId$typeId]))

patitle <- "High-risk patients"


featName <- sapply(strsplit(rownames(riskpats_imp),":"), function(x) {
  
  datatype <- x[1]
  
  if(datatype %in% c("icd_cpt","icd_cm","lab")) {
    
    tmp <- unlist(strsplit(x[2], "_"))
    x1 <- tolower(tmp[2])
    substr(x1,1,1) <- toupper(substr(x1,1,1))
    x1
    
  } else if(datatype=="presc") {
    
    x[2] <- gsub("__","_",x[2])
    
    tmp <- unlist(strsplit(x[2], "_"))
    
    if(length(tmp) > 3) {
      x1 <- tolower(paste(tmp[-c(length(tmp)-1,length(tmp))],collapse = ","))
    } else {
      x1 <- tolower(tmp[1])  
    }
    
    substr(x1,1,1) <- toupper(substr(x1,1,1))
    x1
    
  } else if(datatype == "drg") {
    
    tmp <- unlist(strsplit(x[[2]], "_"))
    
    x1 <- tmp[3]
    x1 <- tolower(x1)
    substr(x1,1,1) <- toupper(substr(x1,1,1))
    x1
    
  } else if(datatype=="notes") {
    x[2]
  }
})

featId <- sapply(strsplit(rownames(riskpats_imp),":"), function(x) {
  
  datatype <- x[1]
  
  if(datatype %in% c("icd_cpt","icd_cm","lab")) {
    
    tmp <- unlist(strsplit(x[2], "_"))[1]
    
    sprintf("(%s)", tmp)
    
  } else if(datatype == "drg") {
    
    tmp <- paste(gsub(" ","",unlist(strsplit(unlist(x)[2], "_"))[1:2]),collapse = ",")
    
    # tmp2 <- paste(tail(unlist(strsplit(unlist(x)[2], "_")),2),collapse = ",")
    
    # sprintf("(%s,%s)", tmp,tmp2)
    
    sprintf("(%s)", tmp)
    
  } else if (datatype=="presc") {
    
    tmp <- paste(tail(unlist(strsplit(x[2], "_")),2),collapse=",")
    
    sprintf("(%s)", tmp)
    
  } else if(datatype=="notes") {
    
    ""
  }
})

rownames(riskpats_obs) <-sprintf("%s %s", featName, featId)

hm_obs <- Heatmap(riskpats_obs, 
                  top_annotation = traitGroupAnnot,
                  name="Term freq",
                  row_names_side = "left",
                  cluster_rows = F,
                  cluster_columns = F,
                  rect_gp = gpar(col = "grey", lwd =0.1),
                  col=c("white", "red"),
                  show_column_names = F,
                  row_title="",
                  column_title=patitle) + groupAnnot

ggout <- "/results/Fig3c_riskpats_topCodes.eps"
cairo_ps(ggout, width=14, height=6)
draw(hm_obs, heatmap_legend_side = "right", padding=unit(c(5,90,5,10), "mm")) # bottom, left, top, right
dev.off()









































