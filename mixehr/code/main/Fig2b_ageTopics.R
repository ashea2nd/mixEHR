library(ggplot2)
library(plyr)
library(ComplexHeatmap)
library(RColorBrewer)
library(reshape2)

rm(list=ls(all=TRUE))

featIdRda <- "../../data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featIdRda)

inpdir <- "/data/precomputed_results/mixmimic"

### load phi ###
K <- 75
niter <- 5497

filein_phi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_phi_normalized.csv", inpdir, K, niter)
filein_eta <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_eta_normalized.csv", inpdir, K, niter)
filein_psi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_psi.csv", inpdir, K, niter)
filein_metaphe <- sprintf("%s/mimic_trainData_mimic_trainData_JCVB0_nmar_K%s_iter%s_metaphe.csv", inpdir, K, niter)
filein_ehrinp <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"

patrda <- "/data/Mimic/version_1_3/processed/PATIENTS.RData"
load(patrda) # age, adm_single

metaphe <- read.csv(filein_metaphe, header=F)

ehrinp <- read.table(filein_ehrinp, header=F)
colnames(ehrinp) <- c("patId","typeId","pheId","stateId","freq")
rownames(metaphe) <- unique(ehrinp$patId)

patage <- age$age[match(rownames(metaphe), age$SUBJECT_ID)]

patage[patage>90] <- 90

colnames(metaphe) <- sub("V","M", colnames(metaphe))

metage <- cor(patage, metaphe, method = 'pearson')

metage <- t(metage)

metage <- metage[order(metage, decreasing = F),,drop=F]

df <- reshape2::melt(metage)

colnames(df) <- c("Metaphe", "Var2", "Cor")

gg1 <- ggplot(df, aes(x=Metaphe, y=Cor)) + theme_bw() +
  geom_bar(stat='identity', fill='grey') + 
  geom_text(aes(label=Metaphe)) + coord_flip()

ggout1 <- "/results/FigS5a_ageToipcCor.eps"
ggsave(ggout1, gg1, width=6.4, height=11)


# topagemeta <- melt(data.frame(metaphe[,head(rownames(metage),3)], age=patage), id.vars = "age")

metaphe_sel <- c(head(rownames(metage),3), tail(rownames(metage),3))

agemeta_sel <- melt(data.frame(metaphe[,metaphe_sel], age=patage), id.vars = "age")

colnames(agemeta_sel)[2:3] <- c("Metaphe","Prob")

gg2 <- ggplot(agemeta_sel, aes(x=age, y=Prob, color=Metaphe)) + 
  geom_smooth(method='lm', fill=NA)

ggout2 <- "/results/FigS5b_ageTopicTrend.eps"
ggsave(ggout2, gg2, width=6.4, height=4.8)


#### get age topic ####
resdir <- "/data/precomputed_results/mixmimic"
resrda <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_combined.RData", resdir, K, niter)
load(resrda)

phisel <- phi[, rev(metaphe_sel)]

for(topW in c(3, 10)) {
    
    topphe <- unique(unlist(lapply(1:ncol(phisel), function(k) {
      head(rownames(phisel)[order(phisel[,k], decreasing = T)], topW)
    })))

    # phisel <- rbind(phisel, psisel)

    topsel <- phisel[topphe,]

    typeColors <- brewer.pal(length(ehrFeatInfo), "Dark2")
    names(typeColors) <- names(ehrFeatInfo)

    df <- melt(as.matrix(topsel))
    colnames(df) <- c("pheId0", "topic", "freq")
    df$pheId0 <- as.character(as.matrix(df$pheId0))
    df$typeId <- unlist(sapply(strsplit(df$pheId0,":"), function(x) x[[1]]))
    pheIds <- sapply(strsplit(df$pheId0,":"), function(x) x[[2]])
    df$pheId <- sapply(strsplit(pheIds,"_"), function(x) tail(x,1))

    groupAnnot <- rowAnnotation(data.frame(category=names(typeColors[df$typeId])),
                                col=list(category=typeColors[df$typeId]))

    rownames(topsel) <- sub("_NA_NA","", rownames(topsel))

    featName <- sapply(strsplit(rownames(topsel),":"), function(x) {

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

    featId <- sapply(strsplit(rownames(topsel),":"), function(x) {

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

    # rownames(topsel) <- sprintf("%s %s", featName, featId)

    rownames(topsel) <- substr(rownames(topsel), 1, 70)

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
        
    if(topW == 3) {
        myw <- 7.5
        myh <- 4
        hmpout <- "/results/Fig2b_phiAgeTopicHmp.eps"
        cairo_ps(hmpout, width=myw, height=myh)        
        draw(hmp, heatmap_legend_side = "right", padding=unit(c(5,53,5,20), "mm")) # bottom, left, top, right        
        dev.off()        
    } else {
        myw <- 11
        myh <- 11
        hmpout <- "/results/FigS5c_phiAgeTopicHmp.eps"        
        cairo_ps(hmpout, width=myw, height=myh)
        draw(hmp, heatmap_legend_side = "right", padding=unit(c(5,140,5,20), "mm")) # bottom, left, top, right        
        dev.off()        
    }
}


#### lab abnormal results under topics ####
etaNrm <- subset(eta, stateId==0)[,-c(1:3)]
etaAbn <- subset(eta, stateId==1)[,-c(1:3)]

colnames(eta) <- sub("V","M",colnames(eta))
colnames(psi) <- sub("V","M",colnames(psi))

etaAbnPsi <- etaAbn[,rev(metaphe_sel)] * psi[,rev(metaphe_sel)]
# etaAbnPsi <- psi[,metaphe_sel]

etaAbnPsi <- scale(etaAbnPsi, F, scale=colSums(etaAbnPsi))

rownames(etaAbnPsi) <- rownames(etaNrm) <- rownames(etaAbn) <- rownames(psi)

topL <- 5

toplab <- unique(as.numeric(sapply(1:ncol(etaAbnPsi), function(k) 
  head(order(etaAbnPsi[,k], decreasing = TRUE), topL))))

etaAbnPsi <- as.matrix(etaAbnPsi[toplab,])

labName <- sapply(strsplit(rownames(etaAbnPsi),"_"), function(x) x[2])
labName <- tolower(labName)
substr(labName,1,1) <- toupper(substr(labName,1,1))
labId <- sapply(strsplit(rownames(etaAbnPsi),"_"), function(x) x[1])

rownames(etaAbnPsi) <- sprintf("%s (%s)", labName, labId)

groupAnnot <- rowAnnotation(data.frame(category=rep("lab",nrow(etaAbnPsi))),
                            col=list(category=rep(typeColors[4],nrow(etaAbnPsi))))

hmp <- Heatmap(etaAbnPsi,
               cluster_columns=FALSE, 
               cluster_rows=FALSE,
               row_names_side = "left",
               row_dend_side = "right",
               rect_gp = gpar(col = "grey", lwd =0.2),
               # col = colorRamp2(c(0, 0.1), c("white", "red")),
               col = c("white", "red"),
               name="Prob",
               column_title="Meta-phenotypes") + groupAnnot

myw <- 8
myh <- 6

hmpout2 <- "/results/FigS5d_labAgeTopicHmp.eps"
cairo_ps(hmpout2, width=myw, height=myh)
draw(hmp, heatmap_legend_side = "right", padding=unit(c(5,60,5,10), "mm"))
# bottom, left, top, right
dev.off()













