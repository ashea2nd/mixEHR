library(ggplot2)
# library(xtable)
library(plyr)
library(scales)
library(gridExtra)
# library(tidytext)
library(ComplexHeatmap)
library(circlize)
library(reshape2)

rm(list=ls(all=TRUE))

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"

load(featId) # ehrFeatInfo

### load phi ###
K <- 75
niter <- 5497

inpdir <- "/data/precomputed_results/mixmimic"

filein_phi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_phi_normalized.csv", inpdir, K, niter)
filein_psi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_psi.csv", inpdir, K, niter)

inpdir <- dirname(filein_phi)

phi <- read.csv(filein_phi, header=F)
psi <- read.csv(filein_psi, header=F)

phi <- rbind(phi, psi)

colnames(phi)[c(1,2)] <- c("typeId","pheId")

colnames(phi)[-c(1,2)] <- sprintf("M%s", 1:(ncol(phi)-2))

# from mimic/prep/icd9groupings.R
icd9rda <- "/data/Mimic/version_1_3/processed/metadata/icd9ref.RData"

load(icd9rda) # icd9ref

philist <- split(phi, phi$typeId)

trainData <- "/data/Mimic/version_1_3/processed/mimic_trainData.txt"
ehrinp <- read.table(trainData, skip=1)
colnames(ehrinp) <- c("patId","typeId","pheId","stateId","freq")

x <- philist[[1]]

phi_normalized <- do.call(rbind, lapply(philist, function(x) {
  
  featId <- x[,c(1,2)]
  
  mytypeid <- x[1,"typeId"]
  
  ehrinp_tar <- subset(ehrinp, typeId==mytypeid)
  
  patcnt <- table(ehrinp_tar$pheId)[as.character(x[,"pheId"])]
  
  ehrFeatInfo_t <- ehrFeatInfo[[mytypeid]]
  
  x <- as.matrix(x[,-c(1,2)])
  
  x <- cbind(type=names(ehrFeatInfo)[ehrFeatInfo_t$typeId],
             phe=as.character(ehrFeatInfo_t$pheName[match(featId[,"pheId"], ehrFeatInfo_t$pheId)]), 
             patcnt=as.numeric(patcnt),
             as.data.frame(x))
  x
}))

thres <- 100

df <- subset(phi_normalized, patcnt > thres & type=="icd_cm")

rownames(df) <- df$phe

x <- cor(t(df[,-c(1:3)]), method='pearson')

idx <- grep("^[VE]", rownames(x), invert = T)

x1 <- x[idx,idx]


#### icd taxonomy info ####
# from icd9groupings.R
icd9refrda <- "/data/Mimic/version_1_3/processed/metadata/icd9ref.RData"
load(icd9refrda) # icd9ref
icd9ref$ICD_KEY <- sprintf("%s_%s", icd9ref$ICD9_CODE, icd9ref$SHORT_TITLE)

icd9grp <- icd9ref$level1[match(rownames(x1), icd9ref$ICD_KEY)]

# groupColors <- rainbow(length(unique(icd9grp)), s = 0.5)

groupColors <- rainbow(length(unique(icd9grp)))

set.seed(123)
groupColors <- sample(groupColors, length(groupColors))

names(groupColors) <- unique(icd9grp)

groupAnnot <- rowAnnotation(data.frame(category=icd9grp),
                            col=list(category=groupColors[icd9grp]))

hmAnnot <- HeatmapAnnotation(data.frame(category=icd9grp),
                             col=list(category=groupColors[icd9grp]), show_legend=FALSE)

hmp <- Heatmap(x1,
               name="icd9Cor",
               top_annotation=hmAnnot,
               column_dend_side="bottom",
               row_dend_side="right",
               row_names_gp = gpar(fontsize = 3),
               rect_gp = gpar(col = "grey", lwd =0.2),
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               column_names_gp = gpar(fontsize = 3),
               show_heatmap_legend=FALSE) + groupAnnot

hmpout1 <- "/results/FigS6_icd9corhmp.pdf"
cairo_pdf(hmpout1, width=23, height=16)
draw(hmp, heatmap_legend_side = "right", padding=unit(c(2,2,2,20), "mm")) # bottom, left, top, right
dev.off()


#### cluster raw ICD ####
tarTypeId <- which(names(ehrFeatInfo)=="icd_cm")

tarpheid <- ehrFeatInfo$icd_cm$pheId[match(rownames(x1), ehrFeatInfo$icd_cm$pheName)]

ehrinp_sel <- subset(ehrinp, typeId == tarTypeId & pheId %in% tarpheid)

# x <- as.matrix(cast_sparse(ehrinp_sel, patId, pheId, freq))

x0 <- acast(ehrinp_sel[,c("patId","pheId","freq")], patId~pheId, sum)

colnames(x0) <- ehrFeatInfo$icd_cm$pheName[match(colnames(x0), ehrFeatInfo$icd_cm$pheId)]

x2 <- cor(x0, method='pearson')

hmp <- Heatmap(x2,
               name="icd9Cor",
               # cluster_rows=FALSE,
               # cluster_columns=FALSE,
               # row_order=unlist(row_order(ht)),
               row_order=unlist(row_order(hmp)),
               column_order=column_order(hmp)$icd9Cor,
               top_annotation=hmAnnot,
               column_dend_side="bottom",
               row_dend_side="right",
               row_names_gp = gpar(fontsize = 3),
               rect_gp = gpar(col = "grey", lwd =0.2),
               col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
               column_names_gp = gpar(fontsize = 3),
               show_heatmap_legend=FALSE) + groupAnnot

hmpout2 <- "/results/FigS6_icd9corhmp_raw.pdf"
cairo_pdf(hmpout2, width=23, height=16)
draw(hmp, heatmap_legend_side = "right", padding=unit(c(2,2,2,20), "mm")) # bottom, left, top, right
dev.off()



#### cluster metaphe ####
# x2 <- cor(phi_normalized[,-c(1:3)], method='pearson')

# hmp <- Heatmap(x2,
#                name="metaPheCor",
#                column_dend_side="bottom",
#                row_dend_side="right",
#                row_names_gp = gpar(fontsize = 7),
#                rect_gp = gpar(col = "grey", lwd =0.2),
#                col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
#                column_names_gp = gpar(fontsize = 7))

# hmpout2 <- sprintf("%s/metaphe_corhmp.pdf", inpdir)
# cairo_pdf(hmpout2, width=8, height=8)
# draw(hmp, heatmap_legend_side = "left", padding=unit(c(2,2,2,2), "mm"))
# dev.off()
# system(sprintf("open %s", hmpout2))

library(wordcloud)
library(RColorBrewer)

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId)
typeColors <- brewer.pal(length(ehrFeatInfo), "Dark2")
names(typeColors) <- names(ehrFeatInfo)

cormeth <- "pearson"
trait <- "ptsd"

cormeth <- "pearson"

options(warn=2)

figid <- 7

for(trait in c("alz","bip","scz","ptsd")) {

  traitPat <- switch(trait,
                     alz="3310_A",
                     bip="29680_B",
                     scz="^29590",
                     ptsd="^30981")

  x1 <- t(phi_normalized[grep(traitPat, phi_normalized$phe),-c(1:3)])

  set.seed(1234)
  bgcor <- t(cor(sample(x1,length(x1)), t(phi_normalized[,-c(1:3)]), method=cormeth))

  x1 <- t(cor(x1, t(phi_normalized[,-c(1:3)]), method=cormeth))
  x1 <- data.frame(type=phi_normalized$type, phe=phi_normalized$phe, spmcor=as.numeric(x1))

  thres <- quantile(bgcor, 1 - 1/100)

  wcout <- sprintf("/results/FigS%s_comorbid_%s.pdf", figid, trait)
    
  figid <- figid + 1

  print(sum(x1$spmcor>thres))

  s <- 7
  
  maxwords <- 100
  
  sigidx <- which(x1$spmcor > thres)
  
  tmp1 <- t(phi_normalized[grep(traitPat, phi_normalized$phe),-c(1:3)])
  tmp2 <- t(phi_normalized[,-c(1:3)])
  
  sigscore <- sapply(sigidx, function(i) {
    # -log10(cor.test(tmp1, tmp2[,i])$p.value)
    cor.test(tmp1, tmp2[,i])$statistic
  })
  
  sigscore[!is.finite(sigscore)] <- max(sigscore[is.finite(sigscore)])
  sigscore[sigscore>500] <- 500
  
  x1 <- x1[sigidx,]
  x1$score <- sigscore
  
  x1 <- arrange(x1, score, decreasing=TRUE)
  
  pheName <- sprintf("%s:%s", as.character(x1$type), as.character(x1$phe))
  
  featName <- sapply(strsplit(pheName,":"), function(x) {
    
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
  
  featId <- sapply(strsplit(pheName,":"), function(x) {
    
    datatype <- x[1]
    
    if(datatype %in% c("icd_cpt","icd_cm","lab")) {
      
      tmp <- unlist(strsplit(x[2], "_"))[1]
      
      sprintf("(%s)", tmp)
      
    } else if(datatype == "drg") {
      
      tmp <- paste(gsub(" ","",unlist(strsplit(unlist(x)[2], "_"))[1:2]),collapse = ",")
      
      tmp2 <- paste(tail(unlist(strsplit(unlist(x)[2], "_")),2),collapse = ",")
      
      sub(",NA,NA","",sprintf("(%s,%s)", tmp,tmp2))
      
    } else if (datatype=="presc") {
      
      tmp <- paste(tail(unlist(strsplit(x[2], "_")),2),collapse=",")
      
      sprintf("(%s)", tmp)
      
    } else if(datatype=="notes") {
      
      ""
    }
  })
  
  x1$phe <- sprintf("%s %s", featName, featId)

  cairo_pdf(wcout, width=s, height=s)

  set.seed(1234)

  tmp3 <- try(wordcloud(words=as.character(x1$phe), 
                        freq=x1$score,
                        min.freq=thres,
                        max.words = maxwords,
                        rot.per=0,
                        random.order = FALSE, # rot.per=0,
                        ordered.colors=TRUE,
                        colors=rep(typeColors[x1$type])))

  dev.off()

  while(length(grep("Error", tmp3))>0) {

    s <- s + 1

    cairo_pdf(wcout, width=s, height=s)

    set.seed(1234)

    tmp3 <- try(wordcloud(words=as.character(x1$phe), 
                          freq=x1$score,
                          min.freq=thres,
                          max.words = maxwords,
                          rot.per=0,
                          random.order = FALSE, # rot.per=0,
                          ordered.colors=TRUE,
                          colors=rep(typeColors[x1$type])))
    print(s)

    dev.off()
  }
  
  
  #### barplot ####
  if(trait %in% c("scz","ptsd")) {
      ggout <- sprintf("/results/Fig2c_comorbid_%s_barplot.eps",  trait)      
  }
      
  ggout <- sprintf("/results/FigS%s_comorbid_%s_barplot.eps",  figid, trait)
  figid <- figid + 1
      
  df <- head(x1, 50)
  
  df$phe <- as.character(df$phe)
  
  df$phe <- factor(df$phe, levels=rev(df$phe))
  
  ggbar <- ggplot(df, aes(x=phe, y=score, fill=type)) + 
    geom_bar(stat='identity') + 
    scale_x_discrete(position='top') + 
    scale_y_reverse() + 
    theme(legend.position = "left") +
    xlab("") + ylab("Correlation score") +
    scale_fill_manual(values=typeColors[df$type]) +
    coord_flip()
  
  ggsave(ggout, ggbar, width=8, height=6)  
}












































