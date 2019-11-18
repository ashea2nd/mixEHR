library(Rcpp)
library(gridExtra)
library(plyr)
library(ggplot2)
library(reshape)

rm(list=ls(all=TRUE))

K <- 75
niter <- 5497

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId)

inpdir <- "/data/precomputed_results/mixmimic"

filein_phi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_phi_normalized.csv", inpdir, K, niter)
phi <- read.csv(filein_phi, header=F)
colnames(phi)[c(1,2)] <- c("typeId","pheId")
colnames(phi)[-c(1,2)] <- sprintf("M%s", 1:(ncol(phi)-2))


filein_psi <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_psi.csv", inpdir, K, niter)
inpdir <- dirname(filein_phi)
psi <- read.csv(filein_psi, header=F)
colnames(psi)[c(1:2)] <- c("typeId","pheId")
colnames(psi)[-c(1:2)] <- sprintf("M%s", 1:(ncol(psi)-2))

tmp <- dimnames(psi)
psi <- scale(psi, center=F, scale=colSums(psi))
dimnames(psi) <- tmp

phidf <- melt(phi[,-2], id.vars = "typeId")
philist <- split(phi, phi$typeId)

x <- philist[[1]]

phi <- do.call(rbind, lapply(philist, function(x) {
  
  featId <- as.matrix(x[,c(1,2)])
  
  tpid <- as.character(featId[1,1])
  
  ehrFeatInfo_t <- ehrFeatInfo[[as.numeric(x[1,"typeId"])]]
  
  x <- as.matrix(x[,-c(1,2)])  
  
  rownames(x) <- sprintf("%s:%s",
                         names(ehrFeatInfo)[ehrFeatInfo_t$typeId][1],
                         ehrFeatInfo_t$pheName[match(featId[,"pheId"], ehrFeatInfo_t$pheId)])
  x
}))


featName <- sapply(strsplit(rownames(phi),":"), function(x) {
  
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

featId <- sapply(strsplit(rownames(phi),":"), function(x) {
  
  datatype <- x[1]
  
  if(datatype %in% c("icd_cpt","icd_cm","lab")) {
    
    tmp <- unlist(strsplit(x[2], "_"))[1]
    
    sprintf("(%s)", tmp)
    
  } else if(datatype == "drg") {
    
    tmp <- paste(gsub(" ","",unlist(strsplit(unlist(x)[2], "_"))[1:2]),collapse = ",")
        
    sprintf("(%s)", tmp)
    
  } else if (datatype=="presc") {
    
    tmp <- paste(tail(unlist(strsplit(x[2], "_")),2),collapse=",")
    
    sprintf("(%s)", tmp)
    
  } else if(datatype=="notes") {
    
    ""
  }
})

datatypes <- sapply(strsplit(rownames(phi),":"), function(x) x[[1]])

rownames(phi) <- sprintf("%s:%s %s", datatypes, featName, featId)

topW <- 500

topFeat <- lapply(colnames(phi), function(k) {
  
  x <- phi[,k,drop=F]
  
  x[head(order(x, decreasing = TRUE),topW),,drop=F]
})

names(topFeat) <- colnames(phi)

x <- topFeat[[1]]

library(wordcloud)

typeColors <- brewer.pal(length(ehrFeatInfo), "Dark2")
names(typeColors) <- names(ehrFeatInfo)

k <- 1
K <- length(topFeat)

outdir <- sprintf("%s/disease_wordclouds_composite", dirname(filein_phi))
if(!dir.exists(outdir)) dir.create(outdir)

options(warn=2)

sel <- c("M31", "M35", "M50")
# sel <- colnames(phi) # uncomment for plotting all 75 topic word clouds                    

for(k in sel) {
  
  print(k)
  
  x <- topFeat[[k]]
  typeIds <- sapply(strsplit(rownames(x),":"), function(x) x[[1]])
  pheIds <- sapply(strsplit(rownames(x),":"), function(x) x[[2]])
  # pheIds <- sapply(strsplit(pheIds,"_"), function(x) tail(x,1))
  
  tmp <- split(pheIds, typeIds)
  tmp2 <- split(x, typeIds)
  
  wcout <- sprintf("/results/Fig3b_%s_phi.eps", k)
  
  thres <- quantile(phi[,k],1-1e-3)

  s <- 7
  
  cairo_ps(wcout, width=s, height=s)
  
  set.seed(1234)
  
  tmp3 <- try(wordcloud(words=pheIds, 
                        freq=as.numeric(x),
                        min.freq=thres,
                        random.order = FALSE,
                        rot.per=0,
                        ordered.colors=TRUE,
                        # use.r.layout=TRUE,
                        colors=rep(typeColors[typeIds])))
  
  dev.off()
  
  while(length(grep("Error", tmp3))>0) {
    
    s <- s + 1
    
    cairo_ps(wcout, width=s, height=s)
    
    tmp3 <- try(wordcloud(words=pheIds, 
                          freq=as.numeric(x),
                          min.freq=thres,
                          random.order = FALSE, 
                          rot.per=0,
                          ordered.colors=TRUE,
                          # use.r.layout=TRUE,
                          colors=rep(typeColors[typeIds])))
    print(s)
    
    dev.off()
  }  
}

























