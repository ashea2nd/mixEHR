library(ggplot2)
library(plyr)
library(gridExtra)
library(glmnet)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize) # colorRamp2
library(reshape)

rm(list=ls(all=TRUE))

featId <- "/data/Mimic/version_1_3/processed/ehrFeatId.RData"
load(featId)

source('/code/main/eval_fun.R')

K <- 75

inpdir <- "/data/precomputed_results/mortality"
deathrda <- sprintf("%s/predictDeath_mixehr_K%s.RData", inpdir, K)
load(deathrda)

effectSizes <- as.matrix(coef(fitFull, s="lambda.min"))[-1]

resdir <- "/data/results/mixmimic"

#### death coefficients ####
df <- data.frame(topic=1:K, effectSize=effectSizes)
df <- subset(df, effectSize != 0)
df <- arrange(df, effectSize, decreasing=FALSE)
df$topic <- factor(df$topic, levels=df$topic)

df$type <- ""
df$type[df$topic %in% tail(df$topic,3)] <- "top4"
df$type[df$topic %in% head(df$topic,1)] <- "bottom1"

gg <- ggplot(df, aes(x=topic, y=effectSize, fill=type)) + 
  geom_bar(stat='identity') + coord_flip() + 
  geom_text(aes(label=topic)) +
  xlab("Mortality topic coefficients") + 
  ylab("Linear coefficients for mortality") +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  scale_fill_manual(values=c("black","green","cyan")) + theme(legend.position = "none")

ggout <- "/results/Fig5b_deathCoef.eps"
ggsave(ggout, gg, width=3.6, height=2.4)



#### phi,psi,eta from diseaseToipcsExport.R ####
K <- 75
niter <- 5497
resdir <- "/data/precomputed_results/mixmimic"
resrda <- sprintf("%s/mimic_trainData_JCVB0_nmar_K%s_iter%s_combined.RData", resdir, K, niter)
load(resrda)

df <- arrange(df, effectSize, decreasing=TRUE)
df$topic <- as.numeric(as.character(df$topic))

metaphe_sel <- c(head(df$topic, 3), tail(df$topic, 1))

phisel <- phi[, metaphe_sel]

topW <- 5

topphe <- unique(unlist(lapply(1:ncol(phisel), function(k) {
  head(rownames(phisel)[order(phisel[,k], decreasing = T)], topW)
})))


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
               column_title="Meta-phenotypes") + groupAnnot

myw <- 10
myh <- 4

hmpout <- sprintf("/results/Fig5b_phiDeathTopicHmp.eps")
cairo_ps(hmpout, width=myw, height=myh)
draw(hmp, heatmap_legend_side = "right", padding=unit(c(5,110,5,20), "mm")) # bottom, left, top, right
dev.off()



#### lab abnormal results under topics ####
etaNrm <- subset(eta, stateId==0)[,-c(1:3)]
etaAbn <- subset(eta, stateId==1)[,-c(1:3)]

etaAbnPsi <- etaAbn[,metaphe_sel] * psi[,metaphe_sel]
# etaAbnPsi <- psi[,metaphe_sel]

etaAbnPsi <- scale(etaAbnPsi, F, scale=colSums(etaAbnPsi))

rownames(etaAbnPsi) <- rownames(etaNrm) <- rownames(etaAbn) <- rownames(psi)

topL <- 5

toplab <- unique(as.numeric(sapply(1:ncol(etaAbnPsi), function(k) 
  head(order(etaAbnPsi[,k], decreasing = TRUE), topL))))

etaAbnPsi <- as.matrix(etaAbnPsi[toplab,])

hmp <- Heatmap(etaAbnPsi,
               cluster_columns=FALSE, 
               cluster_rows=FALSE,
               row_names_side = "left",
               row_dend_side = "right",
               rect_gp = gpar(col = "grey", lwd =0.2),
               # col = colorRamp2(c(0, 0.1), c("white", "red")),
               col = c("white", "red"),
               name="Prob",
               column_title="Meta-phenotypes")

myw <- 7.2
myh <- 4

# hmpout2 <- sprintf("/results/FigS_labDeathTopicHmp.eps")
# cairo_ps(hmpout2, width=myw, height=myh)
# draw(hmp, heatmap_legend_side = "right", padding=unit(c(5,80,5,10), "mm")) # bottom, left, top, right
# dev.off()





































