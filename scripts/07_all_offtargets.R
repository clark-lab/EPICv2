library(ggplot2)
library(GenomicRanges)
library(effsize)

#
load("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/Mvalues_NODROPOUT.RData")
rm(EPICv1_Ms)
offtargets <- read.csv("~/Projects/EPIC_V2/Tim_Analysis/Offtargets_FINAL.csv")
otgr <- GRanges(paste(offtargets$CHR, offtargets$MAPINFO, sep=":"))
names(otgr) <- offtargets$ProbeID
values(otgr) <- offtargets[,4:6]
otgr <- otgr[otgr %over% CpGs.hg38]
require(parallel)
otgr <- otgr[-which(!names(otgr) %in% rownames(EPICv2_Ms)),]
EPICnumvals <- mclapply(names(otgr), function (x) sum(!is.na(EPICv2_Ms[x,])), mc.cores = 16)
otgr$EPICnonas <- unlist(EPICnumvals)
save(otgr, file="~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/OTGR.RData")
names(CpGs.hg38) <- paste(CpGs.hg38)
olot <- findOverlaps(otgr, CpGs.hg38)
otgr <- otgr[queryHits(olot),]
values(otgr) <- cbind(values(otgr), values(CpGs.hg38)[subjectHits(olot),])
nas <- apply(values(otgr)[,5:22], 1, function (x) sum(!is.na(x)))
keep <- nas >= 12 & otgr$EPICnonas >=12
otgr <- otgr[keep]
save(otgr, file="~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/OTGR.RData")

#Now do the target CpGs from WGBS
otmanifest <- new_manifest[names(otgr),]
manifestgr <- GRanges(paste(otmanifest$CHR, otmanifest$MAPINFO, sep=":"))
names(manifestgr) <- names(otgr)
olom <- findOverlaps(manifestgr, CpGs.hg38)
manifestgr <- manifestgr[queryHits(olom)]
values(manifestgr) <- values(CpGs.hg38)[subjectHits(olom),]
nas <- apply(values(manifestgr), 1, function (x) sum(!is.na(x)))
keep <- nas >=12
manifestgr <- manifestgr[keep]
manifestgr <- manifestgr[!duplicated(names(manifestgr))]
otgr <- otgr[names(otgr) %in% names(manifestgr)]
save(manifestgr, otgr, EPICv2_Ms, file="~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/Readyfor_offtargets.RData")

getmses <- function(probe, plot=FALSE){
  message(probe)
  epic <- as.numeric(as.matrix(EPICv2_Ms[probe,]))
  wgbson <- as.numeric(as.matrix(values(manifestgr[probe,])))
  wgbsoffs <- as.matrix(values(otgr[names(otgr) %in% probe, 5:22]))
  offsgr <- otgr[names(otgr) %in% probe]
  values(offsgr) <- NULL
  #if(any(is.na(c(epic, wgbson, wgbsoffs)))){
  #  nas <- union(which(is.na(epic)), which(is.na(wgbson)))
  #  nas <- union(nas, which(apply(wgbsoffs, 2, function (x) any(is.na(x)))))
  #  epic <- epic[-nas]
  #  wgbson <- wgbson[-nas]
  #  wgbsoffs <- wgbsoffs[,-nas,drop=FALSE]
  #}
  stopifnot(length(epic)==length(wgbson))
  stopifnot(length(epic)==ncol(wgbsoffs))
  
  onmse <- mean((wgbson-epic)^2, na.rm=T)
  offmses <- apply(wgbsoffs, 1, function (x) mean((x-epic)^2, na.rm=T))
  offsgr$offmses <- offmses
  
  if(plot){
    ylimit <- range(na.omit(c(epic, wgbson, wgbsoffs)))
    plot(epic, wgbson, ylim=ylimit, pch=16, col="deepskyblue", xlab="M-values, EPICv2", ylab="M-values, WGBS", main=probe)
    abline(0, 1, lwd=2)
    abline(lm(wgbson~epic), col="deepskyblue", lwd=2)
    for (i in 1:nrow(wgbsoffs)){
      points(epic, wgbsoffs[i,], pch=16, col="red")
      abline(lm(wgbsoffs[i,]~epic), col="red", lwd=2)
    }
    legend("topleft", c("Target CpG", "Off-target CpGs", "Expected value"), text.col=c("deepskyblue", "red", "black"), bty="n")
  }
  
  list(ONMSE=onmse, OFFMSE=offsgr)
}

otmses <- mclapply(names(manifestgr), getmses, mc.cores = 16)
manifestgr$ONRMSE <- sqrt(unlist(lapply(otmses, function (x) x$ONMSE)))
manifestgr$MINOFFMSE <- sqrt(unlist(lapply(otmses, function (x) min(x$OFFMSE$offmses))))
manifestgr$lowestofftarget <- unlist(lapply(1:length(otmses), function (x) paste(otmses[[x]]$OFFMSE[which.min(otmses[[x]]$OFFMSE$offmses)])))

offtargetres <- manifestgr
values(offtargetres) <- values(offtargetres)[,19:21]
offtargetres$CH_wgbs_evidence <- offtargetres$ONRMSE > offtargetres$MINOFFMSE
wgbsev <- names(offtargetres)[offtargetres$ONRMSE > offtargetres$MINOFFMSE]
likelyCH <- offtargetres[wgbsev,"lowestofftarget"]

save(offtargetres, otmses, file="~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/Offtarget_results.RData")

#Figure 6d
smoothScatter(sqrt(offtargetres$ONRMSE), sqrt(offtargetres$MINOFFMSE), xlab="RMSE with WGBS, target CpG site",
              ylab="Minimum RMSE with WGBS, off-target CpG sites", colramp = colorRampPalette(brewer.pal(9, "Greens")))
title("Cross-hybridising probe preference, 17928 probes")
abline(0, 1, lwd=2, lty=2)
text(1.3, 3, "34% of probes have a higher MSE\nwhen regressed against the off-target")
text(2.5, 0.8, "66% of probes have a higher MSE when\nregressed against the manifest target")

new_manifest$CH_BLAT <- ifelse(new_manifest$IlmnID %in% offtargets$ProbeID, "Y", "")
new_manifest$CH_WGBS_evidence <- ifelse(new_manifest$IlmnID %in% wgbsev, "Y", "")
new_manifest$Likely_offtarget <- ""
m <- match(names(likelyCH), rownames(new_manifest))
new_manifest$Likely_offtarget[m] <- likelyCH$lowestofftarget
new_manifest$Num_offtargets <- unlist(mclapply(new_manifest$IlmnID, function (x) sum(offtargets$ProbeID == x), mc.cores = 16))
new_manifest$Suggested_offtarget <- new_manifest$Likely_offtarget
new_manifest <- new_manifest[,-75]
save(new_manifest, file="new_manifest_with_CH.RData")

stopifnot(all(rownames(EPICv2_Ms)==rownames(new_manifest)))
manifestgr <- GRanges(paste(new_manifest$CHR, new_manifest$MAPINFO, sep=":"))
names(manifestgr) <- paste(manifestgr)
olon <- findOverlaps(manifestgr, CpGs.hg38)
keepprobename <- rownames(EPICv2_Ms)[queryHits(olon)]
manifestgr <- manifestgr[queryHits(olon),]
EPICv2_forcor <- EPICv2_Ms[queryHits(olon),]
values(manifestgr) <- values(CpGs.hg38)[subjectHits(olon),]
epicnas <- apply(EPICv2_forcor, 1, function (x) sum(!is.na(x)))
wgbsnas <- apply(values(manifestgr), 1, function (x) sum(!is.na(x)))
keep <- epicnas >= 12 & wgbsnas >= 12
EPICv2_forcor <- EPICv2_forcor[keep,]
manifestgr <- manifestgr[keep]

stopifnot(all(paste(new_manifest[rownames(EPICv2_forcor), "CHR"], new_manifest[rownames(EPICv2_forcor), "MAPINFO"], sep=":")==names(manifestgr)))


getontarget <- function(idx){
  epic <- as.numeric(EPICv2_forcor[idx,])
  wgbs <- as.numeric(as.matrix(values(manifestgr[idx])))
  mean((wgbs-epic)^2, na.rm=T)
}

mses <- unlist(mclapply(1:nrow(EPICv2_forcor), getontarget, mc.cores = 16))
rmses <- sqrt(mses)
names(rmses) <- rownames(EPICv2_forcor)
new_manifest$RMSE_W_WGBS <- NA 
new_manifest[names(rmses), "RMSE_W_WGBS"] <- rmses

CH_type <- ifelse(new_manifest$CH_WGBS_evidence=="Y", "WGBS + BLAT", ifelse(new_manifest$CH_BLAT=="Y", "BLAT only", "No evidence"))
#Anova to determine p-values
data <- data.frame(Probe=new_manifest$IlmnID, RMSE=new_manifest$RMSE_W_WGBS, Type=factor(CH_type))
data <- data[!is.na(data$RMSE),]
model <- aov(RMSE~Type, data=data)
c1 <- c(0.5, -1, 0.5)
c2 <- c(-1, 0, 1)
mat <- cbind(c1, c2)
contrasts(data$Type) <- mat
model1 <- aov(RMSE~Type, data=data)
j <- summary.aov(model1, split=list(Type=list("CHvsNCH"=1, "WGBSvsBLATonly"=2)))
#P-values are extremely low, try Cohen d
cohen.d(data$RMSE[data$Type!="No evidence"], data$RMSE[data$Type=="No evidence"])
#Cohen's d

#d estimate: 2.663204 (large)
#95 percent confidence interval:
#   lower    upper 
#2.650413 2.675996 

cohen.d(data$RMSE[data$Type=="WGBS + BLAT"], data$RMSE[data$Type=="BLAT only"])

#Cohen's d

#d estimate: 0.9314771 (large)
#95 percent confidence interval:
#    lower     upper 
#0.9061261 0.9568282 
data$Type <- relevel(data$Type, "No evidence")

#Figure 6e
vp1 <-  ggplot(data, aes(x=Type, y=RMSE, fill=Type)) +
  geom_violin() + labs(y="RMSE, EPICv2 vs. WGBS", x="Evidence for cross-hybridisation") + 
  ggtitle("Effect of cross-hybridisation on similarity to WGBS") +
  scale_y_continuous(trans = "log", labels=function(x) sprintf("%.2f", x)) +
  stat_summary(mapping = aes(x=Type, y=RMSE, fill=Type),
               fun.min = function(z) { quantile(z,0.25) },
               fun.max = function(z) { quantile(z,0.75) },
               fun = median) +
  theme(legend.position = "none")
vp1

#Plot a few bad ones
data <- data[data$Type=="WGBS + BLAT",]
names(otmses) <- unlist(lapply(otmses, function (x) names(x$OFFMSE)[1]))
data$minofftarget <- unlist(lapply(otmses[data$Probe], function (x) sqrt(min(x$OFFMSE$offmses))))
data$diff <- data$RMSE - data$minofftarget
data <- data[order(data$diff, decreasing = T),]
load("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/Readyfor_offtargets.RData")
m <- match(data$Probe, rownames(new_manifest))

data$numoff <- new_manifest$Num_offtargets[m]
dataunder10 <- data[data$numoff < 10,]
load("~/hg38.WGBS.Rdata")
#Plot some bad probes

#Figures 6f - 6j
pdf("Plots/badofftargets_regression.pdf", 14, 10)
par(mfrow=c(3, 3))
sapply(1:180, function (x) getmses(dataunder10$Probe[x], plot=TRUE))
dev.off()

#Add "Is_CpG to offtargets"
length(otgr)

#Pie charts
#How many brought back from 450K
chmf <- new_manifest[new_manifest$CH_BLAT=="Y",]
chmf$EPICv1 <- !is.na(chmf$EPICv1locmatch)
chmf$K450 <- !is.na(chmf$K450locmatch) 
chmf$K27 <- !is.na(chmf$K27locmatch)

chmf$history <- ""
chmf$history[!chmf$EPICv1 & !chmf$K450 & !chmf$K27] <- "New probe"
chmf$history[chmf$EPICv1] <- "On EPICv1"
chmf$history[!chmf$EPICv1 & chmf$K450] <- "On 450K"
chmf$history[!chmf$EPICv1 & !chmf$K450 & chmf$K27] <- "On 27K"

tb <- table(chmf$history)
tbpc <- round(table(chmf$history)*100/nrow(chmf))
tbpc

#New probe    On 27K   On 450K On EPICv1 
#       81         0         5        14 

#Fig 6a
library(plotrix)
library(RColorBrewer)
pie3D(tb, labels=c("New probe (24929)", "On 27K (25)", "On 450K (1437)", "On EPICV1 (4236)"),
      labelcex=1.2, explode=0.2, col=brewer.pal(4, "Accent"), mar=c(2, 2, 2, 2), theta=pi/3,
      main="Cross-hybridising probe history")

#Fig 6b
nttable <- table(offtargets$base)
pie3D(nttable, labels=c("A (3%)", "C (70%)", "G\n(3%)", "N\n(0%)", "T (24%)"),
      labelcex=1.2, explode=0.2, col=c("lightyellow", "cyan", "lightcoral", "gray", "lightgreen"), mar=c(2, 2, 2, 2), theta=pi/3,
      main="Off-target nucleotide")


#Fig 6c
cpgtable <- c(1629882, 2908181-1629882)
names(cpgtable) <- c("CpG (56%)", "CpH (44%)")
pie3D(cpgtable, labels=names(cpgtable),
      labelcex=1.2, explode=0.2, col=c("lightcoral", "grey70"), mar=c(2, 2, 2, 2), theta=pi/3,
      main="Off-target cytosine")