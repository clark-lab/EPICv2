library(sesame)
library(minfi)
library(limma)
library(data.table)
library(consensus)
library(rtracklayer)
library(plyr)

setwd("~/hippo/Gagri_Immunogenomics/Tim_Peters/EPICv2/")

#Crossplatform data

load("~/hippo/Gagri_Cancer-Epigenetics/Google_Drive/Ruth_Projects/EPIC_V2/Tim_Analysis/Crossplatform/Mvalues_NODROPOUT.RData")

#Load manifest
new_manifest <- read.csv("~/hippo/Gagri_Immunogenomics/Tim_Peters/EPICv2/AdditionalFile4.csv", header = T)
rownames(new_manifest) <- new_manifest$IlmnID

names(CpGs.hg38) <- paste(CpGs.hg38)
namedups <- unique(new_manifest$Name[duplicated(new_manifest$Name)])
namemf <- new_manifest[new_manifest$Name %in% namedups,]
namemf$namegroups <- sapply(namemf$Name, function (x) paste(sort(namemf$IlmnID[namemf$Name==x]), collapse=";"))

seqmf <- new_manifest[new_manifest$seqrep_IlmnIDs!="",]
seqmf$namegroups <- sapply(1:nrow(seqmf), function (x) paste(sort(unlist(strsplit(seqmf$seqrep_IlmnIDs[x], ";"))), collapse=";"))

posmf <- new_manifest[new_manifest$posrep_IlmnIDs!="",]
posmf$namegroups <- sapply(1:nrow(posmf), function (x) paste(sort(unlist(strsplit(posmf$posrep_IlmnIDs[x], ";"))), collapse=";"))


#Name - consensus

namemf$namevector <- ""
namemfcons <- namemf[!is.na(namemf$EPICv1probeID),]

runconsensusmean <- function(probe){
  message(probe)
  epic1 <- EPICv1_Ms[probe,]
  epic1na <- which(is.na(epic1))
  epic2idxs <- new_manifest$IlmnID[new_manifest$Name==probe]
  epic2 <- EPICv2_Ms[epic2idxs,]
  epic2 <- rbind(epic2, mean=colMeans(epic2, na.rm = T))
  epic2na <- which(apply(epic2, 2, function (x) any(is.na(x))))
  submf <- new_manifest[epic2idxs,]
  coords <- unique(paste(submf$CHR, submf$MAPINFO, sep=":"))
  wgbs <- as.numeric(as.matrix(values(CpGs.hg38)[coords,]))
  wgbsna <- which(is.na(wgbs))
  
  remove <- sort(union(wgbsna, union(epic1na , epic2na)))
  if(length(remove) > 0){
    epic1 <- epic1[-remove]
    wgbs <- wgbs[-remove]
    epic2 <- epic2[,-remove]
  }
  
  if(length(epic1) < 12){
    return(NA)
  }
  
  ###Replicate
  epic1 <- matrix(rep(epic1, nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
  rownames(epic1) <- rownames(epic2)
  colnames(epic1) <- colnames(epic2)
  wgbs <- matrix(rep(wgbs, nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
  rownames(wgbs) <- rownames(epic2)
  colnames(wgbs) <- colnames(epic2)
  
  
  
  pfnames=c("EPICv1", "EPICv2", "WGBS")
  EPICmm <- MultiMeasure(names=pfnames, data=list(epic1, epic2, wgbs))
  fit <- fitConsensus(EPICmm)
  fit
}

namefits <- mclapply(unique(namemfcons$Name), runconsensusmean, mc.cores = 16)
names(namefits) <- unique(namemfcons$Name)
nas <- unlist(lapply(namefits, is.na))
nas <- names(nas)[nas]
namemf$namevector[namemf$Name %in% nas] <- "Insufficient evidence"
namefits <- namefits[!names(namefits) %in% nas]

namedfres <- data.frame(EPICv2name=unlist(lapply(namefits, function (x) rownames(x@a_i))),
                        v1sens=unlist(lapply(namefits, function (x) x@b_i[,1])),
                        v2sens=unlist(lapply(namefits, function (x) x@b_i[,2])),
                        v1prec=unlist(lapply(namefits, function (x) x@d_i[,1])),
                        v2prec=unlist(lapply(namefits, function (x) x@d_i[,2])))
namedfres$EPICv1name <- substr(rownames(namedfres), 1, 10)
namedfres$namegroup <- sapply(namedfres$EPICv1name, function (x) paste(namedfres$EPICv2name[gsub("_.*", "", namedfres$EPICv2name)==x & namedfres$EPICv2name!="mean"], collapse=";"))
stopifnot(all((namedfres$namegroup %in% namemf$namegroups)))
namedfres$namegroup <- paste(namedfres$namegroup, "mean", collapse=";")
apply(namedfres[,2:5], 2, function (x) any(is.na(x)))
#v1sens v2sens v1prec v2prec 
#FALSE  FALSE  FALSE  FALSE 

namedfres$namegroup <- paste(namedfres$namegroup, "mean", sep=";")
namemf$namegroups <- paste(namemf$namegroups, "mean", sep=";")

getrank <- function(df, probegroup){
  sub <- df[df$namegroup==probegroup,]
  rownames(sub) <- sub$EPICv2name
  sub <- sub[sort(unlist(strsplit(sub$namegroup[1], ";"))),]
  sub$sens <- rank(-sub$v2sens)
  sub$prec <- rank(sub$v2prec)
  sub
}

nameranks <- lapply(unique(namedfres$namegroup), function (x) getrank (namedfres, x))
finalname <- rbind.fill(nameranks)
meannameresults <- finalname[finalname$EPICv2name=="mean",]
bestmeanname <- meannameresults[meannameresults$prec==1,]
nrow(bestmeanname)
#1369
meannameprobes <- unlist(sapply(bestmeanname$namegroup, function (x) strsplit(x, ";")))
meannameprobes <- meannameprobes[meannameprobes!="mean"]
length(meannameprobes)
#3078

new_manifest$Best_precision_as_group_mean_NAME <- new_manifest$IlmnID %in% meannameprobes
new_manifest$Best_precision_as_group_mean_NAME[!new_manifest$IlmnID %in% namemfcons$IlmnID] <- NA

table(new_manifest$Best_precision_as_group_mean_NAME)/nrow(namemfcons)
#FALSE     TRUE 
#0.634572 0.365428

#Then WGBS
#Then use just WGBS to do the walkoff
load("~/hippo/Gagri_Cancer-Epigenetics/Google_Drive/Ruth_Projects/EPIC_V2/Tim_Analysis/Crossplatform/Readyfor_offtargets.RData")
namewgbs <- namemf[grep("WGBS", namemf$Rep_results_by_NAME),]

getmeanrmse <- function(namegroup){
  sub <- namewgbs[namewgbs$namegroups==namegroup, c("IlmnID", "Name", "RMSE_with_WGBS", "namegroups", "namevector")]
  epic <- colMeans(as.matrix(EPICv2_Ms[sub$IlmnID,]))
  coords <- unique(sapply(sub$IlmnID, function (x) paste(new_manifest[x, c("CHR", "MAPINFO")], collapse=":")))
  wgbs <-  as.matrix(values(CpGs.hg38[coords,]))
  rmses <- apply(wgbs, 1, function (x) mean((x-epic)^2, na.rm=T))
  if(length(rmses) > 1){return(rmses)
  } else {
    sub <- rbind(sub, c(IlmnID="mean", Name=sub$Name[1], RMSE_with_WGBS=rmses, namegroups=sub$namegroups[1], namevector=sub$namevector[1]))
    return(sub) 
    }
 }

wgbsnameres <- mclapply(unique(namewgbs$namegroups), getmeanrmse, mc.cores = 10)
wgbsnameres <- rbind.fill(wgbsnameres)

getrank <- function(df, probegroup){
  sub <- df[df$namegroup==probegroup,]
  sub$rmserank <- rank(sub$RMSE_with_WGBS)
  sub
}

rankres <- lapply(unique(wgbsnameres$namegroups), function (x) getrank(wgbsnameres, x))
wgbsnameres <- rbind.fill(rankres)
wgbsmeanres <- wgbsnameres[wgbsnameres$IlmnID=="mean",]
wgbsmeanbest <- wgbsmeanres[wgbsmeanres$rmserank==1,]
wgbsmeanprobes <- unlist(strsplit(wgbsmeanbest$namegroups, ";"))
wgbsmeanprobes <- wgbsmeanprobes[wgbsmeanprobes!="mean"]
new_manifest$Rep_results_by_NAME[new_manifest$IlmnID %in% wgbsmeanprobes] <- "Superior group mean (by WGBS)"




#sequence - consensus

#seq - consensus

seqmf$seqvector <- ""
seqmfcons <- seqmf[!is.na(seqmf$EPICv1probeID),]

runconsensusmean <- function(namegroup){
  message(namegroup)
  epic1 <- EPICv1_Ms[unique(seqmfcons$Name[seqmfcons$namegroups==namegroup]),]
  epic1na <- which(is.na(epic1))
  epic2idxs <- unlist(strsplit(namegroup, ";"))
  epic2 <- EPICv2_Ms[epic2idxs,]
  epic2 <- rbind(epic2, mean=colMeans(epic2, na.rm = T))
  epic2na <- which(apply(epic2, 2, function (x) any(is.na(x))))
  submf <- new_manifest[epic2idxs,]
  coords <- unique(paste(submf$CHR, submf$MAPINFO, sep=":"))
  wgbs <- as.numeric(as.matrix(values(CpGs.hg38)[coords,]))
  wgbsna <- which(is.na(wgbs))
  
  remove <- sort(union(wgbsna, union(epic1na , epic2na)))
  if(length(remove) > 0){
    epic1 <- epic1[-remove]
    wgbs <- wgbs[-remove]
    epic2 <- epic2[,-remove]
  }
  
  if(length(epic1) < 12){
    return(NA)
  }
  
  ###Replicate
  epic1 <- matrix(rep(epic1, nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
  rownames(epic1) <- rownames(epic2)
  colnames(epic1) <- colnames(epic2)
  wgbs <- matrix(rep(wgbs, nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
  rownames(wgbs) <- rownames(epic2)
  colnames(wgbs) <- colnames(epic2)
  
  pfnames=c("EPICv1", "EPICv2", "WGBS")
  EPICmm <- MultiMeasure(names=pfnames, data=list(epic1, epic2, wgbs))
  fit <- fitConsensus(EPICmm)
  fit
}

seqfits <- mclapply(unique(seqmfcons$namegroups), runconsensusmean, mc.cores = 16)
names(seqfits) <- unique(seqmfcons$namegroups)
nas <- unlist(lapply(seqfits, is.na))
nas <- names(nas)[nas]
seqmf$seqvector[seqmf$namegroups %in% nas] <- "Insufficient evidence"
seqfits <- seqfits[!names(seqfits) %in% nas]

seqdfres <- data.frame(Probeseq=unlist(lapply(seqfits, function (x) rownames(x@a_i))),
                       v1sens=unlist(lapply(seqfits, function (x) x@b_i[,1])),
                       v2sens=unlist(lapply(seqfits, function (x) x@b_i[,2])),
                       v1prec=unlist(lapply(seqfits, function (x) x@d_i[,1])),
                       v2prec=unlist(lapply(seqfits, function (x) x@d_i[,2])))
seqdfres$seqgroup <- substr(rownames(seqdfres), 1, nchar(rownames(seqdfres))-1)
seqdfres$seqgroup[grep("cg06373096_TC11", seqdfres$seqgroup)] <- "cg06373096_TC11;cg06373096_TC110;cg06373096_TC12;cg06373096_TC13;cg06373096_TC14;cg06373096_TC15;cg06373096_TC16;cg06373096_TC17;cg06373096_TC18;cg06373096_TC19"
stopifnot(all((seqdfres$seqgroup %in% seqmf$namegroups)))
apply(seqdfres[,2:5], 2, function (x) any(is.na(x)))
#v1sens v2sens v1prec v2prec 
#FALSE  FALSE  FALSE  FALSE 

seqdfres$seqgroup <- paste(seqdfres$seqgroup, "mean", sep=";")
seqmf$seqgroups <- paste(seqmf$seqrep_IlmnIDs, "mean", sep=";")

getrank <- function(df, probegroup){
  sub <- df[df$seqgroup==probegroup,]
  rownames(sub) <- sub$Probeseq
  sub <- sub[sort(unlist(strsplit(sub$seqgroup[1], ";"))),]
  sub$sens <- rank(-sub$v2sens)
  sub$prec <- rank(sub$v2prec)
  sub
}


seqranks <- lapply(unique(seqdfres$seqgroup), function (x) getrank (seqdfres, x))
finalseq <- rbind.fill(seqranks)
meanseqresults <- finalseq[finalseq$Probeseq=="mean",]
bestmeanseq <- meanseqresults[meanseqresults$prec==1,]
nrow(bestmeanseq)
#1337
meanseqprobes <- unlist(sapply(bestmeanseq$seqgroup, function (x) strsplit(x, ";")))
meanseqprobes <- meanseqprobes[meanseqprobes!="mean"]
length(meanseqprobes)
#3007

new_manifest$Best_precision_as_group_mean_SEQUENCE <- new_manifest$IlmnID %in% meanseqprobes
new_manifest$Best_precision_as_group_mean_SEQUENCE[!new_manifest$IlmnID %in% seqmfcons$IlmnID] <- NA

table(new_manifest$Best_precision_as_group_mean_SEQUENCE)/nrow(seqmfcons)
#FALSE      TRUE 
#0.6361766 0.3638234

#Sequence

seqwgbs <- seqmf[grep("WGBS", seqmf$Rep_results_by_SEQUENCE),]

getmeanrmse <- function(seqgroup){
  sub <- seqwgbs[seqwgbs$seqgroups==seqgroup, c("IlmnID", "Name", "RMSE_with_WGBS", "seqgroups", "seqvector")]
  epic <- colMeans(as.matrix(EPICv2_Ms[sub$IlmnID,]))
  coords <- unique(sapply(sub$IlmnID, function (x) paste(new_manifest[x, c("CHR", "MAPINFO")], collapse=":")))
  wgbs <-  as.matrix(values(CpGs.hg38[coords,]))
  rmses <- apply(wgbs, 1, function (x) mean((x-epic)^2, na.rm=T))
  if(length(rmses) > 1){return(rmses)
  } else {
    sub <- rbind(sub, c(IlmnID="mean", seq=sub$Name[1], RMSE_with_WGBS=rmses, seqgroups=sub$seqgroups[1], seqvector=sub$seqvector[1]))
    return(sub) 
  }
}

wgbsseqres <- mclapply(unique(seqwgbs$seqgroups), getmeanrmse, mc.cores = 10)
wgbsseqres <- rbind.fill(wgbsseqres)

getrank <- function(df, probegroup){
  sub <- df[df$seqgroup==probegroup,]
  sub$rmserank <- rank(sub$RMSE_with_WGBS)
  sub
}

rankres <- lapply(unique(wgbsseqres$seqgroups), function (x) getrank(wgbsseqres, x))
wgbsseqres <- rbind.fill(rankres)
wgbsmeanres <- wgbsseqres[wgbsseqres$IlmnID=="mean",]
wgbsmeanbest <- wgbsmeanres[wgbsmeanres$rmserank==1,]
wgbsmeanprobes <- unlist(strsplit(wgbsmeanbest$seqgroups, ";"))
wgbsmeanprobes <- wgbsmeanprobes[wgbsmeanprobes!="mean"]
new_manifest$Rep_results_by_SEQUENCE[new_manifest$IlmnID %in% wgbsmeanprobes] <- "Superior group mean (by WGBS)"


#Position - consensus

#pos - consensus

posmf$posvector <- ""
posmfcons <- posmf[!is.na(posmf$EPICv1probeID),]

runconsensusmean <- function(namegroup){
  message(namegroup)
  epic1 <- EPICv1_Ms[unique(posmfcons$Name[posmfcons$namegroups==namegroup]),]
  epic1na <- which(is.na(epic1))
  epic2idxs <- unlist(strsplit(namegroup, ";"))
  epic2 <- EPICv2_Ms[epic2idxs,]
  epic2 <- rbind(epic2, mean=colMeans(epic2, na.rm = T))
  epic2na <- which(apply(epic2, 2, function (x) any(is.na(x))))
  submf <- new_manifest[epic2idxs,]
  coords <- unique(paste(submf$CHR, submf$MAPINFO, sep=":"))
  wgbs <- as.numeric(as.matrix(values(CpGs.hg38)[coords,]))
  wgbsna <- which(is.na(wgbs))
  
  remove <- sort(union(wgbsna, union(epic1na , epic2na)))
  if(length(remove) > 0){
    epic1 <- epic1[-remove]
    wgbs <- wgbs[-remove]
    epic2 <- epic2[,-remove]
  }
  
  if(length(epic1) < 12){
    return(NA)
  }
  
  ###Replicate
  epic1 <- matrix(rep(epic1, nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
  rownames(epic1) <- rownames(epic2)
  colnames(epic1) <- colnames(epic2)
  wgbs <- matrix(rep(wgbs, nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
  rownames(wgbs) <- rownames(epic2)
  colnames(wgbs) <- colnames(epic2)
  
  pfnames=c("EPICv1", "EPICv2", "WGBS")
  EPICmm <- MultiMeasure(names=pfnames, data=list(epic1, epic2, wgbs))
  fit <- fitConsensus(EPICmm)
  fit
}

posfits <- mclapply(unique(posmfcons$namegroups), runconsensusmean, mc.cores = 16)
names(posfits) <- unique(posmfcons$namegroups)
nas <- unlist(lapply(posfits, is.na))
nas <- names(nas)[nas]
posmf$posvector[posmf$namegroups %in% nas] <- "Insufficient evidence"
posfits <- posfits[!names(posfits) %in% nas]

posdfres <- data.frame(Probepos=unlist(lapply(posfits, function (x) rownames(x@a_i))),
                       v1sens=unlist(lapply(posfits, function (x) x@b_i[,1])),
                       v2sens=unlist(lapply(posfits, function (x) x@b_i[,2])),
                       v1prec=unlist(lapply(posfits, function (x) x@d_i[,1])),
                       v2prec=unlist(lapply(posfits, function (x) x@d_i[,2])))
posdfres$posgroup <- substr(rownames(posdfres), 1, nchar(rownames(posdfres))-1)
posdfres$posgroup[grep("cg06373096_TC11", posdfres$posgroup)] <- "cg06373096_TC11;cg06373096_TC110;cg06373096_TC12;cg06373096_TC13;cg06373096_TC14;cg06373096_TC15;cg06373096_TC16;cg06373096_TC17;cg06373096_TC18;cg06373096_TC19"
stopifnot(all((posdfres$posgroup %in% posmf$namegroups)))
apply(posdfres[,2:5], 2, function (x) any(is.na(x)))
#v1sens v2sens v1prec v2prec 
#FALSE  FALSE  FALSE  FALSE 

posdfres$posgroup <- paste(posdfres$posgroup, "mean", sep=";")
posmf$posgroups <- paste(posmf$posrep_IlmnIDs, "mean", sep=";")

getrank <- function(df, probegroup){
  sub <- df[df$posgroup==probegroup,]
  rownames(sub) <- sub$Probepos
  sub <- sub[sort(unlist(strsplit(sub$posgroup[1], ";"))),]
  sub$sens <- rank(-sub$v2sens)
  sub$prec <- rank(sub$v2prec)
  sub
}


posranks <- lapply(unique(posdfres$posgroup), function (x) getrank (posdfres, x))
finalpos <- rbind.fill(posranks)
meanposresults <- finalpos[finalpos$Probepos=="mean",]
bestmeanpos <- meanposresults[meanposresults$prec==1,]
nrow(bestmeanpos)
#1369
meanposprobes <- unlist(sapply(bestmeanpos$posgroup, function (x) strsplit(x, ";")))
meanposprobes <- meanposprobes[meanposprobes!="mean"]
length(meanposprobes)
#3078

new_manifest$Best_precision_as_group_mean_LOCATION <- new_manifest$IlmnID %in% meanposprobes
new_manifest$Best_precision_as_group_mean_LOCATION[!new_manifest$IlmnID %in% posmfcons$IlmnID] <- NA

table(new_manifest$Best_precision_as_group_mean_LOCATION)/nrow(posmfcons)
#FALSE      TRUE 
#0.6361766 0.3638234

#Location

poswgbs <- posmf[grep("WGBS", posmf$Rep_results_by_LOCATION),]

getmeanrmse <- function(posgroup){
  sub <- poswgbs[poswgbs$posgroups==posgroup, c("IlmnID", "Name", "RMSE_with_WGBS", "posgroups", "posvector")]
  epic <- colMeans(as.matrix(EPICv2_Ms[sub$IlmnID,]))
  coords <- unique(sapply(sub$IlmnID, function (x) paste(new_manifest[x, c("CHR", "MAPINFO")], collapse=":")))
  wgbs <-  as.matrix(values(CpGs.hg38[coords,]))
  rmses <- apply(wgbs, 1, function (x) mean((x-epic)^2, na.rm=T))
  if(length(rmses) > 1){return(rmses)
  } else {
    sub <- rbind(sub, c(IlmnID="mean", pos=sub$Name[1], RMSE_with_WGBS=rmses, posgroups=sub$posgroups[1], posvector=sub$posvector[1]))
    return(sub) 
  }
}

wgbsposres <- mclapply(unique(poswgbs$posgroups), getmeanrmse, mc.cores = 10)
wgbsposres <- rbind.fill(wgbsposres)

getrank <- function(df, probegroup){
  sub <- df[df$posgroup==probegroup,]
  sub$rmserank <- rank(sub$RMSE_with_WGBS)
  sub
}

rankres <- lapply(unique(wgbsposres$posgroups), function (x) getrank(wgbsposres, x))
wgbsposres <- rbind.fill(rankres)
wgbsmeanres <- wgbsposres[wgbsposres$IlmnID=="mean",]
wgbsmeanbest <- wgbsmeanres[wgbsmeanres$rmserank==1,]
wgbsmeanprobes <- unlist(strsplit(wgbsmeanbest$posgroups, ";"))
wgbsmeanprobes <- wgbsmeanprobes[wgbsmeanprobes!="mean"]
new_manifest$Rep_results_by_LOCATION[new_manifest$IlmnID %in% wgbsmeanprobes] <- "Superior group mean (by WGBS)"

#####################################

new_manifest$Rep_results_by_NAME[new_manifest$Rep_results_by_NAME=="Superior probe" & new_manifest$Best_precision_as_group_mean_NAME] <- "Best sensitivity"
new_manifest$Rep_results_by_NAME[new_manifest$Rep_results_by_NAME=="Inferior probe" & new_manifest$Best_precision_as_group_mean_NAME] <- "Best precision by group mean"
new_manifest$Rep_results_by_NAME[new_manifest$Rep_results_by_NAME=="Best precision" & new_manifest$Best_precision_as_group_mean_NAME] <- "Best precision by group mean"

new_manifest$Rep_results_by_SEQUENCE[new_manifest$Rep_results_by_SEQUENCE=="Superior probe" & new_manifest$Best_precision_as_group_mean_SEQUENCE] <- "Best sensitivity"
new_manifest$Rep_results_by_SEQUENCE[new_manifest$Rep_results_by_SEQUENCE=="Inferior probe" & new_manifest$Best_precision_as_group_mean_SEQUENCE] <- "Best precision by group mean"
new_manifest$Rep_results_by_SEQUENCE[new_manifest$Rep_results_by_SEQUENCE=="Best precision" & new_manifest$Best_precision_as_group_mean_SEQUENCE] <- "Best precision by group mean"

new_manifest$Rep_results_by_LOCATION[new_manifest$Rep_results_by_LOCATION=="Superior probe" & new_manifest$Best_precision_as_group_mean_LOCATION] <- "Best sensitivity"
new_manifest$Rep_results_by_LOCATION[new_manifest$Rep_results_by_LOCATION=="Inferior probe" & new_manifest$Best_precision_as_group_mean_LOCATION] <- "Best precision by group mean"
new_manifest$Rep_results_by_LOCATION[new_manifest$Rep_results_by_LOCATION=="Best precision" & new_manifest$Best_precision_as_group_mean_LOCATION] <- "Best precision by group mean"

#Figure S14a
nametable <- table(new_manifest$Rep_results_by_NAME)[-1]
nametable <- nametable[c("Superior probe", "Inferior probe", "Best sensitivity", "Best precision", "Best precision by group mean", "Superior by WGBS", "Inferior by WGBS", "Superior group mean (by WGBS)", "Insufficient evidence")]
names(nametable) <- sub(" ", "\n", names(nametable))
names(nametable)[5] <- "Best precision\nby group mean"
names(nametable)[8] <- "Group mean\n(by WGBS)"
barplot(nametable, xlab="Classification", ylab="Number of probes", col="white")


#Figure S14b
nametable <- table(new_manifest$Rep_results_by_SEQUENCE)[-1]
nametable <- nametable[c("Superior probe", "Inferior probe", "Best sensitivity", "Best precision", "Best precision by group mean", "Superior by WGBS", "Inferior by WGBS", "Superior group mean (by WGBS)", "Insufficient evidence")]
names(nametable) <- sub(" ", "\n", names(nametable))
names(nametable)[5] <- "Best precision\nby group mean"
names(nametable)[8] <- "Group mean\n(by WGBS)"
barplot(nametable, xlab="Classification", ylab="Number of probes", col="white")

#Figure 7a
nametable <- table(new_manifest$Rep_results_by_LOCATION)[-1]
nametable <- nametable[c("Superior probe", "Inferior probe", "Best sensitivity", "Best precision", "Best precision by group mean", "Superior by WGBS", "Inferior by WGBS", "Superior group mean (by WGBS)", "Insufficient evidence")]
names(nametable) <- sub(" ", "\n", names(nametable))
names(nametable)[5] <- "Best precision\nby group mean"
names(nametable)[8] <- "Group mean\n(by WGBS)"
barplot(nametable, xlab="Classification", ylab="Number of probes", col="white")

new_manifest <- new_manifest[,1:80]
save(new_manifest, file="../new_manifest.RData")
write.csv(new_manifest, file="~/hippo/Gagri_Immunogenomics/Tim_Peters/EPICv2/AdditionalFile4.csv", row.names = F)


#Find if there are trends
dups_manifest <- new_manifest[new_manifest$Rep_results_by_LOCATION!="",]

#By Infinium type
mixed <- sapply(unique(dups_manifest$Name), function (x) length(unique(dups_manifest$Infinium_Design[dups_manifest$Name==x])))
probetype <- names(mixed)[mixed==2]
sub_manifest <- dups_manifest[dups_manifest$Name %in% probetype,]
res <- table(sub_manifest$Rep_results_by_LOCATION, sub_manifest$Infinium_Design_Type)
res <- res[c(9, 5, 7, 4, 8, 3, 1, 2, 6),]

#Figure 7b
par(mar=c(6.1, 4.1, 4.1, 2.1))
x <- barplot(t(res), beside=T, col=c("cadetblue1", "dodgerblue1"), xaxt="n", ylab="# of probes",
             main="Replicate probe classification,\nheterogeneous by Infinium Design type")
labs <- c("Superior\nprobe", "Inferior\nprobe", "Superior\nby WGBS", "Inferior\nby WGBS", "Superior\ngroup mean\n(by WGBS)",
          "Best\nsensitivity", "Best\nprecision", "Best precision\nby group mean", "Insufficient\nevidence")
text(cex=1, x=x[1,]+0.5, y=-4, labs, xpd=TRUE, srt=45)
legend("topleft", c("Type I", "Type II"), fill=c("cadetblue1", "dodgerblue1"), bty="n")

