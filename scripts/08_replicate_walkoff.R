library(consensus)
library(ggplot2)
library(GenomicRanges)
library(effsize)
library(Metrics)
library(plyr)

###Walk-off between namereps

load("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/Mvalues_NODROPOUT.RData")
load("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/new_manifest_with_CH.RData")   
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
namemfcons <- namemf[!is.na(namemf$vecEPICv1probeID),]

#Check if all duplicates are also duplicates for position

checknamesameid <- function(probe){
  epic2idxs <- new_manifest$IlmnID[new_manifest$Name==probe]
  submf <- new_manifest[epic2idxs,]
  coords <- unique(paste(submf$CHR, submf$MAPINFO, sep=":"))
  length(coords)
}

lengths <- mclapply(namemfcons$Name, checknamesameid, mc.cores=16)
stopifnot(all(unlist(lengths)==1))

checkepicv1id <- function(probe){
  epic2idxs <- new_manifest$IlmnID[new_manifest$Name==probe]
  submf <- new_manifest[epic2idxs,]
  ids <- unique(submf$vecEPICv1probeID)
  length(ids)
}

lengths <- mclapply(namemfcons$Name, checkepicv1id, mc.cores=16)
stopifnot(all(unlist(lengths)==1))

runconsensus <- function(probe){
  message(probe)
  epic1 <- EPICv1_Ms[probe,]
  epic1na <- which(is.na(epic1))
  epic2idxs <- new_manifest$IlmnID[new_manifest$Name==probe]
  epic2 <- EPICv2_Ms[epic2idxs,]
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

require(parallel)



namefits <- mclapply(unique(namemfcons$Name), runconsensus, mc.cores = 16)
names(namefits) <- unique(namemfcons$Name)
nas <- unlist(lapply(namefits, is.na))
nas <- names(nas)[nas]
namemf$namevector[namemf$Name %in% nas] <- "Insufficient evidence"
namefits <- namefits[!names(namefits) %in% nas]

namedfres <- data.frame(Probename=unlist(lapply(namefits, function (x) rownames(x@a_i))),
                        v1sens=unlist(lapply(namefits, function (x) x@b_i[,1])),
                        v2sens=unlist(lapply(namefits, function (x) x@b_i[,2])),
                        v1prec=unlist(lapply(namefits, function (x) x@d_i[,1])),
                        v2prec=unlist(lapply(namefits, function (x) x@d_i[,2])))
epicname <- substr(namedfres$Probename, 1, 10)
namedfres$namegroup <- sapply(epicname, function (x) paste(namedfres$Probename[gsub("_.*", "", namedfres$Probename)==x], collapse=";"))
stopifnot(all((namedfres$namegroup %in% namemf$namegroups)))
apply(namedfres[,2:5], 2, function (x) any(is.na(x)))
#v1sens v2sens v1prec v2prec 
#FALSE  FALSE  FALSE  FALSE 

getrank <- function(df, probegroup){
  sub <- df[df$namegroup==probegroup,]
  rownames(sub) <- sub$Probename
  sub <- sub[sort(unlist(strsplit(sub$namegroup[1], ";"))),]
  sub$sens <- rank(-sub$v2sens)
  sub$prec <- rank(sub$v2prec)
  sub
}

nameranks <- lapply(unique(namedfres$namegroup), function (x) getrank (namedfres, x))
finalname <- rbind.fill(nameranks)
superior <- unique(finalname$Probename[finalname$sens==1 & finalname$prec==1])
inferior <- unique(finalname$Probename[finalname$sens!=1 & finalname$prec!=1])
bettersens <- unique(finalname$Probename[finalname$sens==1 & finalname$prec!=1])
betterprec <- unique(finalname$Probename[finalname$sens!=1 & finalname$prec==1])

namemf[superior, "namevector"] <- "Superior probe"
namemf[inferior, "namevector"] <- "Inferior probe"
namemf[bettersens, "namevector"] <- "Best sensitivity"
namemf[betterprec, "namevector"] <- "Best precision"

#Then use just WGBS to do the walkoff
stopifnot(all(is.na(namemf$vecEPICv1probeID[namemf$namevector==""])))

#Proof that there are no more EPICv1 probes to compare

namegroups <- unique(namemf$namegroup[namemf$namevector==""])

getminwgbs <- function(namegroup){
  sub <- namemf[namemf$namegroups==namegroup, c("IlmnID", "Name", "RMSE_W_WGBS", "namegroups", "namevector")]
  rmses <- sub$RMSE_W_WGBS
  names(rmses) <- rownames(sub)
  if(length(na.omit(rmses)) < 2){
    sub$namevector <- "Insufficient evidence"
  } else {
    sub$namevector[sub$RMSE_W_WGBS==min(sub$RMSE_W_WGBS, na.rm=T)] <- "Superior by WGBS"
    sub$namevector[sub$RMSE_W_WGBS!=min(sub$RMSE_W_WGBS, na.rm=T)] <- "Inferior by WGBS"
  }
  sub$namevector[is.na(sub$RMSE_W_WGBS)] <- "Insufficient evidence"
  sub
}

wgbsres <- lapply(namegroups, getminwgbs)
wgbsres <- rbind.fill(wgbsres)
rownames(wgbsres) <- wgbsres$IlmnID
namemf[rownames(wgbsres),"namevector"] <- wgbsres$namevector

new_manifest$Dup_results_by_NAME <- ""
new_manifest[rownames(namemf),"Dup_results_by_NAME"] <- namemf$namevector

######################################################
#Sequence reps
seqmf$seqvector <- ""
seqgroups <- unique(seqmf$namegroups)
seqmfcons <- seqmf[!is.na(seqmf$vecEPICv1probeID),]

#Check if all duplicates are also duplicates for position

checkseqsameid <- function(seqgroup){
  epic2idxs <- unlist(strsplit(seqgroup, ";"))
  submf <- new_manifest[epic2idxs,]
  coords <- unique(paste(submf$CHR, submf$MAPINFO, sep=":"))
  length(coords)
}

lengths <- mclapply(unique(seqmf$namegroups), checkseqsameid)
stopifnot(all(unlist(lengths)==1))

checkepicv1id <- function(namegroup){
  epic2idxs <- seqmfcons$IlmnID[seqmfcons$namegroups==namegroup]
  submf <- new_manifest[epic2idxs,]
  ids <- unique(submf$vecEPICv1probeID)
  length(ids)
}

lengths <- mclapply(unique(seqmfcons$namegroups), checkepicv1id, mc.cores=16)
stopifnot(all(unlist(lengths)==1))

runconsensus <- function(namegroup){
  message(namegroup)
  epic1 <- EPICv1_Ms[unique(seqmfcons$Name[seqmfcons$namegroups==namegroup]),]
  epic1na <- which(is.na(epic1))
  epic2idxs <- unlist(strsplit(namegroup, ";"))
  epic2 <- EPICv2_Ms[epic2idxs,]
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

require(parallel)

seqfits <- mclapply(unique(seqmfcons$namegroups), runconsensus, mc.cores = 16)
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
seqdfres$seqgroup <- sapply(seqdfres$Probeseq, function (x) unique(unlist(seqmfcons$namegroups[grep(x, seqmfcons$namegroups)])))
stopifnot(all((seqdfres$seqgroup %in% seqmf$namegroups)))
apply(seqdfres[,2:5], 2, function (x) any(is.na(x)))
#v1sens v2sens v1prec v2prec 
#FALSE  FALSE  FALSE  FALSE 

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
superior <- unique(finalseq$Probeseq[finalseq$sens==1 & finalseq$prec==1])
inferior <- unique(finalseq$Probeseq[finalseq$sens!=1 & finalseq$prec!=1])
bettersens <- unique(finalseq$Probeseq[finalseq$sens==1 & finalseq$prec!=1])
betterprec <- unique(finalseq$Probeseq[finalseq$sens!=1 & finalseq$prec==1])

seqmf[superior, "seqvector"] <- "Superior probe"
seqmf[inferior, "seqvector"] <- "Inferior probe"
seqmf[bettersens, "seqvector"] <- "Best sensitivity"
seqmf[betterprec, "seqvector"] <- "Best precision"

#Then use just WGBS to do the walkoff
stopifnot(all(is.na(seqmf$vecEPICv1probeID[seqmf$seqvector==""])))

#Proof that there are no more EPICv1 probes to compare

seqgroups <- unique(seqmf$namegroups[seqmf$seqvector==""])

getminwgbs <- function(seqgroup){
  sub <- seqmf[seqmf$namegroups==seqgroup, c("IlmnID", "Name", "RMSE_W_WGBS", "namegroups", "seqvector")]
  rmses <- sub$RMSE_W_WGBS
  names(rmses) <- rownames(sub)
  if(length(na.omit(rmses)) < 2){
    sub$seqvector <- "Insufficient evidence"
  } else {
    sub$seqvector[sub$RMSE_W_WGBS==min(sub$RMSE_W_WGBS, na.rm=T)] <- "Superior by WGBS"
    sub$seqvector[sub$RMSE_W_WGBS!=min(sub$RMSE_W_WGBS, na.rm=T)] <- "Inferior by WGBS"
  }
  sub$seqvector[is.na(sub$RMSE_W_WGBS)] <- "Insufficient evidence"
  sub
}

wgbsres <- lapply(seqgroups, getminwgbs)
wgbsres <- rbind.fill(wgbsres)
rownames(wgbsres) <- wgbsres$IlmnID
seqmf[rownames(wgbsres),"seqvector"] <- wgbsres$seqvector

new_manifest$Dup_results_by_SEQUENCE <- ""
new_manifest[rownames(seqmf),"Dup_results_by_SEQUENCE"] <- seqmf$seqvector

##############################################

######################################################
#position reps
posmf$posvector <- ""
posgroups <- unique(posmf$namegroups)
posmfcons <- posmf[!is.na(posmf$vecEPICv1probeID),]

#Check if all duplicates are also duplicates for position

checkpossameid <- function(posgroup){
  epic2idxs <- unlist(strsplit(posgroup, ";"))
  submf <- new_manifest[epic2idxs,]
  coords <- unique(paste(submf$CHR, submf$MAPINFO, sep=":"))
  length(coords)
}

lengths <- mclapply(unique(posmf$namegroups), checkpossameid)
stopifnot(all(unlist(lengths)==1))

checkepicv1id <- function(namegroup){
  epic2idxs <- posmfcons$IlmnID[posmfcons$namegroups==namegroup]
  submf <- new_manifest[epic2idxs,]
  ids <- unique(submf$vecEPICv1probeID)
  length(ids)
}

lengths <- mclapply(unique(posmfcons$namegroups), checkepicv1id, mc.cores=16)
stopifnot(all(unlist(lengths)==1))

runconsensus <- function(namegroup){
  message(namegroup)
  epic1 <- EPICv1_Ms[unique(posmfcons$Name[posmfcons$namegroups==namegroup]),]
  epic1na <- which(is.na(epic1))
  epic2idxs <- unlist(strsplit(namegroup, ";"))
  epic2 <- EPICv2_Ms[epic2idxs,]
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

require(parallel)

posfits <- mclapply(unique(posmfcons$namegroups), runconsensus, mc.cores = 16)
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
posdfres$posgroup <- sapply(posdfres$Probepos, function (x) unique(unlist(posmfcons$namegroups[grep(x, posmfcons$namegroups)])))
stopifnot(all((posdfres$posgroup %in% posmf$namegroups)))
apply(posdfres[,2:5], 2, function (x) any(is.na(x)))
#v1sens v2sens v1prec v2prec 
#FALSE  FALSE  FALSE  FALSE 

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
superior <- unique(finalpos$Probepos[finalpos$sens==1 & finalpos$prec==1])
inferior <- unique(finalpos$Probepos[finalpos$sens!=1 & finalpos$prec!=1])
bettersens <- unique(finalpos$Probepos[finalpos$sens==1 & finalpos$prec!=1])
betterprec <- unique(finalpos$Probepos[finalpos$sens!=1 & finalpos$prec==1])

posmf[superior, "posvector"] <- "Superior probe"
posmf[inferior, "posvector"] <- "Inferior probe"
posmf[bettersens, "posvector"] <- "Best sensitivity"
posmf[betterprec, "posvector"] <- "Best precision"

#Then use just WGBS to do the walkoff
stopifnot(all(is.na(posmf$vecEPICv1probeID[posmf$posvector==""])))

#Proof that there are no more EPICv1 probes to compare

posgroups <- unique(posmf$namegroups[posmf$posvector==""])

getminwgbs <- function(posgroup){
  sub <- posmf[posmf$namegroups==posgroup, c("IlmnID", "Name", "RMSE_W_WGBS", "namegroups", "posvector")]
  rmses <- sub$RMSE_W_WGBS
  names(rmses) <- rownames(sub)
  if(length(na.omit(rmses)) < 2){
    sub$posvector <- "Insufficient evidence"
  } else {
    sub$posvector[sub$RMSE_W_WGBS==min(sub$RMSE_W_WGBS, na.rm=T)] <- "Superior by WGBS"
    sub$posvector[sub$RMSE_W_WGBS!=min(sub$RMSE_W_WGBS, na.rm=T)] <- "Inferior by WGBS"
  }
  sub$posvector[is.na(sub$RMSE_W_WGBS)] <- "Insufficient evidence"
  sub
}

wgbsres <- lapply(posgroups, getminwgbs)
wgbsres <- rbind.fill(wgbsres)
rownames(wgbsres) <- wgbsres$IlmnID
posmf[rownames(wgbsres),"posvector"] <- wgbsres$posvector

new_manifest$Dup_results_by_POSITION <- ""
new_manifest[rownames(posmf),"Dup_results_by_POSITION"] <- posmf$posvector
save(new_manifest, file="Manifest_FINAL.RData")

#Figure S14a
nametable <- table(new_manifest$Dup_results_by_NAME)[-1]
nametable <- nametable[c("Superior probe", "Inferior probe", "Best sensitivity", "Best precision", "Superior by WGBS", "Inferior by WGBS", "Insufficient evidence")]
names(nametable) <- sub(" ", "\n", names(nametable))
barplot(nametable, col=c("deepskyblue", "coral", "lightblue", "lightblue", "lightgreen", "green4", "grey"),
        xlab="Type", ylab="Count", main="Duplicate probe classification by name")

#Figure S14b
seqtable <- table(new_manifest$Dup_results_by_SEQUENCE)[-1]
seqtable <- seqtable[c("Superior probe", "Inferior probe", "Best sensitivity", "Best precision", "Superior by WGBS", "Inferior by WGBS", "Insufficient evidence")]
names(seqtable) <- sub(" ", "\n", names(seqtable))
barplot(seqtable, col=c("deepskyblue", "coral", "lightblue", "lightblue", "lightgreen", "green4", "grey"),
        xlab="Type", ylab="Count", main="Duplicate probe classification by probe sequence")

#Figure 7a
postable <- table(new_manifest$Dup_results_by_POSITION)[-1]
postable <- postable[c("Superior probe", "Inferior probe", "Best sensitivity", "Best precision", "Superior by WGBS", "Inferior by WGBS", "Insufficient evidence")]
names(postable) <- sub(" ", "\n", names(postable))
barplot(postable, col=c("deepskyblue", "coral", "lightblue", "lightblue", "lightgreen", "green4", "grey"),
        xlab="Type", ylab="Count", main="Duplicate probe classification by genomic position")

#Plot cg04853151
cg04853151_BC11 <- EPICv2_Ms["cg04853151_BC11",]
cg04853151_BC21 <- EPICv2_Ms["cg04853151_BC21",]
cg04853151 <- EPICv1_Ms["cg04853151",]
wgbs <- values(CpGs.hg38)["chr1:54053571",]

epic2 <- rbind(cg04853151_BC11, cg04853151_BC21)

epic1 <- matrix(rep(cg04853151, nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
rownames(epic1) <- rownames(epic2)
colnames(epic1) <- colnames(epic2)
wgbs <- matrix(rep(as.numeric(as.matrix(wgbs)), nrow(epic2)), nrow=nrow(epic2), byrow = TRUE)
rownames(wgbs) <- rownames(epic2)
colnames(wgbs) <- colnames(epic2)
pfnames=c("EPICv1", "EPICv2", "WGBS")
EPICmm <- MultiMeasure(names=pfnames, data=list(epic1, epic2, wgbs))
fit <- fitConsensus(EPICmm)

#Figure 7c
par(mfrow=c(1, 2))
plotOneFit(EPICmm, "cg04853151_BC11", pal=c("deepskyblue", "midnightblue", "forestgreen"), x="bottomright", ncol=1)
plotOneFit(EPICmm, "cg04853151_BC21", x="bottomright", ncol=1)


#Find if there are trends
dups_manifest <- new_manifest[new_manifest$Dup_results_by_POSITION!="",]

#By Infinium type
mixed <- sapply(unique(dups_manifest$Name), function (x) length(unique(dups_manifest$Infinium_Design[dups_manifest$Name==x])))
probetype <- names(mixed)[mixed==2]
sub_manifest <- dups_manifest[dups_manifest$Name %in% probetype,]
res <- table(sub_manifest$Dup_results_by_POSITION, sub_manifest$Infinium_Design_Type)
res <- res[c(7, 4, 6, 3, 1, 2, 5),]

#Figure 7b
par(mar=c(6.1, 4.1, 4.1, 2.1))
x <- barplot(t(res), beside=T, col=c("cadetblue1", "dodgerblue1"), xaxt="n", ylab="# of probes",
        main="Probe performance by Infinium Design type")
labs <- c("Superior\nprobe", "Inferior\nprobe", "Superior\nby WGBS", "Inferior\nby WGBS",
          "Best\nprecision", "Best\nsensitivity", "Insufficient\nevidence")
text(cex=1, x=x[1,]+0.5, y=-6, labs, xpd=TRUE, srt=45)
legend("top", c("Type I", "Type II"), fill=c("cadetblue1", "dodgerblue1"), bty="n")

