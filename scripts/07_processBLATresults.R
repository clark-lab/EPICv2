library(GenomicRanges)
library(plyr)
library(Biostrings)
library(stringi)
library(rtracklayer)
require(parallel)


setwd("~/scratch/EPIC_V2/filtered/")
allmatches <- NULL

for (i in list.files(pattern = ".*psl$")){
  if (file.size(i) > 0){
    thisblat <- read.table(i)
    thisblat$sourcefile <- i
    scorefile <- read.table(gsub(".psl", ".score", i))
    thisblat$score <- scorefile$V5
    allmatches <- rbind(allmatches, thisblat)
  }
  print(i)
}

colnames(allmatches) <- c('match', 'mismatch', 'rep.match', 'Ns', 'Qgapcount', 'Qgapbases', 'Tgapcount', 'T gap bases',	'strand',	'Qname', 'Qsize',
                          'Qstart', 'Qend', 'chr', 'chrsize', 'start', 'end', 'blockcount', 'blocksizes', 'Qstarts',  'Tstarts',  'sourcefile', 'score')

#We filtered out gaps before, but not this time since we now have Blat score
allmatches$probe <- substr(allmatches$Qname, 1, nchar(allmatches$Qname) - 2)
allmatches$chr <- paste0("chr", allmatches$chr)

forwardmatches <- allmatches[grep("forward", allmatches$sourcefile),]
reversematches <- allmatches[grep("reverse", allmatches$sourcefile),]

getuniqueforward <- function(probe){
  reference <- forwardmatches
  matches <- reference[reference$probe %in% probe,]
  tmpranges <- GRanges(matches$chr, IRanges(matches$start, matches$end-1))
  tmpranges$genome <- matches$sourcefile
  tmpranges$score <- matches$score
  ol <- as.data.frame(findOverlaps(tmpranges))
  ol$startdiff <- abs(start(tmpranges[ol$queryHits]) - start(tmpranges[ol$subjectHits]))
  ol <- ol[ol$startdiff <= 3,]
  groups <- unique(tapply(ol$subjectHits, ol$queryHits, function(x) {unique(paste(sort(x), collapse='-'))}))
  groups <- lapply(strsplit(groups, split = '-'), as.numeric)
  final <- NULL
  for (i in 1:length(groups)){
    if(length(unlist(groups[i])) == 1){
      final <- c(final, tmpranges[unlist(groups[i])])
    } else {
      thisgroup <- tmpranges[unlist(groups[i])]
      thisgroup <- thisgroup[thisgroup$score == max(thisgroup$score)]
      if (length(thisgroup) > 1){
        thisgroup <- thisgroup[1]
        thisgroup$genome <- "Forward_meth_and_unmeth"
      }
      final <- c(final, thisgroup)
    }
  }
  print(probe)
  unlist(GRangesList(final))
}



getuniquereverse <- function(probe){
  reference <- reversematches
  matches <- reference[reference$probe %in% probe,]
  tmpranges <- GRanges(matches$chr, IRanges(matches$start, matches$end-1))
  tmpranges$genome <- matches$sourcefile
  tmpranges$score <- matches$score
  ol <- as.data.frame(findOverlaps(tmpranges))
  ol$startdiff <- abs(start(tmpranges[ol$queryHits]) - start(tmpranges[ol$subjectHits]))
  ol <- ol[ol$startdiff <= 3,]
  groups <- unique(tapply(ol$subjectHits, ol$queryHits, function(x) {unique(paste(sort(x), collapse='-'))}))
  groups <- lapply(strsplit(groups, split = '-'), as.numeric)
  final <- NULL
  for (i in 1:length(groups)){
    if(length(unlist(groups[i])) == 1){
      final <- c(final, tmpranges[unlist(groups[i])])
    } else {
      thisgroup <- tmpranges[unlist(groups[i])]
      thisgroup <- thisgroup[thisgroup$score == max(thisgroup$score)]
      if (length(thisgroup) > 1){
        thisgroup <- thisgroup[1]
        thisgroup$genome <- "Reverse_meth_and_unmeth"
      }
      final <- c(final, thisgroup)
    }
  }
  print(probe)
  unlist(GRangesList(final))
}


forwardresults <- mclapply(unique(as.character(forwardmatches$probe)), getuniqueforward, mc.cores = 16)
reverseresults <- mclapply(unique(as.character(reversematches$probe)), getuniquereverse, mc.cores = 16)

names(forwardresults) <- unique(forwardmatches$probe)
names(reverseresults) <- unique(reversematches$probe)

forwardresults <- mclapply(1:length(forwardresults), function (x){print(x);
                                                                  tmp <- forwardresults[x];
                                                                  tmp$probe <- names(forwardresults)[x];
                                                                  tmp}, mc.cores = 50)
reverseresults <- mclapply(1:length(reverseresults), function (x){print(x);
                                                                  tmp <- reverseresults[x];
                                                                  tmp$probe <- names(reverseresults)[x];
                                                                  tmp}, mc.cores = 50)

forwardresults <- unlist(GRangesList(forwardresults))
reverseresults <- unlist(GRangesList(reverseresults))

# Change reverse to forward strand

chr.sizes <- read.table("~/scratch/EPIC_V2/hg38.chrom.sizes", sep=' ')
chr.sizes$V4 <- as.numeric(gsub("REF\t", "", chr.sizes$V4))
chr.sizes$V1 <- paste0("chr", chr.sizes$V1)

reversebychr <- function(chr){
  chrranges <- reverseresults[seqnames(reverseresults) %in% chr]
  chrsize <- chr.sizes$V4[chr.sizes$V1 %in% chr]
  newchrranges <- GRanges(chr, IRanges(chrsize - end(chrranges), chrsize - start(chrranges)))
  mcols(newchrranges) <- mcols(chrranges)
  newchrranges
}

reverse_made_forward <- sapply(unique(seqnames(reverseresults)), reversebychr)
reverse_made_forward <- unlist(GRangesList(reverse_made_forward))

allresults <- c(forwardresults, reverse_made_forward)
#Remove trailing underscore
trailund <- grep("_$|_1$", allresults$probe)
allresults$probe[trailund] <- substr(allresults$probe[trailund], 1, 15)
#Change chrMT to chrM
seqlevels(allresults)[seqlevels(allresults)=="chrMT"] <- "chrM"
seqnames(allresults)[seqnames(allresults)=="chrM"] <- "chrM"

manifest <- read.csv("EPIC-8v2-0_EA.csv", skip = 7)

#Remove controls
manifest <- manifest[manifest$Probe_Type %in% c("cg", "ch", "nv", "rs"),]

manifestGR <- GRanges(paste0(manifest$CHR, ":", manifest$MAPINFO))
names(manifestGR) <- manifest$IlmnID
manifestGR$type <- paste0(manifest$Probe_Type, manifest$Infinium_Design_Type, manifest$Strand_FR, manifest$Strand_TB, manifest$Strand_CO)


annotatetarget <- function(probe){
  print(probe)
  manifestinfo <- manifestGR[probe]
  if (as.character(seqnames(manifestinfo))=="chr0"){
    return(data.frame(probe=probe, from=NA, move=NA, targetstrand=NA, type="No given target coordinate"))
  }
  matches <- allresults[allresults$probe %in% probe,]
  strand(matches) <- ifelse(grepl("orward", matches$genome), "+", "-")
  if (length(matches)== 0){
    return(data.frame(probe=probe, from=NA, move=NA, targetstrand=NA, type="No homology to hg38"))
  }
  #Find pattern
  ontarget <- matches[as.character(seqnames(matches))==as.character(seqnames(manifestinfo))]
  if (length(ontarget) == 0){
    return(data.frame(probe=probe, from=NA, move=NA, targetstrand=NA, type="No on target matches"))
  }
  if (length(ontarget) > 1){
    dists <- t(sapply(1:length(ontarget), function (x) start(manifestinfo) - c(start(ontarget[x]), end(ontarget[x]))))
    otidx <- apply(dists, 1, function (x) any(abs(x) <= 1)) 
    ontarget <- ontarget[otidx]
  }
  if (length(ontarget) == 0){
    return(data.frame(probe=probe, from=NA, move=NA, targetstrand=NA, type="Cis-target"))
  }
  ontarget <- ontarget[!duplicated(ontarget)]
  dists <- start(manifestinfo) - c(start(ontarget), end(ontarget))
  from <- ifelse(abs(dists[1] <= 1), "start", "end")
  move <- ifelse(from=="start", dists[1], dists[2])
  targetstrand <- strand(ontarget)
  data.frame(probe=probe, from=from, move=move, targetstrand=targetstrand, type=manifestinfo$type)
}

targetres <- mclapply(names(manifestGR), annotatetarget, mc.cores=16)
targetres <- rbind.fill(targetres)

#3 probes on cusp of repeat region - correct for ambiguity
targetres$from[targetres$probe=="cg15884558_BC21"] <- "end"
targetres$move[targetres$probe=="cg15884558_BC21"] <- 0
targetres$from[targetres$probe=="cg16365276_BC21"] <- "start"
targetres$move[targetres$probe=="cg16365276_BC21"] <- 0
targetres$from[targetres$probe=="cg20934014_TC21"] <- "start"
targetres$move[targetres$probe=="cg20934014_TC21"] <- 0
targetres <- targetres[!duplicated(targetres$probe),]
rownames(targetres) <- targetres$probe

nohomology <- targetres$probe[targetres$type=="No homology to hg38"]
write.table(nohomology, sep=",", file="No_homology_probes.csv", row.names = F, col.names = F)

notarget <- targetres$probe[targetres$type=="No given target coordinate"]
write.table(notarget, sep=",", file="No_given_target_probes.csv", row.names = F, col.names = F)

#Correct for chrM/chrMT
targetresM <- lapply(names(manifestGR)[as.logical(seqnames(manifestGR)=="chrM")], annotatetarget)
targetresM <- rbind.fill(targetresM)
targetresM$targetstrand <- as.character(targetresM$targetstrand)
rownames(targetresM) <- targetresM$probe
targetres[rownames(targetresM),] <- targetresM

#Off target alignments with no ontarget
offnoon <- targetres$probe[targetres$type %in% c("No on target matches", "Cis-target")]
offnoon <- c(offnoon, na.omit(targetres$probe[abs(targetres$move) > 1]))

assembleoffnoon <- function(probe){
  print(probe)
  manifestinfo <- manifestGR[probe]
  matches <- allresults[allresults$probe %in% probe,]
  matches$targetCHR <- seqnames(manifestinfo)
  matches$targetPOS <- start(manifestinfo)
  matches$genome[grep("forward_meth.*psl", matches$genome)] <- "Forward_methylated"
  matches$genome[grep("forward_unmeth.*psl", matches$genome)] <- "Forward_unmethylated"
  matches$genome[grep("reverse_meth.*psl", matches$genome)] <- "Reverse_methylated"
  matches$genome[grep("reverse_unmeth.*psl", matches$genome)] <- "Reverse_unmethylated"
  colnames(values(matches))[2] <- "BLAT_score"
  values(matches) <- values(matches)[,c("probe", "targetCHR", "targetPOS", "BLAT_score", "genome")]
  names(matches) <- paste0("H", 1:length(matches))
  as.data.frame(matches)
}

offnoonres <- lapply(offnoon, assembleoffnoon)
offnoonres <- rbind.fill(offnoonres)
write.table(offnoonres, sep=",", file="BLAT_hits_no_homology_to_target.csv", row.names = F)


#Extend 2 base pairs each direction
allresultsplustwo <- allresults
seqlevels(allresultsplustwo) <- gsub("chr", "", seqlevels(allresultsplustwo))
seqnames(allresultsplustwo) <- gsub("chr", "", seqnames(allresultsplustwo))
seqlevels(allresultsplustwo)[seqlevels(allresultsplustwo)=="M"] <- "MT"
seqnames(allresultsplustwo)[seqnames(allresultsplustwo)=="M"] <- "MT"
start(allresultsplustwo) <- start(allresultsplustwo) - 2
end(allresultsplustwo) <- end(allresultsplustwo) + 2
start(allresultsplustwo[start(allresultsplustwo) < 0]) <- 1
export.bed(allresultsplustwo, "allresults.bed")
seqs <- system(paste0("bedtools getfasta -fi ../hg38.fa -bed allresults.bed"), intern = T)
allresults$seqs <- seqs[1:length(seqs) %%2 == 0]
allresults$seqs <- DNAStringSet(allresults$seqs)


remove <- c(nohomology, notarget, offnoon)
targetres <- targetres[!targetres$probe %in% remove,]
manifestGR <- manifestGR[!names(manifestGR) %in% remove]

all(names(manifestGR) %in% targetres$probe)
#TRUE
all(targetres$probe %in% names(manifestGR))
#TRUE

table(targetres$type, paste0(targetres$from, targetres$move))

###end-1   end0   end1 start-1 start0 start1
#cgIFBC   29230      0      0       0      0     11
#cgIFBO       0      0    134       0      0      0
#cgIFTC   33978      0      0       0      0     11
#cgIFTO       0      0    131       0      0      0
#cgIIFBC      0 224200      0       0   2412      0
#cgIIFBO      0      1     27       0      0      0
#cgIIFTC      0 170143      0       0   2992      0
#cgIIFTO      0      0    140       0      0      0
#cgIIRBC      0   2402      0       0 224216      1
#cgIIRTC      0   2852      0       0 170235      1
#cgIRBC      15      0      0       0      0  29283
#cgIRBO       0      0      0     150      0      0
#cgIRTC      11      0      0       0      0  33648
#cgIRTO       0      0      0     106      0      0
#chIIFBC      0      0   1145       0      0      0
#chIIFBO      0      0     13       0      0      0
#chIIFTC      0      0    308       0      0      0
#chIIFTO      0      0      3       0      0      0
#chIIRBC      0      0      0       0   1389      0
#chIIRBO      0      0      0       0     12      0
#chIIRTC      0      0      0       0     34      0
#chIIRTO      0      0      0       0      1      0
#nvIFBC       0    244      0       0      0      0
#nvIFBO       0      2     10       0      0      0
#nvIFTC       0     95      0       0      0      0
#nvIFTO       0      1     26       0      0      0
#nvIIFBC      0      0     26       0      0      0
#nvIIFBO      0      0      1       0      0      0
#nvIIFTC      0      0      5       0      0      0
#nvIIRBC      0      0      0       0     24      0
#nvIIRTC      0      0      0       0      3      0
#nvIIRTO      0      0      0       2      0      0
#nvIRBC       0      0      0       0      0    230
#nvIRBO       0      0      0       0      6      0
#nvIRTC       0      0      0       0      0    109
#nvIRTO       0      0      0       0     20      0
#rsIFBC       0     12      0       0      0      0
#rsIFBO       0      1      0       0      0      0
#rsIFTC       0      1      0       0      0      0
#rsIIFBC      0      0     17       0      0      0
#rsIIFBO      0      0      1       0      0      0
#rsIIRBC      0      0      0       0     21      0
#rsIIRTC      0      0      0       0      2      0
#rsIRBC       0      0      0       0      0     10


findofftargets <- function(probe){
  print(probe)
  targets <- allresults[allresults$probe %in% probe]
  targets <- targets[!duplicated(targets)]
  targets <- targets[!duplicated(start(targets))]
  targets <- targets[!duplicated(end(targets))]
  mf <- manifestGR[probe]
  tr <- targetres[probe,]
  fromtr <- paste0(tr$from, tr$move)
  rel <- mapvalues(fromtr, c("start-1", "start0", "start1", "end-1", "end0", "end1"),
                     c("startminus", "start", "startplus", "endminus", "end", "endplus"), warn_missing = F)
  nttargets <- GRanges(seqnames(targets), IRanges(switch(rel, "startminus" = start(targets) - 1,
                                        "start" = start(targets),
                                        "startplus" = start(targets) + 1,
                                        "endminus" = end(targets) - 1,
                                        "end" = end(targets),
                                        "endplus" = end(targets) + 1)))
  seqlevels(mf) <- seqlevels(nttargets)
  targets$target <- nttargets==mf
  probeseq <- targets$seqs[targets$target]
  posn <- switch(rel, "startminus" = 2, "start" = 3, "startplus" = 4, "endminus" = 51, "end" = 52, "endplus" = 53)
  
  if(length(targets)==1){
    targets$target=TRUE
    targets$base=substr(targets$seqs, posn, posn)
    targets$basecoord <- start(mf)
    values(targets) <- values(targets)[,-4]
    return(targets)
  }
  
  if(grepl("orward", targets$genome[targets$target])){
    print("forward probe")
    fwds <- targets[grepl("orward", targets$genome)]
    if(length(fwds) > 0){
      alns <- sapply(fwds$seqs, function(x) pairwiseAlignment(x, probeseq))
      baseposns <- sapply(alns, function (x) nchar(gsub("-", "", substr(as.character(x), 1, posn))))
      fwds$base <- sapply(1:length(fwds), function (x) substr(fwds$seqs[x], baseposns[x], baseposns[x]))
      fwds$basecoord <- start(fwds) + baseposns - 3
    } else {
      fwds <- GRanges()
    }
    revs <- targets[grepl("everse", targets$genome)]
    if(length(revs) > 0){
      alns <- sapply(reverseComplement(revs$seqs), function(x) pairwiseAlignment(x, probeseq))
      baseposns <- sapply(alns, function (x) nchar(gsub("-", "", substr(as.character(x), 1, posn)))) + 1
      revs$base <- sapply(1:length(revs), function (x) substr(stri_reverse(revs$seqs[x]), baseposns[x], baseposns[x]))
      revs$basecoord <- end(revs) - baseposns + 3
    } else {
      revs <- GRanges()
    }
    
  } else {
    print("reverse probe")
    fwds <- targets[grepl("orward", targets$genome)]
    if(length(fwds) > 0){
      alns <- sapply(reverseComplement(fwds$seqs), function(x) pairwiseAlignment(x, probeseq))
      baseposns <- sapply(alns, function (x) nchar(gsub("-", "", substr(as.character(x), 1, posn)))) + 1
      fwds$base <- sapply(1:length(fwds), function (x) substr(stri_reverse(fwds$seqs[x]), baseposns[x], baseposns[x]))
      fwds$basecoord <- end(fwds) - baseposns + 3
      } else {
      fwds <- GRanges()
    }
    revs <- targets[grepl("everse", targets$genome)]
    if(length(revs) > 0){
      alns <- sapply(revs$seqs, function(x) pairwiseAlignment(x, probeseq))
      baseposns <- sapply(alns, function (x) nchar(gsub("-", "", substr(as.character(x), 1, posn))))
      revs$base <- sapply(1:length(revs), function (x) substr(revs$seqs[x], baseposns[x], baseposns[x]))
      revs$basecoord <- start(revs) + baseposns - 3
     } else {
      revs <- GRanges()
    }
  }
  targetsfinal <- c(fwds, revs)
  values(targetsfinal) <- values(targetsfinal)[,-4]
  targetsfinal
}

offtargets <- mclapply(names(manifestGR), findofftargets, mc.cores = 16)

#Check that the probe coordinate and manifest coordinate line up

probes <- unlist(mclapply(offtargets, function(x) x$probe[1], mc.cores = 10))
all(probes==names(manifestGR))
#[1] TRUE
offtargets <- unlist(GRangesList(offtargets))
targets <- offtargets[offtargets$target]
all(targets$probe==names(manifestGR))
#TRUE
all(targets$basecoord==start(manifestGR))
#FALSE
idxs <- which(targets$basecoord!=start(manifestGR))
offbyone <- data.frame(targets[idxs])
offbyoneidx <- which(offtargets$probe %in% offbyone$probe)
offbyone <- offtargets[offbyoneidx]
oboprobes <- unique(offbyone$probe)
offbyone$basecoord[offbyone$probe=="cg11115555_BC21"] <- end(offbyone[offbyone$probe=="cg11115555_BC21"])
offbyone$basecoord[offbyone$probe!="cg11115555_BC21"] <- end(offbyone[offbyone$probe!="cg11115555_BC21"]) + 1
oboGR <- GRanges(gsub("chr", "", seqnames(offbyone)), IRanges(offbyone$basecoord, offbyone$basecoord))
export.bed(oboGR, "offbyone.bed")
seqs <- system(paste0("bedtools getfasta -fi ../hg38.fa -bed offbyone.bed"), intern = T)
offbyone$base <- seqs[1:length(seqs) %%2 == 0]

offtargets[offbyoneidx] <- offbyone
targets <- offtargets[offtargets$target]
all(targets$probe==names(manifestGR))
#TRUE
all(targets$basecoord==start(manifestGR))
#TRUE

#Now check that the base lines up with the base coordinate
alltargetsGR <- GRanges(gsub("chr", "", seqnames(offtargets)), IRanges(offtargets$basecoord, offtargets$basecoord))
seqlevels(alltargetsGR)[seqlevels(alltargetsGR)=="M"] <- "MT"
seqnames(alltargetsGR)[seqnames(alltargetsGR)=="M"] <- "MT"
export.bed(alltargetsGR, "alltargets.bed")
seqs <- system(paste0("bedtools getfasta -fi ../hg38.fa -bed alltargets.bed"), intern = T)
bases <- seqs[1:length(seqs) %%2 == 0]
all(offtargets$base==bases)
#FALSE
nmidx <- offtargets$base!=bases
nonmatches <- offtargets[nmidx]
nmGR <- GRanges(gsub("chr", "", seqnames(nonmatches)), IRanges(nonmatches$basecoord, nonmatches$basecoord))
export.bed(nmGR, "nonmatches.bed")
seqs <- system(paste0("bedtools getfasta -fi ../hg38.fa -bed nonmatches.bed"), intern = T)
nonmatches$base <- seqs[1:length(seqs) %%2 == 0]
offtargets[nmidx] <- nonmatches
all(offtargets$base==bases)
#TRUE


offtargets <- offtargets[!offtargets$target]
offtargetsGR <- GRanges(seqnames(offtargets), IRanges(offtargets$basecoord, offtargets$basecoord))
offtargetsGR$probe <- offtargets$probe
offtargetsGR$genome <- offtargets$genome
offtargetsGR$BLAT_score <- offtargets$score
offtargetsGR$base <- offtargets$base

offtargetsGR <- data.frame(offtargetsGR)
offtargetsGR <- offtargetsGR[,c("probe", "seqnames", "start", "genome", "BLAT_score", "base")]
colnames(offtargetsGR)[1:3] <- c("ProbeID", "CHR", "MAPINFO")
offtargetsGR$genome[grep("_forward_meth", offtargetsGR$genome)] <- "Forward_methylated"
offtargetsGR$genome[grep("_forward_unmeth", offtargetsGR$genome)] <- "Forward_unmethylated"
offtargetsGR$genome[grep("_reverse_meth", offtargetsGR$genome)] <- "Reverse_methylated"
offtargetsGR$genome[grep("_reverse_unmeth", offtargetsGR$genome)] <- "Reverse_unmethylated"
otgr <- GRanges(paste(offtargets$CHR, offtargets$MAPINFO, sep=":"))

#Load in hg38 CpG annotation
load("hg38.WGBS.Rdata")
cs <- CpG[strand(CpG)=="+"]
gs <- CpG[strand(CpG)=="-"]
offtargets$is_cpg <- otgr %over% cs
write.csv(offtargets, file="Offtargets_FINAL.csv", row.names = F)

