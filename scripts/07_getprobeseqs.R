library(stringr)
library(gtools)
library(Biostrings)
require(parallel)
setwd("~/scratch/EPIC_V2/")
manifest <- read.csv("~/hippo/Gagri_Cancer-Epigenetics/Projects/EPIC_V2/EPIC-8v2-0_EA.csv", skip = 7)

#Remove controls
manifest <- manifest[substr(manifest$IlmnID, 1, 2) %in% c("cg", "ch", "nv", "rs"),]

#Type I probe sequences

typeIAseqs <- DNAStringSet(manifest$AlleleA_ProbeSeq[manifest$Infinium_Design_Type=="I"])
names(typeIAseqs) <- paste0(manifest$IlmnID[manifest$Infinium_Design_Type=="I"], "_A")

final <- typeIAseqs
typeIBseqs <- DNAStringSet(manifest$AlleleB_ProbeSeq[manifest$Infinium_Design_Type=="I"])
names(typeIBseqs) <- paste0(manifest$IlmnID[manifest$Infinium_Design_Type=="I"], "_B")

final <- c(final, typeIBseqs)

#Type II probe sequences
#No Rs

manifest$numRs <- str_count(manifest$AlleleA_ProbeSeq, "R")
typeIInoRs <- DNAStringSet(manifest$AlleleA_ProbeSeq[manifest$Infinium_Design_Type=="II" & manifest$numRs == 0])
names(typeIInoRs) <- paste0(manifest$IlmnID[manifest$Infinium_Design_Type=="II" & manifest$numRs == 0], "_1")

final <- c(final, typeIInoRs)

#1 or more Rs

putativeAGseqs <- function(seq, cpgid){
  rseq <- BString(seq)
  rpos <- unlist(gregexpr("R", rseq))
  numRs <- length(rpos)
  replaceAGs <- permutations(2, numRs, v=c("A", "G"), repeats.allowed = T)
  at <- IRanges(rpos, rpos)
  AGseqs <- NULL
  for (i in 1:nrow(replaceAGs)){
    newseq <- replaceAt(rseq, at, replaceAGs[i,])
    AGseqs <- c(AGseqs, newseq)
  }
  AGseqs <- DNAStringSet(AGseqs)
  names(AGseqs) <- paste0(cpgid, "_", 1:length(AGseqs))
  AGseqs
}

idswith1ormoreRs <- manifest$IlmnID[manifest$numRs > 0]

rownames(manifest) <- manifest$IlmnID

seqs_1ormoreRs <- mclapply(idswith1ormoreRs, function (x) putativeAGseqs(manifest[x, "AlleleA_ProbeSeq"], x), mc.cores=16)
lengths <- unlist(lapply(seqs_1ormoreRs, length))

seqs_1ormoreRs <- unlist(DNAStringSetList(seqs_1ormoreRs))

final <- c(final, seqs_1ormoreRs)

writeXStringSet(final, file="allEPICseqs.fasta")




