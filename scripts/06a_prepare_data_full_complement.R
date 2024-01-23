library(sesame)
library(minfi)
library(limma)
library(data.table)
library(consensus)
library(rtracklayer)

#Only need to do once:
#sesameDataCache()

#Load Ruth's custom manifest
load("NewManifest.RData")
rownames(new_manifest) <- new_manifest$IlmnID

#EPICv2 first

setwd("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/EPICv2/idats/")
addr <- sesameAnno_buildAddressFile("../../EPICv2.hg38.manifest.tsv")
targets <- read.metharray.sheet("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/EPICv2/Samplesheet/")
targets$sesameID <- paste(targets$Slide, targets$Array, sep="_")

idats <- searchIDATprefixes(getwd())
#Read in
idatlist <- lapply(idats, function (x) readIDATpair(x, manifest = addr))
#Find TypeI probes needed for dye bias normalisation
idatlist <- lapply(idatlist, inferInfiniumIChannel)
#Dye bias detection
idatlist <- lapply(idatlist, dyeBiasNL)
#Detection P step
#Stricter than minfi::detectionP!!! So use 0.05
idatlist <- lapply(idatlist, function (x) pOOBAH(x, pval.threshold = 0.05))
#Find which of the technical replicates has fewer masked probes
maskedP <- do.call(cbind, lapply(idatlist, function (x) x$mask)) 
maskedP <- maskedP[,targets$sesameID]
idatlist <- idatlist[targets$sesameID]
colnames(maskedP) <- targets$Sample_Name
colSums(maskedP)
#PREC_500_1     LNCAP_500_1       SYN_12543       SYN_10651      PREC_500_2     LNCAP_500_2        SYN_5656       SYN_16599 MCF7_Cardiff_P7 
#42990           62465           43002           45542           49091           89416           50146           48191           54297 
#TAMR_P8      PREC_500_3     LNCAP_500_3       SYN-11452       SYN_15917      FD07743210      FD07745869      FD07745961      FD07745969 
#52156           48910           71552           37371           47784           76573           75897           76585           83232 
#FD07746819      FD07746841      FD07746850      FD07747708 
#70352           51644           44799           67121 

#PREC_500_1 and LNCAP_500_1 have the least dropouts, so we retain them

rm <- which(colnames(maskedP) %in% c("PREC_500_2", "LNCAP_500_2", "PREC_500_3", "LNCAP_500_3"))

maskedP <- maskedP[,-rm]
idatlist <- idatlist[-rm]
targets <- targets[-rm,]
stopifnot(all(colnames(maskedP)==targets$Sample_Name))
stopifnot(all(names(idatlist)==targets$sesameID))

Pmask <- apply(maskedP, 1, any)
sum(Pmask)
#139691

#Out-Of-Band array hybridisation
idatlist <- lapply(idatlist, noob)

#Extract betas. mask=FALSE is important bc noob() filters out cross-hybridisers
EPICv2_betas <- do.call(cbind, lapply(idatlist, function (x) getBetas(x, mask=FALSE)))
#Then post-hoc, filter out the P threshold
EPICv2_betas <- EPICv2_betas[!Pmask,]

#Organise by sample sheet

stopifnot(all(colnames(EPICv2_betas)==targets$sesameID))
colnames(EPICv2_betas) <- targets$Sample_Name

#Beta to M
EPICv2_Ms <- BetaValueToMValue(EPICv2_betas)
#Remove control and rs probes
EPICv2_Ms <- EPICv2_Ms[grep("^ctl|^rs", rownames(EPICv2_Ms), invert = T),]
#Align new manifest with Ilmn_ID
new_manifest <- new_manifest[rownames(EPICv2_Ms),]


rmvecdisc <- new_manifest$vecdisc=="Y"
EPICv2_Ms <- EPICv2_Ms[!rmvecdisc,]
new_manifest <- new_manifest[!rmvecdisc,]
stopifnot(all(rownames(EPICv2_Ms)==rownames(new_manifest)))

#Get annotation, downloaded from http://zwdzwd.github.io/InfiniumAnnotation
manifest <- sesameAnno_buildManifestGRanges("../../EPICv2.hg38.manifest.tsv", genome="hg38")
manifest <- manifest[rownames(EPICv2_Ms)]
stopifnot(all(names(manifest)==rownames(EPICv2_Ms)))

#Name target by genome coordinate
rownames(EPICv2_Ms) <- paste(new_manifest$CHR, new_manifest$MAPINFO, sep=":")

rm(EPICv2_betas, idatlist, manifest, maskedP, Pmask, idats, rm, targets)
gc()

#########################################################################

#EPICv1
#Start with minfi sample sheet as phenotype data
setwd("~/hippo/Gagri_Cancer-Epigenetics/Projects/EPIC_V2/Tim_Analysis/Crossplatform/EPICv1/idats/")
targets <- read.metharray.sheet("~/hippo/Gagri_Cancer-Epigenetics/Projects/EPIC_V2/Tim_Analysis/Crossplatform/EPICv1/")
targets$sesameID <- paste(targets$Slide, targets$Array, sep="_")

idats <- searchIDATprefixes(getwd())
#Read in
idatlist <- lapply(idats, readIDATpair)
#Find TypeI probes needed for dye bias normalisation
idatlist <- lapply(idatlist, inferInfiniumIChannel)
#Dye bias detection
idatlist <- lapply(idatlist, dyeBiasNL)
#Detection P step
#Stricter than minfi::detectionP!!! So use 0.05
idatlist <- lapply(idatlist, function (x) pOOBAH(x, pval.threshold = 0.05))
#Find which of the technical replicates has fewer masked probes
maskedP <- do.call(cbind, lapply(idatlist, function (x) x$mask)) 
maskedP <- maskedP[,targets$sesameID]
idatlist <- idatlist[targets$sesameID]
colnames(maskedP) <- targets$Sample_Name
colSums(maskedP)
#PrEC P5_1       PrEC P5_2     LNCaP P73_1     LNCaP P73_2       SYN_12543       SYN_10651        SYN_5656       SYN_16599       SYN_15917 
#26893           26081           28508           28331           30421           36514           37669           36450           40480 
#SYN-11452 MCF7_Cardiff_P7         TAMR_P8      FD07743210      FD07745869      FD07745961      FD07745969      FD07746819      FD07746841 
#40412           57407           48675           46840           48965           45040           50931           61100           40444 
#FD07746850      FD07747708 
#42285           52280 

#PrEC P5_1 and LNCaP P73_1 have more dropouts, so we remove them

rm <- which(colnames(maskedP) %in% c("PrEC P5_1", "LNCaP P73_1"))

maskedP <- maskedP[,-rm]
idatlist <- idatlist[-rm]
targets <- targets[-rm,]
stopifnot(all(colnames(maskedP)==targets$Sample_Name))
stopifnot(all(names(idatlist)==targets$sesameID))

Pmask <- apply(maskedP, 1, any)
sum(Pmask)

#115178
#Out-Of-Band array hybridisation
idatlist <- lapply(idatlist, noob)

#Extract betas. mask=FALSE is important bc noob() filters out cross-hybridisers
EPICv1_betas <- do.call(cbind, lapply(idatlist, function (x) getBetas(x, mask=FALSE)))
#Then post-hoc, filter out the P threshold
EPICv1_betas <- EPICv1_betas[!Pmask,]

#Organise by sample sheet

stopifnot(all(colnames(EPICv1_betas)==targets$sesameID))
colnames(EPICv1_betas) <- targets$Sample_Name

#Beta to M
EPICv1_Ms <- BetaValueToMValue(EPICv1_betas)
#Remove control and rs probes
EPICv1_Ms <- EPICv1_Ms[grep("^ctl|^rs", rownames(EPICv1_Ms), invert = T),]

#Get annotation, downloaded from http://zwdzwd.github.io/InfiniumAnnotation
manifest <- sesameAnno_buildManifestGRanges("../../EPIC.hg38.manifest.tsv", genome="hg38")
manifest <- manifest[rownames(EPICv1_Ms)]
stopifnot(all(names(manifest)==rownames(EPICv1_Ms)))

#Filter out probes with no target site
keep <- as.logical(seqnames(manifest)!="*")
EPICv1_Ms <- EPICv1_Ms[keep,]
manifest <- manifest[keep]
stopifnot(all(names(manifest)==rownames(EPICv1_Ms)))

#Name target by genome coordinate
rownames(EPICv1_Ms) <- paste(seqnames(manifest), start(manifest), sep=":")

rm(EPICv1_betas, idatlist, maskedP, Pmask, keep, idats, rm, targets)
gc()

#####################################################################################



#Conform EPICv1 and EPICv2 matrices
colnames(EPICv1_Ms) <- gsub("-", "_", gsub(" .*", "", colnames(EPICv1_Ms)))
colnames(EPICv2_Ms) <- gsub("-", "_", gsub("_500.*", "", colnames(EPICv2_Ms)))
colnames(EPICv2_Ms)[1:2] <- colnames(EPICv1_Ms)[1:2]

EPICv2_Ms <- EPICv2_Ms[,colnames(EPICv1_Ms)]

save(EPICv1_Ms, EPICv2_Ms, new_manifest, manifest, file="../../EPICdata.RData")

#Then WGBS
#LNCaP and PrEC liftover
setwd("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/")
lncapprectab <- fread("~/Bis-Seq_Level_3/darlo/Prostate/bigTable/bigTable.tsv.gz")
#All CpGs must have coverage
lncapprectab <- lncapprectab[lncapprectab$PrEC.cov > 0 & lncapprectab$LNCaP.cov > 0,]
CpGs <- GRanges(lncapprectab$chr, IRanges(lncapprectab$position, width=1))
values(CpGs) <- data.frame(PrEC=logit2((lncapprectab$PrEC.C + 0.5)/(lncapprectab$PrEC.cov + 1)), 
                           LNCaP=logit2((lncapprectab$LNCaP.C + 0.5)/(lncapprectab$LNCaP.cov + 1)))
#Liftover
ch <- import.chain("hg19ToHg38.over.chain")
CpGs.hg38 <- unlist(liftOver(CpGs, ch))
rm(CpGs)

#Synergistic and TAMR are already in hg38
syntab <- fread("~/Bis-Seq_Level_3/darlo_luu/Synergistic/bigTable/bigTable.tsv")
syntab <- syntab[syntab$G8_12543.cov > 0 & syntab$G8_10651.cov > 0 & syntab$G6_5656.cov > 0 & syntab$G8_16599.cov > 0 & syntab$G8_15917.cov > 0 & syntab$G6_11452.cov > 0,]
CpGs <- GRanges(syntab$chr, IRanges(syntab$position, width=1))
values(CpGs) <- data.frame(SYN_12543=logit2((syntab$G8_12543.C + 0.5)/(syntab$G8_12543.cov + 1)), 
                           SYN_10651=logit2((syntab$G8_10651.C + 0.5)/(syntab$G8_10651.cov + 1)),
                           SYN_5656=logit2((syntab$G6_5656.C + 0.5)/(syntab$G6_5656.cov + 1)),
                           SYN_16599=logit2((syntab$G8_16599.C + 0.5)/(syntab$G8_16599.cov + 1)),
                           SYN_15917=logit2((syntab$G8_15917.C + 0.5)/(syntab$G8_15917.cov + 1)),
                           SYN_11452=logit2((syntab$G6_11452.C + 0.5)/(syntab$G6_11452.cov + 1)))

CpGs.hg38 <- CpGs.hg38[CpGs.hg38 %over% CpGs]

ol <- findOverlaps(CpGs.hg38, CpGs)
stopifnot(all(1:length(ol)==queryHits(ol)))
values(CpGs.hg38) <- cbind(values(CpGs.hg38), values(CpGs)[subjectHits(ol),])

#TAMR samples

tamrtab <- fread("~/joaach/bigtable_for_tim/bigTable.tsv.gz")
tamrtab <- tamrtab[tamrtab$MCF7_merged.cov > 0 & tamrtab$TAMR_merged.cov > 0,]
CpGs <- GRanges(tamrtab$`#chr`, IRanges(tamrtab$position, width=1))
values(CpGs) <- data.frame(MCF7_Cardiff_P7=logit2((tamrtab$MCF7_merged.C + 0.5)/(tamrtab$MCF7_merged.cov + 1)), 
                           TAMR_P8=logit2((tamrtab$TAMR_merged.C + 0.5)/(tamrtab$TAMR_merged.cov + 1)))

CpGs.hg38 <- CpGs.hg38[CpGs.hg38 %over% CpGs]

ol <- findOverlaps(CpGs.hg38, CpGs)
stopifnot(all(1:length(ol)==queryHits(ol)))
values(CpGs.hg38) <- cbind(values(CpGs.hg38), values(CpGs)[subjectHits(ol),])

#PDX samples
pdxtab <- fread("~/Projects/EPIC_V2/Tim_Analysis/Crossplatform/bigTable_fix.tsv")
pdxtab <- pdxtab[pdxtab$FD07747708.cov > 0 & pdxtab$FD07746819.cov > 0 & pdxtab$FD07746841.cov > 0 & pdxtab$FD07746850.cov > 0 & 
                 pdxtab$FD07743210.cov > 0 & pdxtab$FD07745961.cov > 0 & pdxtab$FD07745969.cov > 0 & pdxtab$FD07745869.cov > 0,]
CpGs <- GRanges(pdxtab$`#chr`, IRanges(pdxtab$position, width=1))
values(CpGs) <- data.frame(FD07743210=logit2((pdxtab$FD07743210.C + 0.5)/(pdxtab$FD07743210.cov + 1)), 
                           FD07745869=logit2((pdxtab$FD07745869.C + 0.5)/(pdxtab$FD07745869.cov + 1)),
                           FD07745961=logit2((pdxtab$FD07745961.C + 0.5)/(pdxtab$FD07745961.cov + 1)),
                           FD07745969=logit2((pdxtab$FD07745969.C + 0.5)/(pdxtab$FD07745969.cov + 1)),
                           FD07746819=logit2((pdxtab$FD07746819.C + 0.5)/(pdxtab$FD07746819.cov + 1)),
                           FD07746841=logit2((pdxtab$FD07746841.C + 0.5)/(pdxtab$FD07746841.cov + 1)),
                           FD07746850=logit2((pdxtab$FD07746850.C + 0.5)/(pdxtab$FD07746850.cov + 1)),
                           FD07747708=logit2((pdxtab$FD07747708.C + 0.5)/(pdxtab$FD07747708.cov + 1)))

CpGs.hg38 <- CpGs.hg38[CpGs.hg38 %over% CpGs]

ol <- findOverlaps(CpGs.hg38, CpGs)
stopifnot(all(1:length(ol)==queryHits(ol)))
values(CpGs.hg38) <- cbind(values(CpGs.hg38), values(CpGs)[subjectHits(ol),])

save(CpGs.hg38, file="../../WGBS_allCpGs.RData")

#Then conform to EPIC

common <- intersect(rownames(EPICv1_Ms), rownames(EPICv2_Ms))
manifest <- manifest[rownames(EPICv1_Ms) %in% common]
EPICv1_Ms <- EPICv1_Ms[rownames(EPICv1_Ms) %in% common,]
stopifnot(all(rownames(EPICv1_Ms)==paste(seqnames(manifest), start(manifest), sep=":")))
new_manifest <- new_manifest[rownames(EPICv2_Ms) %in% common,]
EPICv2_Ms <- EPICv2_Ms[rownames(EPICv2_Ms) %in% common,]
stopifnot(all(rownames(EPICv2_Ms)==paste(new_manifest$CHR, new_manifest$MAPINFO, sep=":")))
stopifnot(all(colnames(EPICv1_Ms)==colnames(EPICv2_Ms)))

stopifnot(all(colnames(values(CpGs.hg38))==colnames(EPICv2_Ms)))
stopifnot(all(colnames(values(CpGs.hg38))==colnames(EPICv1_Ms)))
WGBS_Ms <- values(CpGs.hg38)
rownames(WGBS_Ms) <- paste(CpGs.hg38)

common <- intersect(rownames(WGBS_Ms), rownames(EPICv2_Ms))
manifest <- manifest[rownames(EPICv1_Ms) %in% common]
EPICv1_Ms <- EPICv1_Ms[rownames(EPICv1_Ms) %in% common,]
stopifnot(all(rownames(EPICv1_Ms)==paste(seqnames(manifest), start(manifest), sep=":")))
new_manifest <- new_manifest[rownames(EPICv2_Ms) %in% common,]
EPICv2_Ms <- EPICv2_Ms[rownames(EPICv2_Ms) %in% common,]
stopifnot(all(rownames(EPICv2_Ms)==paste(new_manifest$CHR, new_manifest$MAPINFO, sep=":")))
stopifnot(all(colnames(EPICv1_Ms)==colnames(EPICv2_Ms)))

WGBS_Ms <- data.matrix(WGBS_Ms[rownames(WGBS_Ms) %in% common,])

#Then cut down new manifest

save(EPICv1_Ms, EPICv2_Ms, WGBS_Ms, manifest, new_manifest, file="Mvalues_FINAL.RData")

