#Install packages and data
#install.packages("remotes")
#remotes::install_github("zwdzwd/sesame")
#remotes::install_github("zwdzwd/sesameData")
#remotes::install_github("zwdzwd/sesameData")
#sesameDataCache()

#load in packages
library(sesame)
library(minfi)
library(ggplot2)
library(ggfortify)
library(Hmisc)
library(corrplot)
library(reshape2)
library(RColorBrewer)

filepath<-"../your/file/path/"

#read in EPICv2 manifest (sesame) and order by addresses
addr <- sesameAnno_buildAddressFile(paste(filepath,"metadata/EPICv2.hg38.manifest.tsv",sep=""))
manifest_full <- read.table(paste(filepath,"metadata/EPICv2.hg38.manifest.tsv",sep=""), sep="\t", header=T)

#read in EPICv2 manifest (Illumina) to cross-check next base annotation
manifest_EPICv2_Ilmn <- read.csv(paste(filepath,"metadata/EPIC-8v2-0_A1.csv",sep=""), skip = 7)

#re-order Illumina manifest to match sesame
m2<-manifest_EPICv2_Ilmn[match(manifest_full$Probe_ID,manifest_EPICv2_Ilmn$IlmnID),]

#compare color channel for sesame and Illumina manifest
col_compare<-(m2$col==manifest_full$channel)

table(col_compare)
#col_compare
#FALSE   TRUE 
#478 127059 

mismatchcol<-m2$Strand_CO[which(col_compare==FALSE)]
table(mismatchcol)
#mismatchcol
#C   O 
#45 433 

#Some color channel differences between manifests, a disproportionate number on opposite strand
#Accounts for error observed in cg03254865_BO11

m2[which(m2$IlmnID=="cg03254865_BO11"),1:12]
# IlmnID       Name AddressA_ID
# 6368 cg03254865_BO11 cg03254865    25724474
# AlleleA_ProbeSeq AddressB_ID
# 6368 AAAAGAAAAAGGTATGGTTAAATGTGAATATGTGTTGAAATGGTTAGTGT    24600432
# AlleleB_ProbeSeq Next_Base Color_Channel
# 6368 AAAAGAAAAAGGTATGGTTAAATGTGAATATGTGTTGAAATGGTTAGTGC         T           Red
# col Probe_Type Strand_FR Strand_TB
# 6368   R         cg         F         B

manifest_full[which(manifest_full$Probe_ID=="cg03254865_BO11"),1:12]
# CpG_chrm  CpG_beg  CpG_end address_A address_B target nextBase channel
# 118888     chr5 71494300 71494302  25724474  24600432     CG        G       G
# Probe_ID mapFlag_A mapChrm_A mapPos_A
# 118888 cg03254865_BO11         0      chr5 71494252

#read in sample sheet. #make sample sheet corresponding to EPICv2 data in GEO GSE240482
targets <- read.metharray.sheet("sampledata/")
#extract barcode and position information
targets$sesameID <- paste(targets$Slide, targets$Array, sep="_")

#find IDAT files in specified folder # download EPICv2 data from GEO GSE240482
idats <- searchIDATprefixes("sampledata/")

#Read in IDATs
idatlist <- lapply(idats, function (x) readIDATpair(x, manifest = addr))

#Find Type I probes needed for dye bias normalisation
idatlist <- lapply(idatlist, inferInfiniumIChannel)
#Dye bias detection
idatlist <- lapply(idatlist, dyeBiasNL)
#Detection P step
#Stricter than minfi::detectionP!!! So use 0.05
idatlist <- lapply(idatlist, function (x) pOOBAH(x, pval.threshold = 0.05))

#Find which of the samples has the most masked probes based on pOOBAH
maskedP <- do.call(cbind, lapply(idatlist, function (x) x$mask)) 

#order detection p-value by order of samples in sample sheet
maskedP <- maskedP[,targets$sesameID]

#order list of methylation data by order of samples in sample sheet
idatlist <- idatlist[targets$sesameID]

colnames(maskedP) <- targets$Sample_Name
colSums(maskedP)


#Barplot of number of probes that failed detection p-value for each sample
pdf(paste(filepath,"plots/Barplot_DetecP.pdf",sep=""),width=12)
par(mar=c(12,4,2,2))
barplot(colSums(maskedP),las=2)
abline(h=(nrow(EPICv2_betas)/100)*10,lty=2)
dev.off()


#Write out number of probes that failed detection p-value for each sample
write.csv(colSums(maskedP),paste(filepath,"results/DetecP.csv",sep=""))

Pmask <- apply(maskedP, 1, any)
sum(Pmask)
#331252

#Out-Of-Band array hybridisation
idatlist <- lapply(idatlist, noob)

#Extract betas. mask=TRUE to remove suspicious probes such as cross-hybridisers
EPICv2_betas <- do.call(cbind, lapply(idatlist, function (x) getBetas(x, mask=TRUE)))
dim(EPICv2_betas)#[1] 936866     40

#prepare for data visualisation
rs<-grep("rs",rownames(EPICv2_betas))#65
nv<-grep("nv",rownames(EPICv2_betas))#0
ctl<-grep("ctl",rownames(EPICv2_betas))#635

EPICv2_betas2<-EPICv2_betas[-c(rs,ctl),]

#not going to filter on detection p-value as set to NA, and samples are expected to be compromised (e.g. lower input) and don't want them causing other samples' data to be excluded

#Organise by sample sheet
stopifnot(all(colnames(EPICv2_betas)==targets$sesameID))
colnames(EPICv2_betas) <- targets$Sample_Name

#extract methylation levels at known SNPs
SNPmat<-EPICv2_betas[grep("rs",rownames(EPICv2_betas)),]
nba.m <- melt(SNPmat)

# Cluster SNP methylation based on euclidean distance
clust <- hclust(dist(t(SNPmat)))

#Plot heatmap of SNP methylation levels
pdf(paste(filepath,"plots/SNP_heatmap.pdf",sep=""),width=10,height=10)
ggplot(nba.m, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value), colour = "white") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  scale_x_discrete(limits = colnames(SNPmat)[clust$order])
dev.off()

#read in additional sample data for visualisation
extrasamp<-read.csv("sampledata/SampleSheet_V2_additional.csv")
all(extrasamp$Sample_Name==targets$Sample_Name)#TRUE

source(paste(filepath,"functions/densityPlot2.R",sep=""))
pdf(paste(filepath,"plots/DensityPlot.pdf",sep=""),width=7,height=10)
par(mfrow=c(3,2))
#Plot LNCap and PrEC, 500ng only
LP500<-which((extrasamp$Sample=="PrEC" | extrasamp$Sample=="LNCaP") & extrasamp$Input==500)
densityPlot2(EPICv2_betas[,LP500],sampGroups=extrasamp$Sample[LP500], main="LNCaP and PrEC (500ng)")

#Plot PrEC only, all concentrations
Ponly<-which(extrasamp$Sample=="PrEC")
densityPlot2(EPICv2_betas[,Ponly],sampGroups=extrasamp$Input[Ponly], main="PrEC - varied DNA input")

#Plot LNCaP only, all concentrations
Lonly<-which(extrasamp$Sample=="LNCaP")
densityPlot2(EPICv2_betas[,Lonly],sampGroups=extrasamp$Input[Lonly], main="LNCaP - varied DNA input")

#Plot TAMR only
Tonly<-which(extrasamp$Group=="BreastCell")
densityPlot2(EPICv2_betas[,Tonly],sampGroups=extrasamp$Sample[Tonly], main="Breast cancer cell lines")

#Plot PDX only 
PDXonly<-which(extrasamp$Group=="BreastPDX")
densityPlot2(EPICv2_betas[,PDXonly],sampGroups=extrasamp$Sample[PDXonly], main="Breast cancer PDX")

#Plot Prostate Tumour only 
PrCaonly<-which(extrasamp$Group=="ProstateTumour")
#densityPlot(EPICv2_betas[,PrCaonly])
densityPlot2(EPICv2_betas[,PrCaonly],sampGroups=extrasamp$Sample_Name[PrCaonly], main="Prostate tumour tissue")
dev.off()

pdf(paste(filepath,"plots/MDSPlot_all.pdf",sep=""),width=9)
par(mfrow=c(1,1))
mdsPlot(EPICv2_betas2,sampGroup=extrasamp$Sample2,pch=16,numPositions=10000)
#mdsPlot(EPICv2_betas2,sampGroup=extrasamp$Sentrix_ID,pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2,sampNames=extrasamp$Sample_Name,sampGroup=extrasamp$Sample2,pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2,sampGroup=extrasamp$Tissue,pch=16,numPositions=10000)
dev.off()

#Compare B value densities of Type I and II probes for all normal samples
#PrEC, TAMR, MCF7, PrCa, PDX (no drug)
manifest_full_2<-manifest_full[match(rownames(EPICv2_betas2),manifest_full$Probe_ID),]
t2<-which(manifest_full_2$type=="II")
t1<-which(manifest_full_2$type=="I")

pdf(paste(filepath,"plots/DensityPlot_TypeI_II.pdf",sep=""),width=7,height=10)
par(mfrow=c(3,2))

#Plot PrEC only, separated by type I and II
P500only<-which(extrasamp$Sample=="PrEC" & extrasamp$Input==500)
par(lty=1)
densityPlot2(EPICv2_betas[t1,P500only],sampGroups=as.factor(extrasamp$Sample_Name[P500only]),main="PrEC (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[t2,P500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[P500only]),legend=F)

#Plot LNCaP only, separated by type I and II
L500only<-which(extrasamp$Sample=="LNCaP" & extrasamp$Input==500)
par(lty=1)
densityPlot2(EPICv2_betas[t1,L500only],sampGroups=as.factor(extrasamp$Sample_Name[L500only]),main="LNCaP (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[t2,L500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[L500only]),legend=F)

#Plot TAMR and MCF7 only
TMonly<-which(extrasamp$Sample=="TAMR" | extrasamp$Sample=="MCF7")
par(lty=1)
densityPlot2(EPICv2_betas[t1,TMonly],lty=1,sampGroups=as.factor(extrasamp$Sample[TMonly]),ylim=c(0,8),main="Breast Cancer Cell lines")
par(lty=2)
densityPlot2(EPICv2_betas[t2,TMonly],add=FALSE,lty=2,sampGroups=as.factor(extrasamp$Sample[TMonly]),lty=2,legend=F)

#Plot PDX no drug only
PDXnoonly<-which(extrasamp$Group=="BreastPDX" & extrasamp$Drug=="No")
par(lty=1)
densityPlot2(EPICv2_betas[t1,PDXnoonly],lty=1,sampGroups=as.factor(extrasamp$Sample[PDXnoonly]),main="Breast Cancer PDX")
par(lty=2)
densityPlot2(EPICv2_betas[t2,PDXnoonly],add=FALSE,lty=2,sampGroups=as.factor(extrasamp$Sample[PDXnoonly]),legend=F)

#Plot PrCa only, separated by type I and II
par(lty=1)
densityPlot2(EPICv2_betas[t1,PrCaonly],sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),main="Prostate Tumour Tissue")
par(lty=2)
densityPlot2(EPICv2_betas[t2,PrCaonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),legend=F)

plot(1:10,1:10,axes=F,col="white",xlab="",ylab="")
legend("center", legend=c("Type I","Type II"),lty=c(1,2),cex=2)

dev.off()

#Compare B value densities of retained, reinstated and new probes for all normal samples
#PrEC, TAMR, MCF7, PrCa, PDX (no drug)
#read in new_manifest
load(paste(filepath,"results/NewManifest.RData",sep=""))

#create index of probes
new_manifest2<-new_manifest[match(rownames(EPICv2_betas2),new_manifest$IlmnID),]
retained<-which(!is.na(new_manifest2$epic1seqmatch))
reinstated<-which(is.na(new_manifest2$epic1seqmatch) & (!is.na(new_manifest2$K450seqmatch) | !is.na(new_manifest2$K27seqmatch)) & !is.na(new_manifest2$IlmnID))
new<-which(is.na(new_manifest2$epic1seqmatch) & is.na(new_manifest2$K450seqmatch) & is.na(new_manifest2$K27seqmatch) & !is.na(new_manifest2$IlmnID))

retainedI<-which(!is.na(new_manifest2$epic1seqmatch) & new_manifest2$Infinium_Design_Type=="I")
reinstatedI<-which(is.na(new_manifest2$epic1seqmatch) & (!is.na(new_manifest2$K450seqmatch) | !is.na(new_manifest2$K27seqmatch)) & !is.na(new_manifest2$IlmnID) & new_manifest2$Infinium_Design_Type=="I")
newI<-which(is.na(new_manifest2$epic1seqmatch) & is.na(new_manifest2$K450seqmatch) & is.na(new_manifest2$K27seqmatch) & !is.na(new_manifest2$IlmnID) & new_manifest2$Infinium_Design_Type=="I")

retainedII<-which(!is.na(new_manifest2$epic1seqmatch) & new_manifest2$Infinium_Design_Type=="II")
reinstatedII<-which(is.na(new_manifest2$epic1seqmatch) & (!is.na(new_manifest2$K450seqmatch) | !is.na(new_manifest2$K27seqmatch)) & !is.na(new_manifest2$IlmnID) & new_manifest2$Infinium_Design_Type=="II")
newII<-which(is.na(new_manifest2$epic1seqmatch) & is.na(new_manifest2$K450seqmatch) & is.na(new_manifest2$K27seqmatch) & !is.na(new_manifest2$IlmnID) & new_manifest2$Infinium_Design_Type=="II")


pdf(paste(filepath,"plots/DensityPlot_ProbeGroup.pdf",sep=""),width=7,height=10)
par(mfrow=c(3,2))

#Plot PrEC only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retained,P500only],sampGroups=as.factor(extrasamp$Sample_Name[P500only]),main="PrEC (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[reinstated,P500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[P500only]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[new,P500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[P500only]),legend=F)

#Plot LNCaP only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retained,L500only],sampGroups=as.factor(extrasamp$Sample_Name[L500only]),main="LNCaP (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[reinstated,L500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[L500only]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[new,L500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[L500only]),legend=F)

#Plot TAMR and MCF7 only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retained,TMonly],sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),main="Breast Cancer Cell lines")
par(lty=2)
densityPlot2(EPICv2_betas[reinstated,TMonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[new,TMonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),legend=F)


#Plot PDX no drug only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retained,PDXnoonly],sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),main="Breast Cancer PDX")
par(lty=2)
densityPlot2(EPICv2_betas[reinstated,PDXnoonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[new,PDXnoonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),legend=F)

#Plot PrCa only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retained,PrCaonly],sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),main="Prostate Tumour Tissue")
par(lty=2)
densityPlot2(EPICv2_betas[reinstated,PrCaonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[new,PrCaonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),legend=F)

plot(1:10,1:10,axes=F,col="white",xlab="",ylab="")
legend("center", legend=c("Retained","Reinstated","New"),lty=c(1,2,3),cex=2)
dev.off()

#For Type I and II probes separately
pdf(paste(filepath,"plots/DensityPlot_ProbeGroup_TypeI.pdf",sep=""),width=7,height=10)
par(mfrow=c(3,2))

#Plot PrEC only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedI,P500only],sampGroups=as.factor(extrasamp$Sample_Name[P500only]),main="PrEC (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedI,P500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[P500only]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newI,P500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[P500only]),legend=F)

#Plot LNCaP only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedI,L500only],sampGroups=as.factor(extrasamp$Sample_Name[L500only]),main="LNCaP (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedI,L500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[L500only]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newI,L500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[L500only]),legend=F)

#Plot TAMR and MCF7 only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedI,TMonly],sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),main="Breast Cancer Cell lines")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedI,TMonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newI,TMonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),legend=F)


#Plot PDX no drug only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedI,PDXnoonly],sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),main="Breast Cancer PDX")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedI,PDXnoonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newI,PDXnoonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),legend=F)

#Plot PrCa only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedI,PrCaonly],sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),main="Prostate Tumour Tissue")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedI,PrCaonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newI,PrCaonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),legend=F)

plot(1:10,1:10,axes=F,col="white",xlab="",ylab="")
legend("center", legend=c("Retained","Reinstated","New"),lty=c(1,2,3),cex=2)
dev.off()


pdf(paste(filepath,"plots/DensityPlot_ProbeGroup_TypeII.pdf",sep=""),width=7,height=10)
par(mfrow=c(3,2))

#Plot PrEC only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedII,P500only],sampGroups=as.factor(extrasamp$Sample_Name[P500only]),main="PrEC (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedII,P500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[P500only]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newII,P500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[P500only]),legend=F)

#Plot LNCaP only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedII,L500only],sampGroups=as.factor(extrasamp$Sample_Name[L500only]),main="LNCaP (500ng)")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedII,L500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[L500only]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newII,L500only],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[L500only]),legend=F)

#Plot TAMR and MCF7 only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedII,TMonly],sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),main="Breast Cancer Cell lines")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedII,TMonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newII,TMonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[TMonly]),legend=F)


#Plot PDX no drug only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedII,PDXnoonly],sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),main="Breast Cancer PDX")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedII,PDXnoonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newII,PDXnoonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PDXnoonly]),legend=F)

#Plot PrCa only, separated by retained, reinstated, new
par(lty=1)
densityPlot2(EPICv2_betas[retainedII,PrCaonly],sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),main="Prostate Tumour Tissue")
par(lty=2)
densityPlot2(EPICv2_betas[reinstatedII,PrCaonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),legend=F)
par(lty=3)
densityPlot2(EPICv2_betas[newII,PrCaonly],add=FALSE,sampGroups=as.factor(extrasamp$Sample_Name[PrCaonly]),legend=F)

plot(1:10,1:10,axes=F,col="white",xlab="",ylab="")
legend("center", legend=c("Retained","Reinstated","New"),lty=c(1,2,3),cex=2)
dev.off()
###Comparing technical replicates

#Plot PrEC only
PrEC<-which(extrasamp$Sample=="PrEC")
par(lty=1)
densityPlot(EPICv2_betas[,PrEC],sampGroups=as.factor(extrasamp$Input[PrEC]))

LNCaP<-which(extrasamp$Sample=="LNCaP")
par(lty=1)
densityPlot(EPICv2_betas[,LNCaP],sampGroups=as.factor(extrasamp$Input[LNCaP]))

mdsPlot(EPICv2_betas2[,PrEC],sampGroup=as.factor(extrasamp$Input[PrEC]),pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2[,PrEC],sampGroup=as.factor(extrasamp$Sentrix_ID[PrEC]),pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2[,PrEC],sampGroup=as.factor(extrasamp$Sentrix_ID[PrEC]),sampNames=as.factor(extrasamp$Input[PrEC]),pch=16,numPositions=10000)

mdsPlot(EPICv2_betas2[,LNCaP],sampGroup=as.factor(extrasamp$Input[LNCaP]),pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2[,LNCaP],sampGroup=as.factor(extrasamp$Sentrix_ID[LNCaP]),pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2[,LNCaP],sampGroup=as.factor(extrasamp$Sentrix_ID[LNCaP]),sampNames=as.factor(extrasamp$Input[LNCaP]),pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2[,LNCaP[-2]],sampGroup=as.factor(extrasamp$Sentrix_ID[LNCaP[-2]]),sampNames=as.factor(extrasamp$Input[LNCaP[-2]]),pch=16,numPositions=10000)

LP<-c(LNCaP,PrEC)
mdsPlot(EPICv2_betas2[,LP],sampGroup=as.factor(extrasamp$Sentrix_ID[LP]),sampNames=as.factor(extrasamp$Input[LP]),pch=16,numPositions=10000)
mdsPlot(EPICv2_betas2[,LP[-2]],sampGroup=as.factor(extrasamp$Sentrix_ID[LP[-2]]),sampNames=as.factor(extrasamp$Input[LP[-2]]),pch=16,numPositions=10000)

detec<-(colSums(maskedP))[LNCaP]
Lbind<-cbind(extrasamp[LNCaP,],detec)
Lbind2<-Lbind[,c("Sentrix_ID","Input","detec")]
Lbind2$Sentrix_ID<-as.factor(Lbind2$Sentrix_ID)
Lbind2$Input<-factor(Lbind2$Input,levels=c(500,250,125))
pdf(paste(filepath, "plots/DetecP_LNCapPrEC.pdf",sep=""),width=7,height=7)
ggplot(Lbind2, aes(fill=Input, y=detec, x=Sentrix_ID)) + 
         geom_bar(position="dodge", stat="identity") +
        ggtitle("LNCaP") +
        labs(y="# failed probes")
ggplot(Pbind2, aes(fill=Input, y=detec, x=Sentrix_ID)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("PrEC")+
  labs(y="# failed probes")
dev.off()
detec<-(colSums(maskedP))[PrEC]
Pbind<-cbind(extrasamp[PrEC,],detec)
Pbind2<-Pbind[,c("Sentrix_ID","Input","detec")]
Pbind2$Sentrix_ID<-as.factor(Pbind2$Sentrix_ID)
Pbind2$Input<-factor(Pbind2$Input,levels=c(500,250,125))
ggplot(Pbind2, aes(fill=Input, y=detec, x=Sentrix_ID)) + 
  geom_bar(position="dodge", stat="identity") +
  ggtitle("PrEC")

L500only<-which(extrasamp$Sample=="LNCaP" & extrasamp$Input==500)  
mdsPlot(EPICv2_betas2[,c(P500only,L500only)],sampGroup=as.factor(extrasamp$Sample[c(P500only,L500only)]),pch=16,numPositions=1000)


#Correlation matrix of technical replicates - restricted to LNCaP or PrEC

res <- cor(EPICv2_betas2[,c(P500only,L500only)],use="pairwise.complete.obs")
pdf(paste(filepath,"plots/CorrPlot_LNCaP_PrEC500.pdf",sep=""),width=7,height=7)
corrplot(res, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,method="number",number.digits=3)
dev.off()

sscat<-function(a,b,c){
  est<-cor.test(EPICv2_betas2[,c[a]],EPICv2_betas2[,c[b]])$est
  p<-cor.test(EPICv2_betas2[,c[a]],EPICv2_betas2[,c[b]])$p.
  smoothScatter(EPICv2_betas2[,c[a]],EPICv2_betas2[,c[b]],xlab=extrasamp$Sample_Name[c[a]],ylab=extrasamp$Sample_Name[c[b]])
  legend("bottomright", legend=paste("cor=",round(est,3)," p=",p,sep=""))}

pdf(paste(filepath,"plots/SmoothScatter_LNCaP_PrEC500.pdf",sep=""),width=10,height=7)
par(mfrow=c(2,3))
sscat(1,2,P500only)
sscat(1,3,P500only)
sscat(2,3,P500only)
sscat(1,2,L500only)
sscat(1,3,L500only)
sscat(2,3,L500only)
dev.off()
#comparing input levels
res2 <- cor(EPICv2_betas2[,c(LNCaP,PrEC)],use="pairwise.complete.obs",method="pearson")
pdf(paste(filepath,"plots/CorrPlot_LNCaP_PrEC_all.pdf",sep=""),width=7,height=7)
corrplot(res2, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45,method="number",number.digits = 3,number.cex=0.5)
dev.off()
#LNCaP 500 only
res2_temp<-res2[c(1,4,7),c(1,4,7)]
min(res2_temp) #0.99649

#LNCaP 125 only
res2_temp<-res2[c(3,6,9),c(3,6,9)]
min(res2_temp) #0.9955612

#PrEC 500 only
res2_temp<-res2[c(10,13,16),c(10,13,16)]
min(res2_temp) #0.9969758

#PrEC 125 only
res2_temp<-res2[c(12,15,18),c(12,15,18)]
min(res2_temp) #0.9961904


sscat2<-function(a,b,c){
  est<-cor.test(EPICv2_betas2[,c[a]],EPICv2_betas2[,c[b]])$est
  p<-cor.test(EPICv2_betas2[,c[a]],EPICv2_betas2[,c[b]])$p.
  smoothScatter(EPICv2_betas2[,c[a]],EPICv2_betas2[,c[b]],xlab=extrasamp$Sample_Name[c[a]],ylab=extrasamp$Sample_Name[c[b]])
  legend("bottomright", legend=paste("cor=",round(est,3)," p=",p,sep=""))}

pdf(paste(filepath,"plots/SmoothScatter_LNCaP_PrEC_all.pdf",sep=""),width=7,height=8)
par(mfrow=c(3,3))
sscat2(1,4,PrEC)
sscat2(1,7,PrEC)
sscat2(4,7,PrEC)
sscat2(2,5,PrEC)
sscat2(2,8,PrEC)
sscat2(5,8,PrEC)
sscat2(3,6,PrEC)
sscat2(3,9,PrEC)
sscat2(6,9,PrEC)
sscat2(1,4,LNCaP)
sscat2(1,7,LNCaP)
sscat2(4,7,LNCaP)
sscat2(2,5,LNCaP)
sscat2(2,8,LNCaP)
sscat2(5,8,LNCaP)
sscat2(3,6,LNCaP)
sscat2(3,9,LNCaP)
sscat2(6,9,LNCaP)
dev.off()


#write out bed files of L500, P500 and PrCas
EPICsmall<-EPICv2_betas2[,c(L500only,P500only,PrCaonly)]
colnames(EPICsmall)<-targets$Sample_Name[c(L500only,P500only,PrCaonly)]
bedco<-cbind(new_manifest2$CHR,start=new_manifest2$MAPINFO,end=new_manifest2$MAPINFO)
for(i in 1:ncol(EPICsmall)){
  temp<-cbind(bedco,EPICsmall[,i])
  nas<-which(is.na(temp[,4])==TRUE)
  temp<-temp[-nas,]
  write.table(temp,quote=F,sep="\t",row.names=F,col.names=F,file=paste(filepath,colnames(EPICsmall)[i],".bedGraph",sep=""))
}

paste(extrasamp$Sample2, extrasamp$Input)
temp<-c(10,11,13,14)
EPICsmall2<-EPICv2_betas2[,temp]
colnames(EPICsmall2)<-targets$Sample_Name[temp]
bedco<-cbind(new_manifest2$CHR,start=new_manifest2$MAPINFO,end=new_manifest2$MAPINFO)
for(i in 1:ncol(EPICsmall2)){
  temp2<-cbind(bedco,EPICsmall2[,i])
  nas<-which(is.na(temp2[,4])==TRUE)
  temp2<-temp2[-nas,]
  write.table(temp2,quote=F,sep="\t",row.names=F,col.names=F,file=paste(filepath,colnames(EPICsmall2)[i],".bedGraph",sep=""))
}

#RLM plots
#Extract out beta values of input replicates and convert from beta to M-values
EPICv2LP_Ms <- BetaValueToMValue(EPICv2_betas[,LP])

##RLE plots of different input levels (on M-values)
# Substract medians for RLE plots
medLP = apply(EPICv2LP_Ms, 1, median)
RLM = EPICv2LP_Ms - medLP 

#re-doing without outlying sample
# Substract medians for RLE plots
medLP_B = apply(EPICv2LP_Ms[,-2], 1, median)
RLM_B = EPICv2LP_Ms[,-2] - medLP_B 
RLM_B<-cbind(RLM_B,rep(0,nrow(RLM_B)))
colnames(RLM_B)[18]<-colnames(EPICv2LP_Ms)[2]

#Repeat LNCap only (without outlying sample)
medLP_BL = apply(EPICv2LP_Ms[,-c(2,10,13,16,11,14,17,12,15,18)], 1, median)
RLM_BL = EPICv2LP_Ms[,-c(2,10,13,16,11,14,17,12,15,18)] - medLP_BL 
RLM_BL<-cbind(RLM_BL,rep(0,nrow(RLM_BL)))
colnames(RLM_BL)[9]<-colnames(EPICv2LP_Ms)[2]

#Repeat PrEC only (without outlying sample)

medLP_BP = apply(EPICv2LP_Ms[,-c(1,4,7,2,5,8,3,6,9)], 1, median)
RLM_BP = EPICv2LP_Ms[,-c(1,4,7,2,5,8,3,6,9)] - medLP_BP 

pdf(paste(filepath,"plots/RLM_plots_Nov23.pdf",sep=""),width=6)
par(mar=c(9, 4, 3, 2) + 0.1)
boxplot(RLM[,c(1,4,7,2,5,8,3,6,9,10,13,16,11,14,17,12,15,18)],outline=FALSE,ylim=c(-3,3), ylab="Relative Log Methylation Value",las=2, 
        col=c(rep("dark red",3),rep("red",3),rep("pink",3),rep("dark blue",3),rep("blue",3),rep("light blue",3))) 
boxplot(RLM_B[,c(1,3,6,18,4,7,2,5,8,9,12,15,10,13,16,11,14,17)],outline=FALSE,ylim=c(-3,3), ylab="Relative Log Methylation Value",las=2, 
        col=c(rep("dark red",3),rep("red",3),rep("pink",3),rep("dark blue",3),rep("blue",3),rep("light blue",3))) 
boxplot(RLM_BL[,c(1,3,6,9,4,7,2,5,8)],outline=FALSE,ylim=c(-1,1), ylab="Relative Log Methylation Value",las=2, 
        col=c(rep("dark red",3),rep("red",3),rep("pink",3))) 
boxplot(RLM_BP[,c(1,4,7,2,5,8,3,6,9)],outline=FALSE,ylim=c(-1,1), ylab="Relative Log Methylation Value",las=2, 
        col=c(rep("dark blue",3),rep("blue",3),rep("light blue",3))) 
dev.off()

#DMP analysis of different input levels (on M-values)
library(limma)

design<-model.matrix(~factor(c("Chip1","Chip2","Chip3","Chip1","Chip2","Chip3"))+factor(c("LNCaP","LNCaP","LNCaP","PreC","PreC","PreC")))
fit <- lmFit(EPICv2LP_Ms[,c(1,4,7,10,13,16)], design)
fit2 <- eBayes(fit)
tt_nocov<- topTable(fit2, coef = 4, number = nrow(EPICv2LP_Ms),sort.by="none")
EPICv2_betasLP<-EPICv2_betas[,LP]
meandiff<-rowMeans(EPICv2_betasLP[,c(1,4,7)]-EPICv2_betasLP[,c(10,13,16)])
sigs<-which(tt_nocov$adj.P.Val<0.05 & abs(meandiff)>0.05)
length(sigs) #[1] 464220

fit_125 <- lmFit(EPICv2LP_Ms[,c(3,6,9,12,15,18)], design)
fit_125_2 <- eBayes(fit_125)
tt_nocov_125<- topTable(fit_125_2, coef = 4, number = nrow(EPICv2LP_Ms),sort.by="none")
meandiff_125<-rowMeans(EPICv2_betasLP[,c(3,6,9)]-EPICv2_betasLP[,c(12,15,18)])
sigs_125<-which(tt_nocov_125$adj.P.Val<0.05 & abs(meandiff_125)>0.05)
length(sigs_125)#[1] 425892

int<-intersect(sigs,sigs_125)
length(int) # 409576

#CNV calling


#CNV calling in LNCaP vs PrEC
segtest_500<-cnSegmentation(idatlist[[12]],sdfs.normal=idatlist[c(1,9,25)])
segtest_250<-cnSegmentation(idatlist[[13]],sdfs.normal=idatlist[c(2,10,26)])
segtest_125<-cnSegmentation(idatlist[[14]],sdfs.normal=idatlist[c(3,11,27)])

pdf(paste(filepath,"plots/LNCaP_PrEC_CNVtest.pdf",sep==""))
visualizeSegments(segtest_500)
visualizeSegments(segtest_250)
visualizeSegments(segtest_125)
dev.off()

#agrees with https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5112954/ Fig 4
