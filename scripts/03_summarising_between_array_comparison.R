#load libraries
library(gplots)
library(ggVennDiagram)
library(ggplot2)
library(plyr)
library(beeswarm)

#define file path
filepath<-".../your/file/path/"

#Load in new manifest
load(paste(filepath,"results/NewManifest.RData",sep=""))
#loads object: "new_manifest"

###########Comparing probe content from previous versions############

#For EPICv2 probes how many were on older platforms.
#Using columns
length(which(new_manifest$vecEPICv1probeID!=""))#726597
length(which(new_manifest$epic1seqmatch!=""))#727520
length(which(new_manifest$epic1locmatch!=""))#727605

length(which(new_manifest$vec450probeID!=""))#399178
length(which(new_manifest$K450seqmatch!=""))#399717
length(which(new_manifest$K450locmatch!=""))#399758

length(which(new_manifest$vec27probeID!=""))# 24530
length(which(new_manifest$K27seqmatch!=""))#2178
length(which(new_manifest$K27locmatch!=""))#24571

#unique # Minus 1 in each instance for NAs
length(unique(new_manifest$vecEPICv1probeID))-1#721802
length(unique(new_manifest$epic1seqmatch))-1#722794
length(unique(new_manifest$epic1locmatch))-1#722815
length(unique(new_manifest$EPICv1_Loci))-1#722444

length(unique(new_manifest$vec450probeID))-1#394442
length(unique(new_manifest$K450seqmatch))-1#395052
length(unique(new_manifest$K450locmatch))-1#395037
length(unique(new_manifest$Methyl450_Loci))-1#391241

length(unique(new_manifest$vec27probeID))-1# 24348
length(unique(new_manifest$K27seqmatch))-1#2173
length(unique(new_manifest$K27locmatch))-1#24389
length(unique(new_manifest$Methyl27_Loci))-1#24311

#Read in other manifests
#Read in other Sesame manifests (all hg38)
manifest_EPICv1<-read.table(paste(filepath,"metadata/EPIC.hg38.manifest.tsv",sep=""),sep="\t",header=T)
manifest_450K<-read.table(paste(filepath,,"metadata/HM450.hg38.manifest.tsv",sep=""),sep="\t",header=T)
manifest_27K<-read.table(paste(filepath,,"metadata/HM27.hg38.manifest.tsv",sep=""),sep="\t",header=T)

##Probe name venn diagram
#Make vector of unique probe names across all 4 arrays
fullnames<-unique(c(new_manifest$Name,manifest_EPICv1$Probe_ID,manifest_450K$Probe_ID,manifest_27K$probeID))
length(fullnames)
#[1] 1084471

#make list of unique names in each manifest
fullnameslist<-list(unique(new_manifest$Name),unique(manifest_EPICv1$Probe_ID),unique(manifest_450K$Probe_ID),unique(manifest_27K$probeID))

#extract rs probes in list to make into venn diagram of rs probe overlaps
fullnameslist_rs<-list()
for(i in 1:length(fullnameslist)){
  fullnameslist_rs[[i]]<-(fullnameslist[[i]][grep("rs",fullnameslist[[i]])])
}

#export venn diagram of rs probe overlaps
pdf(paste(filepath,"plots/Venn_SNP_versioncomparison.pdf",sep=""))
ggVennDiagram(fullnameslist_rs,label="count",color="grey", category.names=c("EPICv2","EPICv1" ,'450K', '27K')) + scale_fill_gradient(low="white",high = "red") + scale_color_manual(values = c("EPICv2" = "grey","EPICv1" ="grey",'450K' = 'grey', '27K' = 'grey')) 
dev.off()

lapply(fullnameslist_rs,length)
#[[1]] [1] 62
#[[2]] [1] 59
#[[3]] [1] 65
#[[4]] [1] 0

#identify rs duplicates
new_manifest$Name[grep("rs",(new_manifest$Name))][duplicated(new_manifest$Name[grep("rs",(new_manifest$Name))])]
#[1] "rs11249206" "rs133860"   "rs966367" 

#extract ch probes in list to make into venn diagram of ch probe overlaps
fullnameslist_ch<-list()
for(i in 1:length(fullnameslist)){
  fullnameslist_ch[[i]]<-(fullnameslist[[i]][grep("ch",fullnameslist[[i]])])
}

#export venn diagram of ch probe overlaps
pdf(paste(filepath,"plots/Venn_CpH_versioncomparison.pdf",sep=""))
ggVennDiagram(fullnameslist_ch,label="count",color="grey",category.names=c("EPICv2","EPICv1" ,'450K', '27K')) + scale_fill_gradient(low="white",high = "red") + scale_color_manual(values = c("EPICv2" = "grey","EPICv1" ="grey",'450K' = 'grey', '27K' = 'grey'))
dev.off()

#checking that no replicate locations in the ch probes?
ch_man<-new_manifest[which(new_manifest$Probe_Type=="ch"),]
table(ch_man$posrep)
#        Y 
#2906    8 
ch_man$chr0vec[which(ch_man$posrep=="Y")]
#[1] "Y" "Y" "Y" "Y" "Y" "Y" "Y" "Y"
#Not actually replicate locations, just in unmapped probe section

length(unique(ch_man$AlleleA_ProbeSeq))
#2914

length(unique(ch_man$Name))
#2914
#all ch have unique probe sequences and names, so consider unique for venn

#Extract unique genomic locations for each array (excluding rs probes) #Remove replicate locations from EPICv2
Pos27<-unique(paste(manifest_27K$CpG_chrm,manifest_27K$CpG_beg,manifest_27K$CpG_end))
rs450<-grep("rs",manifest_450K$Probe_ID)
Pos450<-unique(paste(manifest_450K$CpG_chrm[-rs450],manifest_450K$CpG_beg[-rs450],manifest_450K$CpG_end[-rs450]))
rsEPIC<-grep("rs",manifest_EPICv1$Probe_ID)
PosEPICv1<-unique(paste(manifest_EPICv1$CpG_chrm[-rsEPIC],manifest_EPICv1$CpG_beg[-rsEPIC],manifest_EPICv1$CpG_end[-rsEPIC]))
rsEPIC2<-grep("rs",new_manifest$Name)
PosEPICv2<-unique(paste(new_manifest$CpG_chrm[-rsEPIC2],new_manifest$CpG_beg[-rsEPIC2],new_manifest$CpG_end[-rsEPIC2]))

length(Pos27)#27566
length(Pos450)#485474
length(PosEPICv1)#865806
length(PosEPICv2)#929646

#Make vector and matrix of locations across all arrays
fullnames_loc<-unique(c(Pos27,Pos450,PosEPICv1,PosEPICv2))
length(fullnames_loc)
#[1] 1082977
locvenntab<-matrix(data=NA,ncol=4,nrow=length(fullnames_loc))
m1<-match(fullnames_loc,PosEPICv2)
locvenntab[which(is.na(m1)==TRUE),1]<-FALSE
locvenntab[which(is.na(m1)!=TRUE),1]<-TRUE
m2<-match(fullnames_loc,PosEPICv1)
locvenntab[which(is.na(m2)==TRUE),2]<-FALSE
locvenntab[which(is.na(m2)!=TRUE),2]<-TRUE
m3<-match(fullnames_loc,Pos450)
locvenntab[which(is.na(m3)==TRUE),3]<-FALSE
locvenntab[which(is.na(m3)!=TRUE),3]<-TRUE
m4<-match(fullnames_loc,Pos27)
locvenntab[which(is.na(m4)==TRUE),4]<-FALSE
locvenntab[which(is.na(m4)!=TRUE),4]<-TRUE

colnames(locvenntab)<-c("EPICv2","EPICv1","450K","27K")

fullloclist<-list(PosEPICv2,PosEPICv1,Pos450,Pos27)

#export venn diagram showing overlap in position between arrays
pdf(paste(filepath,"plots/Venn_Pos_versioncomparison.pdf",sep=""))
ggVennDiagram(fullloclist,label="count",color="grey",category.names=c("EPICv2","EPICv1" ,'450K', '27K')) + scale_fill_gradient(low="white",high = "red") + scale_color_manual(values = c("EPICv2" = "grey","EPICv1" ="grey",'450K' = 'grey', '27K' = 'grey'))
dev.off()

#define and count number locations targeted by EPICv2 that never targeted by any of the previous arrays 
new<-which(locvenntab[,1]==TRUE & locvenntab[,2]==FALSE & locvenntab[,3]==FALSE & locvenntab[,4]==FALSE )
length(new)#182139
length(new)/length(PosEPICv2)#0.196

#number locations targeted by EPICv2 that were also targeted by previous arrays 
length(PosEPICv2)-length(new)##747507
(length(PosEPICv2)-length(new))/length(PosEPICv2)#0.804

#define and count number locations targeted by EPICv2 that also targeted by EPICv1
retained<-which(locvenntab[,1]==TRUE & locvenntab[,2]==TRUE )
length(retained)#722758

#define and count number locations targeted by EPICv2 that not targeted by EPICv1, but were targeted by 450K and/or 27K
reinstated<-which(locvenntab[,1]==TRUE & locvenntab[,2]==FALSE & (locvenntab[,3]==TRUE | locvenntab[,4]==TRUE ))
length(reinstated)#24749

#define and count number locations targeted by EPICv1 tbut removed for EPICv2
excluded<-which(locvenntab[,1]==FALSE & locvenntab[,2]==TRUE)
length(excluded)#143048

allepic2<-which(locvenntab[,1]==TRUE)
length(allepic2)#929646

#Export each as bed files for later analysis 
a<-t(sapply(fullnames_loc,function(x){
  strsplit(x," ")[[1]]}))

new_bed<-a[new,]
write.table(new_bed,file=paste(filepath,"results/Locations_new.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)

retained_bed<-a[retained,]
write.table(retained_bed,file=paste(filepath,"results/Locations_retained.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)

reinstated_bed<-a[reinstated,]
write.table(reinstated_bed,file=paste(filepath,"results/Locations_reinstated.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)

excluded_bed<-a[excluded,]
write.table(excluded_bed,file=paste(filepath,"results/Locations_excluded.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)

all_bed<-a[allepic2,]
write.table(all_bed,file=paste(filepath,"results/Locations_epic2.bed",sep=""), quote=F, sep="\t", row.names=F, col.names=F)


####Exploring overlap between probes excluded/retained between EPICv1 and EPICv2 and Zhou et al. masked list
#read in masked list, downloaded from https://zwdzwd.github.io/InfiniumAnnotation 
maskEPIC<-read.table(paste(filepath,"metadata/EPIC.hg38.mask_sesame.txt",sep=""),sep="\t",header=T)

#get EPICv1 probe IDs for excluded probes
exc<-fullnames_loc[excluded]
E1_coord<-paste(manifest_EPICv1$CpG_chrm,manifest_EPICv1$CpG_beg,manifest_EPICv1$CpG_end)
E1_IDs<-manifest_EPICv1$Probe_ID[match(exc,E1_coord)]

#Calculate percentage of excluded probes overlapping EPICv1 masked probes 
maskEPIC_excluded<-maskEPIC[match(E1_IDs,maskEPIC$probeID),]
excluded_maskgeneral<-table(maskEPIC_excluded[,12])
(excluded_maskgeneral[2]/nrow(maskEPIC_excluded))*100 #72.81053 

#get EPICv1 probe IDs for retained probes
ret<-fullnames_loc[retained]
E1_IDsR<-manifest_EPICv1$Probe_ID[match(ret,E1_coord)]

#Calculate percentage of retained probes overlapping EPICv1 masked probes 
maskEPIC_retained<-maskEPIC[match(E1_IDsR,maskEPIC$probeID),]
retained_maskgeneral<-table(maskEPIC_retained[,12])
(retained_maskgeneral[2]/nrow(maskEPIC_retained))*100 #0.173087


####Exploring chromosomal distribution of probes
##Calculate number of locations per chromosome for each subset of the probes
#extract chromosomal locations for each subset of the data
new_chr<-table(new_bed[,1])
retained_chr<-table(retained_bed[,1])
reinstated_chr<-table(reinstated_bed[,1])
excluded_chr<-table(excluded_bed[,1])
all_chr<-table(all_bed[,1])

#combine as a list
list_chr<-list(new_chr,retained_chr,reinstated_chr,excluded_chr,all_chr)

#does each vector contain a different set of chromosomes
unlist(lapply(list_chr,length))
#[1] 26 25 24 36 27
#yes

#Combine in 1 table
allNames<-lapply(list_chr,function(y){names(y)})
allNames_vec<-unique(unlist(allNames))
tab_chr<-matrix(data=NA,nrow=length(allNames_vec),ncol=length(list_chr))
for(i in 1:length(list_chr)){
  for(j in 1:length(allNames_vec)){
  tab_chr[j,i]<-list_chr[[i]][match(allNames_vec[j],names(list_chr[[i]]))]
  }
}
rownames(tab_chr)<-allNames_vec
colnames(tab_chr)<-c("new","retained","reinstated","excluded","all")
#reorder tab_chr by number
tab_chr2<-tab_chr[c(1,12,17:23,2:11,13:15,25,26,24,36,28:35,16),]


#Plot absolute values for each chromosome for all EPICv2 probes
pdf(paste(filepath,"plots/Barplot_ChromsomeDistribution_number.pdf",sep=""),width=10)
barplot(t(tab_chr2[1:25,5]),legend.text=colnames(tab_chr2)[5],las=2, ylab="Number of sites")
dev.off()

#Plot sites per million base pairs for each chromosome for all EPICv2 probes
#subset all bed by chromosome
na2rm<-which(is.na(all_bed[,1])==TRUE)

all_GR<-GRanges(seqnames=all_bed[-na2rm,1], IRanges(start=all_bed[-na2rm,2],end=all_bed[-na2rm,3] ))

lev<-levels(seqnames(all_GR))

all_GR_list<-list()
for(i in 1:(length(lev)-1)){
  all_GR_list[[i]]<-all_GR[which(seqnames(all_GR)==lev[i])]
}

#split each chromsome into bins
tablelist<-list()

for(i in 1:length(all_GR_list)){
  min<-sort(all_GR_list[[i]])[1,]
  max<-sort(all_GR_list[[i]])[length(all_GR_list[[i]]),]
  
  fakeGR<-GRanges(seqnames=seqnames(min),IRanges(start=min@ranges@start,end=(max@ranges@start+max@ranges@width-1)))
  tiled<-tile(fakeGR, width=1000000)
  
  #estimate overlap per bin
  fo<-findOverlaps(all_GR_list[[i]],tiled[[1]])
  tablelist[[i]]<-(table(fo@to))
}


##scatterplot
names(tablelist)<-lev[1:25]

pdf(paste(filepath,"plots/Beeswarm_ChromsomeDistribution.pdf",sep=""),width=10)
beeswarm(tablelist,cex=0.3,las=2,col="dark grey",ylab="Number of CpG sites per million base pairs")
bxplot(tablelist, add = TRUE , lwd=2,col=c("gold","indian red","gold"))
dev.off()

#Percentage
totals<-colSums(tab_chr2,na.rm=T)
tab_chr3<-list()
for(i in 1:ncol(tab_chr2)){
  tab_chr3[[i]]<-(tab_chr2[,i]/totals[i])*100
}
tab_chr4<-do.call(cbind,tab_chr3)
colnames(tab_chr4)<-colnames(tab_chr2)

#Plot percentage values for each chromosome for each probe set
pdf(paste(filepath,"plots/Barplot_ChromsomeDistribution_percent.pdf",sep=""),width=10)
barplot(t(tab_chr4[1:25,]),beside=T,legend.text=colnames(tab_chr4),las=2, ylab="Percent of sites (%)")
dev.off()

#Summarise all data in tables
 write.csv(tab_chr2,paste(filepath,"results/ChromsomeDistribution_number.csv",sep=""))
 write.csv(tab_chr4,paste(filepath,"results/ChromsomeDistribution_percent.csv",sep=""))
