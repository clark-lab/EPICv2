#Script matching results section: "Categorisation of EPICv2 probe types" and "Replicate probes on EPICv2"

#load in packages
library(GenomicRanges)

#define file path
filepath<-".../your/file/path/"

#Read in Illumina EPICv2 manifest 
manifest <- read.csv(paste(filepath,"metadata/EPIC-8v2-0_A1.csv",sep=""), skip = 7)

###Summarising numbers of different probe types
#How many probes are there in the manifest
length(which(manifest$Name!=""))#937690

#How many of each probe type are there?
typeprobe<-manifest$Probe_Type
typeprobe
#cg     ch     nv     rs 
#636 933252   2914    824     65

#How many of the rs probes are unique? 
length(unique(manifest$Name[which(typeprobe=="rs")]))
#62

table(manifest$Rep_Num[which(typeprobe=="rs")])
#1  2 
#62  3 

#How many of the nv probes target unique loci?
length(unique(manifest$MAPINFO[which(typeprobe=="nv")]))#462

#How many probes are there of each design type?
table(manifest$Infinium_Design_Type)
#     I     II 
#636 128295 808760

twotypes<-table(manifest$Infinium_Design_Type,manifest$Probe_Type)

#write out summary of design type vs probe type numbers
write.csv(twotypes,paste(filepath,"results/TwoTypes_summary.csv",sep=""))

#Check that probes without design type are control probes
temp<-which(manifest$Probe_Type=="")
ctrltypes<-table(manifest$Name[temp])  
#   BISULFITE CONVERSION I BISULFITE CONVERSION II               EXTENSION 
#1                      10                       4                       4 
#HYBRIDIZATION                NEGATIVE         NON-POLYMORPHIC                  NORM_A 
#3                     411                       9                      27 
#NORM_C                  NORM_G                  NORM_T             RESTORATION 
#58                      27                      58                       1 
#SPECIFICITY I          SPECIFICITY II                STAINING          TARGET REMOVAL 
#12                       3                       6                       2 

#write out summary of control probes
write.csv(ctrltypes,paste(filepath,"results/CtrlTypes_summary.csv",sep=""))

#What is the probe with no name?
manifest[temp[which(manifest$Name[temp]=="")],]
#note 937056 is just a header row - not a probe after all

###Comparing nv probes to COSMIC database

#select out nv probes
nvtab<-manifest[which(typeprobe=="nv"),]

# read in cosmic data 
# Downloaded on 10th March 2023 from https://cancer.sanger.ac.uk/census/
census<-read.csv(paste(filepath,"metadata/COSMIC_Census_allFri Mar 10 00_52_33 2023.csv",sep=""))

#Create genomic ranges objects for nv probes and COSMIC census data
#nv probes GRanges
nvGR<-GRanges(seqnames=nvtab$CHR,IRanges(start=nvtab$MAPINFO,end=nvtab$MAPINFO))
length(nvGR)#824
nvGR2<-unique(nvGR)
length(nvGR2)#462

#This is incomplete so instead use coordinates in names
chr<-start<-end<-vector()

for(i in 1:nrow(nvtab)){
  chr[i]<-strsplit(nvtab$Name[i],"-")[[1]][3]
  start[i]<-strsplit(nvtab$Name[i],"-")[[1]][4]
  end[i]<-strsplit(nvtab$Name[i],"-")[[1]][5]
}
nvGR3<-GRanges(seqnames=chr,IRanges(start=as.numeric(start),end=as.numeric(end)))
length(nvGR3)#824
length(unique(nvGR3))#474

#census GRanges
#remove census data without mapping information
rm<-which(is.na(census$Start)==TRUE)
censusrm<-census[-rm,]

#format cooridnates in census data
chrcensus<-paste("chr",censusrm$Chr,sep="")

censusGR<-GRanges(seqnames=chrcensus,IRanges(start=censusrm$Start,end=censusrm$End))

#find overlaps between datasets
fo<-findOverlaps(nvGR3,censusGR)
length(unique(fo@from))#821

nvGR4<-unique(nvGR3)
fo4<-findOverlaps(nvGR4,censusGR)
length(fo4)#472

#create comparison table between nv probes and census data
censusanno<-matrix(data=NA,ncol=ncol(censusrm),nrow=nrow(nvtab))
censuslist<-list()
for(i in 1:length(fo)){
  censuslist[[fo@from[i]]]<-censusrm[fo@to[i],]
}
for(i in 1:length(censuslist)){
    if(length(censuslist[[i]])>1){
  censusanno[i,]<-as.character(censuslist[[i]])
  }
}
colnames(censusanno)<-colnames(censusrm)
nvtab_anno<-cbind(nvtab,censusanno)

#write out table of how nv probes overlap with cosmic data
write.csv(nvtab_anno,paste(filepath,"results/nv_COSMICanno.csv",sep=""))


###Understanding replicates
##Finding 'exact replicates'
#Identify Illumina defined replicates
repnum2<-which(manifest$Rep_Num==2)
length(repnum2)#5141
table(manifest$Rep_Num)
#1      2      3      4      5      6      7      8      9     10 
#929958   5141   1002     70     52      4      1      1      1      1 

#how many unique probe names among Illumina defined replicates?
repnum_all<-which(manifest$Rep_Num>1)
length(unique(manifest$Name[repnum_all]))#5141
#extract unique probe names among replicates?
uniqnamereps<-unique(manifest$Name[repnum_all])

#find all possible IlmnIDs of Illumina defined replicates
#extract from IlmnID the probe name and first 3 letters only of suffix
suffix3<-sapply(manifest$IlmnID,function(x){substr(x,1,nchar(x)-1)})
#extract the truncated IlmnID corresponding to Illumina defined replicates
temp2<-suffix3[which(manifest$Rep_Num==2)]
length(temp2)#5141
length(unique(temp2))#5141

#extract full IlmnIDs for probes for Illumina defined replicates
lista<-list()
for(i in 1:length(temp2)){
  lista[[i]]<-manifest$IlmnID[grep(temp2[i],suffix3)]
}

#summary table of number of probe replicates per probe
table(unlist(lapply(lista,length)))
#2    3    4    5    6   10 
#4139  932   18   48    3    1

sum(table(unlist(lapply(lista,length))))#5141
length(unlist(lista))#11414

#check that all probes are indeed exact matches (in terms of probe sequence)
seqmatch_tf<-vector()
for(i in 1:length(lista)){
    m<-match(lista[[i]],manifest$IlmnID)
    temp<-paste(manifest$AlleleA_ProbeSeq[m],manifest$AlleleB_ProbeSeq[m])
    seqmatch_tf[i]<-length(unique(temp)) == 1  
}
length(which(seqmatch_tf=="FALSE"))#4
lista[which(seqmatch_tf=="FALSE")]

#summary table of number of probe replicates per probe with non-exact replicates excluded
table(unlist(lapply(lista[-which(seqmatch_tf=="FALSE")],length)))
#  2    3    4    5    6   10 
#4135  932   18   48    3    1 
sum(table(unlist(lapply(lista[-which(seqmatch_tf=="FALSE")],length))))
#[1] 5137

#export summary table of number of probe replicates per 'exact replicate' probe
erep_summary<-table(unlist(lapply(lista[-which(seqmatch_tf=="FALSE")],length)))
write.csv(erep_summary,paste(filepath,"results/ExactRep_summary.csv",sep=""))

#export vector of IlmnIDs of true 'exact replicates'
exact_replicates<-unlist(lista[-which(seqmatch_tf=="FALSE")])
length(exact_replicates)#11406
write.csv(exact_replicates,paste(filepath,"results/Exact_replicate_vector.csv",sep=""))

lista_2<-lapply(lista[-which(seqmatch_tf=="FALSE")],function(x){paste(x,collapse=";")})

write.csv(lista_2,paste(filepath,"results/Exact_replicate_vector_collapsed.csv",sep=""))


##Finding location replicates
#Identify those with same name, but exclude those that are exact replicate IlmnID (i.e. same probe sequences)
length(which(duplicated(manifest$Name)==TRUE))#7017

#create table of only duplicated probes
manifest_dups<-manifest[which(duplicated(manifest$Name)==TRUE),]
#remove those already identified as replicates by Illumina
manifest_dups2<-manifest_dups[which(manifest_dups$Rep_Num==1),]
nrow(manifest_dups2)#124

#count how many uniquely duplicated probes are there(location, but not sequence) and extract their probe names
length(unique(manifest_dups2$Name))#113
unames<-unique(manifest_dups2$Name)

#extract IDs of location probes and check that their sequences really are distinct
listb<-list()
for(i in 1:length(unames)){
  listb[[i]]<-manifest$IlmnID[grep(unames[i],manifest$Name)]
}

#count number of location replicates per probe
table(unlist(lapply(listb,length)))
#2  3  4 
#77 31  5 
sum(table(unlist(lapply(listb,length))))#113

#are there any location replicates that also have exact replicates?
length(which(manifest$Rep_Num[match(unlist(listb),manifest$IlmnID)]==2))#29 - yes

#extract IlmnIDs of all location replicate probes, grouped by probe name
listc<-unlist(listb)[which(manifest$Rep_Num[match(unlist(listb),manifest$IlmnID)]==1)]
length(listc)#237

#Check that probe sequences truly are unique
manifest_loc<-manifest[match(listc,manifest$IlmnID),]
length(paste(manifest_loc$AlleleA_ProbeSeq,manifest_loc$AlleleB_ProbeSeq))#237
length(unique(paste(manifest_loc$AlleleA_ProbeSeq,manifest_loc$AlleleB_ProbeSeq)))#236 - one is not!

#which location replicate has sequence duplicates?
manifest_loc$IlmnID[which(duplicated(paste(manifest_loc$AlleleA_ProbeSeq,manifest_loc$AlleleB_ProbeSeq))==TRUE)]#"cg10335333_BC21"

unique(manifest$AlleleA_ProbeSeq[grep("cg10335333",manifest$Name)])
#3 probe sequences all identical. Not truly a position replicate

#remove from this list
nonlocrep<-manifest$IlmnID[grep("cg10335333",manifest$Name)]
listd<-listc[-na.omit(match(nonlocrep,listc))]
length(listd)#235

#write out list of IlmnID location replicates
write.csv(listd,paste(filepath,"results/Location_replicate_vector.csv",sep=""))

#write out list of IlmnID location replicates, concatenated by probe location
listb_2<-lapply(listb[-113],function(x){paste(x,collapse=";")})

write.csv(listb_2,paste(filepath,"results/Location_replicate_vector_collapsed.csv",sep=""))

length(unique(manifest$Name[match(listd,manifest$IlmnID)]))#112

max(table(manifest$Name[match(listd,manifest$IlmnID)]))#4
min(table(manifest$Name[match(listd,manifest$IlmnID)]))#2


#export summary table of number of probe replicates per 'location replicate' probe
lrep<-table(table(manifest$Name[match(listd,manifest$IlmnID)]))
write.csv(lrep,paste(filepath,"results/LocationRep_summary.csv",sep=""))

#Identify probes that have both location and exact replicates - exclude cg10335333 as doesn't quite fit either criteria
length(which(manifest$Rep_Num[match(unlist(listb),manifest$IlmnID)]==2))#29
doublereps<-unlist(listb)[which(manifest$Rep_Num[match(unlist(listb),manifest$IlmnID)]==2)]#includes "cg10335333_TC22", so remove from this count
doublereps2<-doublereps[-match("cg10335333_TC22",doublereps)]
length(doublereps2)#28
#export table of probes with both location and exact replicates
write.csv(doublereps2,paste(filepath,"results/Double_replicates.csv",sep=""))

##Finding sequence-only replicates
#create new vector of combined alleleA and alleleB probe sequences
alleleall<-paste(manifest$AlleleA_ProbeSeq,manifest$AlleleB_ProbeSeq)

#determine how mnay probe sequences are duplicated between probes
dups<-which(duplicated(alleleall)==TRUE)
length(dups)#6836

#create table of sequence duplicated probes only
nm_seq_dups<-manifest[dups,]

#exclude probes that are already known to be 'exact replicates'
sd_dups_b<-setdiff(nm_seq_dups$IlmnID,exact_replicates)
length(sd_dups_b)#567

#check if remainder overlap with location
any<-intersect(sd_dups_b,listd)
length(any)#0 - no

#check if any overlap by name (a feature of exact replicates)
nm_seq_dups2<-nm_seq_dups[match(sd_dups_b,nm_seq_dups$IlmnID),]
length(unique(nm_seq_dups2$IlmnID))#567
length(unique(nm_seq_dups2$Name))#567 - all have unique names

#How many unique sequences in sequence duplicates?
length(unique(alleleall[match(sd_dups_b,manifest$IlmnID)]))#435

#create vector of unique sequences
allall<-unique(alleleall[match(sd_dups_b,manifest$IlmnID)])

#create list of IlmnIDs for each probe sequence
allalllist<-list()
for(i in 1:length(allall)){
  allalllist[[i]]<-manifest$IlmnID[grep(allall[i],alleleall)]
}

#summary table of how many probe replicates per sequence
table(unlist(lapply(allalllist,length)))
# 2   3   4   5   6   7   8  10  11 
#372  37  10   4   7   1   1   2   1 

#create vector of IlmnIDs for sequence-only replicates
seqonlylist<-unlist(allalllist)

#what chromsome are sequence-only replicates on
table(manifest$CHR[match(seqonlylist,manifest$IlmnID)])
#chr0 
#1003 

#check that names don't match within each set of sequence matched probes
namematch_tf<-vector()
for(i in 1:length(allalllist)){
  temp<-manifest$Name[match(allalllist[[i]],manifest$IlmnID)]
  namematch_tf[i]<-length(unique(temp)) < length(temp)  
}

#look at probe where names are matched
rm33<-which(namematch_tf==TRUE)#431
allalllist[rm33]
#[[1]]
#[1] "cg10335333_TC21" "cg10335333_TC22" "cg10335333_BC21"

#This is the discrepant probe from location replicate. Exclude from sequence-only here as doesn't meet criteria of different probe names

#remove name matched probe and export IlmnIDs for sequence-only replicates
allalllist2<-allalllist[-rm33]
seqonlylist<-unlist(allalllist2)
length(seqonlylist)#1000

write.csv(seqonlylist,paste(filepath,"results/Sequenceonly_replicate_vector.csv",sep=""))

allalllist3<-lapply(allalllist2,function(x){paste(x,collapse=";")})
allalllist4<-unlist(allalllist3)
length(allalllist4)#434
write.csv(allalllist4,paste(filepath,"results/Sequenceonly_replicate_vector_collapsed.csv",sep=""))


#export summary table of number of probe replicates per 'sequence-only replicate' probe
srep<-table(unlist(lapply(allalllist2,length)))
write.csv(srep,paste(filepath,"results/SequenceRep_summary.csv",sep=""))


#Summarising 'other' replicate probes, that don't fit the criteria for exact, location or sequence
manifest$IlmnID[grep("cg12495610",manifest$Name)]#"cg12495610_BC21" "cg12495610_BC22"
manifest$IlmnID[grep("cg22823189",manifest$Name)]#"cg22823189_TC21" "cg22823189_TC22"
manifest$IlmnID[grep("cg13305632",manifest$Name)]#"cg13305632_BC21" "cg13305632_BC22"
manifest$IlmnID[grep("cg18801637",manifest$Name)]#"cg18801637_BC21" "cg18801637_BC22"
                     , , 


###Summarising manifest columns for overlap with probes on previous arrays

#Overlap with HM27K
length(which(manifest$Methyl27_Loci!=""))#24490

#Overlap with HM450K
length(which(manifest$Methyl450_Loci!=""))#395936

#Overlap with EPICv1
length(which(manifest$EPICv1_Loci!=""))#727222

#check probe IDs to show changes
E1name<-manifest$EPICv1_Loci[which(manifest$EPICv1_Loci!="")]
E2name<-manifest$Name[which(manifest$EPICv1_Loci!="")]

length(E1name)#727222
length(E2name)#727222
length(intersect(E1name,E2name))#721431
#shows that probe name changes between platforms


