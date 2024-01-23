#load libraries 
library(GenomicRanges)

#define file path
filepath<-".../your/file/path/"

#Read in Illumina EPICv2 manifest 
manifest_EPICv2_Ilmn <- read.csv(paste(filepath,"metadata/EPIC-8v2-0_A1.csv",sep=""), skip = 7)
nrow(manifest_EPICv2_Ilmn)#937691

#Read in Sesame EPICv2 manifest (downloaded from https://zwdzwd.github.io/InfiniumAnnotation)
manifest_EPICv2<-read.table(paste(filepath,"metadata/EPICv2.hg38.manifest.tsv",sep=""),sep="\t",header=T)
nrow(manifest_EPICv2)#936866

#Read in Sesame manifests for older versions of the array (all hg38) (downloaded from https://zwdzwd.github.io/InfiniumAnnotation)
manifest_EPICv1<-read.table(paste(filepath,"metadata/EPIC.hg38.manifest.tsv",sep=""),sep="\t",header=T)
manifest_450K<-read.table(paste(filepath,"metadata/HM450.hg38.manifest.tsv",sep=""),sep="\t",header=T)
manifest_27K<-read.table(paste(filepath,"metadata/HM27.hg38.manifest.tsv",sep=""),sep="\t",header=T)

###################################
#Comparison between Illumina and Sesame EPICv2 manifest 

#check Illumina and Sesame manifest contain same probes based on probe ID - match the two EPICv2 manifests based on Probe_ID 
match2<-match(manifest_EPICv2_Ilmn$IlmnID,manifest_EPICv2$Probe_ID)
length(which(is.na(match2)==TRUE))#1460

#1460 probes don't match between EPICv2 manifests - what are the probe types of these probes, and how many of each probe type 
table(manifest_EPICv2_Ilmn$Probe_Type[which(is.na(match2)==TRUE)])
#     nv 
#636 824 

#Discrepant 636 are the control probes

length(grep("nv",manifest_EPICv2$Probe_ID))#0
#The nv probes have been excluded from sesame manifest

#Create version of Illumina manifest_EPICv2_Ilmn with non-matching probes removed:
manifest_EPICv2_IlmnTrim<-manifest_EPICv2_Ilmn[-which(is.na(match2)==TRUE),]
nrow(manifest_EPICv2_IlmnTrim)# 936231

#Illumina manifest_EPICv2_IlmnTrim now has less probes than the Sesame manifest, why is this?
m2<-match(manifest_EPICv2$Probe_ID,manifest_EPICv2_Ilmn$IlmnID)
head(manifest_EPICv2$Probe_ID[(which(is.na(m2)==TRUE))])
#[1] "ctl_10609447_NEGATIVE" "ctl_10627500_NEGATIVE" "ctl_10676356_NEGATIVE" "ctl_10714330_NORM_T"   "ctl_10731326_NEGATIVE" "ctl_10732387_NEGATIVE"

#Sesame manifest has 'ctl' added to control probe probe IDs. So the control probes ARE in the Sesame manifest, but with a different name from Illumina manifest
#They are not needed for our data analysis, so we will remove control probe from both versions of the manifest for now

#match the two EPICv2 manifest based on Probe_ID (check completely unique)
match2<-match(manifest_EPICv2_IlmnTrim$IlmnID,manifest_EPICv2$Probe_ID)
length(which(is.na(match2)==TRUE))#0
#All Illumina manifest probes are found in sesame manifest

#Reorder Sesame manifest to match probes and probe order in the trimmed Illumina Manifest 
manifest_EPICv2_ord<-manifest_EPICv2[match2,]

#Check re-ordering has been successful
identical(manifest_EPICv2_IlmnTrim$IlmnID,manifest_EPICv2_ord$Probe_ID)#TRUE

#Confirm probe sequences are also identical between EPICv2 manifests
identical(manifest_EPICv2_IlmnTrim$AlleleA_ProbeSeq,manifest_EPICv2_ord$AlleleA_ProbeSeq)#TRUE
identical(manifest_EPICv2_IlmnTrim$AlleleB_ProbeSeq,manifest_EPICv2_ord$AlleleB_ProbeSeq)#FALSE
identical(manifest_EPICv2_IlmnTrim$AlleleB_ProbeSeq[-which(manifest_EPICv2_IlmnTrim$AlleleB_ProbeSeq=="")],manifest_EPICv2_ord$AlleleB_ProbeSeq[-which(is.na(manifest_EPICv2_ord$AlleleB_ProbeSeq)==TRUE)])#TRUE
#Sequences are identical between manifests, just didn't seem like it at first as Sesame sets Type II probes AlleleB_ProbeSeq to NA, whereas Illumina uses ""

#Determine if genomic target for each probe is identical between Illumina and Sesame manifest
#Sesame uses a different mapping system giving a range for the target CpG site, rather than a single position, so will convert coordinates from both manifests to GenomicRanges objects for comparison

#Create GenomicRanges object from Sesame manifest
manifest_EPICv2_ordGR<-GRanges(seqnames=manifest_EPICv2_ord$CpG_chrm,IRanges(start=manifest_EPICv2_ord$CpG_beg,end=manifest_EPICv2_ord$CpG_end))
#Error in .normarg_seqnames1(seqnames) : 'seqnames' cannot contain NAs

#Receive error message
#How many of the probes have an NA for chromosome?
length(which(is.na(manifest_EPICv2_ord$CpG_chrm)==TRUE))#1
manifest_EPICv2_ord$Probe_ID[which(is.na(manifest_EPICv2_ord$CpG_chrm)==TRUE)]#"cg03136958_TC21"
which(is.na(manifest_EPICv2_ord$CpG_chrm)==TRUE)# 168124

#Does this probe have mapping info in the Illumina file?
manifest_EPICv2_IlmnTrim$CHR[grep("cg03136958_TC21",manifest_EPICv2_IlmnTrim$IlmnID)]#chr2
#Yes

#Will remove this probe from the Sesame GR object for now
manifest_EPICv2_ordGR<-GRanges(seqnames=manifest_EPICv2_ord$CpG_chrm[-168124],IRanges(start=manifest_EPICv2_ord$CpG_beg[-168124],end=manifest_EPICv2_ord$CpG_end[-168124]))

#Create GenomicRanges object from Illumina manifest
manifest_EPICv2_IlTrGR<-GRanges(seqnames=manifest_EPICv2_IlmnTrim$CHR,IRanges(start=manifest_EPICv2_IlmnTrim$MAPINFO,end=manifest_EPICv2_IlmnTrim$MAPINFO))

#Now find overlaps between these 2 manifest using GenomicRanges
#Would expect perfect overlap except for probes with missing chromosome mapping data in sesame...
OL<-findOverlaps(manifest_EPICv2_ordGR,manifest_EPICv2_IlTrGR)
length(OL)#944116 
#Overlap length longer than the row number of either manifest, which suggests duplication in target of CpG probes within arrays (replicates)

#Firstly are all probe positions from the Illumina manifest accounted for in the overlap?
length(unique(OL@to))# 928715
length(manifest_EPICv2_IlTrGR)#936231
#No - which are missing?
nonOLIlm<-manifest_EPICv2_IlTrGR[-unique(OL@to),]
length(nonOLIlm)#7516
table(nonOLIlm@seqnames)
#chr19 chr16  chr2  chr1 chr12 chr22  chr9 chr20 chr13 chr10 chr17 chr11  chr5 chr15 chr14  chr8  chr7  chr6  chr4  chr3 chr21 chr18  chrX  chrY  chr0  chrM 
#1     1     4     3     1     0     0     0     0     0     0     1     1     1     2     0     2     4     1     0     0     1   604     0  6889     0 

#The majority of Illumina manifest locations that did not overlap with Sesame are listed as being in 'chr0' - these are 6889 unique probes that are missing genomic position entirely. 
#The next most common chromosome is ChrX
#Do the chr0 probes have mapping information in Sesame?

tmp<-which(manifest_EPICv2_IlmnTrim$CHR=="chr0")
length(which(tmp>900000))#6883 - majority are towards the end of the manifest, which is suspicious, suggests that mapping just didn't complete somehow?
chr0probes<-manifest_EPICv2_IlmnTrim$IlmnID[which(manifest_EPICv2_IlmnTrim$CHR=="chr0")]
table(manifest_EPICv2_ord$CpG_chrm[match(chr0probes,manifest_EPICv2_ord$Probe_ID)])
# chr1                chr10                chr11                chr12                chr13                chr14                chr15 
# 794                  249                  455                  199                  182                  180                  256 
# chr16                chr17                chr18                chr19                 chr2                chr20                chr21 
# 439                  329                   87                  208                  478                  104                  264 
# chr22 chr22_KI270879v1_alt                 chr3                 chr4                 chr5                 chr6                 chr7 
# 193                    1                  249                  315                  232                  240                  460 
# chr8                 chr9                 chrM                 chrX                 chrY 
# 235                  301                    4                  322                  113 
#Yes these have all have mapping information in Sesame, distributed throughout the genome

#For the non "chr0" discrepant locations, where are they lociated according to the Sesame manifest?
discrepantLocations<-manifest_EPICv2_IlmnTrim[-(unique(OL@to)),]
discIDs<-discrepantLocations$IlmnID[-which(discrepantLocations$CHR=="chr0")]
table(manifest_EPICv2_ord$CpG_chrm[match(discIDs,manifest_EPICv2_ord$Probe_ID)])
#chr1 chr11 chr12 chr15 chr16 chr19  chr2 chr20  chr3  chr4  chr6  chr7  chrX  chrY 
#5     1     3     1     1     1     3     1     1     2     2     2     1   602 
#For the non chr0 discrepant CpGs, there is a bias towards Chr X in Illumina and Y in Sesame makes me wonder if there is a difficulty of mapping to XY homology genes

#Also noticed that even within sesame manifest, some mapping locations differ between probe seq A and B
length(which(manifest_EPICv2_ord$mapChrm_A!=manifest_EPICv2_ord$mapChrm_B))#308
mismatchses<-which(manifest_EPICv2_ord$mapChrm_A!=manifest_EPICv2_ord$mapChrm_B)
length(which(is.na(match(discIDs,manifest_EPICv2_ord$Probe_ID[mismatchses]))==FALSE))#28 
#only 28 of the 308 probes with discrepant mapping between probe A and probe B are discrepant between Illumina and Sesame, so not a big factor in the discrepancy

#make vector of probes with discrepant locations between Sesame and Illumina
discID1<-match(discrepantLocations$IlmnID,manifest_EPICv2_IlmnTrim$IlmnID)
vecdisc<-rep("",nrow(manifest_EPICv2_IlmnTrim))
vecdisc[discID1]<-"Y"
table(vecdisc)
#         Y 
#928715   7516

#make vector of probes that are missing chromosome position in the Illumina manifest
vecmiss<-rep("",nrow(manifest_EPICv2_IlmnTrim))
vecmiss[which(manifest_EPICv2_IlmnTrim$CHR=="chr0")]<-"Y"


###################################
#Explore replicate probes within Illumina manifest 

#make vector of probes that have duplicate Names in Illumina manifest
length(which(duplicated(manifest_EPICv2_IlmnTrim$Name)==TRUE))#6397
uniqname_dups<-unique(manifest_EPICv2_IlmnTrim$Name[which(duplicated(manifest_EPICv2_IlmnTrim$Name)==TRUE)])
dups_name_list<-sapply(uniqname_dups,function(x){grep(x,manifest_EPICv2_IlmnTrim$Name)})
dupsname<-rep("",nrow(manifest_EPICv2_IlmnTrim))
dupsname[unlist(dups_name_list)]<-"Y"

#How many duplicated locations in each EPICv2 manifest
length(which(duplicated(OL@from)==TRUE))#15400
length(which(duplicated(OL@to)==TRUE))#15401

#Why would this be? This led us to look at Illumina manifest and see a 'Rep_Num' column - turns out that some CpGs are targeted multiple times on the array
length(which(manifest_EPICv2_IlmnTrim$Rep_Num>1))#6273
table(manifest_EPICv2_IlmnTrim$Rep_Num)
# 1      2      3      4      5      6      7      8      9     10 
#929958   5141   1002     70     52      4      1      1      1      1 

(5141*2) +  (1002*3) + (70*4) + (52*5) + (4*6) + (1*7) + (1*8) + (1*9) + (1*10)
#[1] 13886

#Summing these replicates, they still fall short of the total number of duplicated locations (see below)
#There must also be other probes that target same CpG but without recognition in the Rep_Num position - likely due to not being exact replicates

#Start with matching on Illumina location
pospos<-paste(manifest_EPICv2_IlmnTrim$CHR,manifest_EPICv2_IlmnTrim$MAPINFO)
loc_dup<-which(duplicated(paste(pospos))==TRUE)
length(loc_dup)#[1] 13230
length(which(manifest_EPICv2_IlmnTrim$CHR[loc_dup]=="chr0"))#6888
#remove chr0 as not true replicates
loc_dup2<-loc_dup[-(which(manifest_EPICv2_IlmnTrim$CHR[loc_dup]=="chr0"))]
length(loc_dup2)#6342

#make a manifest just of replicate probes
rep_manifest<-manifest_EPICv2_IlmnTrim[loc_dup2,]

#how many of the name replicate probes are captured in the 'Rep_Num' column (Illumina identified exact replicates)
table(rep_manifest$Rep_Num)
#  1    2    3    4    5    6    7    8    9   10 
#120 5112  982   69   51    4    1    1    1    1 
#So, the majority are already recognised by Illumina as replicates
#120 are not though...

table(manifest_EPICv2_IlmnTrim$Rep_Num)
#     1      2      3      4      5      6      7      8      9     10 
#929958   5141   1002     70     52      4      1      1      1      1
#More in whole array - likely as chr0 probes (which I've excluded from this analysis) have some reps

table(manifest_EPICv2_IlmnTrim$Rep_Num[which(manifest_EPICv2_IlmnTrim$CHR=="chr0")])
#   1    2    3    4    5 
#6838   29   20    1    1 
#Yes, some of the chr0 probes are replicates

#identify probes with duplicated probe sequences
#create vector of sequences for each probe
joinedseq<-paste(manifest_EPICv2_IlmnTrim$AlleleA_ProbeSeq,manifest_EPICv2_IlmnTrim$AlleleB_ProbeSeq)
#identify duplicates within vector of sequences for each probe
ABdup<-which(duplicated(joinedseq)==TRUE)

#create vector of unique sequences for sequences which are duplicated
ABdup_seq<-unique(joinedseq[ABdup])
length(ABdup_seq)#5571

#identify index of all probes that are sequence duplicates
duplist2<-list()
for(i in 1:length(ABdup_seq)){
  duplist2[[i]]<-which(joinedseq==ABdup_seq[i])
}

#summary table of probes per sequence replicates
table(unlist(lapply(duplist2,length)))
#2    3    4    5    6    7    8   10   11 
#4506  969   28   52   10    1    1    3    1 

#create vector of probes with duplicated sequences
vecsupseq<-unlist(duplist2)
length(vecsupseq)#12407
vecrep<-rep("",nrow(manifest_EPICv2_IlmnTrim))
vecrep[vecsupseq]<-"Y"

#Check sequence replicate probe ids against Illumina defined exact replicates 
knownreps<-which(manifest_EPICv2_IlmnTrim$Rep_Num!=1)
mrep<-match(knownreps,vecsupseq)
length(which(is.na(mrep)==TRUE))#4
nonreps<-manifest_EPICv2_IlmnTrim[knownreps[which(is.na(mrep)==TRUE)],]
#4 sets of replicates are not actually exact replicates

#double check if these 4 probes have replicate sequences
checknonrep<-vector()
for(i in 1:nrow(nonreps)){
  tmp<-which(manifest_EPICv2_IlmnTrim$Name==nonreps$Name[i])
  checknonrep[i]<-identical(manifest_EPICv2_IlmnTrim$AlleleA_ProbeSeq[tmp[1]],manifest_EPICv2_IlmnTrim$AlleleA_ProbeSeq[tmp[2]])
}
checknonrep# FALSE FALSE FALSE FALSE
#Confirmed that 4 of Illumina's 'replicates are not actually replicates, slight sequence differences

##matching instead on location
#identify probes with duplicated locations (excluding chr0)

#create vector of genomic position for each probe
posvec<-paste(manifest_EPICv2_IlmnTrim$CHR,manifest_EPICv2_IlmnTrim$MAPINFO)

#identify duplicate probes according to vector of genomic position for each probe
posdup<-which(duplicated(posvec)==TRUE)
length(posdup)#13230

#extract unique positions from those probes with duplicate positions
uniqposdup<-unique(posvec[posdup])
length(uniqposdup)#[1] 5194

#create list of indexes of all position replicate probes
duplist3<-list()
for(i in 1:length(uniqposdup)){
  duplist3[[i]]<-which(posvec==uniqposdup[i])
}

#create vector of all probes with position replicates
posrepall<-unlist(duplist3)
length(posrepall)#18424
posrep<-rep("",nrow(manifest_EPICv2_IlmnTrim))
posrep[posrepall]<-"Y"

#are the Illumina defined replicates all within the position replicates - expect yes
mrep2<-match(knownreps,posrepall)
length(which(is.na(mrep2)==TRUE))#0
#all known Illumina defined replicates (using Rep_Num) are accounted for in this vector

#How many position replicates unique to posdup?
length(setdiff(posrepall,knownreps))#12151

#make manifest of probes that have replicate locations
manifest_rep<-manifest_EPICv2_IlmnTrim[posrepall,]
#Remove chr0 from replicate manifest
manifest_rep2<-manifest_rep[-which(manifest_rep$CHR=="chr0"),]
nrow(manifest_rep2)#[1] 11535
#Remove known Illumina replicates from replicate manifest
manifest_rep2<-manifest_rep2[-which(manifest_rep2$Rep_Num>1),]
nrow(manifest_rep2)#5313
#Remove those with exact replicates  previously identified
exactseq<-manifest_EPICv2_IlmnTrim$IlmnID[vecsupseq]
manifest_rep3<-manifest_rep2[-na.omit(match(exactseq,manifest_rep2$IlmnID)),]
nrow(manifest_rep3)#204

#calculate how many unique locations are targeted by probes with non-identical probe sequences
manifest_rep3pos<-paste(manifest_rep3$CHR,manifest_rep3$MAPINFO)
length(unique(manifest_rep3pos))#112
length(unique(manifest_rep3$Name))#112
length(table(manifest_rep3$Name))#112
length(intersect(unique(manifest_rep3$Name),unique(manifest_EPICv2_IlmnTrim$Name[vecsupseq])))#28
#112 locations are targeted by more than one probe, each with a distinct sequence, some with a distinct name. 28 of these locations also have an exact replicate 

#Check that Name is always the same for probes matched on locations
checkidrep<-vector()
for(i in 1:length(duplist3)){
  tmp<-manifest_EPICv2_IlmnTrim$Name[duplist3[[i]]]
  checkidrep[i]<-all(sapply(as.list(tmp), FUN = identical, tmp[1]))
}
table(checkidrep)
#FALSE  TRUE 
#1  5193 
table(manifest_EPICv2_IlmnTrim$CHR[duplist3[[which(checkidrep==FALSE)]]])#This turns out to be the chr0 IDs. So all other probes targeting the same location have the same cgID
#chr0 
#6889 

#Have created vectors for probes with duplicated sequences and positions 
#vecdisc indicates discrepant probes between Sesame and Illumina manifests
#vecmiss indicates missing probe info on Illumina manifest 
#vecrep indicates probes with exact sequence matches within EPICv2 manifests
#posrep indicates probes with genomic position overlap within EPICv2 (according to Illumina genomic position)

#These probes will need careful consideration in any cross-platform comparison or DMR analysis

###################################
#Compare EPICv2 Sesame manifest with Sesame manifests of older platforms

###Check whether each sesame manifest has replicate probe IDs, locations or sequences

##Are there any probeID replicates in previous versions? - no
length(which(duplicated(manifest_27K$probeID)==TRUE))#0
length(which(duplicated(manifest_450K$probeID)==TRUE))#0
length(which(duplicated(manifest_EPICv1$probeID)==TRUE))#0

##Any there any probe sequence replicates in previous versions?
AB27<-paste(manifest_27K$ProbeSeq_A,manifest_27K$ProbeSeq_B)
AB450<-paste(manifest_450K$AlleleA_ProbeSeq,manifest_450K$AlleleB_ProbeSeq)
ABEPIC<-paste(manifest_EPICv1$AlleleA_ProbeSeq,manifest_EPICv1$AlleleB_ProbeSeq)
length(which(duplicated(AB27)==TRUE))#23
length(which(duplicated(AB450)==TRUE))#849
length(which(duplicated(ABEPIC)==TRUE))#634

#Extract details of 27K array sequence replicates
AB27dup<-unique(AB27[which(duplicated(AB27)==TRUE)])
AB27duplist<-list()
for(i in 1:length(AB27dup)){
  AB27duplist[[i]]<-grep(AB27dup[i],AB27)
}
tmp<-manifest_27K[unlist(AB27duplist),]
#write out manifest trimmed to 27K sequence replicates
write.csv(tmp,paste(filepath,"results/HM27_seqduplicates.csv",sep=""))

#How many 27K probes have sequence duplicates
length(which(duplicated(paste(tmp$ProbeSeq_A,tmp$ProbeSeq_B))==TRUE))#23
#How many unique 27K sequences are duplicated
length(unique(paste(tmp$ProbeSeq_A,tmp$ProbeSeq_B)))#21
length(unique(tmp$probeID))#44
nrow(tmp)#44
#Interesting that many are the same location and sequence, but different CpG IDs - suggests that originally probe ID was for the probe rather than the location

#Extract details of 450K array sequence replicates 
AB450dup<-manifest_450K[which(duplicated(AB450)==TRUE),]
ctl450<-grep("ctl",AB450dup$Probe_ID)
length(ctl450)==nrow(AB450dup)#TRUE
#all sequence replicates are control probes

#Extract details of EPIC array sequence replicates 
ABEPICdup<-manifest_EPICv1[which(duplicated(ABEPIC)==TRUE),]
ctlE1<-grep("ctl",ABEPICdup$Probe_ID)
length(ctlE1)==nrow(ABEPICdup)#TRUE
#all sequence replicates are control probes

##Are there any location replicates in previous versions

#Create GRange objects of each manifest 
GR27<-GRanges(seqnames=manifest_27K$CpG_chrm,IRanges(start=manifest_27K$CpG_beg,end=manifest_27K$CpG_end))
#GR450<-paste(seqnames=manifest_450K$CpG_chrm,IRanges(start=manifest_450K$CpG_beg,end=manifest_450K$CpG_end))#didn't work as NAs present
#GREPIC<-paste(seqnames=manifest_EPICv1$CpG_chrm,IRanges(start=manifest_EPICv1$CpG_beg,end=manifest_EPICv1$CpG_end))#didn't work as NAs present

#Granges object didn't work for 450K and EPIC as some probes do not have genomic position information
#What are these probes without genomic info. Remove and create GRanges
length(which(is.na(manifest_450K$CpG_chrm)==TRUE))#858
IDsNA_450<-manifest_450K$Probe_ID[which(is.na(manifest_450K$CpG_chrm)==TRUE)]
length((grep('ctl',IDsNA_450)))#850
#Predominantly non-mapping as control probes

length(which(is.na(manifest_EPICv1$CpG_chrm)==TRUE))#649
IDsNA_EPICv1<-manifest_EPICv1$Probe_ID[which(is.na(manifest_EPICv1$CpG_chrm)==TRUE)]
length((grep('ctl',IDsNA_EPICv1)))#635
#Predominantly non-mapping as control probes

rm450GR<-which(is.na(manifest_450K$CpG_chrm)==TRUE)
GR450<-GRanges(seqnames=manifest_450K$CpG_chrm[-rm450GR],IRanges(start=manifest_450K$CpG_beg[-rm450GR],end=manifest_450K$CpG_end[-rm450GR]))

rmEPICv1GR<-which(is.na(manifest_EPICv1$CpG_chrm)==TRUE)
GREPIC<-GRanges(seqnames=manifest_EPICv1$CpG_chrm[-rmEPICv1GR],IRanges(start=manifest_EPICv1$CpG_beg[-rmEPICv1GR],end=manifest_EPICv1$CpG_end[-rmEPICv1GR]))

length(which(duplicated(GR27)==TRUE))#12
length(which(duplicated(GR450)==TRUE))#32
length(which(duplicated(GREPIC)==TRUE))#32

#make a list of EPICv1 probes with duplicated locations
EPIC1locdups<-GREPIC[which(duplicated(GREPIC)==TRUE)]
chr<-start<-vector()
for(i in 1:length(EPIC1locdups)){
  chr[i]<-as.character(EPIC1locdups[i]@seqnames@values)
  start[i]<-as.character(EPIC1locdups[i]@ranges@start)
}
list_chr<-list()
for (i in 1:length(chr)){
  list_chr[[i]]<-manifest_EPICv1$Probe_ID[which(manifest_EPICv1$CpG_chrm==chr[i] & manifest_EPICv1$CpG_beg==start[i])]
}
length(unlist(list_chr))#82
length((list_chr))#41

#make a list of 450K probes with duplicated locations
K450locdups<-GR450[which(duplicated(GR450)==TRUE)]
chr450<-start450<-vector()
for(i in 1:length(K450locdups)){
  chr450[i]<-as.character(K450locdups[i]@seqnames@values)
  start450[i]<-as.character(K450locdups[i]@ranges@start)
}
list_chr450<-list()
for (i in 1:length(chr450)){
  list_chr450[[i]]<-manifest_450K$Probe_ID[which(manifest_450K$CpG_chrm==chr450[i] & manifest_450K$CpG_beg==start450[i])]
}
tmp3<-manifest_450K[match(unlist(list_chr450),manifest_450K$probeID),]
write.csv(tmp3,paste(filepath,"results/HM450_locduplicates.csv",sep=""))
nrow(tmp3)
length(list_chr450)#32
length(unlist(list_chr450))#64

#make a list of 27K probes with duplicated locations
K27locdups<-GR27[which(duplicated(GR27)==TRUE)]
chr27<-start27<-vector()
for(i in 1:length(K27locdups)){
  chr27[i]<-as.character(K27locdups[i]@seqnames@values)
  start27[i]<-as.character(K27locdups[i]@ranges@start)
}
list_chr27<-list()
for (i in 1:length(chr27)){
  list_chr27[[i]]<-manifest_27K$probeID[which(manifest_27K$CpG_chrm==chr27[i] & manifest_27K$CpG_beg==start27[i])]
}
#Some overlapping locations
#Are 27K a subset of those with same probe sequence?
vectmp<-manifest_27K$probeID[which(duplicated(GR27)==TRUE)]
length(which(is.na(match(vectmp,manifest_27K$probeID[unlist(AB27duplist)]))==TRUE))#[1] 1 #Yes mostly, just one with a different sequence
vectmp[which(is.na(match(vectmp,manifest_27K$probeID[unlist(AB27duplist)]))==TRUE)]#cg00896220"
manifest_27K[which(duplicated(GR27)==TRUE)[which(is.na(match(vectmp,manifest_27K$probeID[unlist(AB27duplist)]))==TRUE)],]

tmp2<-manifest_27K[match(unlist(list_chr27),manifest_27K$probeID),]
write.csv(tmp2,paste(filepath,"results/HM27_locduplicates.csv",sep=""))
length(list_chr27)#12
length(unlist(list_chr27))#24

#match between Sesame manifests of different array versions based on 1) probe name, 2) location and 3) sequence 
#pay attention to those probes where there is a discrepancy between the EPICv2 Sesame and Illumina manifest, and deal with replicates appropriately

#Check that Sesame and Illumina manifests are in the same order
identical(manifest_EPICv2_ord$Probe_ID,manifest_EPICv2_IlmnTrim$IlmnID)#TRUE

#How many EPICv2 and EPICv1 probe names in common? 
E2E1_PID<-match(manifest_EPICv2_IlmnTrim$Name,manifest_EPICv1$Probe_ID)
length(which(is.na(E2E1_PID)==FALSE))#726597

#How many EPICv2 and EPICv1 probe sequences in common? 
ABEPICv2<-paste(manifest_EPICv2_ord$AlleleA_ProbeSeq,manifest_EPICv2_ord$AlleleB_ProbeSeq)
E2E1_PSeq<-match(ABEPICv2,ABEPIC)
length(which(is.na(E2E1_PSeq)==FALSE))#727520
#Extract probe names of EPICv1  probes with matched sequences in EPICv2
epic1seqmatch<-manifest_EPICv1$Probe_ID[E2E1_PSeq]

#How many EPICv2 and EPICv1 probe locations in common (using Sesame locations)? 
PosEPICv2<-paste(manifest_EPICv2_ord$CpG_chrm,manifest_EPICv2_ord$CpG_beg,manifest_EPICv2_ord$CpG_end)
PosEPICv1<-paste(manifest_EPICv1$CpG_chrm,manifest_EPICv1$CpG_beg,manifest_EPICv1$CpG_end)
E2E1_Pos<-match(PosEPICv2,PosEPICv1)
length(which(is.na(E2E1_Pos)==FALSE))#727605
epic1locmatch<-manifest_EPICv1$Probe_ID[E2E1_Pos]
#do these probes with matched position between EPICv1 and EPICv2 match any of those with duplicated locations on EPICv1 
length(which(is.na(match(unlist(list_chr),epic1locmatch))==FALSE))#0 - no matches
      

#How many EPICv2 and 450K probe names in common? 
E2_450_PID<-match(manifest_EPICv2_IlmnTrim$Name,manifest_450K$Probe_ID)
length(which(is.na(E2_450_PID)==FALSE))#399178

#How many EPICv2 and 450K probe sequences in common? 
E2_450_PSeq<-match(ABEPICv2,AB450)
length(which(is.na(E2_450_PSeq)==FALSE))#399717
#Extract probe names of 450K  probes with matched sequences in EPICv2
K450seqmatch<-manifest_450K$Probe_ID[E2_450_PSeq]

#How many EPICv2 and 450K probe locations in common (using Sesame locations)? 
Pos450v1<-paste(manifest_450K$CpG_chrm,manifest_450K$CpG_beg,manifest_450K$CpG_end)
E2_450_Pos<-match(PosEPICv2,Pos450v1)
length(which(is.na(E2_450_Pos)==FALSE))#399758
K450locmatch<-manifest_450K$Probe_ID[E2_450_Pos]
##do these probes with matched position between 450K and EPICv2 match any of those with duplicated locations on 450K
length(which(is.na(match(unlist(list_chr450),K450locmatch))==FALSE))#3 duplicate locations targeted


#make additional vector with info about the probes with matched position between 450K and EPICv2, that also have duplicated locations on 450K 
dup450_matchEPIC2<-K450locmatch[match(na.omit(K450locmatch[match(unlist(list_chr450),K450locmatch)]),K450locmatch)]
K450locmatch2<-vector(length=length(K450locmatch))
mat<-na.omit(match(unlist(list_chr450),K450locmatch))
for(i in 1:length(dup450_matchEPIC2)){
  K450locmatch2[mat[i]]<-list_chr450[grep(dup450_matchEPIC2[i],list_chr450)][[1]][2]
}
K450locmatch2[which(K450locmatch2=="FALSE")]<-""

#How many EPICv2 and 27K probe names in common? 
E2_27_PID<-match(manifest_EPICv2_IlmnTrim$Name,manifest_27K$probeID)
length(which(is.na(E2_27_PID)==FALSE))#24530

#How many EPICv2 and 27K probe sequences in common? 
E2_27_PSeq<-match(ABEPICv2,AB27)
length(which(is.na(E2_27_PSeq)==FALSE))#2178
K27seqmatch<-manifest_27K$probeID[E2_27_PSeq]

#How many EPICv2 and 27K probe locations in common (using Sesame locations)? 
Pos27<-paste(manifest_27K$CpG_chrm,manifest_27K$CpG_beg,manifest_27K$CpG_end)
E2_27_Pos<-match(PosEPICv2,Pos27)
length(which(is.na(E2_27_Pos)==FALSE))#24571
K27locmatch<-manifest_27K$probeID[E2_27_Pos]

##do these probes with matched position between 27K and EPICv2 match any of those with duplicated locations on 27K
match(unlist(list_chr27),K27locmatch)#1 duplicate location targeted
match(manifest_27K$probeID[unlist(AB27duplist)],K27locmatch)#1 duplicate probe sequence targeted (same as location in previous line)

#make additional vector with info about the probes with matched position between 27K and EPICv2, that also have duplicated locations on 27K 
dup27_matchEPIC2<-K27locmatch[match(na.omit(K27locmatch[match(unlist(list_chr27),K27locmatch)]),K27locmatch)]
K27locmatch2<-vector(length=length(K27locmatch))
mat<-na.omit(match(unlist(list_chr27),K27locmatch))
K27locmatch2[mat]<-list_chr27[grep(dup27_matchEPIC2,list_chr27)][[1]][2]
K27locmatch2[which(K27locmatch2=="FALSE")]<-""


#Between 27K and EPICv2 there are many similar locations and names, but changed sequences

#make vector of probes matched by name between older array versions and EPICv2
#27Knames
vec27probeID<-manifest_27K$probeID[match(manifest_EPICv2_IlmnTrim$Name,manifest_27K$probeID)]
length(which(is.na(vec27probeID)==FALSE))#[1] 24530

#450Knames
vec450probeID<-manifest_450K$Probe_ID[match(manifest_EPICv2_IlmnTrim$Name,manifest_450K$Probe_ID)]
length(which(is.na(vec450probeID)==FALSE))#[1] 399178
#are any of these duplicated
length(which(duplicated(vec450probeID[which(is.na(vec450probeID)==FALSE)])==TRUE))#4736 - yes

#EPICnames
vecEPICv1probeID<-manifest_EPICv1$Probe_ID[match(manifest_EPICv2_IlmnTrim$Name,manifest_EPICv1$Probe_ID)]
length(which(is.na(vecEPICv1probeID)==FALSE))#[1] 726597
#are any of these duplicated
length(which(duplicated(vecEPICv1probeID[which(is.na(vecEPICv1probeID)==FALSE)])==TRUE))#4795 -yes


#There is one which is not in Rep_Num but has same vecrep and posrep, probably an error (cg10335333)
#Would need to check this against probe ID column
#Column indicating those probes with probeID overlap within EPICv2 - posrep (checked and IDs are all consistent for location matched probes)
#will find those beyond rep_num, summarise and comment

#Use vectors of Illumina/sesame discrepant probes and replicate probes to check if these interfere with comparison with older array
#Reminder of different vectors used:
#vecdisc indicates discrepant probes between Sesame and Illumina manifests
#vecmiss indicates missing probe info on Illumina manifest 
#vecrep indicates probes with exact sequence matches within EPICv2 manifests
#posrep indicates probes with genomic position overlap within EPICv2 (according to Illumina genomic position)

#Do any of the Illumina sequence replicates having discrepant mapping between Sesame and Illumina?
which(vecdisc=="Y" & vecmiss=="" & vecrep=="Y")
#integer(0) - no
#Do any of the Illumina location replicates having discrepant mapping between Sesame and Illumina?
which(vecdisc=="Y" & vecmiss=="" & posrep=="Y")
#integer(0)

#Count number of affected probes in each category
table(vecdisc)
#Y 
#928715   7516 

table(vecmiss)
#Y 
#929342   6889 

table(vecrep)
#Y 
#923824  12407 

table(posrep)
#Y 
#917807  18424 

###################################
##Making new version of the manifest including information about sesame and Illumina manifest comparison, replicate probes and comparison with manifests of older arrays 

#Start with Illumina manifest. Then add the following new columns:
#Sesame coordinates (#Noting Sesame coordinates are a range)
#Discrepant mapping between Illumina and Sesame manifests - vecdisc
#Probes with missing location info on Illumina manifest - vecmiss
#Probes with exact probe ID ('Name') matches within EPICv2 -dupsname
#Probes with exact sequence matches within EPICv2 - vecrep
#Probes with genomic position overlap within EPICv2 - posrep

#Between EPICv2 and EPICv1
#Probes with same probe name - vecEPICv1probeID
#Probes with same genomic location (give probe in this column) - epic1locmatch
#Probes with same sequence (give probe ID in this column) - epic1seqmatch
#Probes with same genomic location where replicate probe locations - none for EPICv1

#Between EPICv2 and 450K
#Probes with same probe name - vec450probeID
#Probes with same genomic location (give probe in this column) - K450locmatch
#Probes with same sequence (give probe ID in this column) - K450seqmatch
#Probes with same genomic location where replicate probe locations - K450locmatch2

#Between EPICv2 and 27K
#Probes with same probe name - vec27probeID
#Probes with same genomic location (give probe in this column) - K27locmatch
#Probes with same sequence (give probe ID in this column) - K27seqmatch
#Probes with same genomic location where replicate probe locations - K27locmatch2

#additionally create vectors for columns to indicate how within-EPICv2 replicates are related 

#create empty vectors for probe IDs and replicate numbers of EPICv2 location replicates
posrep_IlmnIDs<-rep("",nrow(manifest_EPICv2_IlmnTrim))
posrep_RepNum<-rep(1,nrow(manifest_EPICv2_IlmnTrim))

#remove probes with missing locations (chr 0) from location-replicate list duplist3
duplist4<-duplist3[-which(unlist(lapply(duplist3,length))>500)]

#create vector of concatenated location-replicate IDs at location-replicate probes, and corresponding Rep_Num
for(i in 1:length(duplist4)){
  tmp<-paste(manifest_EPICv2_IlmnTrim$IlmnID[duplist4[[i]]],collapse=";")
  posrep_IlmnIDs[duplist4[[i]]]<-tmp
  for(j in 1:length(duplist4[[i]])){
    posrep_RepNum[duplist4[[i]][j]]<-j
  }
}

#amend posrep vector (location-replicate) for those with missing chr info
posrep_IlmnIDs[which(vecmiss=="Y")]<-""
posrep_RepNum[which(vecmiss=="Y")]<-NA

#create empty vectors for probe IDs and replicate numbers of EPICv2 sequence replicates
seqrep_IlmnIDs<-rep("",nrow(manifest_EPICv2_IlmnTrim))
seqrep_RepNum<-rep(1,nrow(manifest_EPICv2_IlmnTrim))

#create vector of concatenated sequence-replicate IDs at sequence-replicate probes, and corresponding Rep_Num
for(i in 1:length(duplist2)){
  tmp<-paste(manifest_EPICv2_IlmnTrim$IlmnID[duplist2[[i]]],collapse=";")
  seqrep_IlmnIDs[duplist2[[i]]]<-tmp
  for(j in 1:length(duplist2[[i]])){
    seqrep_RepNum[duplist2[[i]][j]]<-j
  }
}

#comparison of probes per replicated for location and sequence replicates
table(posrep_RepNum[-(which(vecmiss=="Y"))])
#1      2      3      4      5      6      7      8      9     10 
#923000   5193   1016     74     51      4      1      1      1      1 

table(seqrep_RepNum[-(which(vecmiss=="Y"))]) 
#1      2      3      4      5      6      7      8      9     10 
#923123   5109    982     69     51      4      1      1      1      1 

#so, combining all data:
#Illumina data, minus nv and control probes
new_manifest<-cbind(manifest_EPICv2_IlmnTrim,#Illumina manifest minus nv and control probes
                    manifest_EPICv2_ord[,1:3],#sesame coordinates
                    vecdisc,#discrepant mapping between Illumina and sesame
                    vecmiss,#probes with missing location info on Illumina
                    dupsname,#probes with exact probe ID ('Name') matches within EPICv2
                    vecrep,#probes with exact sequence matches within EPICv2 
                    seqrep_IlmnIDs,#IlmnIDs of those with exact sequence matches within EPICv2 
                    seqrep_RepNum,#Rep Num of IlmnID set of those probes with exact sequence matches within EPICv2 
                    posrep,#probes with genomic position overlap within EPICv2 
                    posrep_IlmnIDs,#IlmnIDs of those with exact genomic position matches within EPICv2 
                    posrep_RepNum,#Rep Num of IlmnID set of those probes with exact genomic position matches within EPICv2 
                    vecEPICv1probeID,#Probes with same probe name match in EPICv1  
                    epic1seqmatch, #Probes with sequence match in EPICv1
                    epic1locmatch,#Probes with genomic location match in EPICv1
                    vec450probeID,#Probes with same probe name match in 450K 
                    K450seqmatch, #Probes with sequence match in 450K
                    K450locmatch,#Probes with genomic location match in 450K 
                    K450locmatch2, #Probes with additional location match in 450K
                    vec27probeID,#Probes with same probe name match in 27K
                    K27seqmatch, #Probes with sequence match in 27K
                    K27locmatch,#Probes with genomic location match in 27K
                    K27locmatch2) #Probes with additional location match in 27K
#Note, location matching between array versions based on sesame data
#Advise that if users want to look at probes, they should use the original Illumina manifest
#This manifest is for enabling cross-platform comparisons (in terms of genomic coverage and data integration)

#Summary stats on 'new' manifest

#How many probes have discrepant locations between Illumina and sesame manifest?
length(which(new_manifest$vecdisc=="Y"))#7516

#How many uniquely located probes have discrepant locations between Illumina and sesame manifest?
length(which(new_manifest$vecdisc=="Y" & new_manifest$posrep!="Y"))#627

#How many probes have discrepant locations between Illumina and sesame manifest excluding those missing mapping info in Illumina?
length(which(new_manifest$vecdisc=="Y" & new_manifest$vecmiss==""))#627

#How many probes are missing locations on Illumina manifest?
length(which(new_manifest$missvec=="Y"))#6889

#How many probes on EPICv2 are duplicates based on probe name?
length(which(new_manifest$dupsname=="Y"))#11622

#How many unique probe names are there within EPICv2 duplicates based on probe name?
length(unique(new_manifest$Name[which(new_manifest$dupsname=="Y")]))#5225

#How many unique probe names are there within EPICv2 probe name duplicates, excluding chromosome 0 probes?
length(unique(new_manifest$Name[which(new_manifest$dupsname=="Y" & new_manifest$vecmiss!="Y")]))#5195

#How many unique probe names are there within EPICv2 probe name duplicates, excluding exact sequence matching probes?
length(unique(new_manifest$Name[which(new_manifest$dupsname=="Y" & new_manifest$vecrep!="Y")]))#116

#How many unique probe names are there within EPICv2 probe name duplicates, that are Illumina recognised replicates but not actually exact sequence matching probes?
length(unique(new_manifest$Name[which(new_manifest$dupsname=="Y" & new_manifest$vecrep!="Y" & new_manifest$Rep_Num=="2")]))#4

#How many probes on EPICv2 are unique based on probe name?
length(unique(new_manifest$Name))#929834

#How many probes on EPICv2 are non-duplicates based on probe name?
length(which(new_manifest$dupsname!="Y"))#924609

#How many probes on EPICv2 are non-duplicates based on probe location?
length(which(new_manifest$posrep!="Y"))#917807 (Note, all probes with missing location removed from posrep)

#How many probes on EPICv2 are non-duplicates based on probe sequence?
length(which(new_manifest$vecrep!="Y"))#923824

#Could use these vectors to subset duplicate names/sequences/locations and then perform own summary stats on this

#How many EPICv2 probe names overlap EPICv1 probe names
length(which(new_manifest$vecEPICv1probeID!=""))# 726597

#How many EPICv2 probe seq overlap EPICv1 probe seq
length(which(new_manifest$epic1seqmatch!=""))# 727520

#How many EPICv2 probe location overlap EPICv1 probe location
length(which(new_manifest$epic1locmatch!=""))# 727605

#How many EPICv2 probe location overlap EPICv1 probe location, at probes with discrepant locations between Sesame and Illumina
length(which(new_manifest$epic1locmatch!="" & vecdisc=="Y"))#392

#Then repeat filtering out duplicate name/seq/locations probes...
length(unique(new_manifest$Name[which(new_manifest$vecEPICv1probeID!="")]))# 721802
length(unique(ABEPICv2[which(new_manifest$epic1seqmatch!="")])) #722794
length(unique(PosEPICv2[which(new_manifest$epic1locmatch!="")])) # 722815
#will have to extract these vectors for then making multi-venn

#Noticed that one probe is an exact replicate, but not listed as such in manifest Rep_Num column
length(unique(manifest_EPICv2_IlmnTrim$Name[which(vecrep=="Y" & posrep=="Y")]))#6137
length(which(manifest_EPICv2_IlmnTrim$Rep_Num==1 & vecrep=="Y" & posrep=="Y"))#6138
id2<-manifest_EPICv2_IlmnTrim$Name[which(manifest_EPICv2_IlmnTrim$Rep_Num==1 & vecrep=="Y" & posrep=="Y")]
id2[which(duplicated(id2))]#"cg10335333"

save(new_manifest,file=paste(filepath,"results/NewManifest.RData",sep=""))

#Compare vecdisc with Sesame masking info
length(which(vecdisc=="Y" & vecmiss!="Y"))#627
discrepant<-new_manifest$IlmnID[which(vecdisc=="Y" & vecmiss!="Y")]

#read in Sesame mask list
maskEPIC2<-read.table(paste(filepath,"metadata/EPICv2.hg38.mask_sesame.tsv",sep=""),sep="\t",header=T)
#match up
ft<-table(maskEPIC2$M_general[match(discrepant,maskEPIC2$Probe_ID)])
1-(ft[1]/ft[2])#99.8
