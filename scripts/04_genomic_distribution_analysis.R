
#define file path
filepath<-".../your/file/path"

# Load packages
library(ChIPpeakAnno)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(ggplot2)
library(dplyr)
library(tibble)
library(data.table)
library(tidyr)
library(AnnotationHub)
library(knitr)
library(rtracklayer)
library(IRanges)
library(stringr)
library(reshape)
library(patchwork)
library(tidyverse)
library(hacksaw)

# GENOMIC ANNOTATION ANALYSIS 

# Load EPICv2 probe category bed files and convert to genomic regions objects 
## Reinstated probes
Locations_reinstated <- read.delim(file.path(filepath,"results/Locations_reinstated.bed"), header = FALSE)
Locations_reinstated.gr <- makeGRangesFromDataFrame(Locations_reinstated, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3", na.rm = TRUE)
## Retained probes
Locations_retained <- read.delim(file.path(filepath,"results/Locations_retained.bed"), header = FALSE)
Locations_retained <- na.omit(Locations_retained)
Locations_retained.gr <- makeGRangesFromDataFrame(Locations_retained, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3")
## New probes 
Locations_new <- read.delim(file.path(filepath,"results/Locations_new.bed"), header = FALSE)
Locations_new.gr <- makeGRangesFromDataFrame(Locations_new, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3", na.rm = TRUE)
## Excluded probes 
Locations_excluded <- read.delim(file.path(filepath,"results/Locations_excluded.bed"), header = FALSE)
Locations_excluded.gr <- makeGRangesFromDataFrame(Locations_excluded, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3", na.rm = TRUE)
## All probes 
Locations_all <- read.delim(file.path(filepath,"results/Locations_all.bed"), header = FALSE)
Locations_all <- na.omit(Locations_all)
Locations_all.gr <- makeGRangesFromDataFrame(Locations_all, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3", na.rm = TRUE)

## Create category "EPICv1" To show the probes present in the 850K version of the array for comparison 
Locations_EPICv1.gr <- c(Locations_excluded.gr, Locations_retained.gr)
# Save bed file of EPICv1 locations 
write.table(Locations_EPICv1.gr, file.path(filepath, "results/Locations_EPICv1.bed"), row.names = FALSE, quote = FALSE, col.names = FALSE)

## create genomic ranges list for all 6 genomic objects 
EPICv2_probes.grl <- GRangesList("Reinstated"=Locations_reinstated.gr, "Excluded"=Locations_excluded.gr, "New"=Locations_new.gr, "Retained"=Locations_retained.gr, "EPICv2"=Locations_all.gr, "EPICv1"=Locations_EPICv1.gr)

# Change default script for genomicElementDistribution from ChIPpeakanno package to allow for plot to be formatted correctly
source(file.path(filepath, "functions/genomicElementDistribution_custom.R"))
#reads in mygrl

# # Plot the genomic locations for each probe *gives plot separated by location rather than probe type* 
# EPICv2_probes.out <- mygrl(EPICv2_probes.grl, TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene, promoterRegion = c(upstream = 2000, downstream = 500), geneDownstream = c(upstream = 0, downstream = 2000), promoterLevel = list(breaks = c(-2000, -1000, -500, 0, 500), labels = c("upstream 1-2Kb", "upstream 0.5-1Kb", "upstream <500b", "TSS - 500b"), colors = c("#FFE5CC", "#FFCA99", "#FFAD65", "#FF8E32")))
# 
# # Plot the genomic locations for each probes *gives a better figure than above and is separated by probes*
# dat1<-EPICv2_probes.out$plot$data
# pp <- ggplot(dat1[which(dat1$percentage>0 & dat1$type!="undefined"),], aes_string(x = "type", y = "percentage", 
#                                                                                   fill = "source")) + geom_bar(stat = "identity", position=position_dodge()) + coord_flip() + 
#   facet_wrap(as.formula("~ category"), ncol = 1, scales = "free_y") + 
#   theme_bw()
# pp
# ggsave(file.path(filepath, "plots/EPICv2_Probes_Genomic_Anno.pdf"), plot = last_plot())
# 
# write.csv(dat1, file.path(filepath, "results/Genomic_Anno.csv"))

# Load in bed files of genic regions, derived from gencode: TSS (promoter), Gene Body and Intergenic (downloaded and formatted using bedtools as in DownloadingGencodeData_Command.log)
Locations_TSS <- read.delim(file.path("metadata/tss_2kb.hierarchy.bed"), header = FALSE)
colnames(Locations_TSS) <- c("Chrom", "Start", "End")
Locations_TSS.gr <- makeGRangesFromDataFrame(Locations_TSS, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "Chrom", start.field = "Start", end.field = "End")

Locations_GeneBody <- read.delim(file.path("metadata/transcript_2kb.hierarchy.bed"), header = FALSE)
colnames(Locations_GeneBody) <- c("Chrom", "Start", "End")
Locations_GeneBody.gr <- makeGRangesFromDataFrame(Locations_GeneBody, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "Chrom", start.field = "Start", end.field = "End")

Locations_Intergenic <- read.delim(file.path("metadata/intergenic.hierarchy.bed"), header = FALSE)
colnames(Locations_Intergenic) <- c("Chrom", "Start", "End")
Locations_Intergenic.gr <- makeGRangesFromDataFrame(Locations_Intergenic, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "Chrom", start.field = "Start", end.field = "End")
# Make GRanges into GRangesList
Genomic_Anno.grl <- GRangesList("TSS"=Locations_TSS.gr, "Gene_Body"=Locations_GeneBody.gr, "Intergenic"=Locations_Intergenic.gr)

## Find overlaps between each probe category and each genomic element category 
Genomic_anno_list<-list()

for (i in 1:length(Genomic_Anno.grl)) {
  
  Genomic_Anno_matrix <- matrix(data = NA, ncol = 2, nrow = length(EPICv2_probes.grl))
  rownames(Genomic_Anno_matrix) <- names(EPICv2_probes.grl)
  
  for(j in 1:length(EPICv2_probes.grl)){
    ol <- findOverlaps(Genomic_Anno.grl[[i]], EPICv2_probes.grl[[j]])
    # Find number of unique probes overlapping Genomic_Anno.grl[[i]]
    Count <- length(unique(ol@to))
    # Find percentage of unique probes overlapping Genomic_Anno.grl[[i]]
    Percent <- (length(unique(ol@to))/length(EPICv2_probes.grl[[j]]))*100
    Genomic_Anno_matrix[j,] <- c(Count, Percent)
  }
  
  Genomic_anno_list[[i]]<-Genomic_Anno_matrix
  
}  

Genomic_Anno_matrixB <- do.call(cbind, Genomic_anno_list)
colnames(Genomic_Anno_matrixB) <- c("TSS_Count", "TSS_Percent", "GeneBody_Count", "GeneBody_Percent", "Intergenic_Count", "Intergenic_Percent")
write.csv(Genomic_Anno_matrixB, file.path(filepath, "results/Genomic_Anno_Summary.csv"))

names(Genomic_anno_list) <- c("TSS", "Gene_Body", "Intergenic")

#Find overlaps between each unique genomic element and probe category 

##TSS
TSS_Coverage_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))

for(j in 1:length(EPICv2_probes.grl)){
  # Find number and percentage of individual probes overlapping TSS
  TSS_density_ol <- findOverlaps(Locations_TSS.gr, EPICv2_probes.grl[[j]])
  TSS_uniq_probe_n <- length(unique(TSS_density_ol@to))
  TSS_uniq_probe_perc <- (length(unique(TSS_density_ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique CpGs targeted by probes 
  TSS_uniq_n <- length(unique(TSS_density_ol@from))
  TSS_uniq_perc <- (length(unique(TSS_density_ol@from))/length(Locations_TSS.gr))*100
  # Put CpG data into matrix
  TSS_Coverage_Matrix[j,] <- c(TSS_uniq_probe_n, TSS_uniq_probe_perc, TSS_uniq_n,TSS_uniq_perc)
}

##Body
Body_Coverage_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))

for(j in 1:length(EPICv2_probes.grl)){
  # Find number and percentage of individual probes overlapping Body
  Body_density_ol <- findOverlaps(Locations_GeneBody.gr, EPICv2_probes.grl[[j]])
  Body_uniq_probe_n <- length(unique(Body_density_ol@to))
  Body_uniq_probe_perc <- (length(unique(Body_density_ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique CpGs targeted by probes 
  Body_uniq_n <- length(unique(Body_density_ol@from))
  Body_uniq_perc <- (length(unique(Body_density_ol@from))/length(Locations_GeneBody.gr))*100
  # Put CpG data into matrix
  Body_Coverage_Matrix[j,] <- c(Body_uniq_probe_n, Body_uniq_probe_perc, Body_uniq_n,Body_uniq_perc)
}

##Intergenic
Ig_Coverage_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))

for(j in 1:length(EPICv2_probes.grl)){
  # Find number and percentage of individual probes overlapping Ig
  Ig_density_ol <- findOverlaps(Locations_Intergenic.gr, EPICv2_probes.grl[[j]])
  Ig_uniq_probe_n <- length(unique(Ig_density_ol@to))
  Ig_uniq_probe_perc <- (length(unique(Ig_density_ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique CpGs targeted by probes 
  Ig_uniq_n <- length(unique(Ig_density_ol@from))
  Ig_uniq_perc <- (length(unique(Ig_density_ol@from))/length(Locations_Intergenic.gr))*100
  # Put CpG data into matrix
  Ig_Coverage_Matrix[j,] <- c(Ig_uniq_probe_n, Ig_uniq_probe_perc, Ig_uniq_n,Ig_uniq_perc)
}


# Save out csv of probe density coverage summary for CpG islands 
write.csv(TSS_Coverage_Matrix, file.path(pwd, "results/TSS_Coverage_Summary.csv"))
write.csv(Body_Coverage_Matrix, file.path(pwd, "results/GeneBody_Coverage_Summary.csv"))
write.csv(Ig_Coverage_Matrix, file.path(pwd, "results/Intergenic_Coverage_Summary.csv"))

# Plot Genomic Anno COUNTS for each category 
GA_Counts_All <- as.data.frame(rbind(Genomic_anno_list$TSS[,1], Genomic_anno_list$Gene_Body[,1], Genomic_anno_list$Intergenic[,1]))
rownames(GA_Counts_All) <- c("TSS", "Gene_Body", "Intergenic")
GA_Counts_All <- rownames_to_column(GA_Counts_All)
GA_Counts_All <- pivot_longer(GA_Counts_All, cols = c("Reinstated", "Excluded", "New", "Retained", "EPICv2", "EPICv1"))
colnames(GA_Counts_All) <- c("Genomic_Element", "Probe_Category", "Count")

options(scipen=999)
GA_Counts_All_Plot <- ggplot(GA_Counts_All, aes(x=Genomic_Element, y=Count, fill=Probe_Category)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label = Count), vjust=-0.5, position=position_dodge(width=0.9), size = 2.7) + 
  ggtitle("Genomic Anno Counts")
GA_Counts_All_Plot

ggsave(file.path(filepath, "plots/Genomic_Anno_Counts.pdf"), plot = last_plot())

# Plot Genomic Anno PERCENTAGES for each probe category
GA_Percent_All <- as.data.frame(rbind(Genomic_anno_list$TSS[,2], Genomic_anno_list$Gene_Body[,2], Genomic_anno_list$Intergenic[,2]))
rownames(GA_Percent_All) <- c("TSS", "Gene_Body", "Intergenic")
GA_Percent_All <- rownames_to_column(GA_Percent_All)
GA_Percent_All <- pivot_longer(GA_Percent_All, cols = c("Reinstated", "Excluded", "New", "Retained", "EPICv2", "EPICv1"))
colnames(GA_Percent_All) <- c("Genomic_Element", "Probe_Category", "Percentage")

GA_Percent_All_Plot <- ggplot(GA_Percent_All, aes(x=Genomic_Element, y=Percentage, fill=Probe_Category)) + 
  geom_bar(position="dodge", stat="identity") + 
  geom_text(aes(label = round(Percentage, digits = 2)), vjust=-0.5, position=position_dodge(width=0.9), size = 2.7) + 
  ggtitle("Genomic Anno Percentage")

GA_Percent_All_Plot

ggsave(file.path(filepath, "plots/Genomic_Anno_Percent.pdf"), plot = last_plot())

#****#

# CpG ELEMENT ANNOTATION ANALYSIS 

# Obtain CpG location info from UCSC 
session <- browserSession()
genome(session) <- "hg38"
query <- ucscTableQuery(session, "cpgIslandExt")
CpGislands <- track(query)
genome(CpGislands) <- NA

# Create CpG island shores from CpG Island locations 
CpGshores <- GenomicRanges::setdiff(resize(CpGislands, width = width(CpGislands)+4000, fix="center"), CpGislands)

# Save out CpG Island and Shore bed files *This is so we can read in as a table to convert to GRanges later 
export(CpGislands, (file.path("metadata/CpGislands.bed")))
export(CpGshores, (file.path("metadata/CpGshores.bed")))

# Creates a matrix for raw counts of probes in each category
CpGMatrix <- matrix(nrow = length(EPICv2_probes.grl), ncol = 3, dimnames = list(names(EPICv2_probes.grl), c("CpG Island", "CpG Shore", "Non Cpg Island")))

# Load CpG Island and Shore bed files and convert to Genomic Ranges object
CpGislands <- read.delim(file.path(filepath, "metadata/CpGislands.bed"), header = FALSE)
CpGislands.gr <- makeGRangesFromDataFrame(CpGislands, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3")
CpGshores <- read.delim(file.path(filepath, "metadata/CpGshores.bed"), header = FALSE)
CpGshores.gr <- makeGRangesFromDataFrame(CpGshores, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# Find overlaps each probe category: CpG Islands, CpG Shores, Non-Island 
CpGprobetallylist <- list()
CpG_Coverage_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))
rownames(CpG_Coverage_Matrix) <- names(EPICv2_probes.grl)
colnames(CpG_Coverage_Matrix) <- c("CpG_uniq_probe_n", "CpG_uniq_probe_perc", "CpG_uniq_n", "CpG_uniq_perc")

for (j in 1:length(EPICv2_probes.grl)) {
  ## Find ol for CpGIslands
  island_ol <- subsetByOverlaps(EPICv2_probes.grl[[j]], CpGislands.gr, invert = FALSE)
  # Find number of probes overlapping CpG Islands 
  CpGIsland_ol <- island_ol@strand@lengths
  ## Find ol for CpGShores
  shore_ol <- subsetByOverlaps(EPICv2_probes.grl[[j]], CpGshores.gr, type = "within", invert = FALSE)
  # Find number of probes overlapping CpG Islands 
  CpGShore_ol <- shore_ol@strand@lengths
  # Find ol for Non-CpG elements
  non_ol <- subsetByOverlaps(EPICv2_probes.grl[[j]], c(island_ol, shore_ol), invert = TRUE)
  # Find number of probes overlapping Non-CpG elements 
  NonCpG_ol <- non_ol@strand@lengths
  CpGMatrix[j,] <- c(CpGIsland_ol, CpGShore_ol, NonCpG_ol)
  # Summarise the number of probes per SE 
  island_density_ol <- findOverlaps(CpGislands.gr, EPICv2_probes.grl[[j]])
  CpGprobetally <- table(table(island_density_ol@from))
  CpGprobetallylist[[j]] <- CpGprobetally
  # Find number and percentage of individual probes overlapping CpG Islands 
  CpG_uniq_probe_n <- length(unique(island_density_ol@to))
  CpG_uniq_probe_perc <- (length(unique(island_density_ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique CpGs targeted by probes 
  CpG_uniq_n <- length(unique(island_density_ol@from))
  CpG_uniq_perc <- (length(unique(island_density_ol@from))/length(CpGislands.gr))*100
  # Put CpG data into matrix
  CpG_Coverage_Matrix[j,] <- c(CpG_uniq_probe_n, CpG_uniq_probe_perc, CpG_uniq_n, CpG_uniq_perc)
}

names(CpGprobetallylist) <- names(EPICv2_probes.grl)
CpGMatrix.df <- as.data.frame(CpGMatrix)
write.csv(CpGMatrix, file.path(filepath, "results/CpG_Raw.csv"))

# Save out csv of probe density coverage summary for CpG islands 
write.csv(CpG_Coverage_Matrix, file.path(filepath, "results/CpG_Coverage_Summary.csv"))

# Plot summary of coverage 
CpG_Coverage_Matrix <- rownames_to_column(as.data.frame(CpG_Coverage_Matrix))

# Plot unique CpGs targeted by probes as %
Uniq_CpG <- ggplot(CpG_Coverage_Matrix, aes(x = rowname, y = CpG_uniq_perc, fill = rowname)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label=CpG_uniq_n), vjust=-0.3, size=4) +
  ggtitle("Unique CpG Islands Targeted by Probes") +
  ylab("Percentage of CpG Islands") 
Uniq_CpG
ggsave(file.path(filepath, "plots/Unique_CpG_AllProbes.pdf"))

# Plot percent of individual probes that overlap with CpG
Uniq_Probe_CpG_OL <- ggplot(CpG_Coverage_Matrix, aes(x = rowname, y = CpG_uniq_probe_perc, fill = rowname)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label=CpG_uniq_probe_n), vjust=-0.3, size=4) +
  ggtitle("Number EPICv2 Probes Overlapping CpG Islands") +
  ylab("Percentage of Probes")
Uniq_Probe_CpG_OL
ggsave(file.path(filepath, "plots/Unique_CpG_Probes_OL.pdf"))

# Create new matrix for probe count as % of probes in each category 
CpGMatrix_per <- round((CpGMatrix/rowSums(CpGMatrix))*100,2)

# Save matrix
write.csv(CpGMatrix_per, file.path(filepath, "results/CpG_Percentage.csv"))

# Convert to data frame and change df so it can be plotted
CpG_df <- as.data.frame(CpGMatrix_per)
CpG_df <- rownames_to_column(CpG_df, var = "Probes")
CpG_df <- pivot_longer(CpG_df, cols = 2:4, names_to = "CpGCategory", values_to = "Percentage")

# Plot above data frame - % of probes within each category 
options(scipen = 999)
CpG_plot <- ggplot(data=CpG_df, aes(x=CpGCategory, y=Percentage, fill=Probes)) + geom_bar(stat = "identity", position = position_dodge())
print(CpG_plot)
ggsave(file.path(filepath, "plots/EPICv2_Probes_CpG_Anno.pdf"), plot = last_plot())

# Plot CpG matrix - raw counts of probes within each category 
CpG_Counts <- as.data.frame(CpGMatrix.df)
CpG_Counts <- rownames_to_column(CpG_Counts, var = "Probes")
CpG_Counts <- pivot_longer(CpG_Counts, cols = 2:4, names_to = "CpGCategory", values_to = "Counts")

options(scipen = 999)
CpG_count_plot <- ggplot(data=CpG_Counts, aes(x=CpGCategory, y=Counts, fill=Probes)) + geom_bar(stat = "identity", position = position_dodge())
print(CpG_count_plot)
ggsave(file.path(filepath, "plots/EPICv2_Probes_CpG_Anno_Counts.pdf"), plot = last_plot())

# plot only v1 & v2 as counts 
pdf(file = file.path(filepath, "plots/CpG_ProbeTally_v1v2_Counts.pdf"), width = 10)
par(mfrow=c(1,2))

CpG_Probetally_All <- barplot((CpGprobetallylist[[5]])[1:20], col = "pink", main = "EPICv2 Probes", xlab = "Probe Number per CpG Island", ylab = "Number Probe Overlap", ylim = c(0,5000))
text(CpG_Probetally_All, CpGprobetallylist[[5]][1:20], labels = CpGprobetallylist[[5]][1:20], cex = 0.6, pos = 3, xpd = TRUE)
CpG_Probetally_EPICv1 <- barplot((CpGprobetallylist[[6]])[1:20], col = "pink", main = "EPICv1 Probes", xlab = "Probe Number per CpG Island", ylab = "Number Probe Overlap", ylim = c(0,5000))
text(CpG_Probetally_EPICv1, CpGprobetallylist[[6]][1:20], labels = CpGprobetallylist[[6]][1:20], cex = 0.6, pos = 3, xpd = TRUE)
dev.off()

#***#

# SUPER ENHANCER ANNOTATION ANALYSIS 

## Read in csv file of super-enhancer locations, downloaded from http://www.licpathway.net/sedb/download.php on 29 March 2023
SE_Package.csv <- read.csv(file.path(filepath, "metadata/SE_package.csv"), header = TRUE)

## SE Package sample IDs are wrong by one "0". Fixing this so they can be matched to SEdb sample list later.
SE_Package.csv$sample.ID <- sub("(.{10})(.*)", "\\10\\2", SE_Package.csv$sample.ID)

## Make GRanges object from SE data 
SE_Package.gr <- makeGRangesFromDataFrame(SE_Package.csv, na.rm = TRUE)

## Make GRanges List of SE data split into groups by Sample.id 
SE_Package.grl <- makeGRangesListFromDataFrame(SE_Package.csv, split.field = "sample.ID", na.rm = TRUE)

# Compare probes in each category with ALL SE locations 

## Make empty matrix for results 
SE_All_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))
rownames(SE_All_Matrix) <- names(EPICv2_probes.grl)
colnames(SE_All_Matrix) <- c("uniq_probe_n", "uniq_probe_perc", "uniq_SE_n", "uniq_SE_perc")

## Make empty list for probe tally 
probetallylist <- list()
uniq_SE_locations_list <- list()

for (j in 1:length(EPICv2_probes.grl)) {
  ol <- findOverlaps(SE_Package.gr, EPICv2_probes.grl[[j]])
  # Find number of unique probes overlapping SE locations 
  uniq_probe_n <- length(unique(ol@to))
  # Find % of unique probes overlapping SE locations 
  uniq_probe_perc <- (length(unique(ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique SEs targeted by probes 
  uniq_SE_n <- length(unique(ol@from))
  # Find % of unique SEs targeted by probes 
  uniq_SE_perc <- (length(unique(ol@from))/length(SE_Package.gr))*100
  # Find unique SE LOCATIONS that are targeted by the probes 
  uniq_SE_locations <- unique(ol@from)
  uniq_SE_locations_list[[j]] <- uniq_SE_locations
  # Summarise the number of probes per SE 
  probetally <- table(table(ol@from))
  SE_All_Matrix[j,] <- c(uniq_probe_n, uniq_probe_perc, uniq_SE_n, uniq_SE_perc)
  probetallylist[[j]] <- probetally
}

names(probetallylist) <- names(EPICv2_probes.grl)
names(uniq_SE_locations_list) <- names(EPICv2_probes.grl)
SE_All.df <- as.data.frame(SE_All_Matrix)

# Change first row to be actual column to allow for plotting 
SE_All.df <- rownames_to_column(SE_All.df, var = "Probes")

# Write csv for counts and percentages for all samples as summary 
write.csv(SE_All.df, file.path(filepath, "results/SE_Overlap_Summary.csv"))

# # Plot unique SEs targeted by probes as %
# All_Uniq_SE <- ggplot(SE_All.df, aes(x = Probes, y = uniq_SE_perc, fill = Probes)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   geom_text(aes(label=uniq_SE_n), vjust=-0.3, size=4) +
#   ggtitle("Unique Super-Enhancers Targeted by Probes") +
#   ylab("Percentage of Super-Enhancers") 
# All_Uniq_SE
# ggsave(file.path(filepath, "plots/Unique_SE_AllProbes.pdf"))
# 
# # Plot percent of unique probes that overlap with SEs
# All_Uniq_ProbeOL <- ggplot(SE_All.df, aes(x = Probes, y = uniq_probe_perc, fill = Probes)) +
#   geom_bar(stat = "identity", show.legend = FALSE) +
#   geom_text(aes(label=uniq_probe_n), vjust=-0.3, size=4) +
#   ggtitle("Unique EPICv2 Probes Overlapping Super-Enhancers") +
#   ylab("Percentage of Probes")
# All_Uniq_ProbeOL
# ggsave(file.path(filepath, "plots/Unique_Probes_SEOL.pdf"))


####subset super-enhancers by TSS, Genebody and Intergenic regions, then repeat overlap analysis using this
INTER_TSS <- intersect(SE_Package.gr, Locations_TSS.gr)

INTER_G <- intersect(SE_Package.gr, Locations_GeneBody.gr)

INTER_Int <- intersect(SE_Package.gr, Locations_Intergenic.gr)

# Find number of unique probes overlapping SE TSS locations 
## Make empty list for probe tally 
## Make empty matrix for results 
SE_TSS_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))
rownames(SE_TSS_Matrix) <- names(EPICv2_probes.grl)
colnames(SE_TSS_Matrix) <- c("uniq_probe_n", "uniq_probe_perc", "uniq_SE_n", "uniq_SE_perc")

probetallylistTSS <- list()
uniq_SE_locations_listTSS <- list()

for (j in 1:length(EPICv2_probes.grl)) {
  ol <- findOverlaps(INTER_TSS, EPICv2_probes.grl[[j]])
  # Find number of unique probes overlapping SE locations 
  uniq_probe_n <- length(unique(ol@to))
  # Find % of unique probes overlapping SE locations 
  uniq_probe_perc <- (length(unique(ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique SEs targeted by probes 
  uniq_SE_n <- length(unique(ol@from))
  # Find % of unique SEs targeted by probes 
  uniq_SE_perc <- (length(unique(ol@from))/length(INTER_TSS ))*100
  # Find unique SE LOCATIONS that are targeted by the probes 
  uniq_SE_locations <- unique(ol@from)
  uniq_SE_locations_listTSS[[j]] <- uniq_SE_locations
  # Summarise the number of probes per SE 
  probetally <- table(table(ol@from))
  SE_TSS_Matrix[j,] <- c(uniq_probe_n, uniq_probe_perc, uniq_SE_n, uniq_SE_perc)
  probetallylistTSS[[j]] <- probetally
}

# Find number of unique probes overlapping SE GB locations 
## Make empty list for probe tally 
## Make empty matrix for results 
SE_GB_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))
rownames(SE_GB_Matrix) <- names(EPICv2_probes.grl)
colnames(SE_GB_Matrix) <- c("uniq_probe_n", "uniq_probe_perc", "uniq_SE_n", "uniq_SE_perc")

probetallylistGB <- list()
uniq_SE_locations_listGB <- list()

for (j in 1:length(EPICv2_probes.grl)) {
  ol <- findOverlaps(INTER_G, EPICv2_probes.grl[[j]])
  # Find number of unique probes overlapping SE locations 
  uniq_probe_n <- length(unique(ol@to))
  # Find % of unique probes overlapping SE locations 
  uniq_probe_perc <- (length(unique(ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique SEs targeted by probes 
  uniq_SE_n <- length(unique(ol@from))
  # Find % of unique SEs targeted by probes 
  uniq_SE_perc <- (length(unique(ol@from))/length(INTER_G))*100
  # Find unique SE LOCATIONS that are targeted by the probes 
  uniq_SE_locations <- unique(ol@from)
  uniq_SE_locations_listGB[[j]] <- uniq_SE_locations
  # Summarise the number of probes per SE 
  probetally <- table(table(ol@from))
  SE_GB_Matrix[j,] <- c(uniq_probe_n, uniq_probe_perc, uniq_SE_n, uniq_SE_perc)
  probetallylistGB[[j]] <- probetally
}

# Find number of unique probes overlapping SE Intergenic locations 
## Make empty list for probe tally 
## Make empty matrix for results 
SE_Int_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))
rownames(SE_Int_Matrix) <- names(EPICv2_probes.grl)
colnames(SE_Int_Matrix) <- c("uniq_probe_n", "uniq_probe_perc", "uniq_SE_n", "uniq_SE_perc")

probetallylistInt <- list()
uniq_SE_locations_listInt <- list()

for (j in 1:length(EPICv2_probes.grl)) {
  ol <- findOverlaps(INTER_Int, EPICv2_probes.grl[[j]])
  # Find number of unique probes overlapping SE locations 
  uniq_probe_n <- length(unique(ol@to))
  # Find % of unique probes overlapping SE locations 
  uniq_probe_perc <- (length(unique(ol@to))/length(EPICv2_probes.grl[[j]]))*100
  # Find number of unique SEs targeted by probes 
  uniq_SE_n <- length(unique(ol@from))
  # Find % of unique SEs targeted by probes 
  uniq_SE_perc <- (length(unique(ol@from))/length(INTER_G))*100
  # Find unique SE LOCATIONS that are targeted by the probes 
  uniq_SE_locations <- unique(ol@from)
  uniq_SE_locations_listInt[[j]] <- uniq_SE_locations
  # Summarise the number of probes per SE 
  probetally <- table(table(ol@from))
  SE_Int_Matrix[j,] <- c(uniq_probe_n, uniq_probe_perc, uniq_SE_n, uniq_SE_perc)
  probetallylistInt[[j]] <- probetally
}

write.csv(SE_TSS_Matrix, file.path(pwd, "results/SE_TSS_Overlap_Summary.csv"))
write.csv(SE_GB_Matrix, file.path(pwd, "results/SE_GeneBody_Overlap_Summary.csv"))
write.csv(SE_Int_Matrix, file.path(pwd, "results/SE_Intergenic_Overlap_Summary.csv"))

#***# 

# ENHANCER ANNOTATION ANALYSIS 

## Read in bed file of enhancer locations, downloaded from https://fantom.gsc.riken.jp/5/datafiles/reprocessed/hg38_latest/extra/enhancer/ on 12 April 2023
Enhancer_locations <- read.delim(file.path(filepath, "metadata/F5_enhancers.bed"), header = FALSE)

## Make GRanges object from enhancer data 
Enhancer_locations.gr <- makeGRangesFromDataFrame(Enhancer_locations, keep.extra.columns=TRUE, ignore.strand=FALSE, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# Compare probes in each category with ALL Enhancer locations 

## Make empty matrix for results 
Enhancer_Matrix <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))
rownames(Enhancer_Matrix) <- names(EPICv2_probes.grl)
colnames(Enhancer_Matrix) <- c("Unique_Probe_No", "Unique_Probe_Percent", "Unique_Enh_No", "Unique_Enh_Perc")

## Make empty list for probe tally 
Enh_probetallylist <- list()
uniq_Enh_locations_list <- list()

for (i in 1:length(EPICv2_probes.grl)) {
  ol <- findOverlaps(Enhancer_locations.gr, EPICv2_probes.grl[[i]])
  # Find number of unique probes overlapping SE locations 
  Unique_Probe_No <- length(unique(ol@to))
  # Find % of unique probes overlapping SE locations 
  Unique_Probe_Percent <- (length(unique(ol@to))/length(EPICv2_probes.grl[[i]]))*100
  # Find number of unique SEs targeted by probes 
  Unique_Enh_No <- length(unique(ol@from))
  # Find % of unique SEs targeted by probes 
  Unique_Enh_Perc <- (length(unique(ol@from))/length(Enhancer_locations.gr))*100
  # Find unique enhancer LOCATIONS that are targeted by the probes 
  uniq_Enh_locations <- unique(ol@from)
  uniq_Enh_locations_list[[i]] <- uniq_Enh_locations
  # Summarise the number of probes per SE 
  Enh_probetally <- table(table(ol@from))
  Enhancer_Matrix[i,] <- c(Unique_Probe_No, Unique_Probe_Percent, Unique_Enh_No, Unique_Enh_Perc)
  Enh_probetallylist[[i]] <- Enh_probetally
}

names(Enh_probetallylist) <- names(EPICv2_probes.grl)
names(uniq_Enh_locations_list) <- names(EPICv2_probes.grl)
Enhancer_Sum.df <- as.data.frame(Enhancer_Matrix)

# Change first row to be actual column to allow for plotting 
Enhancer_Sum.df <- rownames_to_column(Enhancer_Sum.df, var = "Probes")

# Save out counts & percentages for Enhancer overlap summary 
write.csv(Enhancer_Sum.df, file.path(filepath, "results/Enhancer_Overlap_Summary.csv"))

# Plot unique Enhancers targeted by probes as %
Unique_Enh <- ggplot(Enhancer_Sum.df, aes(x = Probes, y = Unique_Enh_Perc, fill = Probes)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label=Unique_Enh_No), vjust=-0.3, size=3.5) +
  ggtitle("Unique Enhancers Targeted by Probes") +
  ylab("Percentage of Enhancers") 
Unique_Enh
ggsave(file.path(filepath, "plots/Unique_Enhancers_Probes.pdf"))

# Plot percent of unique probes that overlap with Enhancers
Unique_Probe_EnhOL <- ggplot(Enhancer_Sum.df, aes(x = Probes, y = Unique_Probe_Percent, fill = Probes)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label=Unique_Probe_No), vjust=-0.3, size=3.5) +
  ggtitle("Unique EPICv2 Probes Overlapping Enhancers") +
  ylab("Percentage of Probes")
Unique_Probe_EnhOL
ggsave(file.path(filepath, "plots/Unique_Probes_EnhOL.pdf"))


#subset super-enhancers by enhancer regions, then repeat overlap analysis using this
INTER_Enh <- intersect(SE_Package.gr, Enhancer_locations.gr)

## Make empty matrix for results 
Enhancer_MatrixSE <- matrix(data = NA, ncol = 4, nrow = length(EPICv2_probes.grl))
rownames(Enhancer_MatrixSE) <- names(EPICv2_probes.grl)
colnames(Enhancer_MatrixSE) <- c("Unique_Probe_No", "Unique_Probe_Percent", "Unique_Enh_No", "Unique_Enh_Perc")

## Make empty list for probe tally 
Enh_probetallylistSE <- list()
uniq_Enh_locations_listSE <- list()

for (i in 1:length(EPICv2_probes.grl)) {
  ol <- findOverlaps(INTER_Enh, EPICv2_probes.grl[[i]])
  # Find number of unique probes overlapping SE locations 
  Unique_Probe_No <- length(unique(ol@to))
  # Find % of unique probes overlapping SE locations 
  Unique_Probe_Percent <- (length(unique(ol@to))/length(EPICv2_probes.grl[[i]]))*100
  # Find number of unique SEs targeted by probes 
  Unique_Enh_No <- length(unique(ol@from))
  # Find % of unique SEs targeted by probes 
  Unique_Enh_Perc <- (length(unique(ol@from))/length(INTER_Enh))*100
  # Find unique enhancer LOCATIONS that are targeted by the probes 
  uniq_Enh_locations <- unique(ol@from)
  uniq_Enh_locations_listSE[[i]] <- uniq_Enh_locations
  # Summarise the number of probes per SE 
  Enh_probetally <- table(table(ol@from))
  Enhancer_MatrixSE[i,] <- c(Unique_Probe_No, Unique_Probe_Percent, Unique_Enh_No, Unique_Enh_Perc)
  Enh_probetallylistSE[[i]] <- Enh_probetally
}

write.csv(Enhancer_MatrixSE, file.path(pwd, "results/SE_Enh_Overlap_Summary.csv"))


## LOR between percentages of probe tally density 
# create a blank list 
blah <- sort(as.numeric(unique(c(names(Enh_probetallylist[[5]]), names(Enh_probetallylist[[6]])))))

v1v2_enh_ratio <- matrix(data = NA, nrow = length(blah), ncol = 2)

rownames(v1v2_enh_ratio) <- blah

for (i in 1:length(blah)) {
  k <- which(names(Enh_probetallylist[[5]]) == blah[i])
  if(length(k)>0){
    v1v2_enh_ratio[i,1] <- Enh_probetallylist[[5]][k]
  } else {
    v1v2_enh_ratio[i,1] <- 0
  }
  l <- which(names(Enh_probetallylist[[6]]) == blah[i])
  if(length(l)>0){
    v1v2_enh_ratio[i,2] <- Enh_probetallylist[[6]][l] 
  } else {
    v1v2_enh_ratio[i,2] <- 0
  }
  
}

colnames(v1v2_enh_ratio) <- c("EPICv2", "EPICv1")

diff <- v1v2_enh_ratio[,1]-v1v2_enh_ratio[,2]

pdf(file = file.path(filepath, "plots/Enh_Ratio_ProbeTallyCounts_Density.pdf"))
diff_plot <- barplot(diff, col = "pink", main = "Ratio of Enh Probe Density Coverage between v2 and v1", xlab = "Probe Number per Enhancer", ylab = "Ratio between probe coverage counts v2 vs v1", ylim = c(-100, 1000))
dev.off()

# Convert to % in terms of probe category total 
v1v2_enh_ratio[,1] <- (v1v2_enh_ratio[,1]/sum(v1v2_enh_ratio[,1]))*100

v1v2_enh_ratio[,2] <- (v1v2_enh_ratio[,2]/sum(v1v2_enh_ratio[,2]))*100

diff <- v1v2_enh_ratio[,1]-v1v2_enh_ratio[,2]

pdf(file = file.path(filepath, "plots/Enh_Ratio_ProbeTallyPercent_Density.pdf"))
diff_plot <- barplot(diff, col = "pink", main = "Ratio of Enh Probe Density Coverage between v2 and v1", xlab = "Probe Number per Enhancer", ylab = "Ratio between % probe coverage v2 vs v1")
dev.off()

# Plot probetally list to get an idea of density of coverage (proportionate to the number of probes in each category)
# Set format of plots and set up so it saves out automatically 

pdf(file = file.path(filepath,"plots/Enh_ProbeTally_Percent.pdf"))
par(mfrow=c(3,2))

# Only plotting the ten top values for number of probes in each SE. This gives a snapshot of the density and if there has been a significant increase in the number of probes added that cover SEs. 
# bars are percentage but text is raw counts 
Enh_Tally_Reinstated <- barplot((Enh_probetallylist[[1]]/sum(Enh_probetallylist[[1]])*100)[1:10], col = "pink", main = "Reinstated Probes", xlab = "Probe Number per Enhancer", ylab = "Percentage Probe Overlap", ylim = c(0,100))
text(Enh_Tally_Reinstated, (Enh_probetallylist[[1]]/sum(Enh_probetallylist[[1]])*100)[1:10], labels = Enh_probetallylist[[1]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_Excluded <- barplot((Enh_probetallylist[[2]]/sum(Enh_probetallylist[[2]])*100)[1:10], col = "pink", main = "Excluded Probes", xlab = "Probe Number per Enhancer", ylab = "Percentage Probe Overlap", ylim = c(0,100))
text(Enh_Tally_Excluded, (Enh_probetallylist[[2]]/sum(Enh_probetallylist[[2]])*100)[1:10], labels = Enh_probetallylist[[2]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_New <- barplot((Enh_probetallylist[[3]]/sum(Enh_probetallylist[[3]])*100)[1:10], col = "pink", main = "New Probes", xlab = "Probe Number per Enhancer", ylab = "Percentage Probe Overlap", ylim = c(0,100))
text(Enh_Tally_New, (Enh_probetallylist[[3]]/sum(Enh_probetallylist[[3]])*100)[1:10], labels = Enh_probetallylist[[3]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_Retained <- barplot((Enh_probetallylist[[4]]/sum(Enh_probetallylist[[4]])*100)[1:10], col = "pink", main = "Retained Probes", xlab = "Probe Number per Enhancer", ylab = "Percentage Probe Overlap", ylim = c(0,100))
text(Enh_Tally_Retained, (Enh_probetallylist[[4]]/sum(Enh_probetallylist[[4]])*100)[1:10], labels = Enh_probetallylist[[4]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_All <- barplot((Enh_probetallylist[[5]]/sum(Enh_probetallylist[[5]])*100)[1:10], col = "pink", main = "EPICv2 Probes", xlab = "Probe Number per Enhancer", ylab = "Percentage Probe Overlap", ylim = c(0,100))
text(Enh_Tally_All, (Enh_probetallylist[[5]]/sum(Enh_probetallylist[[5]])*100)[1:10], labels = Enh_probetallylist[[5]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_EPICv1 <- barplot((Enh_probetallylist[[6]]/sum(Enh_probetallylist[[6]])*100)[1:10], col = "pink", main = "EPICv1 Probes", xlab = "Probe Number per Enhancer", ylab = "Percentage Probe Overlap", ylim = c(0,100))
text(Enh_Tally_EPICv1, (Enh_probetallylist[[6]]/sum(Enh_probetallylist[[6]])*100)[1:10], labels = Enh_probetallylist[[6]][1:10], cex = 0.8, pos = 3, xpd = TRUE)

dev.off()

# Set format of plots and set up so it saves out automatically 

pdf(file = file.path(filepath, "plots/Enh_ProbeTally_Counts.pdf"))
par(mfrow=c(3,2))

# Only plotting the ten top values for number of probes in each SE. This gives a snapshot of the density and if there has been a significant increase in the number of probes added that cover SEs. 
Enh_Tally_Reinstated <- barplot((Enh_probetallylist[[1]])[1:10], col = "pink", main = "Reinstated Probes", xlab = "Probe Number per Enhancer", ylab = "Number Probe Overlap", ylim = c(0,25000))
text(Enh_Tally_Reinstated, Enh_probetallylist[[1]][1:10], labels = Enh_probetallylist[[1]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_Excluded <- barplot((Enh_probetallylist[[2]])[1:10], col = "pink", main = "Excluded Probes", xlab = "Probe Number per Enhancer", ylab = "Number Probe Overlap", ylim = c(0,25000))
text(Enh_Tally_Excluded, Enh_probetallylist[[2]][1:10], labels = Enh_probetallylist[[2]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_New <- barplot((Enh_probetallylist[[3]])[1:10], col = "pink", main = "New Probes", xlab = "Probe Number per Enhancer", ylab = "Number Probe Overlap", ylim = c(0,25000))
text(Enh_Tally_New, Enh_probetallylist[[3]][1:10], labels = Enh_probetallylist[[3]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_Retained <- barplot((Enh_probetallylist[[4]])[1:10], col = "pink", main = "Retained Probes", xlab = "Probe Number per Enhancer", ylab = "Number Probe Overlap", ylim = c(0,25000))
text(Enh_Tally_Retained, Enh_probetallylist[[4]][1:10], labels = Enh_probetallylist[[4]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_All <- barplot((Enh_probetallylist[[5]])[1:10], col = "pink", main = "EPICv2 Probes", xlab = "Probe Number per Enhancer", ylab = "Number Probe Overlap", ylim = c(0,25000))
text(Enh_Tally_All, Enh_probetallylist[[5]][1:10], labels = Enh_probetallylist[[5]][1:10], cex = 0.8, pos = 3, xpd = TRUE)
Enh_Tally_EPICv1 <- barplot((Enh_probetallylist[[6]])[1:10], col = "pink", main = "EPICv1 Probes", xlab = "Probe Number per Enhancer", ylab = "Number Probe Overlap", ylim = c(0,25000))
text(Enh_Tally_EPICv1, Enh_probetallylist[[6]][1:10], labels = Enh_probetallylist[[6]][1:10], cex = 0.8, pos = 3, xpd = TRUE)

dev.off()



#Find regions for plotting prostate cancer enhancer example

#Enh_Usage.matrix downloaded from https://fantom.gsc.riken.jp/5/datafiles/latest/extra/Enhancers/  

Enh_usage <- read.delim(file.path(filepath, "metadata/Enh_Usage.matrix"), header = TRUE)

# Subset out Prostate Epithelial Cell Samples (based upon enhancer IDs in 'Enh_Samples')
Prostate_Epi_Enh <- subset(Enh_usage, select = c("CNhs10882", "CNhs11972", "CNhs12014"))

# Remove enhancer locations that arent present in any of the three prostate epi cell samples 
Prostate_Epi_Enh <- Prostate_Epi_Enh[rowSums(Prostate_Epi_Enh) > 0,]

# Make enhancer location a column instead of rowname for downstream analysis 
Prostate_Epi_Enh <- rownames_to_column(Prostate_Epi_Enh, var = "seqnames")

# Reformat the enhancer location format so that we can split columns later in order to make a GRanges 
Prostate_Epi_Enh <- as.data.frame(gsub("-", " ", Prostate_Epi_Enh$seqnames))
Prostate_Epi_Enh <- as.data.frame(gsub(":", " ", Prostate_Epi_Enh$`gsub("-", " ", Prostate_Epi_Enh$seqnames)`))
colnames(Prostate_Epi_Enh) <- "seqnames"

# Split column into three columns for chr, start and end and make GRanges 
Prostate_Epi_Enh[c('chr', 'start', "end")] <- str_split_fixed(Prostate_Epi_Enh$seqnames, ' ', 3)
write.table(Prostate_Epi_Enh[,2:4],paste(filepath,"results/PrEnh_all3.bed",sep""), row.names = FALSE, quote = FALSE, col.names = FALSE)
Prostate_Epi_Enh.gr <- makeGRangesFromDataFrame(Prostate_Epi_Enh, seqnames.field = "chr", start.field = "start", end.field = "end")

# Count overlaps between prostate epi enhancers and EPICv1 probes 
Prostate_Enh_v1 <- as.data.frame(countOverlaps(Prostate_Epi_Enh.gr, EPICv2_probes.grl$EPICv1))
length(which(Prostate_Enh_v1>0))#2144/3145
#Prostate_Enh_v1 <- as.data.frame(Prostate_Enh_v1[!(Prostate_Enh_v1$`countOverlaps(Prostate_Epi_Enh.gr, Locations_EPICv1.gr)`==0),])

# Count overlaps between prostate epi enhancers and EPICv1 probes 
Prostate_Enh_v2 <- as.data.frame(countOverlaps(Prostate_Epi_Enh.gr, EPICv2_probes.grl$EPICv2))
length(which(Prostate_Enh_v2>0))#2205/3145
#Prostate_Enh_v2 <- as.data.frame(Prostate_Enh_v2[!(Prostate_Enh_v2$`countOverlaps(Prostate_Epi_Enh.gr, Locations_all.gr)`==0),])

PE3<-cbind(Prostate_Enh_v1,Prostate_Enh_v2)

unique_enh_v2_test<-which(PE3[,1]==0 & PE3[,2]>0)
length(unique_enh_v2_test)#266

# Make GRanges for probes that are new and reinstated into EPICv2 
Probes_New_Reinstated <- rbind(Locations_new, Locations_reinstated)
Probes_New_Reinstated.gr <- makeGRangesFromDataFrame(Probes_New_Reinstated, seqnames.field = "seqnames", start.field = "start", end.field = "end")

# Count overlaps between prostate epi enhancers and new and reinstated probes to EPICv2. This should be probes that are only in v2 and NOT in v1 
Prostate_Epi_Enh_only_v2 <- as.data.frame(countOverlaps(Prostate_Epi_Enh.gr, Probes_New_Reinstated.gr))
Prostate_Epi_Enh_only_v2 <- as.data.frame(Prostate_Epi_Enh_only_v2[!(Prostate_Epi_Enh_only_v2$`countOverlaps(Prostate_Epi_Enh.gr, Probes_New_Reinstated.gr)`==0),])

# Count overlaps between prostate epi enhancers and retained probes 
Prostate_Epi_Enh_Retained <- as.data.frame(countOverlaps(Prostate_Epi_Enh.gr, Locations_retained.gr))
Prostate_Epi_Enh_Retained <- as.data.frame(Prostate_Epi_Enh_Retained[!(Prostate_Epi_Enh_Retained$`countOverlaps(Prostate_Epi_Enh.gr, Locations_retained.gr)`==0),])

# Count overlaps between prostate epi enhancers and new EPICv2 probes 
Prostate_Epi_Enh_New <- as.data.frame(countOverlaps(Prostate_Epi_Enh.gr, Locations_new.gr))
Prostate_Epi_Enh_New <- as.data.frame(Prostate_Epi_Enh_New[!(Prostate_Epi_Enh_New$`countOverlaps(Prostate_Epi_Enh.gr, Locations_new.gr)`==0),])


