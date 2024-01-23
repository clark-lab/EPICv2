library(GenomicFeatures)
library(GenomicRanges)
library(data.table)

options(stringsAsFactors=FALSE, scipen=99)

#---- Set Working Directory
setwd("your/file/path")

# #---- Creat TxDb object from GENCODE GTF file
# chrominfo <- read.table("../../Sequences/hg38.chrom.sizes", sep="\t", header=FALSE)
# colnames(chrominfo) <- c("chrom", "length"); chrominfo$is_circular <- FALSE
# 
# txdb <- makeTxDbFromGFF(file=gzfile("./gencode.v25.annotation.gtf.gz"),
#                         format="gtf", chrominfo=chrominfo,
#                         dataSource="GENCODE V25",
#                         organism="Homo Sapiens")
#                                 
# saveDb(txdb, "gencode.v25.txdb.sqlite")

txdb <- loadDb("gencode.v25.txdb.sqlite")

#--- Load transcript metadate
metadata <- fread("gunzip -c transcript_metadata.txt.gz", sep="\t", header=TRUE)
setkey(metadata, transcript_type, level)

#--- TSS centered tables
txs <- transcripts(txdb, columns=c("tx_name"))
txs_protein <- txs[txs$tx_name %in% metadata[.("protein_coding", unique(level))]$transcript_id]
txs_lincrna <- txs[txs$tx_name %in% metadata[.("lincRNA", unique(level))]$transcript_id]

# txs_manual <- txs[txs$tx_name %in% metadata[.(unique(transcript_type), c(1,2))]$transcript_id]
# txs_manual_protein <- txs[txs$tx_name %in% metadata[.("protein_coding", c(1,2))]$transcript_id]
# txs_manual_lincrna <- txs[txs$tx_name %in% metadata[.("lincRNA", c(1,2))]$transcript_id]

# Transcripts
regions <- sort(txs, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions))-1,
             as.numeric(end(regions)),
             as.character(regions$tx_name),
             0,
             as.character(strand(regions)))
connection <- gzfile("transcript.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)


# Transcripts extended +/-2kb
regions <- trim(resize(resize(txs, width=width(txs)+2000, fix="start"), width=width(txs)+4000, fix="end"))
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions))-1,
             as.numeric(end(regions)),
             as.character(regions$tx_name),
             0,
             as.character(strand(regions)))
connection <- gzfile("transcript_2kb.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

# Intergenic regions genome minus transcripts +/-2kb
regions_genic <- trim(reduce(resize(resize(txs, width=width(txs)+2000, fix="start"), width=width(txs)+4000, fix="end")))
regions <- gaps(regions_genic); 
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions))-1,
             as.numeric(end(regions)),
             ".",
             0,
             as.character(strand(regions)))
connection <- gzfile("intergenic.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

# TSSs extended +/-2kb
regions <- trim(resize(resize(txs, 2000, fix="start"), 4000, fix="end"))
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions))-1,
             as.numeric(end(regions)),
             as.character(regions$tx_name),
             0,
             as.character(strand(regions)))
connection <- gzfile("tss_2kb.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

# TSSs extended +/-2kb protein coding
regions <- trim(resize(resize(txs_protein, 2000, fix="start"), 4000, fix="end"))
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions))-1,
             as.numeric(end(regions)),
             as.character(regions$tx_name),
             0,
             as.character(strand(regions)))
connection <- gzfile("tss_2kb.protein.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

# TSSs extended +/-2kb lincRNAs
regions <- trim(resize(resize(txs_lincrna, 2000, fix="start"), 4000, fix="end"))
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions))-1,
             as.numeric(end(regions)),
             as.character(regions$tx_name),
             0,
             as.character(strand(regions)))
connection <- gzfile("tss_2kb.lincrna.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

# TSSs 
# +1 Nucleosome (0bp, 250bp).
regions <- trim(resize(txs, 250, fix="start"))
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions))-1,
             as.numeric(end(regions)),
             as.character(regions$tx_name),
             0,
             as.character(strand(regions)))
connection <- gzfile("nucleosome_plusone.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

#--- Intron centered tables
introns <- intronsByTranscript(txdb, use.names=TRUE)
regions <- unlist(introns, use.names=TRUE)
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions)),
             as.numeric(end(regions)),
             names(regions),
             0,
             as.character(strand(regions)))
connection <- gzfile("introns.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

regions <- reduce(regions, drop.empty.ranges=TRUE)
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions)),
             as.numeric(end(regions)),
             ".",
             0,
             as.character(strand(regions)))
connection <- gzfile("introns_merged.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

#--- Exon centered tables
exons <- exonsBy(txdb, by="tx", use.names=TRUE)
exons <- unlist(exons, use.names=TRUE)
regions <- unlist(exons, use.names=TRUE)
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions)),
             as.numeric(end(regions)),
             names(regions),
             0,
             as.character(strand(regions)))
connection <- gzfile("exons.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)

regions <- reduce(regions, drop.empty.ranges=TRUE)
regions <- sort(regions, ignore.strand=TRUE)
out <- cbind(as.character(seqnames(regions)),
             as.numeric(start(regions)),
             as.numeric(end(regions)),
             ".",
             0,
             as.character(strand(regions)))
connection <- gzfile("exons_merged.bed.gz", "w")
write.table(out, file=connection, sep="\t", col.names=FALSE, row.names=FALSE, quote=FALSE)
close(connection)
