library(data.table)

options(stringsAsFactors=FALSE, scipen=99)

#---- Set Working Directory that gtf file downloaded into
setwd("your/file/path")

#--- Read in GTF file
data <- fread("gunzip -c gencode.v25.annotation.gtf.gz", sep="\t", header=FALSE)
colnames(data) <- c("chr", "source", "feature_type", "start", "end", "score", "strand", "phase", "metadata")
setkey(data, feature_type)

#--- Parse out metadata associated with transcripts
text <- data["transcript", metadata]
# Transcript ID
pattern <-  "transcript_id \"([[:alnum:]_\\.]+)\""
temp <- regmatches(text, regexec(pattern, text))
transcript_id <- sapply(temp, function(x) x[2])
# Transcript type
pattern <-  "transcript_type \"([[:alnum:]_]+)\""
temp <- regmatches(text, regexec(pattern, text))
transcript_type <- sapply(temp, function(x) x[2])
# Gene ID
pattern <-  "gene_id \"([[:alnum:]_\\.]+)\""
temp <- regmatches(text, regexec(pattern, text))
gene_id <- sapply(temp, function(x) x[2])
# Gene type
pattern <-  "gene_type \"([[:alnum:]_]+)\""
temp <- regmatches(text, regexec(pattern, text))
gene_type <- sapply(temp, function(x) x[2])
# Gene name
pattern <-  "gene_name \"([[:alnum:]_\\.-]+)\""
temp <- regmatches(text, regexec(pattern, text))
gene_name <- sapply(temp, function(x) x[2])
# Level
pattern <-  "level ([0-9])"
temp <- regmatches(text, regexec(pattern, text))
level <- sapply(temp, function(x) x[2])

#--- Write out metadata file
rv <- data.frame(transcript_id=transcript_id,
                 transcript_type=transcript_type,
                 gene_id=gene_id,
                 gene_name=gene_name,
                 gene_type=gene_type,
                 level=level)
connection <- gzfile("transcript_metadata.txt.gz", "w")
write.table(rv, file=connection, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
close(connection)
