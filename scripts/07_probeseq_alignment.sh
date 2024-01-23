#!/bin/bash -e

cd ~/scratch
mkdir EPIC_V2
cd EPIC_V2
curl -o hg38.fa.gz 'http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
gunzip hg38.fa.gz

#All CpGs temporarily become ZpGs for methylated genomes

biosed hg38.fa -t CG -r ZG tmp1.fasta
biosed tmp1.fasta -t C -r T tmp2.fasta
biosed tmp2.fasta -t Z -r C forward_meth.fa
rm tmp*.fasta

biosed hg38.fa -t C -r T forward_unmeth.fa

revseq hg38.fa -reverse -complement -outseq hg38rev.fa

biosed hg38rev.fa -t CG -r ZG tmp1.fasta
biosed tmp1.fasta -t C -r T tmp2.fasta
biosed tmp2.fasta -t Z -r C reverse_meth.fa
rm tmp*.fasta

biosed hg38rev.fa -t C -r T reverse_unmeth.fa

#Chromosome sizes

cat hg38.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > hg38.chrom.sizes

#Generate probe sequences from Manifest


Rscript getprobeseqs.R

#BLAT

#In serial
blat forward_meth.fa allEPICseqs.fasta -stepSize=5 -repMatch=2147483647 -minScore=0 -minIdentity=0 -maxIntron=0 blat_results_forward_meth.psl &
blat reverse_meth.fa allEPICseqs.fasta -stepSize=5 -repMatch=2147483647 -minScore=0 -minIdentity=0 -maxIntron=0 blat_results_reverse_meth.psl &
blat forward_unmeth.fa allEPICseqs.fasta -stepSize=5 -repMatch=2147483647 -minScore=0 -minIdentity=0 -maxIntron=0 blat_results_forward_unmeth.psl &
blat reverse_unmeth.fa allEPICseqs.fasta -stepSize=5 -repMatch=2147483647 -minScore=0 -minIdentity=0 -maxIntron=0 blat_results_reverse_unmeth.psl &

mkdir ~/scratch/EPIC_V2/filtered

# Filter for at least 47 bp matching

for file in *.psl; do
  sample="${file%.psl}"
  awk '(NR>5) && ($1 > 46) && ($13 == 50)' $sample.psl > filtered/$sample.psl
  echo "$sample"
done

#Generate alignment scores

for file in filtered/*.psl; do
  sample="${file%.psl}"
  echo "$sample"
  pslScore $sample.psl > $sample.score
done

cd ..

Rscript processBLATresults.R












