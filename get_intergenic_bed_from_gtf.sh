#!/bin/bash
set -o errexit

# This script requires seqtk and bedtools to be installed
# To run 
# get_intergenic_bed_from_gtf.sh <genome_fasta> <gtf_annotation_file>
GENOME=$1
RAWGTF=$2
GTF=$(echo ${RAWGTF} | sed 's/.gtf//')
echo GTF=${GTF}

echo Fetching chromosome sizes
seqtk  comp ${GENOME} | cut -f 1,2 > chromSizes.txt
cat chromSizes.txt | sort -k1,1 -k2,2n > chromSizes_sorted.txt
awk 'OFS="\t" {print $1, "0", $2}' chromSizes.txt | sort -k1,1 -k2,2n > chromSizes.bed

echo Sorting ${GTF}.gtf
cat ${GTF}.gtf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k4,4n -k5,5n"}' > ${GTF}_sorted.gtf

echo Building intergenic bed file: ${GTF}_intergenic_sorted.bed
bedtools complement -i ${GTF}_sorted.gtf -g chromSizes_sorted.txt > ${GTF}_intergenic_sorted.bed

echo Building exon bed file: ${GTF}_exon_sorted.bed
awk -F "\t" 'OFS="\t" $1 ~ /^#/ {print $0;next} {if ($3 == "exon") print $1"\t"$4-1"\t"$5}' ${GTF}_sorted.gtf > ${GTF}_exon_sorted.bed

echo Building intron bed file: ${GTF}_intron_sorted.bed
bedtools complement -i <(cat ${GTF}_exon_sorted.bed ${GTF}_intergenic_sorted.bed | sort -k1,1 -k2,2n) -g chromSizes_sorted.txt > ${GTF}_intron_sorted.bed

echo Done!!
