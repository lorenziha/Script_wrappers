#!/usr/bin/bash

VCF=${1:-"NEW_BAMS_filtered_snps.PASS.vcf"}
PREFIX=${2:-"test"}

module load gatk/4.2.0.0-Java-1.8.0_92
gatk VariantsToTable -V ${VCF} -F CHROM -F POS -GF GT -O ${PREFIX}.table
module purge
