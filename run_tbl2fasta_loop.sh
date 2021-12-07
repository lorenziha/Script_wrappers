#!/usr/bin/bash

RED='\033[0;31m'
GREEN='\033[1;32m'
NC='\033[0m' # No color

for s in `cat optional_samples.prefix`
do
	echo -e "Sample ${RED}${s}${NC}"
	cp core_samples.prefix core_samples_plus_${s}.prefix
	echo ${s} >> core_samples_plus_${s}.prefix
	for i in {5..4} ; do echo -e "Coverage = ${GREEN}${i}${NC}"
	qsub  -cwd -N tb2fa -m be -M lorenziha@nih.gov ~/bin/tbl2fasta.py -t NEW_REF_filtered_snps.PASS_COV${i}.table -o  NEW_REF_filtered_snps.PASS_COV${i}.${s}.fasta -s core_samples_plus_${s}.prefix -m N -S 0.75 --skip_heterozygote_positions
	done
done
