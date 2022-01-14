# Script_wrappers

1) fumiaki_snp_pipeline.sh
  Bash script to call SNPs for Mike Grigg's group. It was adapted from a script provided by Fumiaki
  
2) vcf2tbls.sh
 Bash script to convert vcf files into a table with genotypes (GT field only)
 
3) run_tbl2fasta_loop.sh
 Bash wrapper script to run python script tbl2fasta.py in a loop to convert the output from (2) above into a fasta file that can be used as input for Split-Tree.

4) tbl2fasta.py
 Python script to convert the output from (2) above into a fasta file that can be used as input for Split-Tree.
