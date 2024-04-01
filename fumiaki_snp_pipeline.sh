#!/usr/bin/bash

Help()
{
   # Display Help
   echo "This script runs Fumiaki's SNP pipeline."
   echo "Syntax: ${0} [-r|-R|-p|-o|-m|-g|-P|-w|-h]"
   echo "options:"
   echo "-r     Full path to reference fasta file."
   echo "-R     Path to directory with fastq reads."
   echo "-p     Prefix file, one prefix per line."
   echo "-o     output prefix."
   echo "-w     path to working directory [default: current dir]."
   echo "-m     perform mapping flag (default: skip mapping)."
   echo "-g     Run only gvcf generation."
   echo "-P     input Prefix (only with -g flag)."
   echo "-h     show this help."
}

while getopts "gmhr:R:p:o:w:P:" option; do
   case $option in
        h) # display Help
                Help
                exit;;
        r) REFERENCE=${OPTARG};;
        R) READS=${OPTARG};;
        p) PREFIX_FILE=${OPTARG};;
        o) OUTPUT=${OPTARG};;
	m) DO_MAPPING=true;;
	g) GVCF_ONLY=true;;
	w) WORKDIR=${OPTARG};;
	P) PREFIX=${OPTARG};;
        \?) # incorrect option
                echo
                echo "Error, Invalid option"
                echo
                Help
                exit;;
   esac
done

#set -o errexit

#!/bin/sh

# map and variant calling using bwa and gatk
#$ -N SNP_call

# Execute the script from the Current Working Directory
#$ -cwd

# Merge the output of the script, and any error messages generated to one file
#$ -j y

# Send the output of the script to a directory called 'UGE-output' in the current working directory (cwd)
if [ ! -d "UGE-output" ] 
then #Create output directory in case it does NOT exist
     mkdir UGE-output
fi
#$ -o UGE-output/

# Tell the job your memory requirements
#  #$ -l mem_free=100G,h_vmem=100G

# Send mail when the job is submitted, and when the job completes
#$ -m be

#  Specify an email address to use
#$ -M hernan.lorenzi@nih.gov
# -pe threaded 16  # alternatively use for mapping only

####Have to update version before run

module load BWA/0.7.5a-goolf-1.7.20
module load GATK/3.7-0-Java-1.8.0_92  
module load picard/2.22.7-Java-1.8.0_92
module load SAMtools/1.3-goolf-1.7.20-HTSlib-1.3

# These are all of the module paths to the pieces of GATK that I use to map and quality control the genome
# They're all variables so that I can easily update the script if the paths in Locus changes
CD="${EBROOTPICARD}/picard.jar CreateSequenceDictionary"
PC="${EBROOTPICARD}/picard.jar SortSam"
MD="${EBROOTPICARD}/picard.jar MarkDuplicates"
MR="${EBROOTPICARD}/picard.jar MergeSamFiles"
BI="${EBROOTPICARD}/picard.jar BuildBamIndex"
GA="$EBROOTGATK/GenomeAnalysisTK.jar"

# These are all of the variables you should need to change in a standard mapping BEFORE running the script
# Tgindex is the Toxo reference that I use to map against. I have multiple so I comment out all except the one I'm currently using
# tempdir is the folder with a lot of available space I leave open so java has enough space for all the temporary files it makes while doing calculations
# workingdir is the folder set that I want to do my work in and generate output to
# fastqdir is the folder where my fastq files live I do this so I don't always have to have them in the same folder
# strain_text_file is a text file with each file name created from your sequencing without the extension 'sam'.  The text file should have one ID per line.  This file will indicate with files to run through the processing and variant calling pipeline.

if [ -z ${WORKDIR+y} ]
then
        # -g flag off
        WORKDIR=`pwd -P`
fi

Tgindex=${REFERENCE}
Tgindex_out=/hpcdata/lpd/Fumiaki/Toxo_genomes/TgME49/ToxoDB-51_TgondiiME49_Genome
tempdir=/hpcdata/scratch/lorenziha/TMP
workingdir=${WORKDIR}
fastqdir=${READS}
strain_text_file=${PREFIX_FILE}

if [ ${GVCF_ONLY} ]; then
	if [ ${PREFIX} ]; then
		echo ${PREFIX} > ${PREFIX}.prefix
		strain_text_file=${PREFIX}.prefix	
	else
		echo ERROR, you must specify option -P with flag -g; echo
		exit 1
	fi
fi


## Step 0
## create index of genome for BWA and GATK
if [ $DO_MAPPING ]; then
	if [ ! -f ${Tgindex} ]; then
		echo; echo "Crating bwa index files for ${Tgindex}"; echo
		bwa index -a bwtsw $Tgindex
		samtools faidx $Tgindex
		echo "Done!!"
		echo
	fi

	## Check for GATK dict file for reference
	DICT=${gindex/fasta/dict}
	if [ ! -f ${DICT} ]; then
		echo; echo "Creating dict file for referece ${Tgindex}"; echo
		java -Xmx1G -jar $CD R=$Tgindex O=$Tgindex_out.dict ##picard CCreateSequenceDictionary
		echo "Done!"
		echo
	fi
fi

## Step 1
## map to ME49
## create a text file with each file name without the extension 'sam'.  The text file should have one ID per line.  This file will indicate with files to run through the processing and variant calling pipeline.
## You will need to ensure that the fastq files that go into your mapping are named such that this code can find them using the fastqdir and line names given in your text file combined with the extension before .fastq
cd $workingdir/

while read prefix 
do
	prefix_list+=(${prefix}) # stores prefix names for later
	if [ $DO_MAPPING ]; then
		echo; echo "Running bwa mem on ${prefix}_R1.fastq.gz and ${prefix}_R2.fastq.gz"; echo 
		bwa mem -R "@RG\tID:L\tSM:"$line"\tPL:illumina\tLB:lib1\tPU:unit1" \
		-t 16 \
		-M $Tgindex $fastqdir/${prefix}_R1.fastq.gz \
		$fastqdir/${prefix}_R2.fastq.gz > $workingdir/${prefix}.sam
		echo ${prefix}_BWA
	fi
done < ${strain_text_file}

## Step 2
## run pipeline to process sam file and call variants. 
## the parameters for calling variants are basically the default ones recommended by gatk for the current version 3.1.1

#####Caution!!! HAVE TO CHECK RECOMMENDED PARAMETERS FOR THE CURRRENT CERSION

## the pipeline runs in parallel in the HaplotypeCaller step, everything else is not being parallelized.  My recommendation to run this script for large number of samples would be to submit various jobs.  
## The last step includes a ploidy argument. Change it if your genome is not haploid (1)
while read prefix
do
	echo "Calling variants (gvcf) for  ${prefix} ...." 
	if [ ! -f ${prefix}.sorted.bam ]; then
		echo; echo "Running ${PC} on ${prefix}.sam"; echo 
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $PC INPUT=${prefix}.sam OUTPUT=${prefix}.sorted.bam SORT_ORDER=coordinate CREATE_INDEX=true
	fi 
	if [ ! -f ${prefix}.dedup.bam ]; then
		echo; echo "Running ${MD} on ${prefix}.sorted.bam"; echo
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $MD INPUT=${prefix}.sorted.bam OUTPUT=${prefix}.dedup.bam M=${prefix}.metrics.txt REMOVE_DUPLICATES=true
	fi
	if [ ! -f ${prefix}.dedup.bai ]; then
		echo; echo "Running ${BI} on ${prefix}.sorted.bam"; echo	
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $BI INPUT=${prefix}.dedup.bam
	fi
	if [ ! -f ${prefix}.intervals ]; then 
		echo; echo "Running ${GA} -T RealignerTargetCreator on ${prefix}.dedup.bam "; echo
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $GA -T RealignerTargetCreator -R $Tgindex -I ${prefix}.dedup.bam -o ${prefix}.intervals
	fi
	if [ ! -f ${prefix}.realigned.bam ]; then
		echo; echo "Running ${GA} -T IndelRealigner on ${prefix}.dedup.bam"; echo
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $GA -T IndelRealigner -R $Tgindex -I ${prefix}.dedup.bam \
		-targetIntervals ${prefix}.intervals -o ${prefix}.realigned.bam
	fi
	#if [ ! -f ${prefix}.raw.snps.indels.vcf ]; then
		#echo; echo "Running ${GA} -T HaplotypeCaller -ERC GVCF on ${prefix}.realigned.bam"; echo
		#java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $GA -T HaplotypeCaller -ERC GVCF -R $Tgindex \
		#-I ${prefix}.realigned.bam -stand_call_conf 30.0 -nct 10 -o ${prefix}.raw.snps.indels.vcf -ERCIS 100 \
		#--variant_index_type LINEAR --variant_index_parameter 128000
 	#fi
	if [ ! -f ${prefix}.raw.snps.indels.g.vcf ]; then
		echo; echo "Running ${GA} -T HaplotypeCaller -ERC BP_RESOLUTION on ${prefix}.realigned.bam"; echo
		java -Xmx100G -Djava.io.tmpdir=$tempdir -jar $GA -T HaplotypeCaller -ERC BP_RESOLUTION -R $Tgindex \
		-I ${prefix}.realigned.bam -stand_call_conf 30.0  -o ${prefix}.raw.snps.indels.g.vcf --variant_index_type LINEAR \
		--variant_index_parameter 128000 --sample_ploidy 2 
	fi
 	echo Done running ${prefix}!
 done < ${strain_text_file}

if [ $GVCF_ONLY ]; then
	echo; echo Done running!!; echo
	#exit 0
fi

## Step 3
## Combine the GVCF files (created in the last part of the last script into chunks small enough for the Cluster to handle with the available memory
## These will all be generated by the script and named .raw.snps.indels.g.vcf
## Modify the naming so that it makes sense to you for the -o output files
# cd $workingdir

# Loop  => Number of groups to be processed together = 8
min_prefix_per_group=5
num_of_prefix=${#prefix_list[@]}
min_groups=`expr ${num_of_prefix=$} \/ ${min_prefix_per_group=5} + 1`
merge_param=()
for (( loop1=0; loop1<${min_groups}; loop1++ )) 
do
	for (( j=0; j<${min_prefix_per_group}; j++ ))
	do  
		loop2=`expr ${loop1} \* ${min_prefix_per_group} + ${j}`

		if [ ${prefix_list[${loop2}]} ]
		then
			variant_parameter="${variant_parameter} --variant $workingdir/${prefix_list[${loop2}]}.raw.snps.indels.g.vcf"
		fi

	done
	if [ ! -f CombineGVCFs_${loop1}.g.vcf ]; then
		echo; echo "Running CombineGVCFs in loop number ${loop1}"; echo
		java -Xmx13G -Djava.io.tmpdir=$tempdir -jar $GA \
			-R $Tgindex \
 			-T CombineGVCFs \
 			-o CombineGVCFs_${loop1}.g.vcf \
			${variant_parameter}
		echo "Done with ${loop1}!"
	fi
	merge_param+="--variant ${workingdir}/CombineGVCFs_${loop1}.g.vcf "
	variant_parameter=''
done


# ## Step 4
# # Merge all desired g.vcf files from CombineGVCFs into a single VCF file of all strains
# # You will need to change the output of this name as well as making sure the input matches what you created in Step 3
# # This can run for a VERY long time. 96 Toxo genomes ran for 12 days straight before finishing
#cd $workingdir
if [ ! -f ${OUTPUT}_raw_variants.vcf ]; then
echo; echo "Running GenotypeGVCFs on *.g.vcf files" ; echo
java -Xmx13G -Djava.io.tmpdir=$tempdir -jar $GA -R $Tgindex \
	-T GenotypeGVCFs --max_alternate_alleles 4 -nt 6 \
	-stand_call_conf 30 \
	--disable_auto_index_creation_and_locking_when_reading_rods \
	-o ${OUTPUT}_raw_variants.vcf \
	${merge_param}
fi

## Step 4

##FOR this part , Fumiaki told me to follow the site: https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/tutorials/2806-how-to-apply-hard-filters-to-a-call-set

##I did not use this part but should be:

##1. Extract the SNPs from the call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_raw_snps.vcf ]; then
echo; echo "Running SelectVariants on ${OUTPUT}_raw_variants.vcf" ; echo
java -Xmx13G -Djava.io.tmpdir=$tempdir -jar $GA \
	-T SelectVariants \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_variants.vcf \
	-selectType SNP \
	-o ${OUTPUT}_raw_snps.vcf
fi

##This creates a VCF file called raw_snps.vcf, containing just the SNPs from the original file of raw variants.


##2. Apply the filter to the SNP call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_filtered_snps.vcf ]; then
echo; echo "Running VariantFiltration on ${OUTPUT}_raw_snps.vcf" ; echo
java -Xmx13G -Djava.io.tmpdir=$tempdir -jar $GA \
	-T VariantFiltration \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_snps.vcf \
	--filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
	--filterName "my_snp_filter" \
	-o ${OUTPUT}_filtered_snps.vcf
fi

##3.Extract the Indels from the call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_raw_indels.vcf ]; then
echo; echo "Running SelectVariants on ${OUTPUT}_raw_variants.vcf" ; echo
java -Xmx13G -Djava.io.tmpdir=$tempdir -jar $GA \
	-T SelectVariants \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_variants.vcf \
	-selectType INDEL \
	-o ${OUTPUT}_raw_indels.vcf
fi

##6. Apply the filter to the Indel call set

##Action

##Run the following GATK command:
if [ ! -f ${OUTPUT}_filtered_indels.vcf ]; then
echo; echo "Running VariantFiltration on ${OUTPUT}_raw_indels.vcf" ; echo
java -Xmx13G -Djava.io.tmpdir=$tempdir -jar $GA \
	-T VariantFiltration \
	-R ${REFERENCE} \
	-V ${OUTPUT}_raw_indels.vcf \
	--filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
	--filterName "my_indel_filter" \
	-o ${OUTPUT}_filtered_indels.vcf
fi

