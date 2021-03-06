#!/bin/bash

# Activate Enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/rnaseq/bin:$PATH

# Pull FASTQ files out of subdirectories and place in current directory
find . -type f -mindepth 2 -exec mv -i -- {} . \;

# Clean up empty subdirectory folders

find . -depth -mindepth 1 -type d -empty -exec rmdir {} \;

#Assign variables to annotation and assembly file locations

ASSEMBLY='ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq


##################################
## Generate Salmon Genome Index ##
##################################

# Make a directory to store the genome files

mkdir -p genome2

# Download and unpack the genome assembly

curl $ASSEMBLY -o ./genome2/assembly.fasta

#Create a directory to store the index

mkdir -p genome2/index2

# Index the Genome

salmon index -t ./genome2/assembly.fasta -i ./genome2/index2/hg38_index


#################
## Count Reads ##
#################

# Make quant output directory

mkdir -p ./results/quants

# Create array for FASTQ files

FASTQS=($(find ./sequences -name "*\.fastq"))


#Subset FASTQ's in individual arrays for Read 1 and Read 2 

FASTQSONE=(`echo ${FASTQS[@]} | sed 's/ /\n/g' | grep R1_001`)
FASTQSTWO=(`echo ${FASTQS[@]} | sed 's/ /\n/g' | grep R2_001`)


# Quantify

for n in {0..23}; do

 R1=${FASTQSONE[${n}]}
 R2=${FASTQSTWO[${n}]}
 OUTNAME=results/quants/$(basename ${R1} .fastq)_quant

 echo "Processing sample ${R1}"
 salmon quant -i ./genome2/index2/hg38_index -l A \
         -1 ${R1} \
         -2 ${R2} \
	--gcbias \
         -p 8 --validateMappings -o ${OUTNAME}
done

salmon quant -i genome2/index2/hg38_index -l A \
	-1 sequences/ILMN_906_Ohagan_1/DMNEG12R1_S19_R1_001.fastq.gz \
	-2 sequences/ILMN_906_Ohagan_1/DMNEG12R1_S19_R2_001.fastq.gz \
	-p 8 --validateMappings -o quants/ILMN_906_Ohagan_1/DMNEG12R1_S19_R1_001_quants.fastq.gz
