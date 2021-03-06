#!/bin/bash

# Activate Enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/salmon/bin:$PATH

# Assign variables to annotation and assembly file locations


ASSEMBLY='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz'

ANNOTATION='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz'

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

#Create a directory to store the index

mkdir -p genome/index

##################################
## Generate Salmon Genome Index ##
##################################

# Download and unpack the genome assembly

curl $ASSEMBLY -o ./genome/gencode.v37.fa.gz

# Index the Genome

salmon index -t ./genome/assembly.fa.gz -i ./genome/index/hg38_index --gencode



