#!/bin/bash

# Activate Enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/salmon/bin:$PATH

# Assign variables to annotation and assembly file locations

ASSEMBLY='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz'

ANNOTATION='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz'

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

#Create directories to store the index and quant matrix

mkdir -p genome
mkdir -p quants

##################################
## Generate Salmon Genome Index ##
##################################

# Download and unpack the genome assembly and annotation

curl $ASSEMBLY -o ./genome/gencode.v37.transcripts.fa.gz

# Index the Genome

salmon index -t ./genome/gencode.v37.transcripts.fa.gz -i ./genome/gencode_v37_index --gencode


##################################
## Generate Salmon Pseudocounts ##
##################################

# Quantify reads
# Next run, add flags --seqBias --useVBOpt

for fn in sequences/ILMN_906_Ohagan_{1..24};
 do
 samp=`basename ${fn}`
 salmon quant -i ./genome/gencode_v37_index -l A \
        -1 ${fn}/*\_R1_001.fastq.gz \
        -2 ${fn}/*\_R2_001.fastq.gz \
	-o quants/${samp} \
	--validateMappings \
	--seqBias
done

