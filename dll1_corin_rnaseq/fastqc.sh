#!/bin/bash

# Set path to enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/rnaseq/bin:$PATH

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq/

FASTQS=($(find ./sequences -type f -name "*.fastq.gz"))

fastqc \
 	-t 6 \
	${FASTQS[@]}

mkdir -p results/fastqc/




