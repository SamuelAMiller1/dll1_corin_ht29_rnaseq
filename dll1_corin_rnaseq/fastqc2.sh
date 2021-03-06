#!/bin/bash

# Activate enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/rnaseq/bin:$PATH

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

# Collect fastqs into array
FASTQS=($(find ./sequences -type f -name "*.fastq.gz"))

# Run FASTQ

fastq \
 -t 6 \
 ${FASTQS[@]}

# Make results file

mkdir -p results/fastqc/

# Move .html files

HTML=($(find ./sequences -type f -name "*fastqc.html"))

mv ${HTML[@]} results/fastqc/

# Remove .zip files

ZIP=($(find ./sequences -type f -name "*fastqc.zip"))

for files in "${ZIP[@]}"; do
 rm -f "$files"
done

