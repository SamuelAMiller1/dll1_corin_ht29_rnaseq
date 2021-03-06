#!/bin/bash

## Create a single .tsv file that contains gencode v37 transcript and gene ID's

# Copy transcript ID's to text file

zgrep -P -o 'ENST\d{11}' genome/gencode.v37.transcripts.fa.gz > enst.txt


# Copy gene ID's to text file

zgrep -P -o 'ENSG\d{11}' genome/gencode.v37.transcripts.fa.gz > ensg.txt

# Conjoin txt files into one tsv file

paste enst.txt ensg.txt > gene_map.tsv


