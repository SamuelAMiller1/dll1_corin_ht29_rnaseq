#!/bin/bash

#copy data from parent directory into working directory
cp -r /N/project/OHagan_single_cell/ILMN_906_OHagan_BL_mRNAseq24_Dec2020/ /N/slate/millesaa/

#libraries are prepped with SMART-Seq v4 Ultra Low Input RNA Kit for Sequencing

mkdir -p sequences

# Rename fastq directory to sequences

mv ILMN_906_OHagan_BL_mRNAseq24_Dec2020/ sequences/

# gunzip fastqs

gunzip -r ./sequences

# Activate Enviornment

conda activate rnaseq

ANNOTATION='ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.99.gtf.gz'
ASSEMBLY='ftp://ftp.ensembl.org/pub/release-91/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'

FASTQS=($(find ./sequences -type f -name "*.fastq"))

################################
## Generate STAR Genome Index ##
################################

# Make a directory to store the genome files

mkdir -p genome

# Download and unpack the genome assembly

curl $ASSEMBLY | gunzip > ./genome/assembly.fasta

#Download and unpack the genome annotation

curl $ANNOTATION | gunzip > ./genome/annotation.gtf

#Create a directory to store the index

mkdir -p genome/index

# Create the STAR genome index
# Argument values are specified to minimize memory usage

STAR \
--runThreadN 8 \
--genomeDir genome/index \
--genomeFastaFiles genome/assembly.fasta \
--sjdbGTFfile genome/annotation.gtf \
--genomeSAindexNbases 12 \
--genomeSAsparseD 3 \
--genomeChrBinNbits=16 \
--limitGenomeGenerateRAM=119000000000

#############################
## Aligned Reads to Genome ##
#############################

# Create an output directory for aligned reads

mkdir -p results/aligned

# Align the reads with STAR

FASTQS=($(find ./sequences -name "*\.fastq"))

FASTQSONE=($(echo ${FASTQS[@]} | grep "R1_001"))
FASTQSTWO=($(echo ${FASTQS[@]} | grep "R2_001"))

FASTQSONE=(`echo ${FASTQS[@]} | sed 's/ /\n/g' | grep R1_001`)
FASTQSTWO=(`echo ${FASTQS[@]} | sed 's/ /\n/g' | grep R2_001`)

for n in {0..23}; do

 R1=${FASTQSONE[${n}]}
 R2=${FASTQSTWO[${n}]}

 PREFIX=results/aligned/$(basename ${R1} .fastq)_

 STAR \
	--runThreadN 10 \
	--genomeDir genome/index \
	--readFilesIn ${R1} ${R2} \
	--outFileNamePrefix ${PREFIX} \
	--outSAMtype BAM SortedByCoordinate
done

# To check mapping efficiency >> less <filename>.Log.final.out

# Indexing the BAM files

BAMS=($(find ./results/aligned -name "*\.bam"))

for BAM in ${BAMS[@]}; do
 samtools index $BAM
done

####################
## Count Features ##
####################

# Create an output directory for read counts

mkdir -p results/counts

# Count reads

BAMS=$(find ./results/aligned -name "*\.bam")

featureCounts \
	-T 8 \
	-s 1 \
	-a genome/annotation.gtf \
	-o results/counts/counts.tsv \
	-t exon \
	-g gene_name \
	-p \
	--largestOverlap \
	--primary \
	${BAMS}



##Additional test code...

# Other test code for star alignment
for ((n = 0; n < ${#FASTQSONE[23]}; n++)); do 
R1=${FASTQSONE[${n}]} 
R2=${FASTQSTWO[${n}]} 
PREFIX=results/aligned/$(basename ${FASTQSONE[${@}]} .fastq)_ 
STAR 
--runThreadN 8
--runMode alignReads
--genomeLoad
--LoadAndKeep 
--genomeDir genome/index 
--readFilesIn ${R1} ${R2} 
--outFileNamePrefix ${PREFIX} 
--outSAMtype BAM SortedByCoordinate 
done

## new test code under development for STAR alignment
for i in *_R1_001.fastq; do
PREFIX=results/aligned/$(basename ${i%_R1_001.fastq} .fastq)_
STAR --runMode alignReads --genomeLoad  LoadAndKeep --readFilesCommand zcat --outSAMtype BAM SortedByCoordinate --genomeDir genome/index --readFilesIn $i ${i%_R1_001.fastq.}_R2.001.fastq --runThreadN 12 --outFileNamePrefix ${PREFIX}
done
