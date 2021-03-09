
# Remember to activate rnaseq env

# Load R
R

# Load libraries

library(readr)
library(dplyr)
library(magrittr)
library(tximport)
library(DESeq2)
library(ggplot2)
library(apeglm)


############
# tximport #
############

# Read in sample sheet and gene map

sample_table <- read_tsv("rnaseq_samples.txt")
gene_map <- read_tsv("gene_map.tsv")

# Direct tximport to proper files
# Extract the file names from our sample sheet, and append /quant.sf in a vector

sample_files <- paste0(pull(sample_table, `file`))
sample_files <- paste0('quants/', sample_files, '/quant.sf')

# Name the elements of the vector with the "file" names
# The vector elements names will also become column names in the tximport output

names(sample_files) <- pull(sample_table, `file`)

# tx2gene - summarizes multiple transcript expression values to a single gene value
# Convert pseudocounts to counts using tximport

count_data <- tximport(files = sample_files,
				type = "salmon",
				tx2gene = gene_map,
				ignoreTxVersion = TRUE)
			
# Matching count_data order to sample_table
all(rownames(sample_files) %in% colnames(count_data$counts))

rownames(sample_table) <- pull(sample_table, `file`)

all(rownames(sample_files) == colnames(count_data$counts))

# Need to convert tibbles to dataframes for DESeq2
# Need to remove all spaces if present in columns headers or rows

sample_table <- as.data.frame(sample_table)
sample_table$abv <- as.factor(sample_table$abv)
sample_table$batch <- as.factor(sample_table$batch)

# Design = is the expression of a gene dependent on the design of my exp			

dds <- DESeqDataSetFromTximport(txi = count_data,
					colData = sample_table,
					design = ~ batch+abv)	
				
#######
# PCA #
#######

# PCA is a parametric technique - assumes normality
# log transformed RNA-seq data is approximately normally distributed, so PCA is compatible with log transformed counts
# 0's cause issues with log transormation and RNAseq data is heterostochastic (i.e. variance is unequal accross data)
# Variance stabilizaing transformation and regularized log transformation can both be used to address this

vst <- varianceStabilizingTransformation(dds)

# To visualize counts prior to normalization
boxplot(counts(deseq_dataset, normalized=TRUE))

# To visualize normalized counts
boxplot(assay(vst))

# Make PCA plot (can be modified using ggplot2 functions)

plotPCA(vst, intgroup='abv')

###########################
# Hierarchical clustering #
###########################

# Get distance matrix for vst

d = assay(vst)

# Transorm it

d = t(d)

# Calculate distance matrix

d = dist(d)

# Plot clusters

h = hclust(d)


##########
# DESeq2 #
##########

# 3 Steps to DESeq2 analysis

# 1) estimate size factors (normalization)

dds <- estimateSizeFactors(dds)

# To view factors, normalizationFactors(deseq_dataset)
# View normalized counts, counts(deseq_dataset, normalized=TRUE)[1:6, 1:3]
# View unnormalized counts, count_data$counts[1:6, 1:3]

# 2) estimate dispersions

# For low expression genes technical variance dominates
# For high expression genes biological variance dominates
# Variance needs to be estimated
# Know 3 reps isnt sufficient to get true mean
# Uses the variance across whole experiment (observations for each sample for each gene) 
# Empirical bayesian shrinkage - uses measurements for one gene and compares it to 
# measurements of all the other genes and borrows information about the variance of 
# the other genes to shrink the estimated variance for the gene currently being considered

dds <- estimateDispersions(dds)

# Quality control plot
# Expect to see that as mean norm count increases, that variance decreases
# Circled blue points are considered as outliers

plotDispEsts(dds)

# 3) Apply statistics (Wald Test)
# Statistic applied against a negative binomial distribution

dds <- nbinomWaldTest(dds)


# All three steps can be short-cutted by dds <- DESeq(dds)

# Show comparisons made

resultsNames(dds)

# Snap shot of results from comparison made

result_table <- results(dds)
summary(result_table)

# To quickly visualize counts for an individual gene

plotCounts(dds, gene="ENSG00000122859", intgroup="abv")


