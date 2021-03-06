#!/bin/bash

# Activate Enviornment
#ZZ~/miniconda3/envs/salmon/bin/salmon
export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/salmon/bin:$PATH

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

# Quantify reads
# Next run, add flags --seqBias --useVBOpt

for fn in sequences/ILMN_906_Ohagan_{1..24};
 do
 samp=`basename ${fn}`
 salmon quant -i ./genome2/index2/hg38_index -l A \
	-1 ${fn}/*\_R1_001.fastq.gz \
	-2 ${fn}/*\_R2_001.fastq.gz \
	-p 8 --validateMappings -o quants/${samp}_quant
done
