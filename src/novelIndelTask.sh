#!/bin/bash

code_dir=$1; output_dir=$2; indel_cut_interval=$3; indel_genes=$4; neg_cntrl=$5; 


python ${code_dir}/parseCIGAR.py -i ${output_dir}/Reads -o ${output_dir}/VariantCalls -c $indel_cut_interval -g $indel_genes
	
python ${code_dir}/analyzeMultiAlign.py -i ${output_dir}/VariantCalls -c $indel_cut_interval -g $indel_genes

###Compiles table of indel fractions and total read counts [samples,genes]
###Creates a summary of types of indels observed and number of occurances for each sample,gene pair 
Rscript ${code_dir}/indel_quant_analysis.R ${code_dir} $output_dir $indel_cut_interval $neg_cntrl
Rscript ${code_dir}/indelSizes.R ${output_dir}/VariantCalls ${output_dir}/indel_quant/plots_trial