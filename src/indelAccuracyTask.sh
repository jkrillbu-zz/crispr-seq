#!/bin/bash

output_dir=$1; code_dir=$2; indel_cut_site=$3; indel_genes=$4; ref_fasta=$5; indel_cut_interval=$6;

mkdir ${output_dir}/Power

python ${code_dir}/findWildTypeReads.py -i ${output_dir}/Reads -o ${output_dir}/Power -g $indel_genes
Rscript ${code_dir}/createTestReads.R ${output_dir}/Power $indel_cut_site $indel_genes


bwa mem -M -t 6 $ref_fasta ${output_dir}/Power/tests.fastq > ${output_dir}/Power/tests.sam
samtools view -S -h -b ${output_dir}/Power/tests.sam > ${output_dir}/Power/tests.bam
samtools sort ${output_dir}/Power/tests.bam ${output_dir}/Power/tests.sorted
samtools index ${output_dir}/Power/tests.sorted.bam

python ${code_dir}/parseCIGAR.py -i ${output_dir}/Power -o ${output_dir}/Power/VariantCalls -c $indel_cut_interval -g $indel_genes

python ${code_dir}/analyzeMultiAlign.py -i ${output_dir}/Power/VariantCalls -c $indel_cut_interval -g $indel_genes 

Rscript ${code_dir}/indelAccuracy.R ${output_dir}/Power/VariantCalls/tests ${output_dir}/Power/VariantCalls
