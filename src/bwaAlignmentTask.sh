#!/bin/bash

file_string=$1;
hg19_reference=$2;

file_array=(${file_string//,/ });
#readDIR=`dirname ${arrIN[0]}`;

###Creates an aligned,sorted,indexed bam file for each sample
for FILEPATH in "${file_array[@]}"
do
	FNAME=`basename $FILEPATH`
	SAMPLE=${FNAME%.fastq}
	bwa mem -M -t 6 $hg19_reference $FILEPATH > ${SAMPLE}.sam
	samtools view -S -h -b ${SAMPLE}.sam > ${SAMPLE}.bam
	rm ${SAMPLE}.sam
	samtools sort ${SAMPLE}.bam ${SAMPLE}.sorted
	rm ${SAMPLE}.bam
	samtools index ${SAMPLE}.sorted.bam
done