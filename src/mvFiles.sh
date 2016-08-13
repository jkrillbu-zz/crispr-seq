#!/bin/bash

file_string=$1;
target_dir=$2;

file_array=(${file_string//,/ });
#readDIR=`dirname ${arrIN[0]}`;

###Creates an aligned,sorted,indexed bam file for each sample
for FILEPATH in "${file_array[@]}"
do
	mv $FILEPATH ${target_dir}/
done