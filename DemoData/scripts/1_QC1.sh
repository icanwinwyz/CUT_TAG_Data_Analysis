#!/bin/bash

projPath="./DemoData"
mkdir -p ${projPath}/QC/

cores=5
# Trim adaptors with trim_galore. This wraps cutadpat and fastqc
for samp in {H3K27me3,H3K4me3,IgG}_rep{1..2}; do
	trim_galore \
		-j $cores \
		--paired \
		-o ${projPath}/QC/ \
		${projPath}/fastq/${samp}_R1.fastq.gz \
		${projPath}/fastq/${samp}_R2.fastq.gz
done
