#!/bin/bash
projPath="./DemoData"
ref=/mnt/hpc/reference/Bowtie/hg38/hg38
cores=5
Spike=true
SpikeInRef=/mnt/hpc/reference/Bowtie/fastq_screen/Ecoli/Ecoli

mkdir -p ${projPath}/alignment

# Aling to host genome
for samp in {H3K27me3,H3K4me3,IgG}_rep{1..2}; do
	if [ ! -e ${projPath}/alignment/${samp}.sam ]; then
	# use --local instead of --end-to-end for sequencing longer than 25 bp per read
	time bowtie2 \
		--end-to-end \
		--very-sensitive \
		--no-mixed \
		--no-discordant \
		--phred33 \
		-I 10 \
		-X 700 \
		-p $cores \
		-x $ref \
		-1 ${projPath}/fastq/${samp}_R1.fastq.gz \
		-2 ${projPath}/fastq/${samp}_R2.fastq.gz \
		-S ${projPath}/alignment/${samp}.sam &> ${projPath}/QC/${samp}_Summary.txt
	fi
	# Align to E. coli genome if using homegrown Tn5
	if [  "$Spike" == true ] && [ ! -e ${projPath}/alignment/${samp}.spikeIn.sam ]; then
	time bowtie2 \
		--end-to-end \
		--very-sensitive \
		--no-mixed \
		--no-discordant \
		--phred33 \
		-I 10 \
		-X 700 \
		-p $cores \
		-x $SpikeInRef \
		-1 ${projPath}/fastq/${samp}_R1.fastq.gz \
                -2 ${projPath}/fastq/${samp}_R2.fastq.gz \
		-S ${projPath}/alignment/${samp}.spikeIn.sam &> ${projPath}/QC/${samp}.spikeIn_Summary.txt
	fi
done
