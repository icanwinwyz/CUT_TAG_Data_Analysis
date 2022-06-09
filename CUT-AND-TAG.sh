#!/bin/bash
#$ -cwd 
#$ -pe smp 4 
#$ -q all.q
#$ -j y
#$ -l h=csclprd1-c0v
#$ -l h_vmem=5G,mem_free=5G

module load anaconda3 singularity/3.6.0
conda activate nextflow

nextflow run main.nf \
	-ansi-log false


