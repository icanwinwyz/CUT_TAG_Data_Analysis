#!/bin/bash

projPath="./DemoData"
mkdir -p ${projPath}/data/{H3K27me3,H3K4me3,IgG}_rep{1..2}
mkdir -p ${projPath}/fastq/

# Download data from ENA
# H3K27me3 Rep 1
wget -O $projPath/data/H3K27me3_rep1/H3K27me3_rep1_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR122/017/SRR12246717/SRR12246717_1.fastq.gz
wget -O $projPath/data/H3K27me3_rep1/H3K27me3_rep1_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR122/017/SRR12246717/SRR12246717_2.fastq.gz
# H3K27me3 Rep 2
wget -O $projPath/data/H3K27me3_rep2/H3K27me3_rep2_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/040/SRR11074240/SRR11074240_1.fastq.gz
wget -O $projPath/data/H3K27me3_rep2/H3K27me3_rep2_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/040/SRR11074240/SRR11074240_2.fastq.gz

# H3K4me3 Rep 1
wget -O $projPath/data/H3K4me3_rep1/H3K4me3_rep1_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/054/SRR11074254/SRR11074254_1.fastq.gz
wget -O $projPath/data/H3K4me3_rep1/H3K4me3_rep1_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/054/SRR11074254/SRR11074254_2.fastq.gz
# H3K4me3 Rep 2
wget -O $projPath/data/H3K4me3_rep2/H3K4me3_rep2_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/058/SRR11074258/SRR11074258_1.fastq.gz
wget -O $projPath/data/H3K4me3_rep2/H3K4me3_rep2_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR110/058/SRR11074258/SRR11074258_2.fastq.gz

# IgG Rep 1
wget -O $projPath/data/IgG_rep1/IgG_rep1_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/024/SRR11923224/SRR11923224_1.fastq.gz
wget -O $projPath/data/IgG_rep1/IgG_rep1_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR119/024/SRR11923224/SRR11923224_2.fastq.gz
# IgG Rep 2
wget -O $projPath/data/IgG_rep2/IgG_rep2_R1_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_1.fastq.gz
wget -O $projPath/data/IgG_rep2/IgG_rep2_R2_001.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/001/SRR8754611/SRR8754611_2.fastq.gz
wget -O $projPath/data/IgG_rep2/IgG_rep2_R1_002.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_1.fastq.gz
wget -O $projPath/data/IgG_rep2/IgG_rep2_R2_002.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR875/002/SRR8754612/SRR8754612_2.fastq.gz


# Merge technical reps
for samp in {H3K27me3,H3K4me3,IgG}_rep{1..2}; do
	cat ${projPath}/data/${samp}/*_R1_*.fastq.gz >${projPath}/fastq/${samp}_R1.fastq.gz
	cat ${projPath}/data/${samp}/*_R2_*.fastq.gz >${projPath}/fastq/${samp}_R2.fastq.gz
done


