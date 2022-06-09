
#!/bin/bash
projPath="./DemoData"
cores=5
minQualityScore=2
binLen=500

# Consider a Picard markdup step here, though it is likely not necessary
# if that is done, also add an R script to summzarize the duplication rate


for samp in {H3K27me3,H3K4me3,IgG}_rep{1..2}; do
	# Get fragment length
	if [ ! -e $projPath/alignment/${samp}_FragmentLength.txt ]; then
	samtools view \
		-@ $cores \
		-F 0x04 \
		 ${projPath}/alignment/${samp}.sam | \
		awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | \
		sort | \
		uniq -c | \
		awk -v OFS="\t" '{print $2, $1/2}' \
		> $projPath/alignment/${samp}_FragmentLength.txt
	else
		echo "$projpath/alignment/$samp_FragmentLength.txt already present"
	fi

	# Filter for minimum alignment quality (-q) and for mappability (-F 0x04)
	if [ ! -e ${projPath}/alignment/${samp}.Q$minQualityScore.Mapped.bam ];then
	echo "Converting SAM to BAM"
	time samtools view \
		-@ $cores \
		-bS \
		-F 0x04 \
		-q $minQualityScore \
		${projPath}/alignment/${samp}.sam \
		> ${projPath}/alignment/${samp}.Q$minQualityScore.Mapped.bam
	else
		echo "${projPath}/alignment/${samp}.Q$minQualityScore.Mapped.bam already present"
	fi
	# Convert to bed
	if [ ! -e $projPath/alignment/${samp}.Q$minQualityScore.Mapped.bed ];then
	echo "Converting peaks to BED format"
	time bedtools bamtobed \
		-i $projPath/alignment/${samp}.Q$minQualityScore.Mapped.bam \
		-bedpe \
		> $projPath/alignment/${samp}.Q$minQualityScore.Mapped.bed
	else
		echo "$projPath/alignment/${samp}.Q$minQualityScore.Mapped.bed already present"
	fi

	# Clean bed file to include read on same chr with insert <1000bp
	if [ ! -e $projPath/alignment/${samp}.Q$minQualityScore.Mapped.Clean.bed ]; then
	echo "Filtering BED file for peaks on the same chr and with an insert <1kbp"
	time awk '$1==$4 && $6-$2 < 1000 {print $0}' \
		$projPath/alignment/${samp}.Q$minQualityScore.Mapped.bed \
		> $projPath/alignment/${samp}.Q$minQualityScore.Mapped.Clean.bed
	else
		echo "$projPath/alignment/${samp}.Q$minQualityScore.Mapped.Clean.bed already present"
	fi

	# Only extract the fragment related columns
	if [ ! -e  $projPath/alignment/${samp}.Fragments.bed ];then 
	echo "Outputting fragment-related columns to $projPath/alignment/${samp}.Fragments.bed"
	time cut -f 1,2,6 $projPath/alignment/${samp}.Q$minQualityScore.Mapped.Clean.bed \
		| sort -k1,1 -k2,2n -k3,3n  \
		>  $projPath/alignment/${samp}.Fragments.bed
	else
		echo " $projPath/alignment/${samp}.Fragments.bed already present"
	fi 

	# Get counts across bins to assess replicate reproducibility
	if [ ! -e  $projPath/alignment/${samp}.FragmentCount.bin${binLen}.bed ]; then
	time awk -v \
        	w=$binLen \
        	'{print $1, int(($2 + $3)/(2*w))*w + w/2}' \
        	$projPath/alignment/${samp}.Fragments.bed | \
        	sort -k1,1V -k2,2n | \
       		uniq -c | \
        	awk -v OFS="\t" '{print $2, $3, $1}' | \
        	sort -k1,1V -k2,2n \
        	> $projPath/alignment/${samp}.FragmentCount.bin${binLen}.bed
	else
		echo " $projPath/alignment/${samp}.FragmentCount.bin${binLen}.bed already present"
	fi

done
