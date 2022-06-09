#!/bin/bash
projPath="./DemoData"
chromSize=hg38.chrom.sizes

# Get scaling factor based on the spike-in read depth
for samp in {H3K27me3,H3K4me3,IgG}_rep{1..2}; do
	seqDepthDouble=$(samtools view -F 0x04 $projPath/alignment/${samp}.spikeIn.sam | wc -l)
	seqDepth=$((seqDepthDouble/2))
	echo $seqDepth >$projPath/alignment/${samp}.spikeIn.seqDepth
	if [[ "$seqDepth" -gt "1" ]]; then
	scale_factor=`echo "10000 / $seqDepth" | bc -l`
	echo "Scaling factor for $samp is: $scale_factor!"
	if [ ! -e $projPath/alignment/${samp}.Fragments.Norm.bedgraph ]; then
	bedtools genomecov \
		-bg \
		-scale $scale_factor \
		-i $projPath/alignment/${samp}.Fragments.bed \
		-g $chromSize \
		>  $projPath/alignment/${samp}.Fragments.Norm.bedgraph
	fi
	fi
done
