#!/bin/bash

projPath="./DemoData"
mkdir -p $projPath/peaks

# Use SEACR to call peaks
for samp in {H3K27me3,H3K4me3}_rep{1..2}; do
	igg=$(echo $samp |sed 's/\(.*\)_\(.*\)/IgG_\2/g')
	# Call with paired IgG control
	SEACR.sh \
		$projPath/alignment/${samp}.Fragments.Norm.bedgraph \
		$projPath/alignment/${igg}.Fragments.Norm.bedgraph \
		non \
		stringent \
		$projPath/peaks/${samp}_SEACR_Control.peaks
	# Call with 0.01 FDR of peaks
	SEACR.sh \
                $projPath/alignment/${samp}.Fragments.Norm.bedgraph \
		0.01 \
                non \
                stringent \
                $projPath/peaks/${samp}_SEACR_FDR_0.01.peaks
done

