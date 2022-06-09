// Read Sample sheet
Channel
    .fromPath(params.sampleSheet)
    .splitCsv(sep:'\t', header: ['samp','hist','sex','geno','rep','fastq1','fastq2','CTLsamp'])
    .map{ row-> tuple(row.samp, row.hist, row.sex, row.geno, row.rep , file(row.fastq1), file(row.fastq2), row.CTLsamp) }
    .set { samples_ch }

// Get references
Ref_ch = Channel.fromPath(params.Ref).collect()
SpikeInRef_ch = Channel.fromPath(params.SpikeRef).collect()
Channel.fromPath(params.chromSize).collect().into { chromSize_ch1; chromSize_ch2 }
Channel.fromPath(params.GTF).into { GTF_ch1; GTF_ch2 }
RenameQC_ch = Channel.fromPath(params.RenameQC)

// Get Defaults
minQualityScore_ch = Channel.value(params.minQualityScore)
binLen_ch = Channel.value(params.binLen)
corr_methods = ['spearman', 'pearson']
FDR_Thresh_ch = Channel.value(params.FDR_Thresh)

process Trim {
	tag "$samp"
	publishDir 'results/trim/', pattern: "*_val_*.fq.gz", mode: 'copy'
	cpus 4

	input:
	tuple samp, hist, sex, geno, rep, file(fastq1), file(fastq2), CTLsamp from samples_ch

	output:
	tuple samp, hist, sex, geno, rep, file("${samp}_val_1.fq.gz"), file("${samp}_val_2.fq.gz"), CTLsamp into trimmed_reads1, trimmed_reads2
	file "*trimming_report.txt" into trimgalore_results
	file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports

	script:	
	"""
	trim_galore \
                --fastqc \
                -j ${task.cpus} \
                --paired \
		--basename ${samp} \
                ${fastq1} ${fastq2}
	"""
}

process Align {
	tag "$samp"
	cpus 10

	input:
	tuple samp, hist, sex, geno, rep, file(trim_fastq1), file(trim_fastq2), CTLsamp from trimmed_reads1
	path Ref_ch

	output:
	file "${samp}_Summary.txt" into align_summary_ch
	tuple samp, hist, sex, geno, rep, file("${samp}.sam"), CTLsamp into sam_ch1, sam_ch2
	
	script:
	"""
	INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
	bowtie2 \
		--end-to-end \
                --very-sensitive \
                --no-mixed \
                --no-discordant \
                --phred33 \
                -I 10 \
                -X 700 \
                -p ${task.cpus} \
                -x \$INDEX \
                -1 ${trim_fastq1} \
                -2 ${trim_fastq2} \
		-S ${samp}.sam &> ${samp}_Summary.txt
	"""
}

process AlignSpike {
        tag "$samp"
        cpus 10

	input:
        tuple samp, hist, sex, geno, rep, file(trim_fastq1), file(trim_fastq2), CTLsamp from trimmed_reads2
	path SpikeInRef_ch

	output:
	file "${samp}.spikeIn_Summary.txt" into align_spike_summary_ch
        tuple samp, file("${samp}.spikeIn.sam") into sam_spike_ch

	script:
	"""
        INDEX=`find -L ./ -name "*.rev.1.bt2" | sed 's/.rev.1.bt2//'`
	bowtie2 \
		--end-to-end \
                --very-sensitive \
                --no-mixed \
                --no-discordant \
                --phred33 \
                -I 10 \
                -X 700 \
                -p ${task.cpus} \
                -x \$INDEX \
                -1 ${trim_fastq1} \
                -2 ${trim_fastq1} \
                -S ${samp}.spikeIn.sam &> ${samp}.spikeIn_Summary.txt
	"""
}

process GetFragLength {
	tag "$samp"
	cpus 2

	input:
        tuple samp, hist, sex, geno, rep, file(sam), CTLsamp from sam_ch1

	output:
	file "${samp}_FragmentLength.txt" into FragLength_ch

	script:
	"""
	samtools view \
		-@ 1 \
                -F 0x04 \
                ${sam} | \
                awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs(\$9)}' | \
                sort | \
                uniq -c | \
                awk -v OFS="\t" '{print \$2, \$1/2}' \
                > ${samp}_FragmentLength.txt
	"""

}

process FilterSortSam {
	tag "$samp"
	cpus 5

	input:
        tuple samp, hist, sex, geno, rep, file(sam), CTLsamp from sam_ch2
	val minQualityScore from minQualityScore_ch
	val binLen from binLen_ch

	output:
        tuple samp, hist, sex, geno, rep, file("${samp}.Q${minQualityScore}.Sort.bam"), CTLsamp into bam_ch1, bam_ch2

	script:
	"""
	samtools view \
		-@ 4 \
		-bS \
		-F 0x04 \
		-q ${minQualityScore} \
		${sam} | \
	samtools sort \
                -@ 4 \
		> ${samp}.Q${minQualityScore}.Sort.bam
	"""

}

process MarkDup {
	tag "$samp"
	publishDir 'results/align', pattern: "${samp}.Q${minQualityScore}.Sort.MarkDup.bam", mode: 'copy'
	cpus 2

	input:
        tuple samp, hist, sex, geno, rep, file(bam), CTLsamp from bam_ch1
        val minQualityScore from minQualityScore_ch

	output:
        tuple samp, hist, sex, geno, rep, file("${samp}.Q${minQualityScore}.Sort.MarkDup.bam"), CTLsamp into MarkDup_ch1
	file("${samp}.MarkDup") into MarkDup_summary_ch

	script:
	"""
	java -jar /picard/build/libs/picard.jar MarkDuplicates \
                                -I ${bam} \
                                -M ${samp}.MarkDup \
                                -O ${samp}.Q${minQualityScore}.Sort.MarkDup.bam
	"""
}

process IndexBam {
	tag "$samp"
	cpus 5
	publishDir 'results/align', pattern: "*.bai", mode: 'copy'

	input:
        tuple samp, hist, sex, geno, rep, file(MarkDup), CTLsamp from MarkDup_ch1

	output:
	tuple samp, hist, sex, geno, rep, file(MarkDup), file("*.bai"), CTLsamp into Index_ch1, Index_ch2
	tuple file(MarkDup), file("*.bai") into Index_ch3, Index_ch4, Index_ch5, Index_ch6

	script:
	"""
	samtools index  \
		-@ ${task.cpus} \
                ${MarkDup}
	"""
}

process IdxStat {
	tag "$samp"
	cpus 1

	input:
        tuple samp, hist, sex, geno, rep, file(MarkDup), file(idx), CTLsamp from Index_ch1

	output:
	file "*_idxstat" into idxstat_ch

	script:
	"""
	samtools idxstats \
		${MarkDup} > \
                ${samp}_idxstat
	"""
}

process ConvertBED {
	tag "$samp"
	cpus 2

	input:
        tuple samp, hist, sex, geno, rep, file(bam), CTLsamp from bam_ch2

	output:
        tuple samp, hist, sex, geno, rep, file("*.bed"), CTLsamp into bed_ch1

	script:
	"""
	bedtools bamtobed \
		-i ${bam} \
                -bedpe | \
	awk '\$1==\$4 && \$6-\$2 < 1000 {print \$0}' \
                > ${samp}.bed 
	"""
}

process GetBedFrags {
	tag "$samp"
	cpus 2

	input:
        tuple samp, hist, sex, geno, rep, file(bed), CTLsamp from bed_ch1

	output:
        tuple samp, file("*.Fragments.bed") into bedFrag_ch1
	tuple samp, hist, CTLsamp, file("*.Fragments.bed") into bedFrag_ch2, bedFrag_ch3

	script:
	"""
	cut -f 1,2,6 ${bed} | \
	sort -T ./ -k1,1 -k2,2n -k3,3n  >  ${samp}.Fragments.bed
	"""
}

process GetFragCount {
	tag "$samp"
	publishDir 'results/QC', pattern: "${samp}.FragmentCount.bin${binLen}.bed", mode: 'copy'
	cpus 1

	input:
        tuple samp, file(Frag) from bedFrag_ch1
	val binLen from binLen_ch

	script:
	"""
	awk -v \
		w=${binLen} \
                '{print \$1, int((\$2 + \$3)/(2*w))*w + w/2}' \
                ${Frag} | \
        sort -k1,1V -k2,2n | \
        uniq -c | \
        awk -v OFS="\t" '{print \$2, \$3, \$1}' | \
        sort -k1,1V -k2,2n \
                > ${samp}.FragmentCount.bin${binLen}.bed
	"""
}

process ConvertBigWig {
	tag "$samp"
	publishDir 'results/peaks', pattern: "${samp}.bw", mode: 'copy'
	cpus 5

	input:
        tuple samp, hist, sex, geno, rep, file(MarkDup), file(idx), CTLsamp from Index_ch2
	val binLen from binLen_ch

	output:
	file("*.bw") into BW_ch1

	script:
	"""
	bamCoverage \
		-b ${MarkDup} \
                -o ${samp}.bw \
                -p ${task.cpus} \
                --binSize ${binLen} \
		--extendReads
	"""
}

process CalcSpikeDepth {
        tag "$samp"
        cpus 2
	
	when:
	params.Spike

	input:
        tuple samp, file(sam) from sam_spike_ch

	output:
	tuple samp, env(scale_factor) into spike_depth_ch

	script:
	"""
	seqDepthDouble=\$(samtools view -F 0x04 ${sam} | wc -l)
	if [ "\$seqDepthDouble" -eq "0" ]; then
		scale_factor='false'
	else
	        seqDepth=\$((seqDepthDouble/2))
	        scale_factor=`echo "10000 / \$seqDepth" | bc -l`
	fi
	"""
}

// TODO find a way to deal with the spike in E coli if we ever get one
Control_Histone = 'NA'
process ConvertBedgraphTreatment {
        tag "$samp"
        publishDir 'results/peaks', pattern: "${samp}.bedgraph", mode: 'copy'
        cpus 2

	when:
	hist != Control_Histone

        input:
        tuple samp, val(hist), CTLsamp, file(Frag) from bedFrag_ch2
        path chromSize_ch1

        output:
        tuple val(CTLsamp), val(samp), file("${samp}.bedgraph") into TRT_Bedgraph_ch1, TRT_Bedgraph_ch2, TRT_BedGraph_ch3

        script:
        """
        bedtools genomecov \
                -bg \
                 -i ${Frag} \
                -g ${chromSize_ch1} \
                >  ${samp}.bedgraph
        """
}

process ConvertBedgraphControl {
	tag "$samp"
        publishDir 'results/peaks', pattern: "${samp}.bedgraph", mode: 'copy'
	cpus 2

	when:
	hist == Control_Histone

	input:
        tuple samp, val(hist), CTLsamp, file(Frag) from bedFrag_ch3
	path chromSize_ch2

	output:
	tuple val(samp), file("${samp}.bedgraph") into CTRL_Bedgraph_ch1, CTRL_Bedgraph_ch2, CTRL_Bedgraph_ch3

	script:
	"""
	bedtools genomecov \
                -bg \
                 -i ${Frag} \
                -g ${chromSize_ch2} \
		>  ${samp}.bedgraph
	"""
}

process SEACR_Control {
        tag "$samp"
        publishDir 'results/peaks', pattern: "${samp}_Control.peaks.stringent.bed", saveAs: { filename -> "${samp}_Control.SEACR.bed" }, mode: 'copy'
        cpus 2

        input:
        tuple CTRLsamp, samp, file(TRT), file(CTRL) from TRT_Bedgraph_ch1.combine(CTRL_Bedgraph_ch1, by: 0)

        output:
        tuple samp, file("${samp}_Control.peaks.stringent.bed") into SEACR_CTRL_ch

        script:
        """
        SEACR.sh \
                ${TRT} \
                ${CTRL} \
                norm \
                stringent \
                ${samp}_Control.peaks
        """
}

process SEACR_FDR {
        tag "$samp"
        publishDir 'results/peaks', pattern: "${samp}_FDR_${FDR}.peaks.stringent.bed", saveAs: {filename -> "${samp}_FDR_${FDR}.SEACR.bed" }, mode: 'copy'
        cpus 2

        input:
        val FDR from FDR_Thresh_ch
        tuple CTRLsamp, samp, file(TRT), file(CTRL) from TRT_Bedgraph_ch2.combine(CTRL_Bedgraph_ch2, by: 0)

        output:
        tuple samp, file("${samp}_FDR_${FDR}.peaks.stringent.bed") into SEACR_FDR_ch

        script:
        """
        SEACR.sh \
                ${TRT} \
                ${FDR} \
                norm \
                stringent \
                ${samp}_FDR_${FDR}.peaks
        """
}


process multiBamSummary {
        cpus 10
	memory '90 GB'

        input:
        file("*") from Index_ch3.collect()
        val binLen from binLen_ch

	output:
	file("multiBamSummary.npz") into bam_summary_ch1

        script:
        """
	multiBamSummary bins \
		-bs ${binLen} \
		--bamfiles *.bam \
		-o multiBamSummary.npz \
		-p ${task.cpus}
        """
}

process computeMatrix {
	cpus 10
	memory '90 GB'
	queue 'highmem.q'

	input:
	path GTF_ch1
        file("*") from BW_ch1.collect()

	output:
	file("computeMatrix.mat.gz") into computeMatrix_summary_ch

	script:
	"""
	computeMatrix scale-regions \
		-R ${GTF_ch1} \
		-S *.bw \
		-o computeMatrix.mat.gz \
		-b 3000 \
		-a 3000 \
		-m 5000 \
		-p ${task.cpus}
	"""
}

process bamPEFragmentSize {
	cpus 10
	memory '90 GB'
	
	input:
        file("*") from Index_ch4.collect()

	output:
	file("bamPEFragmentSize.table") into bamPEFragmentSize_table_ch
	file("bamPEFragmentSize.raw") into bamPEFragmentSize_raw_ch

	script:
	"""
	bamPEFragmentSize \
		--bamfiles *.bam \
		--table bamPEFragmentSize.table \
		--outRawFragmentLengths bamPEFragmentSize.raw \
		-p ${task.cpus}
	"""
}


process PlotEnrichment {
	cpus 10

        input:
        file("*") from Index_ch5.collect()
	path GTF_ch2

        output:
        file("plotEnrichment.raw") into PlotEnrichment_raw_ch

        script:
        """
	plotEnrichment \
		--bamfiles *.bam \
		--BED ${GTF_ch2} \
		--outRawCounts plotEnrichment.raw \
		-p ${task.cpus}
        """
}

process PlotFingerprint {
        cpus 10

        input:
        file("*") from Index_ch6.collect()

        output:
        file("plotFingerprint.qual") into plotFingerprint_qual_ch
	file("plotFingerprint.raw") into plotFingerprint_raw_ch

        script:
        """
        plotFingerprint \
                --bamfiles *.bam \
		--outQualityMetrics plotFingerprint.qual \
		--outRawCounts plotFingerprint.raw \
                -p ${task.cpus}
        """
}

process plotCorrelation {
	tag "$corr"
        cpus 10
	memory '110 GB'

        input:
	each corr from corr_methods
        file(npz) from bam_summary_ch1

        output:
	file("${corr}_plotCorrelation.matrix") into plotCorrelation_sumary_ch

        script:
        """
        plotCorrelation \
		-in ${npz} \
		-c ${corr} \
		-p heatmap \
		--outFileCorMatrix ${corr}_plotCorrelation.matrix
        """
}

process plotProfile {
        cpus 2
	memory '30GB'

        input:
        file(mat) from computeMatrix_summary_ch

        output:
        file("plotProfile.tab") into plotProfile_summary_ch

        script:
        """
        plotProfile \
		-m ${mat} \
		-o plotProfile.png \
		--outFileNameData plotProfile.tab
        """
}

process MultiQC {
	publishDir 'results/QC', mode: 'copy'
	cpus 2

	input:
	path RenameQC_ch
        file("trim/*") from trimgalore_results.collect()
        file("trim/*}") from trimgalore_fastqc_reports.collect()
        file("align/*") from align_summary_ch.collect()
        file("align/*") from align_spike_summary_ch.ifEmpty('').collect()
        file("align/*") from MarkDup_summary_ch.collect()
        file("align/*") from idxstat_ch.collect()
        file("align/*") from bamPEFragmentSize_table_ch.collect()
        file("align/*") from bamPEFragmentSize_raw_ch.collect()
        file("peaks/*") from PlotEnrichment_raw_ch.collect()
        file("peaks/*") from plotFingerprint_qual_ch.collect()
        file("peaks/*") from plotFingerprint_raw_ch.collect()
        file("peaks/*") from plotCorrelation_sumary_ch.collect()
        file("peaks/") from plotProfile_summary_ch.collect()

	output:
	file("multiqc_report.html") into multiqc_rep_ch
	path("multiqc_data") into multiqc_data_ch

	script:
	"""
	multiqc \
		--cl_config "sample_names_replace_regex: true" \
		--replace-names ${RenameQC_ch} \
		./
	"""	

}

