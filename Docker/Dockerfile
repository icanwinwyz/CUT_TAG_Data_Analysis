FROM bioconductor/bioconductor_docker:RELEASE_3_13
# FROM rocker/tidyverse:4.1.2

# Install bedtools, bowtie2, fastqc, deeptools, and samtools dependencies
RUN apt-get update && \
	apt-get install -y bc nano less curl bedtools bowtie2 fastqc libgsl-dev libncurses5-dev libncursesw5-dev libbz2-dev liblzma-dev pigz pip && \
	 pip install deeptools cutadapt multiqc

# Install trim galore
RUN curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.6.tar.gz -o trim_galore.tar.gz \
	&& tar xvzf trim_galore.tar.gz \
	&& rm  trim_galore.tar.gz
ENV PATH=/TrimGalore-0.6.6:$PATH

# Install samtools
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 \
	&& tar xfj samtools-1.13.tar.bz2 \
	&& cd samtools-1.13 \
	&& ./configure \
	&& make \
	&& make install \
	&& cd / \
	&& rm samtools-1.13.tar.bz2

# Install Picard
RUN git clone https://github.com/broadinstitute/picard.git && \
	cd picard/ && \
	./gradlew shadowJar
ENV PICARD='java -jar /picard/build/libs/picard.jar'

# Download SEACR
RUN git clone https://github.com/FredHutch/SEACR.git && \
	ln -s /SEACR/SEACR_1.3.sh /bin/SEACR.sh && \
	ln -s /SEACR/SEACR_1.3.R /bin/SEACR_1.3.R

# Install CRAN packages
RUN install2.r --error corrplot cowplot FactoMineR matrixStats patchwork rcartocolor tidyverse viridis

# Install Bioconductor Packages
RUN R -e 'BiocManager::install(c("ChIPQC", "ChIPseeker","chromVAR", "clusterProfiler", "DESeq2","limma", "org.Hs.eg.db", "org.Mm.eg.db", "GenomicAlignments", "GenomicRanges", "IRanges", "Rsamtools", "Rsubread", "rtracklayer","TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene"))'

# Install GitHub R packages
RUN R -e 'devtools::install_github("stjude/ChIPseqSpikeInFree")'
