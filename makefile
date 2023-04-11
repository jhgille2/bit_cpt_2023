SHELL := bash

# Paths to input files
genome_path=/share/bitcpt/S23/referenceGenomes/Glycine_max_Lee_v2/glyma.Lee.gnm2.K7BV.genome_main.fna
annotations_path=/share/bitcpt/S23/referenceGenomes/Glycine_max_Lee_v2/glyma.Lee.gnm2.ann1.1FNT.gene_models_main.AGAT.gtf
transcriptome_path=/share/bitcpt/S23/referenceGenomes/Glycine_max_Lee_v2/glyma.Lee.gnm2_transcriptome.fasta

# Paths to input directories
clean_data_dir=/share/bitcpt/S23/cleandata/Glycine_max
raw_data_dir=/share/bitcpt/S23/RawData/Glycine_max

# Paths to output directories
index=starindices
transcriptome_out=AlignedToTranscriptome
fastqc_dir=fastqc
quant_dir=salmon_align_quant

# Paths to programs
STAR=/usr/local/usrapps/bitcpt/star/bin/STAR
FASTQC=/usr/local/usrapps/bitcpt/fastqc/bin/fastqc
SALMON=/usr/local/usrapps/bitcpt/salmon/bin/salmon

# Sample lists
raw_samples = Gm_SA_Rep1 Gm_SA_Rep2 Gm_SA_Rep3 Gm_SA_Rep4 Gm_SA_Rep5
clean_samples = Gm_OldLeaf_Rep1 Gm_OldLeaf_Rep2 Gm_OldLeaf_Rep3 Gm_OldLeaf_Rep4 Gm_OldLeaf_Rep5 Gm_YoungLeaf_Rep1 Gm_YoungLeaf_Rep2 Gm_YoungLeaf_Rep3 Gm_YoungLeaf_Rep4 Gm_YoungLeaf_Rep5

all_samples := $(raw_samples) $(clean_samples)  

# Output bam files
raw_samples_bam := $(foreach wrd,$(raw_samples),$(transcriptome_out)/$(wrd)_Aligned.toTranscriptome.out.bam)
clean_samples_bam := $(foreach wrd,$(clean_samples),$(transcriptome_out)/$(wrd)_Aligned.toTranscriptome.out.bam)
all_samples_bam := $(raw_samples_bam) $(clean_samples_bam)

# Output quantification files
quant_dir := $(foreach wrd,$(all_samples), $(quant_dir)/$(wrd))
quant_files := $(foreach wrd,$(quant_dir), $(wrd)/quant.sf)

# START TARGETS
##################################################################

# All target to run everything
all: starindices/SA $(raw_samples_bam) $(clean_samples_bam) $(quant_files)

fastqc/*_fastqc.html:
	${FASTQC} /share/bitcpt/S23/RawData/Glycine_max/* -t 12 -outdir ${fastqc_dir}

# Target to index the genome
starindices/SA:
	${STAR} --runThreadN 12 \
	 --runMode genomeGenerate \
	 --genomeSAindexNbases 13 \
	 --genomeDir ${index} \
	 --genomeFastaFiles ${genome_path} \
	 --sjdbGTFfile ${annotations_path} \
	 --sjdbOverhang 100


# Align files to transcriptome
$(raw_samples_bam): starindices/SA
	for SAMPLE in $(raw_samples) ; do \
    	$(STAR) --runThreadN 12 --runMode alignReads --genomeDir $(index) --outFileNamePrefix $(transcriptome_out)/$${SAMPLE}_ --readFilesIn $(raw_data_dir)/$${SAMPLE}_1.fq.gz $(raw_data_dir)/$${SAMPLE}_2.fq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --twopassMode Basic --quantMode TranscriptomeSAM ; \
	done

$(clean_samples_bam): starindices/SA
	for SAMPLE in $(clean_samples); do \
    	$(STAR) --runThreadN 12 --runMode alignReads --genomeDir $(index) --outFileNamePrefix $(transcriptome_out)/$${SAMPLE}_ --readFilesIn $(clean_data_dir)/$${SAMPLE}_1.fp.fq.gz $(clean_data_dir)/$${SAMPLE}_2.fp.fq.gz --readFilesCommand zcat --outSAMtype BAM Unsorted --twopassMode Basic --quantMode TranscriptomeSAM ; \
	done


$(quant_files): $(all_samples_bam)
	for i in ${!all_samples_bam[@]} ; do \
		${SALMON} quant -l A \
		-a ${all_samples_bam[i]} \
		--targets ${transcriptome_path} \
		-0 ${quant_dir[i]} \
	done