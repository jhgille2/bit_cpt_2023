# BIT CPT Final portfolio project
![workflow image](https://github.com/jhgille2/bit_cpt_2023/blob/main/quantification%20workflow.svg)
*General diagram of the project workflow*  
  
    
## Overview
This repository holds the code for my BIT CPT final project for the spring 2023 semester. The analysis is organized into a [makefile](https://github.com/jhgille2/bit_cpt_2023/blob/main/makefile) that can be run on the hpc with a [job script](https://github.com/jhgille2/bit_cpt_2023/blob/main/job.sh). The pipeline can't be run completely "out of the box", but I've tried to make it as painless as possible to work with our data and directory structure for the class.  

## Project description  
We have .fastq files from rna seq reads for three soybean tissue types, and five replicates for each tissue type. The overall goal of the pipeline in this repository is to quantify transcripts for each of these samples so that transcript counts can be compared between the different tissue types. This process involves a few steps: 
1. Quality control of input .fastq files with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  
  **Inputs:** 15 samples x 2 reads per sample = 30 .fastq files.  
  **Outputs:** 30 fastqc html reports (one for each .fastq file).  
2. Index a reference genome with [star](https://github.com/alexdobin/STAR).  
  **Inputs:** genome .fna file, annotations .gtf file.  
  **Outputs:** Index helper files to a directory (`starindices`).  
3. Align reads to the reference genome/transcriptome again with star.  
  **Inputs:** genome index directory, the two paired reads for each sample for each alignment.  
  **Outputs:** One .bam file per sample.  
4. Quantify these aligned reads with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)  
  **Inputs:** Alignment .bam file for each sample, transcriptome file.  
  **Outputs:** Quantification `quant.sf` file.  

## Running the pipeline  
The makefile uses input files from the `/share/bitcpt/S23` directory on the NCSU hazel cluster. The first part of the makefile sets up the paths to expected input and output files. The paths to the input files are set already, but some of the output paths assume a strict directory structure. The `make_folders.sh` file can help to set up this structure. In the class, we had seperate folders for Arabidopsis, Soybean analysis to the v2 Lee genome, and a third directory specifically for the portfolio project. This setup assumes you are starting from inside an empty portfolio directory e.g. `/share/bitcpt/S23/{UnityID}/Portfolio/`.  

First, copy `make_folders.sh`, `job.sh`, and `makefile`from this repository to the `Portfolio` directory. The other necessary folders can then be added by running: 
<br>  
`bash make_folders.sh` in the `Portfolio` directory.   
<br>  
Next, paths to the genome, annotations, and transcriptome files have to be changed in the makefile. These are variables that store the paths to genome, annotation, and transcriptome files but because we have multiple genomes along with their associated annotations and transcriptome files, these will have to be changed to they they actually point at the genome I'm using to the portfolio (they're pointing to the data for the v2 Lee genome right now).   

Once the `AlignedToTranscriptome`, `salmon_align_quant`, `transcriptome`, `fastqc`, `starindices`, and `starOutputfiles` folders have been added, the pipeline can be run with:  
<br>  
`bsub < job.sh` Again from the `Portfolio` directory.  

## BIT CPT specifics
For the data we are using in BIT CPT Spring 2023, I know in advance what kind of, and how many files I want to see in each step. I'll just go through them here for completeness/a sanity check. 
1. fastqc: We have three tissue types that each have five samples. For each sample we have two paired reads that have one fastq file per read. So in total we have 30 fastq files that we need to do quality control on. The fastq step outputs data to the `fastqc` directory off the main directory so after the step runs, we should see 30 *\_fastqc.html files here, each paired with a \*_.zip directory.  
2. Indexing the genome: After this step runs, the STAR software outputs some somewhat cryptically named files with some cryptic content to the `starindices` directory. I am not sure what each of these files do individually and a cursory google search around some forums indicated that I'm not alone there. Basically just look that you have [these files](https://www.biostars.org/p/9471213/), that they are not empty, and the job didn't produce any concerning errors/output when generating the index.  
3. Alignment: In this step, we align each of the two reads for each of our samples to the genome. Because we have 30 samples, this will produce 30 *sets* of files in the `AlignedToTranscriptome` directory. Really, we're just worried about the alignment **.bam** files though, and specifically the alignments to the transcriptome. There is one of these alignment files for each sample and they can be identified with {sample name}\_Aligned.toTranscriptome.out.bam. So in total, there will be 30 of these alignment files produced by this step, one for each sample.  
4. Quantification: This step takes the alignment files produced by the previous steps and quantifies the transcripts using some provided transcriptome .fa file. This produces a quantification `quant.sf` file for each sample, and stores this `quant.sf` file in a directory for each sample within the `salmon_align_quant` directory. Because there are 30 samples, there will be 30 of these quantification files, each one within one of 30 sample directories that are in turn nested within the `salmon_align_quant` directory.  
## Output
One quant.sf file is made for each sample which can be found after running the pipeline in directories named after each sample in the `salmon_align_quant` directory. These can be imported into galaxy/R for downstream analysis.  
