# BIT CPT Final portfolio project
![workflow image](https://github.com/jhgille2/bit_cpt_2023/blob/main/quantification%20workflow.svg)
*General diagram of the project workflow*  
  
    
## Overview
This repository holds the code for my BIT CPT final project for the spring 2023 semester. The analysis is organized into a [makefile](https://github.com/jhgille2/bit_cpt_2023/blob/main/makefile) that can be run on the hpc with a [job script](https://github.com/jhgille2/bit_cpt_2023/blob/main/job.sh). The pipeline can't be run completely "out of the box", but I've tried to make it as painless as possible to work with our data and directory structure for the class.  

## Project description  
We have **.fastq** files for three soybean tissue types, and five replicates for each tissue type. The overall goal of the pipeline in this repository is to quantify transcripts for each of these samples so that transcript counts can be compared between the different tissue types. This process involves a few steps: 
1. Quality control of **.fastq** files with [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).  
2. Index a reference genome with [star](https://github.com/alexdobin/STAR).  
3. Align reads to the reference genome/transcriptome again with star.  
4. Quantify these aligned reads with [salmon](https://salmon.readthedocs.io/en/latest/salmon.html)  

## Running the pipeline  
The makefile uses input files from the `/share/bitcpt/S23` directory on the NCSU hazel cluster. The first part of the makefile sets up the paths to expected input and output files. The paths to the input files are set already, but some of the output paths assume a strict directory structure. The `make_folders.sh` file can help to set up this structure. In the class, we had seperate folders for Arabidopsis, Soybean analysis to the v2 Lee genome, and a third directory specifically for the portfolio project. This setup assumes you are starting from inside an empty portfolio directory e.g. `/share/bitcpt/S23/{UnityID}/Portfolio/`.  

First, copy `make_folders.sh`, `job.sh`, and `makefile`from this repository to the `Portfolio` directory. The other necessary folders can then be added by running:  
`bash make_folders.sh` in the `Portfolio` directory.   

Once the `AlignedToTranscriptome`, `salmon_align_quant`, `transcriptome`, `fastqc`, `starindices`, and `starOutputfiles` folders have been added, the pipeline can be run with:  
<br>
`bsub < job.sh`  
<br>
Again from the `Portfolio` directory.
