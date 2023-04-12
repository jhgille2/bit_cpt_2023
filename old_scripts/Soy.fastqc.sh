#!/bin/tcsh
#BSUB -J fastqc_Soy_Meristem #job name
#BSUB -n 20 #number of nodes
#BSUB -W 2:0 #time for job to complete
#BSUB -o fastqc.%J.out #output file
#BSUB -e fastqc.%J.err #error file

# For running fastqc on all my Arabidopsis samples
# Run in working directory /share/bitcpt/S23/UnityID/At
# Must run this in working directory with subdirectory named /fastqc

# -t specifies number of threads

/usr/local/usrapps/bitcpt/fastqc/bin/fastqc /share/bitcpt/S23/RawData/Glycine_max/* -t 20 -outdir ./fastqc
