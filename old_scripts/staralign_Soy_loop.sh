#!/bin/bash
#BSUB -J staralign_Soy #job name
#BSUB -n 12 #number of threads
#BSUB -W 10:0 #time for job to complete
#BSUB -R span[hosts=1] #to keep tasks on one node
#BSUB -R "rusage[mem=0.2]" #to request a node with 20MB of memory
#BSUB -o ./logs/alignment/staralign_Soy_%J.out #output file
#BSUB -e ./logs/alignment/staralign_Soy_%J.err #error file

#to align RNA-seq reads to indexed genome using STAR
#STAR cannot make use of HPC MPI, must have -R options to set 1 node & memory
#set threads under 12 on Henry2
#input of indexed genome path is /share/bitcpt/S23/UNITYID/At/starindices
#input of sequence reads path is /share/bitcpt/S23/cleandata/Arabidopsis_thaliana/
#output of aligned reads will go into AlignedToTranscriptome subdirectory in working directory

STAR=/usr/local/usrapps/bitcpt/star/bin/STAR

# SET IN VARIABLES
IN=/share/bitcpt/S23/RawData/Glycine_max
index=starindices
out=AlignedToTranscriptome

EN=fq.gz

# Loop trough our five samples, align reads 1 and 2 and output to a directory named 
for SAMPLE in Gm_SA_Rep1 Gm_SA_Rep2 Gm_SA_Rep3 Gm_SA_Rep4 Gm_SA_Rep5
do
    echo ${IN}/${SAMPLE}.${EN}
    ${STAR} --runThreadN 12 --runMode alignReads --genomeDir ${index} --outFileNamePrefix ${out}/${SAMPLE}_ --readFilesIn ${IN}/${SAMPLE}_1.${EN} ${IN}/${SAMPLE}_2.${EN} --readFilesCommand zcat --outSAMtype BAM Unsorted --twopassMode Basic --quantMode TranscriptomeSAM
done