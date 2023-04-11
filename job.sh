#!/bin/bash
#BSUB -J quantify_soy #job name
#BSUB -n 12 #number of threads
#BSUB -W 15:0 #time for job to complete
#BSUB -R span[hosts=1] #to keep tasks on one node
#BSUB -R "rusage[mem=0.2]" #to request a node with 20MB of memory
#BSUB -o ./logs/staralign_Soy_%J.out #output file
#BSUB -e ./logs/staralign_Soy_%J.err #error file

# Run in the /share/bitcpt/S23/jhgille2/Portfolio directory

make