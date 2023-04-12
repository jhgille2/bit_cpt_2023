#!/bin/tcsh
#BSUB -J soy_index_genome #job name
#BSUB -n 10 #number of nodes
#BSUB -W 2:0 #time for job to complete
#BSUB -o ./logs/index/starindices.out.%J #output file
#BSUB -e ./logs/index/starindices.err.%J #error file

# For running star to generate genome index
# Run in working directory /share/bitcpt/S23/UnityID/Soy
# Must run this in working directory with subdirectory named starindices/

set STAR=/usr/local/usrapps/bitcpt/star/bin/STAR
set IN=/share/bitcpt/S23/referenceGenomes/Glycine_max_Lee_v2


${STAR} --runThreadN 10 --runMode genomeGenerate --genomeSAindexNbases 13 --genomeDir starindices/ --genomeFastaFiles ${IN}/glyma.Lee.gnm2.K7BV.genome_main.fna --sjdbGTFfile ${IN}/glyma.Lee.gnm2.ann1.1FNT.gene_models_main.AGAT.gtf --sjdbOverhang 100
