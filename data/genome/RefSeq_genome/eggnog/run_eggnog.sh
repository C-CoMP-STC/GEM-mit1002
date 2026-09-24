#!/bin/bash -l
#$ -pe omp 16
#$ -j y  # Merge the error and output streams into a single file
#$ -o genome/RefSeq_genome/eggnog/run_eggnog.out  # Name of the output file

module load miniconda

# Activate the eggnog environment
conda activate /projectnb/cometsfba/hscott/envs/eggnog

# Run EggNOG-mapper
emapper.py -i genome/RefSeq_genome/GCF_001077695.1_ASM107769v1_protein.faa --output genome/RefSeq_genome/eggnog/eggnog_output --cpu 16