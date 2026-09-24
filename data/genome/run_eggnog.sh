#!/bin/bash -l
#$ -pe omp 16
#$ -j y  # Merge the error and output streams into a single file

module load miniconda

# Activate the eggnog environment
conda activate /projectnb/cometsfba/hscott/envs/eggnog

# Run EggNOG-mapper
emapper.py -i genome/MIT1002_anvio_prot_seqs.fa --output genome/eggnog_output --cpu 16