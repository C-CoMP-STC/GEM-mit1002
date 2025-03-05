#!/bin/bash -l
#$ -pe omp 16
#$ -j y  # Merge the error and output streams into a single file

module load miniconda

# Activate the eggnog environment
conda activate /projectnb/cometsfba/hscott/envs/eggnog

# Run EggNOG-mapper
emapper.py -i genome/RefSeq_genome_from_KBase/GCF_001077695.1_assembly.fa --output genome/eggnog/eggnog_output --cpu 16