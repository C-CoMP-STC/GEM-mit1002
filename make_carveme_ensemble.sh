#!/bin/bash -l
#$ -P cometsfba
#$ -N carveme_ensemble
#$ -pe omp 16
#$ -o make_carveme_ensemble.out
#$ -j y  # Merge the error and output streams into a single file

module load python3/3.6.12
module load cplex

# Unload the python module so that I can load my own conda environment
module unload python3/3.6.12

# Activate the CarveMe environment
module load miniconda
conda activate /projectnb/cometsfba/hscott/envs/carveme-dev


# Run CarveMe on the eggNOG annotations
# carve --egg genome/clean_eggnog_output.emapper.annotations \
#     -o carveme_ensemble.xml \
#     -v --debug -n 60

# Run CarveMe with Diamond on the protein sequences directly
carve genome/MIT1002_anvio_prot_seqs.faa \
    -o carveme_ensemble_diamond.xml \
    -v --debug -n 60