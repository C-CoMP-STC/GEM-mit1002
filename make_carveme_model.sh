#!/bin/bash -l
#$ -P cometsfba
#$ -N carveme
#$ -pe omp 16
#$ -o make_carveme_model.out
#$ -j y  # Merge the error and output streams into a single file

module load python3/3.6.12
module load cplex
module load diamond
module load carveme/1.5.1

# Run CarveMe
carve --egg genome/eggnog_output.emapper.annotations \
    -g mbm_glc__D,l1_glc__D,mbm_ac,l1_ac,mbm_ala__L,l1_ala__L,mbm_pro__L,l1_pro__L,l1_glyc3p,l1_tyr__L,l1_glu__L,l1_val__L \
    --mediadb test/test_files/media/media_database.tsv \
    --fbc2 \
    -o carveme_model.xml \
    -v