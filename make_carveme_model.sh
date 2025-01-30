#!/bin/bash -l

module load python3/3.6.12
module load cplex
module load diamond
module load carveme/1.5.1

# Run CarveMe
carve --egg eggnog_output.emapper.annotations -g LB --fbc2