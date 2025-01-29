#!/bin/bash -l

module load python3/3.6.12
module load cplex
module load diamond
module load carveme/1.5.1

carve test.faa -g LB