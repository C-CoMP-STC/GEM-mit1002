#!/bin/bash -l
#$ -m e
#$ -j y
#$ -l h_rt=72:00:00
#$ -pe omp 8
#$ -t 1-20

module load blast+/2.12.0

psiblast -db $SCC_BLAST_DATA/nr \
  -query batched_input_seqs/batch_$SGE_TASK_ID.fa \
  -out batched_output_files/batch_$SGE_TASK_ID.csv \
  -outfmt "10 qseqid saccver evalue score pident staxids" \
  -taxids "28108" \
  -num_threads $NSLOTS 
