# Start a result file with the header line
echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > query_pg1_dut/results/forward_blast_hits.tsv

# Run BLASTp of PG1 dut sequence against the MIT1002 database
blastp -query query_pg1_dut/pg1_dut_seq.fa -db dbs/amac_db  -outfmt 6 -max_target_seqs 10 >> query_pg1_dut/results/forward_blast_hits.tsv