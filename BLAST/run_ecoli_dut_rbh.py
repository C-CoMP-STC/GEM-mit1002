import os

from auto_rbh import find_reciprocal_best_hits

# Assume the script is in GEM-mit1002/BLAST/
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

final_results = find_reciprocal_best_hits(
    query_seq_file=os.path.join(PROJECT_ROOT, "BLAST/ecoli_dut_query/dut_seq.fa"),
    organism_fasta=os.path.join(PROJECT_ROOT, "genome/Michelle_4106_gene_calls/MIT1002_anvio_prot_seqs.fa"),
    organism_db=os.path.join(PROJECT_ROOT, "BLAST/dbs/amac_db"),
    reference_db=os.path.join(PROJECT_ROOT, "BLAST/dbs/ecoli_db"),
    output_dir=os.path.join(PROJECT_ROOT, "BLAST/ecoli_dut_query/results"),
    num_forward_hits=10,
)