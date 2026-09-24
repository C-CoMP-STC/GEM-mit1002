import os
import sys

from auto_rbh import find_reciprocal_best_hits

BLAST_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(BLAST_DIR))  # make `tools` importable
from tools.paths import GENOME_DIR  # noqa: E402

final_results = find_reciprocal_best_hits(
    query_seq_file=os.path.join(BLAST_DIR, "query_ecoli_dut/dut_seq.fa"),
    organism_fasta=os.path.join(
        GENOME_DIR, "Michelle_4106_gene_calls/MIT1002_anvio_prot_seqs.fa"
    ),
    organism_db=os.path.join(BLAST_DIR, "dbs/amac_db"),
    reference_db=os.path.join(BLAST_DIR, "dbs/ecoli_db"),
    output_dir=os.path.join(BLAST_DIR, "query_ecoli_dut/results"),
    num_forward_hits=10,
)
