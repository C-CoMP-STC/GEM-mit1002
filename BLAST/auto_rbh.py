import os
import subprocess
import sys
from io import StringIO

import pandas as pd
from Bio import SeqIO


def find_reciprocal_best_hits(
    query_seq_file: str,
    organism_fasta: str,
    organism_db: str,
    reference_db: str,
    output_dir: str,
    num_forward_hits: int = 10,
) -> pd.DataFrame:
    """
    Performs a forward and reverse BLAST to find reciprocal best hits.
    All file paths should be absolute or relative to the script's execution directory.
    """
    # --- Setup Paths ---
    os.makedirs(output_dir, exist_ok=True)
    forward_blast_output = os.path.join(output_dir, "forward_blast_hits.tsv")
    candidate_seq_output = os.path.join(output_dir, "candidate_sequences.fa")
    reverse_blast_output = os.path.join(output_dir, "reverse_blast_hits.tsv")

    # --- 1. Run the Forward BLAST ---
    print("--- Step 1: Running Forward BLAST ---")
    candidate_hits_df = get_top_blast_hits_local(
        query_file=query_seq_file, db_path=organism_db, num_hits=num_forward_hits
    )
    candidate_hits_df.to_csv(forward_blast_output, sep="\t", index=False)
    print(f"Forward BLAST results saved to {forward_blast_output}")

    if candidate_hits_df.empty:
        print("No hits found in forward BLAST search. Exiting.")
        return pd.DataFrame()

    candidate_ids = candidate_hits_df["sseqid"].astype(str).tolist()

    # --- 2. Extract Candidate Sequences ---
    print("\n--- Step 2: Extracting Candidate Sequences ---")
    try:
        record_dict = SeqIO.index(organism_fasta, "fasta")
        with open(candidate_seq_output, "w") as out_f:
            for tid in candidate_ids:
                if tid in record_dict:
                    SeqIO.write(record_dict[tid], out_f, "fasta")
                else:
                    print(f"Warning: ID {tid} not found in {organism_fasta}. Skipping.")
        print(f"Candidate sequences saved to {candidate_seq_output}")
    except FileNotFoundError:
        print(f"Error: Fasta file not found at {organism_fasta}", file=sys.stderr)
        return pd.DataFrame()

    # --- 3. Run the Reverse BLAST ---
    print("\n--- Step 3: Running Reverse BLAST ---")
    reverse_hits_df = get_top_blast_hits_local(
        query_file=candidate_seq_output,
        db_path=reference_db,
        num_hits=1,
    )
    reverse_hits_df.to_csv(reverse_blast_output, sep="\t", index=False)
    print(f"Reverse BLAST results saved to {reverse_blast_output}")

    # --- 4. Combine and Return Results ---
    print("\n--- Step 4: Combining Results ---")
    original_query_id = SeqIO.read(query_seq_file, "fasta").id
    final_df = pd.merge(
        candidate_hits_df,
        reverse_hits_df,
        left_on="sseqid",
        right_on="qseqid",
        how="left",
        suffixes=("_forward", "_reverse"),
    )
    final_df["is_rbh"] = final_df["sseqid_reverse"] == original_query_id
    print(final_df.to_markdown(index=False))
    return final_df


def get_top_blast_hits_local(
    query_file: str, db_path: str, num_hits: int = 10
) -> pd.DataFrame:
    """
    Runs a local BLAST search and returns the top N hits as a pandas DataFrame.
    """
    outfmt_spec = "6 qseqid sseqid pident length evalue bitscore"
    outfmt_cols = outfmt_spec.split(" ")[1:]
    cmd = [
        "blastp",
        "-query",
        query_file,
        "-db",
        db_path,
        "-outfmt",
        outfmt_spec,
        "-max_target_seqs",
        str(num_hits),
    ]
    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, encoding="utf-8"
        )
        blast_output = result.stdout
        if not blast_output.strip():
            return pd.DataFrame(columns=outfmt_cols)
        return pd.read_csv(
            StringIO(blast_output), sep="\t", header=None, names=outfmt_cols
        )
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"An error occurred during BLAST execution: {e}", file=sys.stderr)
        if hasattr(e, "stderr"):
            print(f"BLAST stderr:\n{e.stderr}", file=sys.stderr)
        return pd.DataFrame(columns=outfmt_cols)
