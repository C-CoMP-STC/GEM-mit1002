import os
import subprocess
import sys
from io import StringIO

import pandas as pd
from Bio import SeqIO

# --- Configuration (Set these file paths and IDs based on your project) ---
FASTA_FILE = "genome/Michelle_4106_gene_calls/MIT1002_anvio_prot_seqs.fa"
ORG_BLAST_DB = "BLAST/dbs/amac_db"
REF_BLAST_DB = "BLAST/dbs/ecoli_db"
# This file contains the single E. coli dUTPase protein sequence
QUERY_SEQ_FILE = "BLAST/rbh_for_dut/dut_seq.fa"
BLAST_OUTPUT = "BLAST/rbh_for_dut/FBH_candidate.tsv"
REVERSE_BLAST_OUTPUT = "BLAST/rbh_for_dut/RBH_check.tsv"
# Output file for the extracted protein sequence of the FBH
TARGET_OUTPUT = "BLAST/rbh_for_dut/AMAC_dut_candidate.fa"
TARGET_ID = None
n_hits = 10  # Number of top hits to retrieve


def main():
    # --- Run the Forward BLAST ---
    print("Running Forward BLAST to find the FBH...")
    candidate_hits_df = get_top_blast_hits_local(
        query_file=QUERY_SEQ_FILE, db_path=ORG_BLAST_DB, num_hits=n_hits
    )
    # Save the forward BLAST results
    candidate_hits_df.to_csv(BLAST_OUTPUT, sep="\t", index=False)
    print(f"Forward BLAST results saved to {BLAST_OUTPUT}")

    # --- Parse the TSV to get the Forward Best Hit (FBH) ID ---
    print("Parsing BLAST output to identify the Locus Tag...")
    if candidate_hits_df.empty:
        print("No hits found in BLAST search. Exiting.")
        sys.exit(1)
    else:
        print(f"{len(candidate_hits_df)} BLAST hits found:")
        print(candidate_hits_df)
        TARGET_ID = candidate_hits_df["sseqid"].astype(str).tolist()
        print(f"Identified candidate IDs: {TARGET_ID}")

    # --- Extract the sequence using Biopython ---
    if TARGET_ID:
        record_dict = SeqIO.index(FASTA_FILE, "fasta")
        with open(TARGET_OUTPUT, "w") as out_f:
            for tid in TARGET_ID:
                print(f"Extracting sequence for {tid}...")
                SeqIO.write(record_dict[tid], out_f, "fasta")
                print(f"Success! Sequence for {tid} saved to {out_f.name}")
    else:
        print("No TARGET_ID found. Exiting without extracting sequences.")
        sys.exit(1)

    # --- Run the Reverse BLAST ---
    print("Running Reverse BLAST to confirm the FBH...")
    reverse_hits_df = get_top_blast_hits_local(
        query_file=TARGET_OUTPUT,
        db_path=REF_BLAST_DB,
        num_hits=1,  # For the reverse, really only want the top hit
    )
    # Save the reverse BLAST results
    reverse_hits_df.to_csv(REVERSE_BLAST_OUTPUT, sep="\t", index=False)
    print(f"Reverse BLAST results saved to {REVERSE_BLAST_OUTPUT}")


# Helper Functions
def get_top_blast_hits_local(
    query_file: str, db_path: str, num_hits: int = 10
) -> pd.DataFrame:
    """
    Runs a local BLAST search and returns the top N hits as a pandas DataFrame.

    Args:
        query_file: Path to the input FASTA file.
        db_path: Path to the BLAST database.
        num_hits: The maximum number of hits to return.

    Returns:
        A pandas DataFrame containing the top BLAST hits, or an empty DataFrame
        if no hits are found.
    """
    print(f"Running local BLAST for {query_file}... requesting top {num_hits} hits.")

    # 1. Define the column headers for the tabular output format.
    # This makes the code easier to read and maintain.'
    outfmt_spec = "6 qseqid sseqid pident length evalue bitscore"
    outfmt_cols = outfmt_spec.split(" ")[1:]  # Skip the '6'

    # 2. Build the BLAST command with the specified number of hits.
    cmd = [
        "blastp",  # or blastn, etc.
        "-query",
        query_file,
        "-db",
        db_path,
        "-outfmt",
        outfmt_spec,
        "-max_target_seqs",
        str(num_hits),  # Use the num_hits parameter here
    ]

    try:
        # Execute the BLAST command
        result = subprocess.run(
            cmd, capture_output=True, text=True, check=True, encoding="utf-8"
        )

        # Parse the multi-line output directly into a pandas DataFrame.
        # Using StringIO is efficient as it treats the string output as a file.
        blast_output = result.stdout
        if not blast_output.strip():
            print("No BLAST hits found.")
            return pd.DataFrame(columns=outfmt_cols)

        df = pd.read_csv(
            StringIO(blast_output), sep="\t", header=None, names=outfmt_cols
        )
        return df

    except subprocess.CalledProcessError as e:
        print(f"An error occurred during BLAST execution: {e}")
        print(f"BLAST stderr:\n{e.stderr}")
        return pd.DataFrame(columns=outfmt_cols)


if __name__ == "__main__":
    main()
