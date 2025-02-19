#!~/opt/miniconda3/envs/med4-hot1a3/bin/python
import os

import pandas

FILE_DIR = os.path.dirname(os.path.abspath(__file__))


def main():
    # Set the input files
    gff_file = os.path.join(FILE_DIR, "KBase_derived_mit1002_rast.gff")
    fasta_file = os.path.join(
        FILE_DIR, "Michelle's 4106 gene calls", "MIT1002_anvio_prot_seqs.fa"
    )

    # Convert the GFF file to a table
    df = convert_gff_to_table(gff_file, fasta_file)

    # Save the DataFrame to a CSV file
    output_file = gff_file.replace(".gff", ".csv")
    df.to_csv(output_file, index=False)


def convert_gff_to_table(gff_file: str, fasta_file: str) -> pandas.DataFrame:
    """
    Convert the GFF file (the RAST output exported from KBase) to a table
    format as a pnadas dataframe with the following columns:
    1. Gene ID (matching the fasta file)
    2. Sequence (from the fasta file)
    3-?: The annotaionts from RAST, split into different columns

    Parameters
    ----------
    gff_file : str
        The path to the GFF file

    Returns
    -------
    pandas.DataFrame
        The DataFrame with the gene IDs, sequences, and annotations
    """

    # Read the GFF file
    with open(gff_file, "r") as f:
        lines = f.readlines()

    # Get the gene IDs and sequences from the fasta file
    gene_ids = []
    sequences = []
    with open(fasta_file, "r") as f:
        current_seq_lines = []
        current_gene_id = None
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_gene_id is not None:
                    # Append the accumulated sequence for the previous gene
                    sequences.append("".join(current_seq_lines))
                # Start a new gene entry
                current_gene_id = line[1:]
                gene_ids.append(current_gene_id)
                current_seq_lines = []
            else:
                current_seq_lines.append(line)
        # Don't forget to add the last sequence
        if current_gene_id is not None:
            sequences.append("".join(current_seq_lines))

    # Get the annotations from the GFF file
    product_dict = {}
    sso_dict = {}
    for line in lines:
        if line.startswith("##FASTA"):
            break
        if not line.startswith("#"):
            parts = line.strip().split("\t")
            # TODO: Remove magic numbers
            # Get the gene ID as the number that comes after "ID"
            gene_id = parts[8].split(";")[0].split("=")[1].split('_')[0]
            product = parts[8].split(";")[1].split('=')[1]
            sso = parts[8].split(";")[2].split('=')[1].split(':')[1]
            # TODO: Remove duplicative code
            # Save the product for the gene
            if gene_id in product_dict:
                product_dict[gene_id].append(product)
            else:
                product_dict[gene_id] = [product]
            # Save the sso for the gene
            if gene_id in sso:
                sso_dict[gene_id].append(sso)
            else:
                sso_dict[gene_id] = [sso]

    # Sort the products and the sso in a list so the genes are in the same order as the fasta file
    product_list = []
    sso_list = []
    for gene_id in gene_ids:
        product_list.append(products.get(gene_id, ""))
        sso_list.append(sso.get(gene_id, ""))

    # Create a DataFrame with the gene IDs, sequences, products, and sso
    data = {"Gene ID": gene_ids, "Sequence": sequences}
    for gene_id, product in products.items():
        data[gene_id] = product
    for gene_id, sso in sso.items():
        data[gene_id] = sso

    df = pandas.DataFrame(data)

    return df


if __name__ == "__main__":
    main()
