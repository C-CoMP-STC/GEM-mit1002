import cobra

def main():
    gff_file = 'genome/RefSeq_genome_from_KBase/KBase_derived_Alteromonas_macleodii.gff'
    model_file = 'model.xml'  # Replace with the actual path to your model file

    gene_to_protein = parse_gff(gff_file)
    model_w_protein_ids = update_model(model_file, gene_to_protein)

    cobra.io.write_sbml_model(model_w_protein_ids, 'model.xml')


def parse_gff(gff_file):
    gene_to_protein = {}
    with open(gff_file, 'r') as file:
        for line in file:
            if line.startswith('#') or line.strip() == '':
                continue
            columns = line.strip().split('\t')
            attributes = [s.strip() for s in columns[-1].split(';')]
            attr_dict = {attr.split('=')[0]: attr.split('=')[1] for attr in attributes}
            gene_id = attr_dict.get('ID')
            protein_id = attr_dict.get('protein_id')
            if gene_id and protein_id:
                gene_to_protein[gene_id] = protein_id
    return gene_to_protein

def update_model(model_file, gene_to_protein):
    model = cobra.io.read_sbml_model(model_file)
    gene_counter = 0
    protein_counter = 0
    for gene in model.genes:
        gene_counter += 1
        if gene.id in gene_to_protein:
            protein_counter += 1
            gene.annotation['refseq'] = gene_to_protein[gene.id]

    print(f'Updated {protein_counter} of {gene_counter} genes with protein IDs')

    return model


if __name__ == '__main__':
    main()