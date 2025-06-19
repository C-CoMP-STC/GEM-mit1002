# parse_pangenome_data.py
'''
Figure out which genes in the pangenome were annotated as encoding enzymes that
catalyze reactions that aren't already in the model and add them to the model
Since the pangenome analysis identified genes using different IDs than the ones
already used in the model, also have to map gene/protein IDs
'''

import pandas as pd
import cobra

# read in all input files
gff_raw = pd.read_csv(
    '/projectnb/cometsfba/hscott/GEM-repos/GEM-mit1002/genome/' +\
    'RefSeq_genome_from_KBase/KBase_derived_Alteromonas_macleodii.gff',
    sep = '\t', header = None, comment = '#'
)
eggnog_raw = pd.read_csv(
    '/projectnb/cometsfba/hscott/GEM-repos/GEM-mit1002/genome/RefSeq_genome/'+\
    'eggnog/eggnog_output.emapper.annotations',
    sep = '\t', skiprows = 4, skipfooter = 3, engine = 'python'
)
pangenome_raw = pd.read_csv('Database_MIT1002GeneCalls.csv').drop(
    [
        'Inferred Presence', 'StepName', 'Notes', 'EnzymeName', 'Name',
        'Definition'
    ], axis = 1
).rename(columns = {
    'KO' : 'kegg_ortholog',
    'KEGG.Reaction' : 'kegg_reaction',
    'KEGG.Enzyme' : 'ecs_pangenome',
    'Gene' : 'gene_symbol'
})
modelseed_rxns_all = pd.read_csv(
    'all_modelseed_reactions.tsv',
    sep = '\t',
    usecols = [
        'id', 'aliases', 'name', 'stoichiometry', 'definition',
        'reversibility', 'ec_numbers', 'status', 'is_obsolete'
    ]
).rename(columns = {
    'id' : 'modelseed_id',
    'stoichiometry' : 'rxn_def',
    'definition' : 'equation',
    'reversibility' : 'rev_dir',
    'ec_numbers' : 'ecs_modelseed'
})
model = cobra.io.load_json_model('model.json')

# use the GFF file to make a dict of protein IDs (which are the primary IDs in
# the eggNOG file) to the gene IDs used in the model
prot_to_gene = {'gene_id' : list(), 'protein_id' : list()}
for row in gff_raw[8].str.split('; '):
    # skip rows that don't have a protein_id
    if any(bit.startswith('protein_id') for bit in row):
        # the first bit is always the gene ID
        gene_id = row[0].split('ID=')[1]
        # also skip rows that have a gene ID that contains "CDS" cuz those
        # always seem to be have duplicates
        if 'CDS' not in gene_id:
            prot_to_gene['gene_id'].append(gene_id)
            for bit in row[1:]:
                if bit.startswith('protein_id'):
                    # strip the label off of the actual ID
                    prot_to_gene['protein_id'].append(
                        bit.split('protein_id=')[1]
                    )

# use the eggNOG file to make a dict of KEGG ortholog IDs to the protein IDs
# used in the GFF file
# skip rows without a KEGG ortholog ID
eggnog_useful = eggnog_raw[eggnog_raw['KEGG_ko'] != '-']
ko_to_prot = pd.DataFrame({
    'protein_id' : eggnog_useful['#query'].to_list(),
    # strip the ko: prefixes
    'kegg_ortholog' : eggnog_useful['KEGG_ko'].str.lstrip('ko:').to_list()
})

# the pangenome file sometimes had multiple KEGG IDs in a single row; split
# those up so that each row has only one KEGG ID so that we can merge on them
pangenome_raw['kegg_reaction'] = pangenome_raw['kegg_reaction'].str.split(';')
pangenome_raw = pangenome_raw.explode('kegg_reaction').drop_duplicates()

# merge on KEGG ortholog IDs and protein_ids
pangenome_useful = pangenome_raw.merge(
    ko_to_prot.merge(pd.DataFrame(prot_to_gene)), how = 'left'
)

# since we currently have one row for each gene-reaction pair, make a separate
# table that just has one row for each unique gene and the information about
# that gene so that we can group pangenome_useful by reaction ID and only have
# to collapse the gene IDs
gene_info = pangenome_useful[[
    'gene_id', 'gene_symbol', 'protein_id', 'kegg_ortholog', 'Gene.Call.ID'
]].drop_duplicates()

# now drop all the columns we just put over in gene_info and drop all rows with
# missing gene_ids since we can't actually add those genes to the model
pangenome_by_rxns = pangenome_useful.dropna(subset = 'gene_id').drop(
    ['protein_id', 'gene_symbol', 'kegg_ortholog', 'Gene.Call.ID'], axis = 1
).groupby('kegg_reaction', as_index = False).agg(
    # only gene_id and ecs_from_pangenome have multiple unique values
    lambda col: ';'.join(col.astype(str).unique())
).rename(columns = {'kegg_reaction' : 'kegg_id'})

# prepare table from ModelSEED for merging with pangenome data
# extract KEGG IDs from the "aliases" column into their own column
modelseed_rxns_less = modelseed_rxns_all.dropna(subset = 'aliases')
modelseed_rxns_less['kegg_id'] = modelseed_rxns_less['aliases'].apply(lambda x:
    x.split('KEGG: ')[1].split('|')[0].split(';') if 'KEGG' in x else ''
)
# ModelSEED sometimes associated one reaction with multiple KEGG IDs
modelseed_rxns_kegg = modelseed_rxns_less.explode('kegg_id').drop_duplicates()
# now drop everything without a KEGG ID or balance issues
modelseed_useful = modelseed_rxns_kegg[
    (modelseed_rxns_kegg['kegg_id'] != '') &
    (modelseed_rxns_kegg['status'] == 'OK') &
    (modelseed_rxns_kegg['is_obsolete'] == 0)
].drop(
    ['status', 'is_obsolete'], axis = 1
# ModelSEED often has multiple identical copies of the same reaction, so group
# by KEGG IDs and only keep the one with the lowest ID number
).sort_values('modelseed_id').groupby('kegg_id', as_index = False).first()

# now merge in the ModelSEED info
full_df = pangenome_by_rxns.merge(modelseed_useful, how = 'left')

# combine the lists of E.C. numbers from ModelSEED and the pangenome analysis
full_df['ec_numbers'] = full_df[[
    # make the NAs in these two columns empty strings so that they just get
    # skipped instead of putting "nan" in the lists of E.C. numbers
    'ecs_modelseed', 'ecs_pangenome'
]].fillna('').apply(
    lambda row: ';'.join(
        # replace |s with ;s in ecs_modelseed before appending
        [row['ecs_pangenome']] + [
            e for e in row['ecs_modelseed'].split('|')
            if e not in row['ecs_pangenome']
        ]
    ), axis = 1
)

# add a column indicating if the reaction is already in the model
in_model = [r.id[:-3] for r in model.reactions]
full_df['already_in_model'] = full_df['modelseed_id'].apply(
    lambda x: 'yes' if x in in_model else 'no'
)

# write to a file so we can manually curate before actually adding anything to
# the model
full_df = full_df[[
    'modelseed_id', 'kegg_id', 'equation', 'gene_id', 'ec_numbers',
    'already_in_model', 'name', 'rxn_def', 'rev_dir'
]]

full_df.to_csv('pangenome_rxns_uncurated.csv', index = False)

# also write the gene_info file from earlier
gene_info.to_csv('pangenome_new_genes.csv', index = False)
