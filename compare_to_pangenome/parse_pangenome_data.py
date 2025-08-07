# parse_pangenome_data.py
'''
Use blast results to map the "gene call IDs" to NCBI protein accession IDs so we
can figure out which of the genes identified in the pangenome analyses were and
weren't already in the model
'''

import pandas as pd
import cobra

# read in all input files
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
id_lookup = pd.read_csv('blast_stuff/blast_output_filtered.csv')
modelseed_rxns_all = pd.read_csv(
    'modelseed_reactions_all.tsv',
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
model = cobra.io.load_json_model('../model.json')

# start by mapping the arbitrary largely meaningless "Gene.Call.ID"s to
# NCBI protein accessions and dropping everything that didn't have a
# Gene.Call.ID or that we couldn't get a protein accession for cuz we can't do
# anything with those genes
pangenome_useful = pangenome_raw.merge(
    id_lookup, how = 'left', left_on = 'Gene.Call.ID', right_on = 'pangenome_id'
).dropna(subset = 'protein_accession').drop('Gene.Call.ID', axis = 1).rename(
    columns = {'protein_accession' : 'gene_id'}
)

# we currently have one row for each gene-reaction pair but want to collapse
# all the genes associated with each reaction but want to keep all the
# information we have about each gene, so split off a table with one row per
# gene and everything we have on it then group by KEGG reaction ID and only
# collapse the protein accessions into lists of associated proteins for each
# reaction
gene_info = pangenome_useful[[
    'gene_id', 'gene_symbol', 'kegg_ortholog'
]].drop_duplicates()
pangenome_by_rxns = pangenome_useful.drop(
    ['gene_symbol', 'kegg_ortholog'], axis = 1
).groupby('kegg_reaction', as_index = False).agg(
    lambda col: ';'.join(col.astype(str).unique())
).rename(columns = {'kegg_reaction' : 'kegg_id'})

# now get the ModelSEED table ready to merge on KEGG IDs
# extract KEGG IDs from the "aliases" column into their own column
modelseed_rxns_less = modelseed_rxns_all.dropna(subset = 'aliases').copy()
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
gene_info.to_csv('pangenome_gene_info.csv', index = False)
