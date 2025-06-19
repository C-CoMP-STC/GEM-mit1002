# add_pangenome_rxns.py
'''
Read in pangenome_rxns_curated.csv, filter down to reactions that passed manual
curation, then add them to the model (along with all necessary metabolites)
'''

import pandas as pd
import cobra

def make_new_metabolite(met_id, all_met_info):
    this_met = all_met_info[all_met_info['id'] == met_id]
    met = cobra.Metabolite(
        met_id + '_c0',
        name = this_met['name'].iloc[0] + ' [c0]',
        compartment = 'c0',
        formula = this_met['formula'].iloc[0],
        charge = int(this_met['charge'].iloc[0])
    )
    met.annotation = {
        'sbo' : 'SBO:0000247',
        'seed.compound' : met_id,
    }
    return(met)

def make_new_reaction(
    rxn_id, rxn_name, direction, definition, gpr, annotations,
    model, all_met_info
):
    rxn = cobra.Reaction(rxn_id, name = rxn_name)
    if direction in ['=', '?']:
        rxn.lower_bound = -1000
        rxn.upper_bound = 1000
    elif direction in ['>', '<']:
        # going to flip all reactions that can only go backwards later
        rxn.lower_bound = 0
        rxn.upper_bound = 1000
    met_dict = dict()
    for bit in definition.split(';'):
        # see if this metabolite is already in the model or if we have to add it
        (coeff, met_id, _, _, met_name) = bit.split(':')
        try:
            met_obj = model.metabolites.get_by_id(met_id + '_c0')
        except KeyError:
            met_obj = make_new_metabolite(met_id, all_met_info)
            print(f'{met_obj.id} | {met_obj.name} | {met_obj.formula} | {met_obj.charge}')
        if direction == '<':
            met_dict[met_obj] = -float(coeff)
        else:
            met_dict[met_obj] = float(coeff)
    rxn.add_metabolites(met_dict)
    rxn.gene_reaction_rule = gpr
    rxn.annotation = annotations
    return(rxn)

all_rxns = pd.read_csv('pangenome_rxns_curated.csv')
all_met_info = pd.read_csv('all_modelseed_metabolites.tsv', sep = '\t')
all_gene_info = pd.read_csv('pangenome_new_genes.csv')
model = cobra.io.load_json_model('../model.json')

rxns_to_make = all_rxns[(
    (all_rxns['manual_curation'] == 'ok') |
    all_rxns['manual_curation'].str.startswith('fixed')
)]
# adding each reaction object as we create it is slower than making them all
# then adding them all in one go
rxns_to_add = list()
print('metabolite ID | name | formula | charge\n---|---|---|---')
for (idx, row) in rxns_to_make.iterrows():
    if ';' in row['ec_numbers']:
        row['ec_numbers'] = row['ec_numbers'].split(';')
    rxns_to_add.append(make_new_reaction(
        row['modelseed_id'] + '_c0',
        row['name'] + ' [c0]',
        row['rev_dir'],
        row['rxn_def'],
        row['gene_id'],
        {
            'sbo' : 'SBO:0000176',
            'seed.reaction' : row['modelseed_id'],
            'kegg.reaction' : row['kegg_id'],
            'ec-code' : row['ec_numbers']
        },
        model,
        all_met_info
    ))
model.add_reactions(rxns_to_add)

# now see if any of the genes associated with those reactions were new (and thus
# only have IDs and no names)
useful_gene_info = all_gene_info[
    # if there were multiple gene symbols for a particular gene, don't bother
    ~all_gene_info['gene_symbol'].str.contains(';')
].dropna(subset = ['gene_id', 'gene_symbol'])
gene_symbol_lookup = dict(zip(
    useful_gene_info['gene_id'], useful_gene_info['gene_symbol']
))
print('\ngene ID | old name | new name\n---|---|---')
for g in model.genes:
    if (g.name == '') or (g.name.startswith('G_')):
        try:
            old_name = g.name
            g.name = gene_symbol_lookup[g.id]
            print(f'{g.id} | {old_name} | {g.name}')
        except KeyError:
            pass
print()

# also fix the metabolite IDs that you accidentally left in the names a while
# ago
for m in model.metabolites:
    if m.name.endswith(f'({m.id})'):
        old_name = m.name
        m.name = m.name.split(f' ({m.id})')[0]

print('\nreaction ID | equation | GPR\n---|---|---')
for r in rxns_to_add:
    print(f'{r.id} | {r.build_reaction_string(True)} | {r.gene_name_reaction_rule} ({r.gene_reaction_rule})')

cobra.io.save_json_model(model, '../model.json')
