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
    new_mets, model, all_met_info
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
        # see if this metabolite is already in the model
        (coeff, met_id, _, _, met_name) = bit.split(':')
        try:
            met_obj = model.metabolites.get_by_id(met_id + '_c0')
        except KeyError:
            # see if we've already made a metabolite object for this and just
            # haven't added it to the model yet (cuz adding them as we make them
            # is much slower than adding them all in one go)
            try:
                met_obj = [m.id for m in new_mets if m.id == met_id][0]
            except IndexError:
                met_obj = make_new_metabolite(met_id, all_met_info)
                new_mets.append(met_obj)
        if direction == '<':
            met_dict[met_obj] = -float(coeff)
        else:
            met_dict[met_obj] = float(coeff)
    rxn.add_metabolites(met_dict)
    rxn.gene_reaction_rule = gpr
    rxn.annotation = annotations
    return((rxn, new_mets))

all_rxns = pd.read_csv('pangenome_rxns_curated.csv')
all_met_info = pd.read_csv(
    'modelseed_metabolites_all.tsv', sep = '\t', low_memory = False
)
all_gene_info = pd.read_csv('pangenome_gene_info.csv')
model = cobra.io.load_json_model('../model.json')

rxns_to_make = all_rxns[all_rxns['manual_curation'] == 'pass']
# adding each reaction object as we create it is slower than making them all
# then adding them all in one go
rxns_to_add = list()
new_mets = list()
for (idx, row) in rxns_to_make.iterrows():
    if ';' in row['ec_numbers']:
        row['ec_numbers'] = row['ec_numbers'].split(';')
    # if a reaction doesn't have a name, just use its ID as its name
    if isinstance(row['name'], float):
        row['name'] = row['modelseed_id']
    # only add the c0 suffix to real ModelSEED IDs and not our ad-hoc IDs
    if row['modelseed_id'].startswith('rxn'):
        row['modelseed_id'] += '_c0'
    (new_rxn, new_mets) = make_new_reaction(
        row['modelseed_id'],
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
        new_mets,
        model,
        all_met_info
    )
    rxns_to_add.append(new_rxn)
model.add_reactions(rxns_to_add)
# all new metabolite objects will be added along with the reactions; we just
# kept that list for our convenience

# also fix the metabolite IDs that you accidentally left in the names a while
# ago
for m in model.metabolites:
    if m.name.endswith(f'({m.id})'):
        old_name = m.name
        m.name = m.name.split(f' ({m.id})')[0]

print('metabolite ID | name | formula | charge\n---|---|---|---')
for m in new_mets:
    print(f'{m.id} | {m.name} | {m.formula} | {m.charge}')

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
print('\nreaction ID | equation | GPR\n---|---|---')
for r in rxns_to_add:
    print(f'{r.id} | {r.build_reaction_string(True)} | {r.gene_name_reaction_rule} ({r.gene_reaction_rule})')

cobra.io.save_json_model(model, '../model.json')
