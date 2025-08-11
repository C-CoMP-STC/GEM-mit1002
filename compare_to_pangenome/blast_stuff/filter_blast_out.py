# filter_blast_out.py
'''
Filter out all of the hits blast gave us with really low sequence identities
'''

import pandas as pd

all_hits = pd.read_csv('blast_output.csv')

# get the highest-scoring hit(s) for each query sequence, then get the hightest
# scoring hit with the highest sequence identity, then the lowest e-value. only
# about a dozen proteins still have multiple best hits after all of that, so
# just drop those
all_best_hits = all_hits.groupby('pangenome_id', as_index = False).apply(
    lambda g: g[g['score'] == g['score'].max()]
).groupby('pangenome_id', as_index = False).apply(
    lambda g: g[g['seq_identity_pct'] == g['seq_identity_pct'].max()]
).groupby('pangenome_id', as_index = False).apply(
    lambda g: g[g['e_value'] == g['e_value'].min()]
).groupby('pangenome_id', as_index = False).filter(
    lambda g: len(g) == 1
).reset_index(drop = True)
# also drop the 24 proteins whose best hit only had a sequence identity of less
# than 90%, because it turns out that the next-highest sequence identity after
# 90% is only 50%, so those are all garbage
good_best_hits = all_best_hits[all_best_hits['seq_identity_pct'] > 90]

good_best_hits[['pangenome_id', 'protein_accession']].to_csv(
    'blast_output_filtered.csv', index = False
)
