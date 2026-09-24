# Alteromonas pangenome

Pangenome data for *Alteromonas* used to add reactions to the MIT1002 model
(see [`code/pangenome/`](../../code/pangenome/)).

The data comes from [(Veseli, 2024)](https://www.nature.com/articles/s41597-024-03778-z):
> Veseli, I., DeMers, M. A., Cooper, Z. S., Schechter, M. S., Miller, S., Weber, L., ... & Braakman, R. (2024). Digital Microbe: a genome-informed data integration framework for team science on emerging model organisms. Scientific Data, 11(1), 967.

All data is available via via https://doi.org/10.5281/zenodo.7430118/

## Not in git

The raw files are not tracked here, due to their size and previous publication above.

To recreate the analyses, create a folder called `digital_microbe` here and in it, download the following files:

What the folder holds:

| File | What it is |
| --- | --- |
| `2738541267_aagenesequences.fa`, `2738541267_genecalls.txt` | MIT1002 protein sequences and gene calls (IMG genome 2738541267) |
| `Genome_Metadata.txt` | The genomes in the pangenome, with genome size and ANI to MIT1002 |
| `Database_MIT1002GeneCalls.csv` | Pathway-step database mapped to MIT1002 gene calls; the input to `code/pangenome/parse_pangenome_data.py` and `make_table.py` |
| `Heatmaps/`, `RawFiles/` | Per-pathway heatmaps with trees, and the underlying per-pathway tables |
| `annotated_gff.xlsx`, `michelle_kbase_comparison.xlsx` | Annotated gene calls, and a comparison against the KBase draft model |

`database_w_*.csv` and `clean_compariosn.*` in the same folder are **not**
raw data: `code/pangenome/make_table.py` and `clean_comparison.py` wrote
them there.
