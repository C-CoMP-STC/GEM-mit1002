# pangenome

Pangenome data for *Alteromonas* used to add reactions to the MIT1002 model
(see [`code/pangenome/`](../../code/pangenome/)).

## Not in git

The raw files came from Michelle and live in `Pangenome from Michelle/` in
this folder, which is gitignored. Anyone reproducing the pangenome step needs
to get them from her.

<!-- TODO(Helen): add Michelle's full name, the date you received the files, and how. File dates suggest Feb–Apr 2024. -->

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
