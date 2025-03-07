[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://hgscott.github.io/mit1002-model)

# MIT1002-model

This repo contains the *Alteromonas macleodii* MIT1002 model, and code associated with its creation, curation, and testing.

The model was generated using GenBank genome for MIT1002 (Accession Number: NZ_JXRW01000001), accessible via KBase.
The narrative for generating the draft model, is available here: https://narrative.kbase.us/narrative/208605

This repo uses GtiHub actions to automatically test the model.
Upon every push, pull request, manual trigger:
1. A new MEMOTE report is generated, and saved as "index.html"
2. The MACAW tests are run. Results are saved in `macaw_results.csv` and `macaw_edge_list.csv`
3. Run custom tests
    * Validate the SBML file
    * Test for growth on no carbon sources
    * Test known growth phenotypes, and regenerate the experimental vs predicted growth heatmap figure
    * Run the MEMOTE test to search for ATP generating cycles
3. The model is exported to JSON and excel formats

## To contribute to the model
1. Make a GitHub account
2. Make a fork/branch of this repo
3. Make your edits
4. Open a pull request
