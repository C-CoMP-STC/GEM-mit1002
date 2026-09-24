[![memote tested](https://img.shields.io/badge/memote-tested-blue.svg?style=plastic)](https://hgscott.github.io/mit1002-model)

# MIT1002-GEM: A manually curated metabolic model for *Alteromonas macleodii* MIT1002

This repo contains the *Alteromonas macleodii* MIT1002 model, and code associated with its creation, curation, and testing.

The model was generated using GenBank genome for MIT1002 (Accession Number: NZ_JXRW01000001), accessible via KBase.
The narrative for generating the draft model, is available here: https://narrative.kbase.us/narrative/208605

This repo uses GtiHub actions to automatically test the model.
Upon every push, pull request, manual trigger:
1. A new MEMOTE report is generated, and saved as "index.html"
2. Run custom tests
    * Validate the SBML file
    * Test for growth on no carbon sources
    * Test known growth phenotypes, and regenerate the experimental vs predicted growth heatmap figure
    * Run the MEMOTE test to search for ATP generating cycles
    * Check that nothing in the deprecated identifier lists is back in the model
3. The model is exported to JSON and excel formats

Note: MACAW is **not** run as part of the action due to the longer run time of the dilution test.
To run MACAW use:
```
python code/scripts/run_macaw.py
```
Its results are written to `code/scripts/results/` with the other generated reports.

## Repository layout

Three directories hold code, and the distinction between them is about *what the
code does*, not what it is about. Please put new code in the matching one.

| Directory | Contains | How it runs |
| --- | --- | --- |
| `code/test/` | Checks that assert something about the model and pass or fail | Automatically, via `pytest` in CI. A failure blocks the PR |
| `code/scripts/` | Code that generates an artifact for a person to look at — a table, a plot, an exported file. No pass/fail | Automatically in CI, writing to `code/scripts/results/` |
| `code/tools/` | Importable functions and definitions, and command-line utilities a curator runs deliberately | By hand, or imported by the above |
| `data/` | Experimental observations, media provenance, and derived tables. See [`data/README.md`](data/README.md) | Read by the above |

Examples of the third kind:

* `code/tools/deprecate.py` — you invoke it yourself when removing a reaction, and
  `code/scripts/export_model.py` and `code/test/test_deprecated.py` both import from it
* `code/tools/media.py` — defines the growth media as importable dictionaries, so
  anything needing a medium does `from tools.media import MEDIA`
* `code/tools/plot_styles.py` — shared colour palettes and figure styling, imported by
  every plotting script so figures stay consistent

Two other directories hold code that is neither of these: `code/curation_process/`
analyses the history of the curation effort itself across past PRs, and
`code/biomass/`, `data/genome/`, `code/escher/` and similar hold the exploratory work behind
particular parts of the model.

Note that [standard-GEM](https://github.com/MetabolicAtlas/standard-GEM), which
yeast-GEM and Human-GEM follow, asks for a single `code/` directory instead. The
split above is a deliberate refinement of that; `code/` is also a poor Python
package name because it shadows a standard-library module.

## Setting Up the Environment
To ensure a smooth setup and avoid system conflicts, follow these steps to create and activate a Python virtual environment before installing dependencies.

1. Check Your Python Version
First, make sure you have Python 3.11 or 3.10 installed.
Run the following command to check your version:
```
python3 --version
```
If it shows Python 3.13, we strongly recommend using Python 3.11 or 3.10, as some dependencies may not yet support Python 3.13.

To install Python 3.11 via Homebrew (if needed), run:
```
brew install python@3.11
```
2. Create a Virtual Environment

Once you have the correct Python version, create a virtual environment:
```
python3.11 -m venv .venv  # Use python3.10 if needed
```
This creates a .venv folder in your project directory, isolating dependencies from the system Python.

3. Activate the Virtual Environment

Before installing packages, activate the virtual environment:

Mac/Linux:
```
source .venv/bin/activate
```
Windows (Command Prompt):
```
.venv\Scripts\activate
```
Windows (PowerShell):
```
.\.venv\Scripts\Activate
```
Your terminal should now show (.venv) at the beginning of the prompt, indicating the environment is active.

4. Upgrade Pip

Inside the virtual environment, upgrade pip to avoid compatibility issues:
```
pip install --upgrade pip setuptools
```
5. Install Dependencies

Now, install all required dependencies:
```
pip install -r requirements.txt
```
If you encounter an error about "externally managed environment", add the following flag:
```
pip install -r requirements.txt --break-system-packages
```
6. Verify Installation

To ensure everything is working correctly, run:
```
python --version  # Should be 3.11 or 3.10
pip list  # Should show installed dependencies
```
If everything looks good, you're ready to start using the project! 🎉
