import sys

import cobra
from cobrakbase.core.kbasefba.newmodeltemplate_builder import NewModelTemplateBuilder
from modelseedpy import FBAHelper, KBaseMediaPkg, MSBuilder, MSGenome, RastClient
from modelseedpy.core.msbuilder import build_biomass, core_atp
from modelseedpy.helpers import get_classifier, get_template

# TODO: Check that this genome file is the correct one to be using
genome = MSGenome.from_fasta(
    "genome/Michelle's 4106 gene calls/MIT1002_anvio_prot_seqs.fa", split=" "
)

# Annotate the genome
# FIXME: Do I need this? Can I use the exact same annotations as Michelle?
rast = RastClient()
rast.annotate_genome(genome)

# Export the annotated genome to be able to compare it to Michelle's?

# # Classify the genome
# # Assuming it is classifying the genome as gram positive or negative
# # N would mean gram negative
# genome_classifier = get_classifier("knn_filter")
# genome_class = genome_classifier.classify(genome)

# # Get the template model
# if genome_class == "G":
#     # Not actually sure if G is gram positive, but that is what copilot says
#     template = get_template("template_gram_pos")
# elif genome_class == "N":
#     template = get_template("template_gram_neg")
# else:
#     # Throw an error if the classifier returns something unexpected
#     raise ValueError(f"Unexpected genome class: {genome_class}")

# Getting a weird error within "get_classifier" that I don't understand
# Traceback (most recent call last):
#   File "/Users/helenscott/Documents/PhD/Segre-lab/GEM-repos/mit1002-model/make_model.py", line 26, in <module>
#     genome_classifier = get_classifier("knn_filter")
#   File "/Users/helenscott/opt/miniconda3/envs/gem-reconstruction/lib/python3.8/site-packages/modelseedpy/helpers.py", line 43, in get_classifier
#     model_filter = pickle.load(fh)
# ModuleNotFoundError: No module named 'sklearn.neighbors._dist_metrics'
# So I'm just going to define that the template is the gram negative one, since I know that
template = get_template("template_gram_neg")

# Convert to Cobra?
# Not sure what this is doing
cobra_template = NewModelTemplateBuilder.from_dict(template, None).build()

# Make the model
model_builder = MSBuilder(genome, cobra_template)
base_model = model_builder.build(
    # TODO: Look up what these options are
    "model",
    "0",
    allow_all_non_grp_reactions=True,
    annotate_with_rast=False,
)
# Add the biomass reaction?
# Not sure what this is doing
# FIXME: Check the template is right (in the demo was using core template)
base_model.add_reactions(
    [build_biomass("bio2", base_model, cobra_template, core_atp, "0")]
)

# Save the model
cobra.io.write_sbml_model(base_model, "modelseedpy_model.xml")
