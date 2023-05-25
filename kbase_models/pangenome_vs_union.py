import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import cobra

# Load the models
gem_union = cobra.io.read_sbml_model('kbase_models/MIT1002_union.sbml')
gem_bayesian_core = cobra.io.load_json_model('kbase_models/alteromonas_bayesian_core.json')
 
# Second way
venn2([set([r.id[0:-3] for r in gem_union.reactions]),
       set([r.id for r in gem_bayesian_core.reactions])])

# Style

# Save
print("foo")