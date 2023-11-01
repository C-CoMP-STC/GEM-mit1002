import cobra
import memote.support.consistency as consistency

# Load the model
model = cobra.io.read_sbml_model("model.xml")

met = "MNXM3"  # ATP

print("Finding ATP generating cycles")
print("----------------------------")
# For a specfied number of times
for i in range(0, 10):
    # Make a copy of the model
    model_copy = model.copy()
    # Run just the function for finding ATP generating cycles
    rxns = consistency.detect_energy_generating_cycles(model_copy, met)
    # Print the results
    print(rxns)
