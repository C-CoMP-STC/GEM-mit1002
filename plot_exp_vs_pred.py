import pandas as pd
import os
import matplotlib.pyplot as plt
import pickle

# Load the results
with open('glc_results.pkl', 'rb') as f:
    glc_experiment = pickle.load(f)

with open('ace_results.pkl', 'rb') as f:
    ace_experiment = pickle.load(f)

# Path to Zac's results
# Assuming you are running from the root of the repository
results_path = '../../CUE/Zac txt data/'

########################################################################
# Experimental and Predicted Biomass
########################################################################
# Load the OD data
od = pd.read_csv(os.path.join(results_path, 'MIT1002_singles_OD600.txt'),
                 sep='\t')

# Get the mean and double standard deviation of the OD for the replicates of
# growth on glucose only
od['glucose_mean'] = od[['MIT1002_glucose',
                         'MIT1002_glucose.1',
                         'MIT1002_glucose.2',
                         'MIT1002_glucose.3',
                         'MIT1002_glucose.4']].mean(axis=1)
od['glucose_double_std'] = od[['MIT1002_glucose',
                               'MIT1002_glucose.1',
                               'MIT1002_glucose.2',
                               'MIT1002_glucose.3',
                               'MIT1002_glucose.4']].std(axis=1) * 2

# Get the mean and double standard deviation of the OD for the replicates of
# growth on acetate only
od['acetate_mean'] = od[['MIT1002_acetate',
                         'MIT1002_acetate.1',
                         'MIT1002_acetate.2',
                         'MIT1002_acetate.3',
                         'MIT1002_acetate.4']].mean(axis=1)
od['acetate_double_std'] = od[['MIT1002_acetate',
                               'MIT1002_acetate.1',
                               'MIT1002_acetate.2',
                               'MIT1002_acetate.3',
                               'MIT1002_acetate.4']].std(axis=1) * 2

fig, ax = plt.subplots(figsize=(20,10)) 

# Plot the OD data on one y axis, with the mean as a scatter plot and the
# double standard deviation as error bars
od.plot(x = 'Time',
        y = 'glucose_mean',
        kind='scatter',
        yerr='glucose_double_std',
        ax = ax,
        label = 'Experimental OD600 (Glucose)')

od.plot(x = 'Time',
        y = 'acetate_mean',
        kind='scatter',
        yerr='acetate_double_std',
        ax = ax,
        label = 'Experimental OD600 (Acetate)',
        color = 'orange')


# Convert the cycles in the total biomass dataframe to hours by dividing
# by 100 and shift it to account for the lag phase
glc_experiment.total_biomass['Time'] = glc_experiment.total_biomass['cycle']/100 + 4

# Plot the FBA restul as Biomass on the secondary y axis
glc_experiment.total_biomass.plot(x = 'Time',
                              y = '', # Because the model doesn't have a name
                              ax = ax,
                              secondary_y = True,
                              label = 'FBA Predicted Biomass (Glucose)')

# Convert the cycles in the total biomass dataframe to hours by dividing
# by 100 and shift it to account for the lag phase
ace_experiment.total_biomass['Time'] = ace_experiment.total_biomass['cycle']/100 + 4

# Plot the FBA restul as Biomass on the secondary y axis
ace_experiment.total_biomass.plot(x = 'Time',
                              y = '', # Because the model doesn't have a name
                              ax = ax,
                              secondary_y = True,
                              label = 'FBA Predicted Biomass (Acetate)')

# Label the axes
ax.set_xlabel('Time (hours)')
ax.set_ylabel('OD600')
ax.right_ax.set_ylabel('Biomass (gr.)')

# Save the plot
plt.tight_layout()
plt.savefig('exp_vs_pred_biomass.png')