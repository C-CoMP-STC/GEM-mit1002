import comets as c
import math
import matplotlib.pyplot as plt

# Load Prochlorococcus Genome-scale model
model = c.model('iSO595v6.xml');
model.initial_pop = [1, 1, 1e-7]
model.obj_style = 'MAX_OBJECTIVE_MIN_TOTAL'

# The ratio of chlorophyll is extracted from the model biomass-function
ci_dvchla = 0.016                 # gr/gDW (Partensky 1993 / Casey 2016)
ci_dvchlb = 0.0013                # gr/gDW (Partensky 1993 / Casey 2016)
absorption_dvchla_680 = 0.0184    # m^2 mg^-1 (Bricaud et al., 2004)
absorption_dvchlb_680 = 0.0018    # m^2 mg^-1 (Bricaud et al., 2004)
absorption_water_680 = 0.465      # m^-1 (Pope and Fry, 1997)
wavelength = 680                  # nm

diameter = 0.6                    # um (Morel et al., 1993)
n_dash = 13.77*1e-3               # imaginary part of refractive index at 675 nm (Stramski et al. 2001)
size_parameter_alpha = diameter*1e3*math.pi/wavelength    # The ratio between the cell size and wavelength
rho_dash = 4*size_parameter_alpha*n_dash
Q_a = 1+(2*math.exp(-rho_dash)/rho_dash)+2*(math.exp(-rho_dash)-1)/rho_dash**2
packaging_effect = 1.5*Q_a/rho_dash

# Calculate the Prochlorococcus specific biomass absorption coefficient in units m2/ g DW biomass
absorption_biomass = packaging_effect*(ci_dvchla*1e3*absorption_dvchla_680+ci_dvchlb*1e3*absorption_dvchlb_680)

model.add_light('LightEX', absorption_biomass, absorption_water_680)

# Make layout with the COMETS toolbox
layout = c.layout(model);

# Define medium
metabs = ['Ammonia[e]', 'HCO3[e]', 'CO2[e]', 'H[e]', 'Orthophosphate[e]', 'H2O[e]', 'Cadmium[e]', 
          'Calcium_cation[e]', 'Chloride_ion[e]', 'Cobalt_ion[e]', 'Copper[e]', 'Fe2[e]', 
          'Magnesium_cation[e]', 'Molybdenum[e]', 'K[e]','Selenate[e]', 'Sodium_cation[e]', 
          'Strontium_cation[e]', 'Sulfate[e]', 'Zn2[e]', 'Hydrogen_sulfide[e]']

for i in metabs:
    layout.set_specific_metabolite(i, 1000)
    layout.set_specific_static(i, 1000)

# Set light conditions by defining parameters
layout.set_global_periodic_media(metabolite='Photon[e]', function='half_sin', amplitude=0.04,
                                 period=24, phase=0, offset=0)

# Set simulation parameters
sim_params = c.params()

sim_params.all_params['maxCycles'] = 480
sim_params.all_params['timeStep'] = 0.1
sim_params.all_params['defaultDiffConst'] = 0

# Runs comets and produce the output files mediaLog.m and biomassLog.m
simulation = c.comets(layout, sim_params)
simulation.JAVA_CLASSPATH = '/home/djordje/Dropbox/COMETS_RUN/lib/jmatio.jar:/home/djordje/Dropbox/COMETS_RUN/lib/jdistlib-0.4.5-bin.jar:/home/djordje/Dropbox/COMETS_RUN/lib/commons-math3-3.6.1.jar:/home/djordje/Dropbox/COMETS_RUN/lib/commons-lang3-3.9.jar:/home/djordje/Dropbox/COMETS_RUN/lib/colt.jar:/home/djordje/Dropbox/COMETS_RUN/lib/concurrent.jar:/home/djordje/Dropbox/COMETS_RUN/bin/comets_2.8.2.jar:/opt/gurobi901/linux64/lib/gurobi.jar'
simulation.run(delete_files=False)

print(simulation.run_output)

# Plot results