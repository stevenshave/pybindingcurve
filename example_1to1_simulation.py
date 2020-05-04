"""Simulation example 1:1 binding"""

import numpy as np
import pybindingcurve as pbc

# Define our system, 1:1 binding has p and l which (usually) relate two protein
# and ligand concentration, although can be any two species which bind.  kdpl
# is the dissociation constant between the two species.

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.
system_parameters = {'p': 1, 'l': 10, 'kdpl': 1}

# Make a pbc BindingCurve defined by the simple 1:1 binding system
my_system = pbc.BindingCurve("1:1kinetic")
print("Simulating 1:1 binding system with these parameters:")
print(system_parameters)
print("pl=", my_system.query(system_parameters))

# Simulate and visualise a binding curve.
# First, we redefine the system parameters so that one variable is changing
# in this case, we choose protein, performing a titration from 0 to 10 uM.
system_parameters = {"p": np.linspace(0, 20), "l": 10, "kdpl": 1}

# We can now add the curve to the plot, name it with an optional name= value.
my_system.add_curve(system_parameters)

# Lets change the KD present in the system parameters into something higher
# affinity (lower KD) and add it to the curve
system_parameters2 = {"p": np.linspace(0, 20), "l": 10, "kdpl": 0.5}
my_system.add_curve(system_parameters2)

# Show the plot
my_system.show_plot()
