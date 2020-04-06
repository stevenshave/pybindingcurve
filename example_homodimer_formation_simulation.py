"""Simulation example homodimer formation"""

import numpy as np
import pybindingcurve as pbc
import time

# Define our system, homodimer formation has only:
# p: protein, or monomer concentration
# kdpp: the dissociation constant between the two species.
# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.


# Define the system
system_parameters = {"p": np.linspace(0, 10), "kdpp": 10}

# Make a pbc BindingCurve defined by the 'homodimer formation' binding system
mySystem = pbc.BindingCurve("homodimer formation")

# We can now add the curve to the plot, name it with an optional name= value.
mySystem.add_curve(system_parameters)
print(mySystem.curves[0].ycoords)
# Show the plot
mySystem.show_plot()
