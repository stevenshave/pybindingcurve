"""Simulation example homodimer breaking"""

import numpy as np
import pybindingcurve as pbc
import time

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

# Define our homodimer breaking system, titrating in inhibitor
system_parameters = {"p": 30, "kdpp": 10, "i": np.linspace(0,60), "kdpi": 1}

# Create the PBC BindingCurve object, expecting a 'homodimer breaking' system.
my_system = pbc.BindingCurve("homodimer breaking")

# Add the system to PBC, generating a plot.
my_system.add_curve(system_parameters)

# Display the plot
my_system.show_plot()
