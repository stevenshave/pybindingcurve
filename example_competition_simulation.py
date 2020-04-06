"""Simulation example 1:1:1 comptition binding"""

import numpy as np
import pybindingcurve as pbc

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

# Create the PBC BindingCurve object, expecting a 'competition' system.
mySystem = pbc.BindingCurve("competition")

# First, lets simulate a curve with no inhibitor present (essentially 1:1)
mySystem.add_curve(
    {"p": np.linspace(0, 40, 20), "l": 10, "i": 0, "kdpi": 1, "kdpl": 10}, "pl"
)
# Add curve with inhibitor (i)
mySystem.add_curve(
    {"p": np.linspace(0, 40, 20), "l": 10, "i": 10, "kdpi": 1, "kdpl": 10}, "pl"
)
# Add curve with more inhibtor (i)
mySystem.add_curve(
    {"p": np.linspace(0, 40, 20), "l": 10, "i": 25, "kdpi": 1, "kdpl": 10}, "pl"
)

# Display the plot
mySystem.show_plot()
