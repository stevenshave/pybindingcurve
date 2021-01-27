"""Simulation example 1:1:1 comptition binding"""

import numpy as np
import pybindingcurve as pbc

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

# Create the PBC BindingCurve object, expecting a 'competition' system.
my_system = pbc.BindingCurve("competition")

# First, lets simulate a curve with no inhibitor present (essentially 1:1)
my_system.add_curve(
    {"p": np.linspace(0, 100), "l": 10, "i": 0, "kdpi": 1, "kdpl": 10}, "No inhibitor"
)
# Add curve with inhibitor (i)
my_system.add_curve(
    {"p": np.linspace(0, 100), "l": 10, "i": 10, "kdpi": 1, "kdpl": 10}, "[i] = 10 µM"
)
# Add curve with more inhibtor (i)
my_system.add_curve(
    {"p": np.linspace(0, 100), "l": 10, "i": 25, "kdpi": 1, "kdpl": 10}, "[i] = 25 uM"
)

# Display the plot
my_system.show_plot()
