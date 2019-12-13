#!/usr/bin/env python
"""Fitting example, determining Kd from 1:1:1 competition data"""

import numpy as np
import pybindingcurve as pbc
import sys

# Experimental data
xcoords = np.array([0.,   4.2,  8.4, 16.8, 21.1, 31.6, 35.8, 40.])
ycoords = np.array([0.15, 0.33, 1.05, 3.08, 4.30, 6.33, 6.49, 6.96])

# Construct the PyBindingCurve object, operating on a 1:1:1 (compeittion) system and add experimental data to the plot
mySystem = pbc.BindingCurve("1:1:1")
mySystem.add_points_to_plot(xcoords, ycoords)

# Known system parameters, kdpl will be added to this by fitting
system_parameters = {'p': xcoords, 'l': 10, 'i': 10, 'kdpl': 10}

# Now we call fit, passing the known parameters, followed by a dict of parameters to be fitted along
# with an initial guess, pass the ycoords, and what the readout (ycoords) is
fitted_system, fit_accuracy = mySystem.fit(
    system_parameters, {'kdpi': 0}, ycoords)

# Print out the fitted parameters
for k, v in fit_accuracy.items():
    print(f"Fit: {k}={fitted_system[k]} +/- {v}")

# Assign more points to 'p' to make a smooth plot
fitted_system['p'] = np.linspace(0, np.max(xcoords))

# Add a new curve, simulated using fitted parameters to our BindingCurve object
mySystem.add_curve(fitted_system)

# Show the plot
mySystem.show_plot()
