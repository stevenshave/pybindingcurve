#!/usr/bin/env python
"""Fitting example, determining Kd for homodimer breaking data"""

import numpy as np
import pybindingcurve as pbc
import sys

# Experimental data
xcoords=np.array([0.,16.7,33.3,50.,66.7,83.3, 100.])
ycoords=np.array([0.,0.004,0.021,0.094,0.312,1.188, 3.854])

# Construct the PyBindingCurve object, operating on a 1:1:1 (compeittion) system and add experimental data to the plot
mySystem = pbc.BindingCurve("homodimerbreaking")
mySystem.add_points_to_plot(xcoords, ycoords)

# Known system parameters, kdpl will be added to this by fitting
system_parameters = {'p': xcoords, 'i': 100, 'kdpp': 10}

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
