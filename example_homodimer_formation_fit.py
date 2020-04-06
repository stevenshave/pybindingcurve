"""Fitting example, determining Kd for homodimer breaking data"""

import numpy as np
import pybindingcurve as pbc
import sys

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

# Experimental data
xcoords = np.array([0.0,2,4,6,8,10])
ycoords = np.array([0., 0.22, 0.71, 1.24,1.88,2.48])

# Construct the PyBindingCurve object, operating on a homodimer formation 
# system and add experimental data to the plot
mySystem = pbc.BindingCurve("homodimer formation")
mySystem.add_scatter(xcoords, ycoords)

# Known system parameters, kdpp will be added to this by fitting
system_parameters = {"p": xcoords}

# Now we call fit, passing the known parameters, followed by a dict of parameters to be fitted along
# with an initial guess, pass the ycoords, and what the readout (ycoords) is
fitted_system, fit_accuracy = mySystem.fit(system_parameters, {"kdpp": 0}, ycoords)

# Print out the fitted parameters
for k, v in fit_accuracy.items():
    print(f"Fit: {k}={fitted_system[k]} +/- {v}")

# Assign more points to 'p' to make a smooth plot
fitted_system["p"] = np.linspace(0, np.max(xcoords))

# Add a new curve, simulated using fitted parameters to our BindingCurve object
mySystem.add_curve(fitted_system)

# Show the plot
mySystem.show_plot()
