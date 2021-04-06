"""Fitting example, determining Kd from 1:1:1 competition data"""

import numpy as np
import pybindingcurve as pbc
import sys

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.

# Experimental data
xcoords = np.array([0.0, 4.2, 8.4, 16.8, 21.1, 31.6, 35.8, 40.0])
ycoords = np.array([150, 330, 1050, 3080, 4300, 6330, 6490, 6960])

# Construct the PyBindingCurve object, operating on a 1:1:1 (compeittion) system and add experimental data to the plot
my_system = pbc.BindingCurve("1:1:1")
my_system.add_scatter(xcoords, ycoords)

# Known system parameters, kdpl will be added to this by fitting
system_parameters = {
    "p": xcoords,
    "l": 10,
    "i": 10,
    "kdpl": 10,
    "ymin": np.min(ycoords),
}

# Now we call fit, passing the known parameters, followed by a dict of parameters to be fitted along
# with an initial guess, pass the ycoords, and what the readout (ycoords) is
fitted_system, fit_accuracy = my_system.fit(
    system_parameters, {"kdpi": 0, "ymax": np.max(ycoords)}, ycoords
)

# Print out the fitted parameters
for k, v in fit_accuracy.items():
    print(f"Fit: {k}={fitted_system[k]} +/- {v}")

# Assign more points to 'p' to make a smooth plot
fitted_system["p"] = np.linspace(0, np.max(xcoords))

# Add a new curve, simulated using fitted parameters to our BindingCurve object
my_system.add_curve(fitted_system)

# Show the plot
my_system.show_plot(ylabel="Signal")
