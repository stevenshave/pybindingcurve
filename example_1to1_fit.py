"""Fitting example, determining Kd from 1:1 binding data"""

import numpy as np
import pybindingcurve as pbc


# Experimental data
xcoords = np.array(
    [0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0]
)
ycoords = np.array(
    [0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225]
)

# Construct the PyBindingCurve object, operating on a simple 1:1 system and add experimental data to the plot
mySystem = pbc.BindingCurve("1:1")
mySystem.add_scatter(xcoords, ycoords)


# Known system parameters, kdpl will be added to this by fitting
system_parameters = {"p": xcoords, "l": 10}

# Now we call fit, passing the known parameters, followed by a dict of parameters to be fitted along
# with an initial guess, pass the ycoords, and what the readout (ycoords) is
fitted_system, fit_accuracy = mySystem.fit(system_parameters, {"kdpl": 0}, ycoords)

# Print out the fitted parameters
for k, v in fit_accuracy.items():
    print(f"Fit: {k}={fitted_system[k]} +/- {v}")

# Assign more points to 'p' to make a smooth plot
fitted_system["p"] = np.linspace(0, np.max(xcoords), num=200)

# Add a new curve, simulated using fitted parameters to our BindingCurve object
mySystem.add_curve(fitted_system)

# Show the plot
mySystem.show_plot()
