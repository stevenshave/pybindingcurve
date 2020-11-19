"""Simulation example 1:1 binding"""

import numpy as np
import pybindingcurve as pbc

system_parameters = {'p': 1, 'l': 10, 'kdpl': 1}

my_system = pbc.BindingCurve("simple")
system_parameters = {"p": np.linspace(0, 20), "l": 10, "kdpl": 1}
my_system.add_curve(system_parameters)
system_parameters2 = {"p": np.linspace(0, 20), "l": 10, "kdpl": 0.5}
my_system.add_curve(system_parameters2)
my_system.show_plot()

"""Fitting example, determining Kd from 1:1 binding data"""

import numpy as np
import pybindingcurve as pbc

xcoords = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0])
ycoords = np.array([0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225])
#my_system = pbc.BindingCurve("1:1")
my_system.add_scatter(xcoords, ycoords)
system_parameters = {"p": xcoords, "l": 10}
fitted_system, fit_accuracy = my_system.fit(system_parameters, {"kdpl": 0}, ycoords)
for k, v in fit_accuracy.items():
    print(f"Fit: {k}={fitted_system[k]} +/- {v}")
fitted_system["p"] = np.linspace(0, np.max(xcoords), num=200)
my_system.add_curve(fitted_system)
my_system.show_plot(ylabel="Signal")