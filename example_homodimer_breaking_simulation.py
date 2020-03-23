"""Simulation example homodimer breaking"""

import numpy as np
import pybindingcurve as pbc
import time

# Simulate a binding curve
system_parameters = {"p": np.linspace(
    0, 100, num=100), "kdpp": 10, "i": 100, "kdpi": 1}
# mySystem = pbc.BindingCurve("homodimer breaking")
mySystem = pbc.BindingCurve(pbc.systems.System_analytical_homodimerbreaking_pp)
mySystem.add_curve(system_parameters)
mySystem.add_curve(
    {"p": np.linspace(0, 10, num=20), "kdpp": 1, "i": 10, "kdpi": 1}, "Blah curve")
mySystem.show_plot()
