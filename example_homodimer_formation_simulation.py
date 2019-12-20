#!/usr/bin/env python
"""Simulation example homodimer formation"""
# XXX Not working
import numpy as np
import pybindingcurve as pbc
import time

# Show the difference in speed between a kinetic solution and the highly efficient direct analytical method for homodimer formation
start = time.time()
print(
    f"Homodimer formation kinetic\t\tpl={pbc.System_kinetic_homodimerformation_pp().query({'p':4, 'kdpp':1})}, \ttook {round(time.time()-start,5)} seconds")
start = time.time()
print(
    f"Homodimer formation analytical \tpl={pbc.System_analytical_homodimerformation_pp().query({'p':4, 'kdpp':1})}, \ttook {round(time.time()-start,5)} seconds")

# Simulate a binding curve
system_parameters = {'p': np.linspace(0,10000), 'kdpp':10}
mySystem = pbc.BindingCurve("homodimer formation")
mySystem.add_curve(system_parameters)
mySystem.show_plot()
