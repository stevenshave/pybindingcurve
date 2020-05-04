"""Simulation example 1:1 binding"""

import numpy as np
import pybindingcurve as pbc

#mySystem = pbc.BindingCurve(pbc.systems.System_analytical_homodimerbreaking_pp)
mySystem = pbc.BindingCurve("homodimerbreaking")
system_parameters={'p':14.7, 'kdpp':3.5, 'i':20, 'kdpi':1.2}

print(mySystem.query(system_parameters))
