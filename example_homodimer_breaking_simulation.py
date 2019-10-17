#!/usr/bin/env python
"""Simulation example homodimer breaking"""

import numpy as np
import pybindingcurve as pbc
import time

# Simulate a binding curve
system_parameters = {'p': np.linspace(0,100, num=7), 'kdpp': 10, 'i':100, 'kdpi':1}
mySystem = pbc.BindingCurve("homodimer breaking")
mySystem.add_curve(system_parameters, 'pp')
print(mySystem.curves[-1].xcoords)
print(mySystem.curves[-1].ycoords)
mySystem.show_plot()
