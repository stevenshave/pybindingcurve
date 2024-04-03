"""
pytest tests for PyBindingCurve


PyBindingCurve source may be tested to ensure internal consistency (agreement)
amongst simulation methods, and externally consistent (agreement with)
literature values. With pytest installed in the local python environment
(pip install pytest), simply run 'pytest' to run the testsuite.
"""

import pybindingcurve as pbc
import numpy as np

##########################
### Test 1:1:1 competition
##########################

# Default
def test_competition_simulation_default(data_res_competition):
	my_system = pbc.BindingCurve("competition")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - data_res_competition)) < 1e-6

# Default
def test_competition_simulation_analytical(data_res_competition):
	my_system = pbc.BindingCurve("competitionanalytical")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - data_res_competition)) < 1e-6

# Minimizer
def test_competition_simulation_minimzer(data_res_competition):
	my_system = pbc.BindingCurve("competitionmin")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - data_res_competition)) < 1e-6

# Kinetic
def test_competition_simulation_kinetic(data_res_competition):
	my_system = pbc.BindingCurve("competitionkinetic")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - data_res_competition)) < 1e-6

# Custom system definition
def test_competition_simulation_custom_definition(data_res_competition):
	my_system = pbc.BindingCurve("p+l<->pl*,p+i<->pi")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kd_p_i_pi": 1, "kd_p_l_pl": 10}) - data_res_competition)) < 1e-6
