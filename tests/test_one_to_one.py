"""
pytest tests for PyBindingCurve


PyBindingCurve source may be tested to ensure internal consistency (agreement)
amongst simulation methods, and externally consistent (agreement with)
literature values. With pytest installed in the local python environment
(pip install pytest), simply run 'pytest' to run the testsuite.
"""

import pybindingcurve as pbc
import numpy as np

########################
### Test 1:1 simulation
########################

# Default
def test_1_to_1_simulation_default_approach(data_res_one_to_one):
	my_system = pbc.BindingCurve("1:1")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - data_res_one_to_one)) < 1e-6

# Analytical
def test_1_to_1_simulation_analytical(data_res_one_to_one):
	my_system = pbc.BindingCurve("1:1analytical")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - data_res_one_to_one)) < 1e-6

# Minimizer
def test_1_to_1_simulation_minimizer(data_res_one_to_one):
	my_system = pbc.BindingCurve("1:1min")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - data_res_one_to_one)) < 1e-6

# Lagrange
def test_1_to_1_simulation_lagrange(data_res_one_to_one):
	my_system = pbc.BindingCurve("1:1lagrange")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - data_res_one_to_one)) < 1e-6

# Kinetic
def test_1_to_1_simulation_kinetic(data_res_one_to_one):
	my_system = pbc.BindingCurve("1:1kinetic")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - data_res_one_to_one)) < 1e-6

# Custom system definition
def test_1_to_1_simulation_custom_definition(data_res_one_to_one):
	my_system = pbc.BindingCurve("P+L<->PL*")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kd_p_l_pl": 1}) - data_res_one_to_one)) < 1e-6
