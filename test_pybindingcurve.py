"""
pytest tests for PyBindingCurve


PyBindingCurve source may be tested to ensure internal consistency (agreement)
amongst simulation methods, and externally consistent (agreement with)
literature values. With pytest installed in the local python environment
(pip install pytest), simply run 'pytest' to run the testsuite.
"""

import pybindingcurve as pbc
import numpy as np


res_one_to_one = np.array(
	[
		0.0, 0.36976681,0.73678745,1.10079404,1.46148871,1.81854052,2.17158225,
		2.52020739,2.86396726,3.20236871,3.5348727,3.8608942,4.17980399,
		4.49093304,4.79358009,5.08702315,5.37053529,5.64340495,5.90496012,
		6.15459543,6.39180018,6.61618461,6.82750179,7.0256622,7.21073915,
		7.38296437,7.5427142,7.69048839,7.82688409,7.95256788,8.06824855,
		8.17465245,8.27250287,8.36250371,8.44532763,8.52160808,8.59193462,
		8.65685092,8.71685466,8.77239888,8.8238942,8.87171161,8.91618555,
		8.95761706,8.99627687,9.03240837,9.06623039,9.09793973,9.12771349,
		9.15571123,
	]
)

res_competition=np.array([
	4.33809621, 4.22751888, 4.11497749, 4.00061453, 3.88460108, 3.76713883,
	3.64846169, 3.52883679, 3.4085647,  3.28797858, 3.16744212, 3.04734594,
	2.92810239, 2.81013866, 2.69388832, 2.57978151, 2.46823425, 2.35963741,
	2.25434621, 2.15267084, 2.05486895, 1.96114064, 1.87162611, 1.78640606,
	1.70550464, 1.62889436, 1.55650261, 1.48821894, 1.42390269, 1.36339043,
	1.30650277, 1.25305051, 1.20283978, 1.1556763,  1.11136872, 1.06973114,
	1.03058497, 0.99376017, 0.95909597, 0.9264413,  0.8956548,  0.86660478,
	0.83916885, 0.81323356, 0.78869392, 0.76545293, 0.74342101, 0.72251552,
	0.70266024, 0.68378487
])


########################
### Test 1:1 simulation
########################

# Default
def test_1_to_1_simulation_default_approach():
	my_system = pbc.BindingCurve("1:1")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - res_one_to_one)) < 1e-6

# Analytical
def test_1_to_1_simulation_analytical():
	my_system = pbc.BindingCurve("1:1analytical")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - res_one_to_one)) < 1e-6

# Minimizer
def test_1_to_1_simulation_minimizer():
	my_system = pbc.BindingCurve("1:1min")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - res_one_to_one)) < 1e-6

# Lagrange
def test_1_to_1_simulation_lagrange():
	my_system = pbc.BindingCurve("1:1lagrange")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - res_one_to_one)) < 1e-6

# Kinetic
def test_1_to_1_simulation_kinetic():
	my_system = pbc.BindingCurve("1:1kinetic")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kdpl": 1}) - res_one_to_one)) < 1e-6

# Custom system definition
def test_1_to_1_simulation_custom_definition():
	my_system = pbc.BindingCurve("P+L<->PL*")
	assert np.sum(np.abs(my_system.query({"p": np.linspace(0, 20), "l": 10, "kd_p_l_pl": 1}) - res_one_to_one)) < 1e-6


##########################
### Test 1:1:1 competition
##########################

# Default
def test_competition_simulation_default():
	my_system = pbc.BindingCurve("competition")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - res_competition)) < 1e-6

# Default
def test_competition_simulation_analytical():
	my_system = pbc.BindingCurve("competitionanalytical")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - res_competition)) < 1e-6

# Minimizer
def test_competition_simulation_minimzer():
	my_system = pbc.BindingCurve("competitionmin")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - res_competition)) < 1e-6

# Kinetic
def test_competition_simulation_kinetic():
	my_system = pbc.BindingCurve("competitionkinetic")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kdpi": 1, "kdpl": 10}) - res_competition)) < 1e-6

# Custom system definition
def test_competition_simulation_custom_definition():
	my_system = pbc.BindingCurve("p+l<->pl*,p+i<->pi")
	assert np.sum(np.abs(my_system.query({"p": 12, "l": 10, "i": np.linspace(0,25), "kd_p_i_pi": 1, "kd_p_l_pl": 10}) - res_competition)) < 1e-6


# # Testing of fits.
# # lmfit returns proper 95 % confidence intervals for where the real KD value is.
# # Taking into account this +/- amount is difficult.
def test_1_to_1_fit_default():
	xcoords = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0])
	ycoords = np.array([0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225])
	my_system = pbc.BindingCurve("1:1")
	system_parameters = {"p": xcoords, "l": 10}
	fitted_system, fit_accuracy = my_system.fit(system_parameters, {"kdpl": 0}, ycoords)
	assert np.abs(fitted_system["kdpl"] - 16.371571968273663) < 1e-6


def test_1_to_1_fit_analytical():
	xcoords = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0])
	ycoords = np.array([0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225])
	my_system = pbc.BindingCurve("1:1analytical")
	system_parameters = {"p": xcoords, "l": 10}
	fitted_system, fit_accuracy = my_system.fit(system_parameters, {"kdpl": 0}, ycoords)
	assert np.abs(fitted_system["kdpl"] - 16.371571968273663) < 1e-6
