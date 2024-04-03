"""
pytest tests for PyBindingCurve

PyBindingCurve source may be tested to ensure internal consistency (agreement)
amongst simulation methods, and externally consistent (agreement with)
literature values. With pytest installed in the local python environment
(pip install pytest), simply run 'pytest' to run the testsuite.
"""
import pytest
import pybindingcurve as pbc
import numpy as np

# Testing of fits.
# lmfit returns proper 95 % confidence intervals for where the real KD value is.
# Taking into account this +/- amount is difficult.

def test_1_to_1_fit_default():
	xcoords = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0])
	ycoords = np.array([0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225])
	my_system = pbc.BindingCurve("1:1")
	system_parameters = {"p": xcoords, "l": 10}
	fitted_system, fit_accuracy = my_system.fit(system_parameters, {"kdpl": 0}, ycoords)
	assert pytest.approx(fitted_system["kdpl"])==16.3715989783099

def test_1_to_1_fit_analytical():
	xcoords = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0])
	ycoords = np.array([0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225])
	my_system = pbc.BindingCurve("1:1analytical")
	system_parameters = {"p": xcoords, "l": 10}
	fitted_system, fit_accuracy = my_system.fit(system_parameters, {"kdpl": 0}, ycoords)
	assert pytest.approx(fitted_system["kdpl"])==16.3715989783099
 