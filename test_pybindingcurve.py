import pybindingcurve as pbc
import numpy as np

def test_1_to_1_simulation():
    system_parameters = {'p': 1, 'l': 10, 'kdpl': 1}
    # Make a pbc BindingCurve defined by the simple 1:1 binding system
    my_system = pbc.BindingCurve("1:1")
    pl=float(my_system.query(system_parameters))
    print(pl)
    assert np.abs(pl-0.9009804864072152)<1e-6

def test_1_to_1_fit():
    xcoords = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0])
    ycoords = np.array([0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225])
    my_system = pbc.BindingCurve("1:1")
    system_parameters = {"p": xcoords, "l": 10}
    fitted_system, fit_accuracy = my_system.fit(system_parameters, {"kdpl": 0}, ycoords)
    assert np.abs(fitted_system["kdpl"]-16.371571968273663)<1e-6

def test_custom_definition_of_1_to_1():
    system_parameters = {'p0': 1, 'l0': 10, 'kd_p_l': 1}
    # Make a pbc BindingCurve defined by the simple 1:1 binding system
    my_system = pbc.BindingCurve("P+L<->PL*")
    pl=float(my_system.query(system_parameters))
    assert np.abs(pl-0.9009804864072152)<1e-6