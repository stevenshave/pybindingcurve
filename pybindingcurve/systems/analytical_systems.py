import numpy as np
from .analytical_equations import *
from mpmath import mp
from .binding_system import BindingSystem


class System_analytical_one_to_one__pl(BindingSystem):
    """
    Analytical 1:1 binding system

    Class defines 1:1 binding, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
    """

    def __init__(self):
        super().__init__(system01_analytical_one_to_one__pl, analytical=True)
        self.default_readout = "pl"

    def query(self, parameters: dict):
        mp.dps = 100
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return (
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_analytical_competition__pl(BindingSystem):
    """
    Analytical 1:1:1 competition binding system

    Class defines 1:1:1 competition, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_competition.html
    """

    def __init__(self):
        super().__init__(system02_analytical_competition__pl, analytical=True)
        self.default_readout = "pl"

    def query(self, parameters: dict):
        mp.dps = 100
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return (
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_analytical_homodimerformation__pp(BindingSystem):
    """
    Analytical homodimer formation system

    Class defines homodimer formation, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
    """

    def __init__(self):
        super().__init__(system03_analytical_homodimer_formation__pp, analytical=True)
        self.default_readout = "pp"

    def query(self, parameters: dict):
        mp.dps = 100
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.nan_to_num(
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / (parameters["p"] / 2.0)
                )
        else:
            return super().query(parameters)


class System_analytical_homodimerbreaking_pp(BindingSystem):
    """
    Analytical homodimer breaking system

    Class defines homodimer breaking, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(
            system04_analytical_homodimer_breaking__pp,
            analytical=True,
        )
        self.default_readout = "pp"
        self.num_solutions = 2

    def query(self, parameters: dict):
        mp.dps = 100
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.nan_to_num(
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / (parameters["p"] / 2.0)
                )
        else:
            return super().query(parameters)
