from scipy.optimize import fsolve
from autograd import grad
from .binding_system import BindingSystem
from .lagrange_binding_system_factory import LagrangeBindingSystemFactory
import numpy as np


# 1:1 binding - see https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
def system01_lagrange(p, l, kdpl):
    def F(X):  # Augmented Lagrange function
        response = (X[0] * X[1]) / kdpl
        return (
            response - X[2] * (p - (response + X[0])) - X[3] * (l - (response + X[1]))
        )

    dfdL = grad(F)  # Gradients of the Lagrange function
    pf, lf, lam1, lam2 = fsolve(dfdL, [p, l, 1.0, 1.0])
    return {"pf": pf, "lf": lf, "pl": (pf * lf) / kdpl}


# 1:1:1 competition - see https://stevenshave.github.io/pybindingcurve/simulate_competition.html
def system02_lagrange(p, l, i, kdpl, kdpi):
    def F(X):  # Augmented Lagrange function
        response = (X[0] * X[1]) / kdpl
        constraint1 = p - (response + X[0] + (X[0] * X[2] / kdpi))
        constraint2 = l - (response + X[1])
        constraint3 = i - (X[2] + (X[0] * X[2] / kdpi))
        return response - X[3] * constraint1 - X[4] * constraint2 - X[5] * constraint3

    dfdL = grad(F, 0)  # Gradients of the Lagrange function
    pf, lf, inhf, lam1, lam2, lam3 = fsolve(dfdL, [p, l, i, 1.0, 1.0, 1.0])
    return {"pf": pf, "lf": lf, "pl": (pf * lf) / kdpl}


# Homodimer formation - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
def system03_lagrange(p, kdpp):
    def F(X):  # Augmented Lagrange function
        response = (X[0] * X[0]) / kdpp
        constraint1 = p - (response * 2 + X[0])
        return response - X[1] * constraint1

    dfdL = grad(F, 0)  # Gradients of the Lagrange function
    pf, lam1 = fsolve(dfdL, [p, 1.0])
    return {"pf": pf, "pp": ((pf * pf) / kdpp)}


# Homodimer breaking - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
def system04_lagrange(p, i, kdpp, kdpi):
    def F(X):  # Augmented Lagrange function
        pp = (X[0] * X[0]) / kdpp
        pi = (X[0] * X[1]) / kdpi
        constraint1 = p - (pp * 2 + pi + X[0])
        constraint2 = i - (pi + X[1])
        return pp - X[2] * constraint1 - X[3] * constraint2

    dfdL = grad(F, 0)  # Gradients of the Lagrange function
    pf, i_f, lam1, lam2 = fsolve(dfdL, [p, i, 1.0, 1.0])
    return {"pf": pf, "if": i_f, "pp": ((pf * pf) / kdpp), "pi": (pf * i_f) / kdpi}


# 1:2 binding
def system12_lagrange(p, l, kdpl1, kdpl2):
    def F(X):  # Augmented Lagrange function
        pf = X[0]
        lf = X[1]
        pl1 = pf * lf / kdpl1
        pl2 = pf * lf / kdpl2
        pl12 = (pl1 * lf + pl2 * lf) / (kdpl1 + kdpl2)
        constraint1 = p - (pf + pl1 + pl2 + pl12)
        constraint2 = l - (lf + pl1 + pl2 + 2 * pl12)
        return pl12 - X[2] * constraint1 - X[3] * constraint2

    dfdL = grad(F, 0)  # Gradients of the Lagrange function

    pf, lf, lam1, lam2 = fsolve(dfdL, [p, l] + [1.0] * 2)
    pl1 = pf * lf / kdpl1
    pl2 = pf * lf / kdpl2
    pl12 = (pl1 * lf + pl2 * lf) / (kdpl1 + kdpl2)
    return {"pf": pf, "lf": lf, "pl1": pl1, "pl2": pl2, "pl12": pl12}


# 1:3 binding
def system13_lagrange(p, l, kdpl1, kdpl2, kdpl3):
    def F(X):  # Augmented Lagrange function
        pf = X[0]
        lf = X[1]
        pl1 = pf * lf / kdpl1
        pl2 = pf * lf / kdpl2
        pl3 = pf * lf / kdpl3
        pl12 = (pl1 * lf + pl2 * lf) / (kdpl1 + kdpl2)
        pl13 = (pl1 * lf + pl3 * lf) / (kdpl1 + kdpl3)
        pl23 = (pl2 * lf + pl3 * lf) / (kdpl2 + kdpl3)
        pl123 = (pl12 * lf + pl13 * lf + pl23 * lf) / (kdpl1 + kdpl2 + kdpl3)
        constraint1 = p - (pf + pl1 + pl2 + pl3 + pl12 + pl13 + pl23 + pl123)
        constraint2 = l - (
            lf + pl1 + pl2 + pl3 + 2 * (pl12 + pl13 + pl23) + 3 * (pl123)
        )
        nonzero_constraint = (constraint1 - abs(constraint1)) - (
            constraint2 - abs(constraint2)
        )
        return (
            pl123 - X[2] * constraint1 - X[3] * constraint2 - X[4] * nonzero_constraint
        )

    dfdL = grad(F, 0)  # Gradients of the Lagrange function
    pf, lf, lam1, lam2, lam3 = fsolve(dfdL, [p, l] + [1.0] * 3)
    pl1 = pf * lf / kdpl1
    pl2 = pf * lf / kdpl2
    pl3 = pf * lf / kdpl3
    pl12 = (pl1 * lf + pl2 * lf) / (kdpl1 + kdpl2)
    pl13 = (pl1 * lf + pl3 * lf) / (kdpl1 + kdpl3)
    pl23 = (pl2 * lf + pl3 * lf) / (kdpl2 + kdpl3)
    pl123 = (pl12 * lf + pl13 * lf + pl23 * lf) / (kdpl1 + kdpl2 + kdpl3)
    return {
        "pf": pf,
        "lf": lf,
        "pl1": pl1,
        "pl2": pl2,
        "pl3": pl3,
        "pl12": pl12,
        "pl13": pl13,
        "pl23": pl23,
        "pl123": pl123,
    }


class System_lagrange_one_to_one__pl(BindingSystem):
    """
    Lagrange 1:1 binding system

    Class defines 1:1 binding, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
    """

    def __init__(self):
        super().__init__(system01_lagrange)
        self.default_readout = "pl"

    def query(self, parameters: dict):
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.nan_to_num(
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_lagrange_competition_pl(BindingSystem):
    """
    Lagrange 1:1:1 competition binding system

    Class defines 1:1:1 competition, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_competition.html
    """

    def __init__(self):
        super().__init__(system02_lagrange)
        self.default_readout = "pl"

    def query(self, parameters: dict):
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.nan_to_num(
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_lagrange_homodimerformation__pp(BindingSystem):
    """
    Lagrange homodimer formation system

    Class defines homodimer formation, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
    """

    def __init__(self):
        super().__init__(system03_lagrange, False)
        self.default_readout = "pp"

    def query(self, parameters: dict):
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


class System_lagrange_homodimerbreaking__pp(BindingSystem):
    """
    Lagrange homodimer breaking system

    Class defines homodimer breaking, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system04_lagrange, False)
        self.default_readout = "pp"

    def query(self, parameters: dict):
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


class System_lagrange_homodimerbreaking__pl(BindingSystem):
    """
    Lagrange homodimer breaking system

    Class defines homodimer breaking, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system04_lagrange, False)
        self.default_readout = "pl"

    def query(self, parameters: dict):
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


class System_lagrange_1_to_2__pl12(BindingSystem):
    """
    Lagrange 1:2 binding system.

    Class defines 1:2 protein:ligand binding, readout is PL12, meaning protein
    with 2 ligands
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system12_lagrange, False)
        self.default_readout = "pl12"

    def query(self, parameters: dict):
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.nan_to_num(
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_lagrange_1_to_3__pl123(BindingSystem):
    """
    Lagrange 1:2 binding system.

    Class defines 1:2 protein:ligand binding, readout is PL123, meaning protein
    with 3 ligands
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system13_lagrange, False)
        self.default_readout = "pl123"

    def query(self, parameters: dict):
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.nan_to_num(
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_lagrange_custom(BindingSystem):
    """
    Lagrange custom binding system

    Class uses LagrangeBindingSystemFactory to make a custom lagrange function.
    """

    def __init__(self, system_string):
        custom_system = LagrangeBindingSystemFactory(system_string)
        super().__init__(custom_system.binding_function)
        self.default_readout = custom_system.default_readout

    def query(self, parameters: dict):
        if self._are_ymin_ymax_present(parameters):
            parameters_no_min_max = self._remove_ymin_ymax_keys_from_dict_return_new(
                parameters
            )
            value = super().query(parameters_no_min_max)
            with np.errstate(divide="ignore", invalid="ignore"):
                return np.nan_to_num(
                    parameters["ymin"]
                    + ((parameters["ymax"] - parameters["ymin"]) * value)
                    / parameters["l"]
                )
        else:
            return super().query(parameters)
