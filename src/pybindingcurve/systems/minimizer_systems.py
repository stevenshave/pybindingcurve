import numpy as np
from .binding_system import BindingSystem
from mpmath import mpf, findroot, mp, almosteq
from .minimizer_binding_system_factory import MinimizerBindingSystemFactory

mpf_zero=mpf(0)
mpf_tol=mpf("1e-10")
max_iters=1e6

# 1:1 binding - see https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
def system01_minimizer(p, l, kdpl):
    p = mpf(p)
    l = mpf(l)
    kdpl = mpf(kdpl)
    if almosteq(kdpl, mpf_zero, mpf_tol):
        kdpl+=mpf_tol
    def f(p_f, l_f):
        pl = p_f * l_f / kdpl
        return p - (p_f + pl), l - (l_f + pl)
    p_f, l_f = findroot(f, [mpf_zero, mpf_zero], tol=mpf_tol, maxsteps=1e6)
    return {"pf": p_f, "lf": l_f, "pl": (p_f * l_f) / kdpl}


# 1:1:1 competition - see https://stevenshave.github.io/pybindingcurve/simulate_competition.html
def system02_minimizer(p, l, i, kdpl, kdpi):
    kdpl = mpf(kdpl)
    kdpi = mpf(kdpi)
    if almosteq(kdpl, mpf_zero, mpf_tol):
        kdpl+=mpf_tol
    if almosteq(kdpi, mpf_zero, mpf_tol):
        kdpi+=mpf_tol
    p = mpf(p)
    l = mpf(l)
    i = mpf(i)
    def f(p_f, l_f, i_f):
        pl = p_f * l_f / kdpl
        pi = p_f * i_f / kdpi
        return p - (p_f + pl + pi), l - (l_f + pl), i - (i_f + pi)
    p_f, l_f, i_f = findroot(f, [mpf_zero, mpf_zero, mpf_zero], tol=mpf_tol, maxsteps=1e6)
    return {
        "pf": p_f,
        "lf": l_f,
        "if": i_f,
        "pl": (p_f * l_f) / kdpl,
        "pi": (p_f * i_f) / kdpi,
    }


# Homodimer formation - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
def system03_minimizer(p, kdpp):
    p = mpf(p)
    kdpp = mpf(kdpp)
    if almosteq(kdpp, mpf_zero, mpf_tol):
        kdpp+=mpf_tol
    def f(p_f):
        pp = p_f * p_f / kdpp
        return p - (p_f + 2 * pp)
    p_f = findroot(f, [mpf_zero], tol=mpf_tol, maxsteps=1e6)
    return {"pf": p_f, "pp": (p_f * p_f) / kdpp}


# Homodimer breaking - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
def system04_minimizer(p, i, kdpp, kdpi):
    p = mpf(p)
    i = mpf(i)
    kdpp = mpf(kdpp)
    kdpi = mpf(kdpi)
    if almosteq(kdpp, mpf_zero, mpf_tol):
        kdpp+=mpf_tol
    if almosteq(kdpi, mpf_zero, mpf_tol):
        kdpi+=mpf_tol
    def f(p_f, i_f):
        pp = p_f * p_f / kdpp
        pi = p_f * i_f / kdpi
        return p - (p_f + pi + 2 * pp), i - (i_f + pi)
    p_f, i_f = findroot(f, (mpf_zero,mpf_zero), tol=mpf_tol, maxsteps=1e6)
    return {"pf": p_f, "if": i_f, "pp": (p_f * p_f) / kdpp, "pi": (p_f * i_f) / kdpi}


# 1:2 binding
def system12_minimizer(p, l, kdpl1, kdpl2):
    p = mpf(p)
    l = mpf(l)
    kdpl1 = mpf(kdpl1)
    kdpl2 = mpf(kdpl2)
    if almosteq(kdpl1, mpf_zero, mpf_tol):
        kdpl1+=mpf_tol
    if almosteq(kdpl2, mpf_zero, mpf_tol):
        kdpl2+=mpf_tol
    def f(p_f, l_f):
        pl1 = p_f * l_f / kdpl1
        pl2 = p_f * l_f / kdpl2
        pl12 = (pl1 * l_f + pl2 * l_f) / (kdpl1 + kdpl2)
        return p - (p_f + pl1 + pl2 + pl12), l - (l_f + pl1 + pl2 + 2 * pl12)
    p_f, l_f = findroot(f, [mpf_zero, mpf_zero], tol=mpf_tol, maxsteps=1e6)
    pl1 = (p_f * l_f) / kdpl1
    pl2 = (p_f * l_f) / kdpl2
    return {
        "pf": p_f,
        "lf": l_f,
        "pl1": pl1,
        "pl2": pl2,
        "pl12": (pl1 * l_f + pl2 * l_f) / (kdpl1 + kdpl2),
    }


# 1:3 binding
def system13_minimizer(p, l, kdpl1, kdpl2, kdpl3):
    p = mpf(p)
    l = mpf(l)
    kdpl1 = mpf(kdpl1)
    kdpl2 = mpf(kdpl2)
    kdpl3 = mpf(kdpl3)
    if almosteq(kdpl1, mpf_zero, mpf_tol):
        kdpl1+=mpf_tol
    if almosteq(kdpl2, mpf_zero, mpf_tol):
        kdpl2+=mpf_tol
    if almosteq(kdpl3, mpf_zero, mpf_tol):
        kdpl3+=mpf_tol
    def f(p_f, l_f):
        pl1 = p_f * l_f / kdpl1
        pl2 = p_f * l_f / kdpl2
        pl3 = p_f * l_f / kdpl3
        pl12 = (pl1 * l_f + pl2 * l_f) / (kdpl1 + kdpl2)
        pl23 = (pl2 * l_f + pl3 * l_f) / (kdpl2 + kdpl3)
        pl13 = (pl1 * l_f + pl3 * l_f) / (kdpl1 + kdpl3)
        pl123 = (pl12 * l_f + pl23 * l_f + pl13 * l_f) / (kdpl1 + kdpl2 + kdpl3)
        return p - (p_f + pl1 + pl2 + pl3 + pl12 + pl13 + pl23), l - (
            l_f + pl1 + pl2 + pl3 + 2 * (pl12 + pl13 + pl23) + 3 * pl123
        )
    p_f, l_f = findroot(f, [mpf_zero, mpf_zero], tol=mpf_tol, maxsteps=1e6)
    pl1 = p_f * l_f / kdpl1
    pl2 = p_f * l_f / kdpl2
    pl3 = p_f * l_f / kdpl3
    pl12 = (pl1 * l_f + pl2 * l_f) / (kdpl1 + kdpl2)
    pl23 = (pl2 * l_f + pl3 * l_f) / (kdpl2 + kdpl3)
    pl13 = (pl1 * l_f + pl3 * l_f) / (kdpl1 + kdpl3)
    return {
        "pf": p_f,
        "lf": l_f,
        "pl1": pl1,
        "pl2": pl2,
        "pl3": pl3,
        "pl12": pl12,
        "pl13": pl13,
        "pl23": pl23,
        "pl123": (pl12 * l_f + pl23 * l_f + pl13 * l_f) / (kdpl1 + kdpl2 + kdpl3),
    }

class System_minimizer_one_to_one__pl(BindingSystem):
    """
    Minimizer-based one to one binding

    Class defines 1:1 binding, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
    """

    def __init__(self):
        super().__init__(system01_minimizer, False)
        self.default_readout = "pl"

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


class System_minimizer_homodimerformation__pp(BindingSystem):
    """
    Minimizer-based homodimer formation system

    Class defines homodimer formation, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
    """

    def __init__(self):
        super().__init__(system03_minimizer, False)
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


class System_minimizer_competition__pl(BindingSystem):
    """
    Minimizer-based competition system

    Class defines 1:1:1 competition, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_competition.html
    """

    def __init__(self):
        super().__init__(system02_minimizer, False)
        self.default_readout = "pl"

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


class System_minimizer_homodimerbreaking__pp(BindingSystem):
    """
    Minimizer-based homodimer breaking system

    Class defines homodimer breaking, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system04_minimizer, False)
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


class System_minimizer_homodimerbreaking__pl(BindingSystem):
    """
    Minimzer-based homodimer breaking system

    Class defines homodimer breaking, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system04_minimizer, False)
        self.default_readout = "pl"

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


class System_minimizer_1_to_2__pl12(BindingSystem):
    """
    Lagrange 1:2 binding system.

    Class defines 1:2 protein:ligand binding, readout is PL12, meaning protein
    with 2 ligands
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system12_minimizer, False)
        self.default_readout = "pl12"

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
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_minimizer_1_to_3__pl123(BindingSystem):
    """
    Lagrange 1:2 binding system.

    Class defines 1:2 protein:ligand binding, readout is PL123, meaning protein
    with 3 ligands
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system13_minimizer, False)
        self.default_readout = "pl123"

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
                    / parameters["l"]
                )
        else:
            return super().query(parameters)


class System_minimizer_custom(BindingSystem):
    """
    Lagrange custom binding system

    Class uses LagrangeBindingSystemFactory to make a custom lagrange function.
    """

    def __init__(self, system_string):
        custom_system = MinimizerBindingSystemFactory(system_string)
        super().__init__(custom_system.binding_function)
        self.all_species=custom_system.all_species
        self.default_readout = custom_system.readout
        self.arguments=custom_system.binding_function_arguments
        
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
                    / parameters["l"]
                )
        else:
            return super().query(parameters)
