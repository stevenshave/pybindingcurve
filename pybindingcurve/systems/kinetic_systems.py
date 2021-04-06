import numpy as np
from .binding_system import BindingSystem
from scipy.integrate import solve_ivp

# 1:1 binding - see https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
def system01_kinetic(p, l, kdpl, interval=(0, 100)):
    def ode(concs, t, kdpl):
        p, l, pl = concs
        r1 = -p * l + kdpl * pl
        dpdt = r1
        dldt = r1
        dpldt = -r1
        return [dpdt, dldt, dpldt]

    ode_result = solve_ivp(
        lambda t, y: ode(y, t, kdpl), interval, [p, l, 0.0], rtol=1e-12, atol=1e-12
    ).y[:, -1]
    return {"p": ode_result[0], "l": ode_result[1], "pl": ode_result[2]}


# 1:1:1 competition - see https://stevenshave.github.io/pybindingcurve/simulate_competition.html
def system02_kinetic(p, l, i, kdpl, kdpi, interval=(0, 100)):
    def ode(concs, t, kdpl, kdpi):
        p, l, i, pl, pi = concs
        r1 = -p * l + kdpl * pl
        r2 = -p * i + kdpi * pi
        dpdt = r1 + r2
        dldt = r1
        didt = r2
        dpldt = -r1
        dpidt = -r2
        return [dpdt, dldt, didt, dpldt, dpidt]

    ode_result = solve_ivp(
        lambda t, y: ode(y, t, kdpl, kdpi),
        interval,
        [p, l, i, 0.0, 0.0],
        rtol=1e-12,
        atol=1e-12,
    ).y[:, -1]
    return {
        "p": ode_result[0],
        "l": ode_result[1],
        "i": ode_result[2],
        "pl": ode_result[3],
        "pi": ode_result[4],
    }


# Homodimer formation - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
def system03_kinetic(p, kdpp, interval=(0, 100)):
    def ode(concs, t, kdpp):
        p, pp = concs
        r1 = -(p * p) + kdpp * pp
        dpdt = 2 * r1
        dppdt = -r1
        return [dpdt, dppdt]

    ode_result = solve_ivp(
        lambda t, y: ode(y, t, kdpp), interval, [p, 0.0], rtol=1e-12, atol=1e-12
    ).y[:, -1]
    return {"p": ode_result[0], "pp": ode_result[1]}


# Homodimer breaking - see https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
def system04_kinetic(p, i, kdpp, kdpi, interval=(0, 100)):
    def ode(concs, t, kdpp, kdpi):
        p, i, pp, pi = concs
        r_pp = -(p * p) + kdpp * pp
        r_pi = -p * i + kdpi * pi
        dpdt = 2 * r_pp + r_pi
        didt = r_pi
        dppdt = -r_pp
        dpldt = -r_pi
        return [dpdt, didt, dppdt, didt]

    ode_result = solve_ivp(
        lambda t, y: ode(y, t, kdpp, kdpi),
        interval,
        [p, i, 0.0, 0.0],
        rtol=1e-12,
        atol=1e-12,
    ).y[:, -1]
    return {
        "p": ode_result[0],
        "i": ode_result[1],
        "pp": ode_result[2],
        "pi": ode_result[3],
    }


class System_kinetic_one_to_one__pl(BindingSystem):
    """
    Kinetic 1:1 binding system

    Class defines 1:1 binding, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_1to1.html
    """

    def __init__(self):
        super().__init__(system01_kinetic)
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


class System_kinetic_competition_pl(BindingSystem):
    """
    Kinetic 1:1:1 competition binding system

    Class defines 1:1:1 competition, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_competition.html
    """

    def __init__(self):
        super().__init__(system02_kinetic)
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


class System_kinetic_homodimerformation__pp(BindingSystem):
    """
    Kinetic homodimer formation system

    Class defines homodimer formation, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerformation.html
    """

    def __init__(self):
        super().__init__(system03_kinetic, False)
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


class System_kinetic_homodimerbreaking__pp(BindingSystem):
    """
    Kinetic homodimer breaking system

    Class defines homodimer breaking, readout is PP
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system04_kinetic, False)
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


class System_kinetic_homodimerbreaking__pl(BindingSystem):
    """
    Kinetic homodimer breaking system

    Class defines homodimer breaking, readout is PL
    See https://stevenshave.github.io/pybindingcurve/simulate_homodimerbreaking.html
    """

    def __init__(self):
        super().__init__(system04_kinetic, False)
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
