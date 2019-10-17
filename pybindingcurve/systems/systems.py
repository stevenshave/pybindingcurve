from pybindingcurve.systems import analyticalsystems, kineticsystems
from inspect import signature
import numpy as np


class BindingSystem():
    system = None
    analytical = False
    arguments = []
    default_readout = None

    def _find_changing_parameters(self, params: dict):
        changing_list = []
        for p in params.keys():
            if type(params[p]) == np.ndarray or type(params[p]) == list:
                changing_list.append(p)
        if len(changing_list) == 0:
            return None
        else:
            return changing_list

    def __init__(self, _system, _analytical=False):
        self.system = _system
        self.analytical = _analytical
        self.arguments = list(signature(_system).parameters.keys())
        if not _analytical and "interval" in self.arguments:
            self.arguments.remove("interval")

    def query(self, parameters, readout):
        # Check that all required parameters are present and abort if not.  Dont worry about all_concs which is used for multiple queries
        missing = sorted(
            list((set(self.arguments) - set(parameters.keys()))-set(['all_concs'])))

        if len(missing) > 0:
            print("The following parameters were missing!", missing)
            return None

        # Are any parameters changing?
        changing_parameters = self._find_changing_parameters(parameters)
        if changing_parameters is None:
            # Querying single point
            if self.analytical:
                return self.system(**parameters)
            else:
                return self.system(**parameters)[readout]

        else:
            # At least 1 changing parameter
            if len(changing_parameters) == 1:
                # 1 changing parameter
                if self.analytical:  # Using an analytical solution
                    return self.system(**parameters)
                else:  # Changing parameter on kinetic solution
                    results_array = np.empty(
                        len(parameters[changing_parameters[0]]))
                    for i in range(results_array.shape[0]):
                        tmp_params = dict(parameters)
                        tmp_params[changing_parameters[0]
                                   ] = parameters[changing_parameters[0]][i]
                        results_array[i] = self.system(**tmp_params)[readout]
                    return results_array

            else:
                print("Only 1 parameter may change, currently changing: ",
                      changing_parameters)
                return None


class System_one_to_one_analytical_pl(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system01_p_l_kd__pl, True)
        self.analytical = True
        self.default_readout="pl"

    def query(self, parameters):
        return super().query(parameters, self.default_readout)


class System_homodimer_formation_analytical_pp(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system03_p_kdpp__pp, True)
        self.analytical = True
        self.default_readout="pp"

    def query(self, parameters):
        return super().query(parameters,self.default_readout)


class System_one_to_one(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system01_p_l_kd__pl, False)

    def query(self, parameters, readout):
        return super().query(parameters, readout)


class System_competition(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system02_p_l_i_kdpl_kdpi__pl, False)

    def query(self, parameters, readout):
        return super().query(parameters, readout)


class System_homodimer_formation(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system03_p_kdpp__pp, False)

    def query(self, parameters, readout):
        return super().query(parameters, readout)


class System_homodimer_breaking(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system04_p_i_kdpp_kdpi__pp, False)

    def query(self, parameters, readout):
        return super().query(parameters, readout)
