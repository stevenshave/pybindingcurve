from pybindingcurve.systems import analyticalsystems, kineticsystems
from inspect import signature
import numpy as np




class BindingSystem():
    system = None
    analytical = False
    arguments = []
    default_readout = None

    def scale_ymin_ymax(self, ymin, ymax, divisor, params_no_ymin_ymax):
        return ymin + ((ymax - ymin) * self.system(**params_no_ymin_ymax)) / divisor

    def _find_changing_parameters(self, params: dict):
        changing_list = []
        for p in params.keys():
            if isinstance(
                    params[p],
                    np.ndarray) or isinstance(
                    params[p],
                    list):
                changing_list.append(p)
        if len(changing_list) == 0:
            return None
        else:
            return changing_list

    def __init__(self, bindingsystem, analytical=False):
        self._system = bindingsystem
        self.analytical = analytical
        self.arguments = list(signature(bindingsystem).parameters.keys())

        if not analytical and "interval" in self.arguments:
            # Make sure interval is not included in arguments
            self.arguments.remove("interval")


    def _remove_ymin_ymax_keys_from_dict(self, d):
        if 'ymin' in d.keys():
            del d['ymin']
        if 'ymax' in d.keys():
            del d['ymax']
        return d




    def query(self, parameters):
        results = None
        # Remove ymin ymax keys if present
        parameters_no_ymin_ymax = self._remove_ymin_ymax_keys_from_dict(dict(parameters))
        
        # Check that all required parameters are present and abort if not.
        # Dont worry about all_concs which is used for multiple queries
        missing = sorted(list((set(self.arguments) -
                               set(parameters_no_ymin_ymax.keys())) -
                              set(['all_concs'])))
        if len(missing) > 0:
            print("The following parameters were missing!", missing)
            return None

        # Are any parameters changing?
        changing_parameters = self._find_changing_parameters(
            parameters_no_ymin_ymax)
        if changing_parameters is None:  # Querying single point
            if self.analytical:
                results = self._system(**parameters_no_ymin_ymax)  # Analytical
            else:
                results = self._system(
                    **parameters_no_ymin_ymax)[self.default_readout]  # Kinetic
        else: 
            # At least 1 changing parameter
            if len(changing_parameters) == 1:
                results = np.empty(len(parameters_no_ymin_ymax[changing_parameters[0]]))
                # 1 changing parameter
                if self.analytical:  # Using an analytical solution
                    for i in range(results.shape[0]):
                        tmp_params = dict(parameters_no_ymin_ymax)
                        tmp_params[changing_parameters[0]] = parameters_no_ymin_ymax[changing_parameters[0]][i]
                        results[i] = self._system(**tmp_params)
                else:  # Changing parameter on kinetic solution
                    for i in range(results.shape[0]):
                        tmp_params = dict(parameters_no_ymin_ymax)
                        tmp_params[changing_parameters[0]] = parameters_no_ymin_ymax[changing_parameters[0]][i]
                        results[i] = self._system(**tmp_params)[self.default_readout]
            else:
                print("Only 1 parameter may change, currently changing: ",
                      changing_parameters)
                return None
        return np.nan_to_num(results)

    def _are_ymin_ymax_present(self, parameters: dict):
        if 'ymin' in parameters.keys():
            if 'ymax' in parameters.keys():
                return True
            else:
                print(
                    "Warning: Only ymin was in parameters, missing ymax, setting ymax to 1.0")
                parameters['ymax'] = 1.0
                return True
        else:
            if 'ymax' in parameters.keys():
                print(
                    "Warning: Only ymax was in parameters, missing ymin, setting ymin to 0.0")
                parameters['ymin'] = 0.0
                return True
        return False


class System_analytical_one_to_one_pl(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system01_one_to_one__p_l_kd__pl, analytical=True)
        self.analytical = True
        self.default_readout = "pl"
    def query(self, parameters: dict):
        return super().query(parameters)

class System_analytical_one_to_one_yl(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system01_one_to_one__p_l_kd__pl, analytical=True)
        self.analytical = True
        self.default_readout = "pl"
    def query(self, parameters: dict):
        parameters_no_min_max=self._remove_ymin_ymax_keys_from_dict(dict(parameters))
        
        return self.super().query(parameters_no_min_max)

class System_analytical_competition_pl(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system02_competition__p_l_i_kdpl_kdpi__pl, analytical=True)
        self.default_readout = "pl"
    def query(self, parameters: dict):
        return super().query(parameters)

class System_analytical_homodimerformation_pp(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system03_homodimer_formation__p_kdpp__pp, analytical=True)
        self.analytical = True
        self.default_readout = "pp"
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_one_to_one_pl(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system01_p_l_kd__pl)
        self.default_readout='pl'
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_one_to_one_p(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system01_p_l_kd__pl)
        self.default_readout='p'
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_one_to_one_l(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system01_p_l_kd__pl)
        self.default_readout='l'
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_competition_pl(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system02_p_l_i_kdpl_kdpi__pl)
        self.default_readout='pl'
    def query(self, parameters: dict):
        return super().query(parameters)

class System_kinetic_competition_pi(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system02_p_l_i_kdpl_kdpi__pl)
        self.default_readout='pi'
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_competition_p(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system02_p_l_i_kdpl_kdpi__pl)
        self.default_readout='p'
    def query(self, parameters: dict):
        return super().query(parameters)

class System_kinetic_competition_l(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system02_p_l_i_kdpl_kdpi__pl)
        self.default_readout='l'
    def query(self, parameters: dict):
        return super().query(parameters)

class System_kinetic_competition_i(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system02_p_l_i_kdpl_kdpi__pl)
        self.default_readout='i'
    def query(self, parameters: dict):
        return super().query(parameters)

class System_kinetic_homodimerformation_pp(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system03_p_kdpp__pp, False)
        self.default_readout='pp'
    def query(self, parameters: dict):
        return super().query(parameters)

class System_kinetic_homodimerformation_p(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system03_p_kdpp__pp, False)
        self.default_readout='p'
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_homodimerbreaking_pp(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system04_p_i_kdpp_kdpi__pp, False)
        self.default_readout='pp'
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_homodimerbreaking_p(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system04_p_i_kdpp_kdpi__pp, False)
        self.default_readout='p'
    def query(self, parameters: dict):
        return super().query(parameters)


class System_kinetic_homodimerbreaking_l(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system04_p_i_kdpp_kdpi__pp, False)
        self.default_readout='l'
    def query(self, parameters: dict):
        return super().query(parameters)

class System_kinetic_homodimerbreaking_pl(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system04_p_i_kdpp_kdpi__pp, False)
        self.default_readout='pl'
    def query(self, parameters: dict):
        return super().query(parameters)
