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

    def query(self, parameters, readout, normalised_to):
        have_ymin_ymax=self._are_ymin_ymax_present(parameters)
        

        results=None      
        # Remove ymin ymax keys if present
        parameters_no_ymin_ymax=dict(parameters)
        if 'ymin' in parameters_no_ymin_ymax.keys():
            del parameters_no_ymin_ymax['ymin']
        if 'ymax' in parameters_no_ymin_ymax.keys():
            del parameters_no_ymin_ymax['ymax']
        # Check that all required parameters are present and abort if not.  Dont worry about all_concs which is used for multiple queries
        missing = sorted(
            list((set(self.arguments) - set(parameters_no_ymin_ymax.keys()))-set(['all_concs'])))

        if len(missing) > 0:
            print("The following parameters were missing!", missing)
            return None

        # Are any parameters changing?
        changing_parameters = self._find_changing_parameters(parameters_no_ymin_ymax)
        if changing_parameters is None:
            # Querying single point
            if self.analytical:
                results=self.system(**parameters_no_ymin_ymax)
            else:
                results=self.system(**parameters_no_ymin_ymax)[readout]
        else:
            # At least 1 changing parameter
            if len(changing_parameters) == 1:
                # 1 changing parameter
                if self.analytical:  # Using an analytical solution
                    results=self.system(**parameters_no_ymin_ymax)
                else:  # Changing parameter on kinetic solution
                    results = np.empty(
                        len(parameters_no_ymin_ymax[changing_parameters[0]]))
                    for i in range(results.shape[0]):
                        tmp_params = dict(parameters_no_ymin_ymax)
                        tmp_params[changing_parameters[0]
                                   ] = parameters_no_ymin_ymax[changing_parameters[0]][i]
                        results[i] = self.system(**tmp_params)[readout]
            else:
                print("Only 1 parameter may change, currently changing: ",
                      changing_parameters)
                return None
        if have_ymin_ymax:
            if type(normalised_to) is str:
                with np.errstate(divide='ignore', invalid='ignore'):
                    results=parameters['ymin']+((parameters['ymax']-parameters['ymin'])*results)/parameters[normalised_to]
                results[results==-np.inf]=0
            else:
                with np.errstate(divide='ignore', invalid='ignore'):
                    results=parameters['ymin']+((parameters['ymax']-parameters['ymin'])*results)/normalised_to
                results[results==-np.inf]=0

        return np.nan_to_num(results)

    def _are_ymin_ymax_present(self, parameters:dict):
        if 'ymin' in parameters.keys():
            if 'ymax' in parameters.keys():
                return True
            else:
                print("Warning: Only ymin was in parameters, missing ymax, setting ymax to 1.0")
                parameters['ymax']=1.0
                return True
        else:
            if 'ymax' in parameters.keys():
                print("Warning: Only ymax was in parameters, missing ymin, setting ymin to 0.0")
                parameters['ymin']=0.0
                return True
        return False
                    


class System_one_to_one_analytical_pl(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system01_p_l_kd__pl, True)
        self.analytical = True
        self.default_readout="pl"
    def query(self, parameters:dict):
        return super().query(parameters, self.default_readout, 'l')

class System_homodimer_formation_analytical_pp(BindingSystem):
    def __init__(self):
        super().__init__(analyticalsystems.system03_p_kdpp__pp, True)
        self.analytical = True
        self.default_readout="pp"
    def query(self, parameters:dict):
        return super().query(parameters,self.default_readout, parameters['p']/2)

class System_one_to_one(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system01_p_l_kd__pl, False)
    def query(self, parameters:dict, readout:str):
        return super().query(parameters, readout, 'l')

class System_competition(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system02_p_l_i_kdpl_kdpi__pl, False)
    def query(self, parameters:dict, readout:str):
        return super().query(parameters, readout, 'l')


class System_homodimer_formation(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system03_p_kdpp__pp, False)
    def query(self, parameters:dict, readout:str):
        return super().query(parameters, readout, 'pp')


class System_homodimer_breaking(BindingSystem):
    def __init__(self):
        super().__init__(kineticsystems.system04_p_i_kdpp_kdpi__pp, False)
    def query(self, parameters:dict, readout:str):
        return super().query(parameters, readout, 'pp')
