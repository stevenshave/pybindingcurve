from inspect import signature
import numpy as np


class BindingSystem:
    """
    BindingSystem class, used to determine the type of binding systems being used

    BindingSystem objects are governed with differnt protein-ligand binding
    system being represented. It also provides methods for visualisation,
    querying and the fitting of system parameters.

    Parameters
    ----------
    bindingsystem : func
        A function deriving the concentration of complex formed. This fuction
        may be varied based on different protein-ligand binding systems.
    analytical: bool
        Perform a analytical analysis (default = False)
    """

    system = None
    analytical = False
    arguments = []
    default_readout = None

    def _find_changing_parameters(self, params: dict):
        """
        Find the changing parameter

        Determine which parameter is changing with titration, including the
        concentration of the protein or ligand. A set of varied concentrations
        of parameter are in the form of array-like or list data type

        Parameters
        ----------
        params : dict
            Parameters defining the binding system to be simulated.

        Returns
        -------
        A list containing the keys of params indicating the name of
            changing parameters.
        """
        changing_list = []
        for p in params.keys():
            if isinstance(params[p], np.ndarray) or isinstance(params[p], list):
                changing_list.append(p)
        if len(changing_list) == 0:
            return None
        else:
            return changing_list

    def __init__(self, bindingsystem: callable, analytical: bool = False):
        """
        Construct BindingSystem objects

        BindingSystem objects are initialised with various subclasses under
        the BindingSystem super class definding different protein-ligand
        binding systems, such as 'System_analytical_one_to_one_pl' or
        'System_analytical_homodimerformation_pp', etc.

        Parameters
        ----------
        bindingsystem : func
            A function deriving the concentration of complex formed. This
            fuction may be varied based on different protein-ligand binding
            systems.
        analytical: bool
            Perform a analytical analysis (default = False)
        """
        self._system = bindingsystem
        self.analytical = analytical
        self.arguments = list(signature(bindingsystem).parameters.keys())

        if not analytical and "interval" in self.arguments:
            # Make sure interval is not included in arguments
            self.arguments.remove("interval")

    def _remove_ymin_ymax_keys_from_dict_in_place(self, d: dict):
        """
        Remove minimum and maximum readout from the orignal system parameters
        (Replace the original dict)

        Output the original dict contains the orignal system parameters with
        removal the minimum and maximum readout signel of the system

        Parameters
        ----------
        d : dict
            Parameters defining the binding system to be simulated

        Returns
        -------
        dict
            Return the original dict contains original system parameters without
            minimum and maximum readout.
        """
        if "ymin" in d.keys():
            del d["ymin"]
        if "ymax" in d.keys():
            del d["ymax"]
        return d

    def _remove_ymin_ymax_keys_from_dict_return_new(self, input: dict):
        """
        Remove minimum and maximum readout from the orignal system parameters
        (Return a new dict)

        Output a new dict contains the orignal system parameters with removal
        the minimum and maximum readout signel of the system

        Parameters
        ----------
        input : dict
            Parameters defining the binding system to be simulated

        Returns
        -------
        dict
            Generate a new dict contains original system parameters without
            minimum and maximum readout.
        """
        d = dict(input)
        if "ymin" in d.keys():
            del d["ymin"]
        if "ymax" in d.keys():
            del d["ymax"]
        return d

    def query(self, parameters: dict):
        """
        Query a binding system

        Get the readout from from a set of system parameters

        Parameters
        ----------
        parameters : dict
            Parameters defining the binding system to be simulated.

        Returns
        -------
        Single floating point of the concentration of the binding complex, or
            array-like Response/signal of the system.
        """
        results = None

        # Check that all required parameters are present and abort if not.
        # Dont worry about all_concs which is used for multiple queries
        missing = sorted(
            list((set(self.arguments) - set(parameters.keys())) - set(["all_concs"]))
        )
        assert len(missing)==0, f"The following parameters were missing: {missing}"
        
        # Are any parameters changing?
        changing_parameters = self._find_changing_parameters(parameters)
        if changing_parameters is None:  # Querying single point
            if self.analytical:
                results = self._system(**parameters)  # Analytical
            else:
                results = self._system(**parameters)[self.default_readout]  # Kinetic
        else:
            # At least 1 changing parameter
            if len(changing_parameters) == 1:
                results = None
                num_solutions = getattr(
                    self, "num_solutions", 1
                )  # Get attribure of num_solutions in BindingSystem class (default = 1)
                if num_solutions == 1:
                    results = np.empty(len(parameters[changing_parameters[0]]))
                else:
                    results = np.empty(
                        (len(parameters[changing_parameters[0]]), num_solutions)
                    )
                # 1 changing parameter
                if self.analytical:  # Using an analytical solution
                    for i in range(results.shape[0]):
                        tmp_params = dict(parameters)
                        tmp_params[changing_parameters[0]] = parameters[
                            changing_parameters[0]
                        ][i]
                        results[i] = self._system(**tmp_params)
                else:  # Changing parameter on kinetic solution
                    for i in range(results.shape[0]):
                        tmp_params = dict(parameters)
                        tmp_params[changing_parameters[0]] = parameters[
                            changing_parameters[0]
                        ][i]
                        results[i] = self._system(**tmp_params)[self.default_readout]
            else:
                print(
                    "Only 1 parameter may change, currently changing: ",
                    changing_parameters,
                )
                return None
        if isinstance(results, (np.ndarray)):
            if results.ndim > 1:
                results = results.T
            return np.nan_to_num(results)
        return results  # We get here if its not a numpy array, but a system with multiple solutions queried at for a single point

    def _are_ymin_ymax_present(self, parameters: dict):
        """
        Check the existance of the minimum or maximum readout

        Check the minimum or maximum readout signal set by user is in
        the binding system parameters

        Parameters
        ----------
        parameters : dict
            Parameters defining the binding system to be simulated.

        Returns
        -------
        Boolean (True or False)
        """
        if "ymin" in parameters.keys():
            if "ymax" in parameters.keys():
                return True
            else:
                print(
                    "Warning: Only ymin was in parameters, missing ymax, setting ymax to 1.0"
                )
                parameters["ymax"] = 1.0
                return True
        else:
            if "ymax" in parameters.keys():
                print(
                    "Warning: Only ymax was in parameters, missing ymin, setting ymin to 0.0"
                )
                parameters["ymin"] = 0.0
                return True
        return False
