"""PyBindingCurve - a package to simulate protein-ligand binding systems

"""


import numpy as np
import matplotlib.pyplot as plt
import lmfit
from pybindingcurve.systems import *
from typing import Union

pbc_plot_style = {
    "axis_label_size": 12,
    "axis_label_font": "DejaVu Sans",
    "title_size": 12,
    "title_font": "DejaVu Sans",
    # 'figure_width': 9,
    # 'figure_height': 9,
    "x_tick_label_font_size": 10,
    "y_tick_label_font_size": 10,
    "legend_font_size": 9,
    "dpi": 300,
    "x_axis_labelpad": None,
    "y_axis_labelpad": None,
    "title_labelpad": None,
    "fig_size": (5 * 1.2, 4 * 1.2),
}


class Readout:
    """
    Class to change the system simulation readouts

    Readout is a container for static functions that change the system
    readout when passed as an argument to add_curve. These functions
    take as arguments the system parameters, and then the calculated
    y-values as raw output from systems equations.  Return types are either
    arrays or singular values.
    """

    @staticmethod
    def fraction_l(system_parameters: dict, y):
        """
        Transform readout to fraction ligand bound

        Readout function to change y values into fraction ligand bound

        Parameters
        ----------
        system_parameters : dict
            Dictionary of system parameters defining the system
        y : float or array-like
            System simulation or query values

        Returns
        -------
        tuple
            Tuple's first value is the y-axis label which should be used with
            this trasformation when plotting.
            The next value is either the single, or array like object
            containing the transformed values.
        """

        return "Fraction l bound", y / system_parameters["l"]

    @staticmethod
    def fraction_possible_dimer(system_parameters: dict, y):
        """
        Transform readout to fraction possible dimer

        Readout function to change y values into fraction possible dimer.
        As maximal achievable dimer concentration is half the total monomer
        concentration, we must take this into account when calculating fraction
        possible dimer

        Parameters
        ----------
        system_parameters : dict
            Dictionary of system parameters defining the system
        y : float or array-like
            System simulation or query values

        Returns
        -------
        tuple
            Tuple's first value is the y-axis label which should be used with
            this trasformation when plotting.
            The next value is either the single, or array like object
            containing the transformed values.
        """
        return "Fraction possible dimer", y / (system_parameters["p"] / 2.0)

    @staticmethod
    def complex_concentration(system_parameters: dict, y):
        """
        Readout as complex concentration, (redundant as NULL achieves the same)

        Redundant readout function, never used, but can be used for
        better understanding Readouts.

        Parameters
        ----------
        system_parameters : dict
            Dictionary of system parameters defining the system
        y : float or array-like
            System simulation or query values

        Returns
        -------
        tuple
            Tuple's first value is the y-axis label which should be used with
            this trasformation when plotting, although in this case the value
            is None to indicate no transformation has been carried out and the
            default which is assigned by the system should be used.
            The next value is either the single, or array like object
            containing the non-transformed values.
        """

        return None, y


class _Curve:
    """
    Curve class, represents a binding curve

    Consists of X and Y coordinates, along with a name and boolean
    flag denoting that the curve is a traced, real physical solution.
    """

    def __init__(self, xcoords: np.array, ycoords: np.array, series_name: str = ""):
        """
        Curve constructor

        Called by PyBindingCurve

        Parameters
        ----------
        xcoors : np.array
            X coordinates of the binding system to be used to present a binding curve
        ycoors : np.array
            Y coordinates of the binding system to be used to present a binding curve
        series_name : str, Optional
            Name of curve to appear in plot legends
        """
        self.xcoords = xcoords
        self.ycoords = ycoords
        self.name = series_name


class BindingCurve:
    """
    BindingCurve class, used to simulate systems

    BindingCurve objects are governed by their underlying system, defining the
    (usually) protein-ligand binding system being represented. It also
    provides the main interface for simulation, visualisation, querying and
    the fitting of system parameters.

    Parameters
    ----------
    binding_system : BindingSystem or str
        Define the binding system which will govern this BindingCurve object.
        Can either be a BindingSystem object, a shortcut string describing a
        system (such as '1:1' or 'competition', etc), or a custom binding
        system definition string.
    """
    arguments=None
    system = None
    _last_custom_readout = None
    curves = []
    fig = None
    axes = None
    plot_solution_colours = list("krgbycmrgbycmrgby") + list(
        np.linspace(0.1, 0.9, num=20)
    )
    _min_x_axis = 0.0
    _max_x_axis = 0.0
    _min_y_axis = 0.0
    _max_y_axis = 0.0
    _num_added_traces = 0
    _num_added_sets_of_points = 0
    _last_known_changing_parameter = "X"


    def query(self, parameters, readout: Readout = None):
        """
        Query a binding system

        Get the readout from from a set of system parameters

        Parameters
        ----------
        parameters : dict
            System parameters defining the system being queried.  Will usually
            contain protein, ligand etc concentrations, and KDs.
        readout : func or None
            Change the readout of the system, can be None for unmodified
            (usually complex concentration), a static member function from
            the pbc.Readout class, or a custom written function following the
            the same defininition as those in pbc.Readout.

        Returns
        -------
        Single floating point, or array-like
            Response/signal of the system
        """

        # Compound return statement, return the query, otherwide apply the
        # readout object to it and return the second returned variable
        return (
            self.system.query(parameters)
            if readout is None
            else readout(parameters, self.system.query(parameters))[1]
        )

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
        changing_list = [
            p
            for p in params.keys()
            if isinstance(params[p], np.ndarray) or isinstance(params[p], list)
        ]
        return changing_list if len(changing_list) > 0 else None

    def __init__(self, binding_system:Union[str, BindingSystem], disable_signal_warnings:bool=False):
        """
        Construct BindingCurve objects

        BindingCurve objects are initialised with different protein-ligand
        binding systems definded by user. This method can also initialises
        the objects as an interface throught the pbc.BindingSystem class,
        providing for simulation, visualisation, querying and the fitting
        of system parameters.

        Parameters
        ----------
        binding_system : str or BindingSystem
            Define the binding system which will govern this BindingCurve
            object. Can either be a BindingSystem object or a human readable
            string shortcut, such as '1:1' or 'competition', etc.
        disable_signal_warnings : bool
            If True, then no warning is given if y values are larger than x
            values. This was put in originally to force people to put in a
            ymin value to their system to specify a signal, rather than a
            concentration is being used.
        """
        self.disable_signal_warning=disable_signal_warnings
        if isinstance(binding_system, str):
            # Check if its a custom defined system - containing <->
            if binding_system.find("<->") != -1:
                # Custom binding system, check if it starts with Lagrange
                if binding_system.lower().startswith("lagrange"):
                    self.system = System_lagrange_custom(binding_system[8:])
                else:
                    self.system = System_minimizer_custom(binding_system)

            binding_system = binding_system.lower().replace(" ", "")
            # 1:1
            if binding_system in ["simple", "1:1", "1:1analytical"]:
                self.system = System_analytical_one_to_one__pl()
            # 1:1 lagrange
            if binding_system in ["simplelagrange", "1:1lagrange"]:
                self.system = System_lagrange_one_to_one__pl()
            # 1:1 kinetic
            if binding_system in ["simplekinetic", "1:1kinetic"]:
                self.system = System_kinetic_one_to_one__pl()
            # 1:1 minimised
            if binding_system in [
                "simplemin",
                "simpleminimized",
                "simpleminimised",
                "1:1min",
                "1:1minimized",
                "1:1minimised",
            ]:
                self.system = System_minimizer_one_to_one__pl()

            # Homodimer formation
            if binding_system in ["homodimerformation", "homodimer"]:
                self.system = System_analytical_homodimerformation__pp()
            # Homodimer formation lagrange
            if binding_system in ["homodimerformationlagrange", "homodimerlagrange"]:
                self.system = System_lagrange_homodimerformation__pp()
            # Homodimer formation kinetic
            if binding_system in ["homodimerformationkinetic", "homodimerkinetic"]:
                self.system = System_kinetic_homodimerformation__pp()
            # Homodimer formation minimizer
            if binding_system in [
                "homodimerformationmin",
                "homodimerformationminimiser",
                "homodimerformationminimizer",
                "homodimermin",
                "homodimerminimiser",
                "homodimerminimizer",
            ]:
                self.system = System_minimizer_homodimerformation__pp()

            # Competition
            if binding_system in ["competition", "1:1:1", "competitionanalytical", "1:1:1analytical"]:
                self.system = System_analytical_competition__pl()
            if binding_system in ["competitionlagrange", "1:1:1lagrange"]:
                self.system = System_lagrange_competition_pl()
            if binding_system in ["competitionkinetic", "1:1:1kinetic"]:
                self.system = System_kinetic_competition_pl()
            if binding_system in [
                "competitionmin",
                "competitionminimiser",
                "competitionminimizer",
                "1:1:1min",
                "1:1:1minimiser",
                "1:1:1minimizer",
            ]:
                self.system = System_minimizer_competition__pl()

            # Homodimer breaking minimizer
            if binding_system in ["homodimerbreaking", "homodimerbreakingmin"]:
                self.system = System_minimizer_homodimerbreaking__pp()
            # Homodimer breaking lagrange
            if binding_system in ["homodimerbreakinglagrange"]:
                self.system = System_lagrange_homodimerbreaking__pp()
            # Homodimer breaking analytical
            if binding_system in ["homodimerbreakinganalytical"]:
                self.system = System_analytical_homodimerbreaking_pp()
            # Homodimer breaking kinetic
            if binding_system in ["homodimerbreakingkinetic"]:
                self.system = System_kinetic_homodimerbreaking__pp()

            # 1:2 minimizer
            if binding_system in ["1:2", "1:2min", "1:2minimiser", "1:2minimizer"]:
                self.system = System_minimizer_1_to_2__pl12()
            # 1:2 lagrange
            if binding_system in ["1:2lagrange"]:
                self.system = System_lagrange_1_to_2__pl12()

            # 1:3 minimizer
            if binding_system in ["1:3", "1:3min", "1:3minimiser", "1:3minimizer"]:
                self.system = System_minimizer_1_to_3__pl123()
            # 1:3 lagrange
            if binding_system in ["1:3lagrange"]:
                self.system = System_lagrange_1_to_3__pl123()
        else:
            if issubclass(binding_system, BindingSystem):
                self.system = binding_system()
        assert self.system is not None, "Invalid system specified, try one of: ['simple', 'homodimer formation', 'competition', 'homdimer breaking'], pass a system object, or define a custom binding system"

    def _initialize_plot(self):
        """
        Initialise setup to being ready for curve plotting

        Control setups needed for plotting a binding plot including layouts
        of subplots, grid lines, and y-axis view limits, etc.
        """

        if self.fig is None:

            # Here we reset a lot as the user may be making multiple
            # plots one after the other
            self._last_custom_readout = None
            self.curves = []
            self.axes = None

            self._min_x_axis = 0.0
            self._max_x_axis = 0.0
            self._min_y_axis = 0.0
            self._max_y_axis = 0.0
            self._num_added_traces = 0
            self._num_added_sets_of_points = 0
            self._last_known_changing_parameter = "X"

            self.fig, self.axes = plt.subplots(
                nrows=1, ncols=1, figsize=pbc_plot_style["fig_size"]
            )
            self.axes.grid(True, which="both")
            self.axes.set_ylim(0, 1)
            plt.tight_layout(rect=(0.05, 0.05, 0.95, 0.92))

    def add_curve(self, parameters: dict, name: str = None, readout: Readout = None):
        """
        Add a curve to the plot

        Add a curve as specified by the system parameters to the
        pbc.BindingSystem's internal plot using the underlying binding system
        specified on intitialisation.

        Parameters
        ----------
        parameters : dict
            Parameters defining the system to be simulated
        name : str or None, optional
            Name of curve to appear in plot legends
        readout : Readout.function, optional
            Change the system readout to one described by a custom readout
            function.  Predefined standard readouts can be found in the static
            pbc.Readout class.
        """
        if self.system is None:
            print("No system defined, could not proceed")
            return None
        self._initialize_plot()
        changing_parameters = self._find_changing_parameters(parameters)
        if changing_parameters is None:
            print("No changing parameters detected. Plotting curves requires something to be changing")
            return
        if not len(changing_parameters) == 1:
            print("Must have 1 changing parameter, no curves added.")
            return

        y_values = self.system.query(parameters)

        if readout is not None:
            self._last_custom_readout, y_values = readout(parameters, y_values)

        # It may be that we have multiple solutions from a direct analytical
        # solution.  If so, then we need to add all curves.
        if y_values.ndim > 1:  # Multiple solutions
            for i in range(y_values.ndim):
                self.curves.append(
                    _Curve(parameters[changing_parameters[0]], y_values[i])
                )
        else:
            self.curves.append(  # Only one solution
                _Curve(parameters[changing_parameters[0]], y_values)
            )
        self._last_known_changing_parameter = changing_parameters[0]

        for curve_it, curve in enumerate(self.curves[self._num_added_traces :]):
            self._num_added_traces += 1
            curve_name_with_number = None
            if name is None:
                curve_name_with_number = f"Curve {self._num_added_traces}"
            else:
                if y_values.ndim == 1:
                    curve_name_with_number = name
                else:
                    curve_name_with_number = name + " " + str(curve_it + 1)
            self.axes.plot(
                curve.xcoords,
                curve.ycoords,
                self.plot_solution_colours[self._num_added_traces] + "-",
                label=curve_name_with_number,
                linewidth=2,
            )
            self._max_x_axis = np.nanmax(
                [self._max_x_axis, parameters[changing_parameters[0]][-1]]
            )
            self._min_x_axis = np.nanmin(
                [self._min_x_axis, parameters[changing_parameters[0]][0]]
            )
            self._min_y_axis = np.nanmin([self._min_y_axis, np.nanmin(curve.ycoords)])
            self._max_y_axis = np.nanmax([self._max_y_axis, np.nanmax(curve.ycoords)])

    def add_scatter(self, xcoords, ycoords, name: str = None):
        """
        Add scatterpoints to a plot, useful to represent real measurement data

        X and Y coordinates may be added to the internal plot, useful when
        fitting to experimental data, and wanting to plot the true experimental
        values alongside a curve generated with fitted parameters.

        Parameters
        ----------
        xcoords : list or array-like
            x-coordinates
        ycoords : list or array-like
            y-coordinates
        name : str or None, optional
            Name of series to appear in plot legends
        """
        self._initialize_plot()
        self._num_added_sets_of_points += 1
        if name is None:
            name = f"Data " + str(self._num_added_sets_of_points)
        self.axes.scatter(xcoords, ycoords, label=name)
        self._min_x_axis = min(self._min_x_axis, np.min(np.real(xcoords)))
        self._max_x_axis = max(self._max_x_axis, np.max(np.real(xcoords)))
        self._min_y_axis = min(self._min_y_axis, np.min(np.real(ycoords)))
        self._max_y_axis = max(self._max_y_axis, np.max(np.real(ycoords)))

    def show_plot(
        self,
        title: str = "System simulation",
        xlabel: str = None,
        ylabel: str = None,
        min_x: float = None,
        max_x: float = None,
        min_y: float = None,
        max_y: float = None,
        log_x_axis: bool = False,
        log_y_axis: bool = False,
        pbc_plot_style: dict = pbc_plot_style,
        png_filename: str = None,
        svg_filename: str = None,
        show_legend: bool = True,
    ):
        """
        Show the PyBindingCurve plot

        Function to display the internal state of the pbc BindingCurve objects
        plot.

        Parameters
        ----------
        title : str
            The title of the plot (default = "System simulation")
        xlabel: str
            X-axis label (default = None)
        ylabel : str
            Y-axis label (default = None, causing label to be "[Complex]")
        min_x : float
            X-axis minimum (default = None)
        max_x : float
            X-axis maximum (default = None)
        min_y : float
            Y-axis minimum (default = None)
        max_y : float
            Y-axis maximum (default = None)
        log_x_axis : bool
            Log scale on X-axis (default = False)
        log_y_axis : bool
            Log scale on Y-axis (default = False)
        ma_style : bool
            Apply MA styling, making plots appear like GraFit plots
        png_filename :  str
            File name/location where png will be written
        svg_filename : str
            File name/location where svg will be written
        """
        assert self.fig is not None, "Nothing to plot"

        if min_x is not None:
            self._min_x_axis = min_x
        if max_x is not None:
            self._max_x_axis = max_x
        if min_y is not None:
            self._min_y_axis = min_y
        if max_y is not None:
            self._max_y_axis = max_y

        if max_y is None:
            self.axes.set_ylim(self._min_y_axis, self._max_y_axis * 1.1)
        else:
            self.axes.set_ylim(self._min_y_axis, self._max_y_axis)
        self.axes.set_xlim(self._min_x_axis, self._max_x_axis)
        if log_x_axis:
            self.axes.set_xscale("log", nonposx="clip")
        if log_y_axis:
            self.axes.set_yscale("log", nonposx="clip")

        if xlabel is None:
            self.axes.set_xlabel(
                "[" + self._last_known_changing_parameter.upper() + "]",
                fontsize=pbc_plot_style["axis_label_size"],
                fontname=pbc_plot_style["axis_label_font"],
                labelpad=pbc_plot_style["x_axis_labelpad"],
            )
        else:
            self.axes.set_xlabel(
                xlabel,
                fontsize=pbc_plot_style["axis_label_size"],
                fontname=pbc_plot_style["axis_label_font"],
                labelpad=pbc_plot_style["x_axis_labelpad"],
            )

        if ylabel is None:
            if self._last_custom_readout is None:
                self.axes.set_ylabel(
                    "[" + self.system.default_readout.upper() + "]",
                    fontsize=pbc_plot_style["axis_label_size"],
                    fontname=pbc_plot_style["axis_label_font"],
                    labelpad=pbc_plot_style["y_axis_labelpad"],
                )
            else:
                self.axes.set_ylabel(
                    self._last_custom_readout,
                    fontsize=pbc_plot_style["axis_label_size"],
                    fontname=pbc_plot_style["axis_label_font"],
                    labelpad=pbc_plot_style["y_axis_labelpad"],
                )
        else:
            self.axes.set_ylabel(
                ylabel,
                fontsize=pbc_plot_style["axis_label_size"],
                fontname=pbc_plot_style["axis_label_font"],
                labelpad=pbc_plot_style["y_axis_labelpad"],
            )

        self.axes.set_title(
            title,
            fontsize=pbc_plot_style["title_size"],
            fontname=pbc_plot_style["title_font"],
            pad=pbc_plot_style["title_labelpad"],
        )

        if show_legend:
            self.axes.legend(prop={"size": pbc_plot_style["legend_font_size"]})

        for tick in self.axes.xaxis.get_major_ticks():
            tick.label.set_fontsize(pbc_plot_style["x_tick_label_font_size"])
        for tick in self.axes.yaxis.get_major_ticks():
            tick.label.set_fontsize(pbc_plot_style["y_tick_label_font_size"])

        if png_filename is not None:
            plt.savefig(
                png_filename,
                dpi=pbc_plot_style["dpi"],
                metadata={"Title": "pyBindingCurve plot"},
            )
        if svg_filename is not None:
            plt.savefig(svg_filename, metadata={"Title": "pyBindingCurve plot"})
        plt.show()
        # Calling show displays, and then consumes the figure, so we should
        # set it to Null so that the next initialisation sets up a new figure.
        self.fig = None
        plt.clf()
        plt.cla()
        plt.close()

    def fit(
        self,
        system_parameters: dict,
        to_fit: dict,
        ycoords: np.array,
        bounds: dict = None,
    ):
        """
        Fit the parameters of a system to a set of data points

        Fit the system to a set of (usually) experimental datapoints.
        The fitted parameters are stored in the system_parameters dict
        which may be accessed after running this function.  It is
        possible to fit multiple parameters at once and define bounds
        for the parameters.  The function returns a dictionary of the
        accuracy of fitted parameters, which may be captured, or not.

        Parameters
        ----------
        system_parameters : dict
            Dictionary containing system parameters, will be used as arguments
            to the systems equations.
        to_fit : dict
            Dictionary containing system parameters to fit.
        xcoords : np.array
            X coordinates of data the system parameters should be fit to
        ycoords : np.array
            Y coordinates of data the system parameters should be fit to
        bounds : dict
            Dictionary of tuples, indexed by system parameters denoting the
            lower and upper bounds of a system parameter being fit, optional,
            default = None

        Returns
        -------
        tuple (dict, dict)
            Tuple containing a dictionary of best fit systems parameters,
            then a dictionary containing the accuracy for fitted variables.
        """

        # Make a copy of system parameters and operate on the copy, so that
        # the original remains untouched.
        system_parameters_copy = dict(system_parameters)
        # Check we have parameters to fit, and nothing is missing
        if len(to_fit.keys()) == 0:
            print("Nothing to fit, insert parameters to fit into to_fit dictionary")
            return None
        missing = sorted(
            list(
                set(self.system.arguments) - set([*system_parameters_copy] + [*to_fit])
            )
        )
        if len(missing) > 0:
            print(
                "Not all system parameters included in system_parameters or to_fit dictionaries, check all variables for the used equation are included"
            )
            print("Missing variables are: ", missing)
            return None

        # Often, people forget to set the readout, or include ymin as a system
        # parameter which indicates we are dealing with fitting signal data.
        # Check if the changing parameter exceeds the matching ycoord at any
        # point. If so, warn.

        changing_parameter = self._find_changing_parameters(system_parameters)
        assert (
            changing_parameter is not None
        ), "Nothing changes in the data in order to fit... cannot fit."
        if any((system_parameters[changing_parameter[0]] - ycoords) < 0):
            if (
                "ymin" not in system_parameters.keys()
                and "ymax" not in system_parameters.keys()
                and not self.disable_signal_warning
            ):
                print(
                    "Warning: Some x-values are greater than y-values, implying you forgot to include a ymin or ymax to indicate your ycoords are signal, not concentrations, please provide a ymin in the fit parameters to indicate it is a signal"
                )

        # Add parameters for lmfit, accounting for bounds
        if bounds is None:
            bounds = {}
        params = lmfit.Parameters()
        for varname in to_fit.keys():
            # Do not be tempted to set bnd_min to 0 to help the minimizer as
            # lmfit (at least of version 1.0.1) then fails to work... no
            # explanation as to why is present, in lmfit documentation,
            # only that little advantage is gained by narrowing bounds and
            # even outlandishly wide bounds has little effect.
            # See https://lmfit.github.io/lmfit-py/bounds.html
            bnd_min = -np.inf
            bnd_max = np.inf
            if varname in bounds.keys():
                bnd_min = bounds[varname][0]
                bnd_max = bounds[varname][1]
            params.add(varname, value=to_fit[varname], min=bnd_min, max=bnd_max)

        lmmini = lmfit.Minimizer(
            self._residual, params, fcn_args=(system_parameters_copy, to_fit, ycoords)
        )
        result = lmmini.minimize()

        # Check that any fitted KDs are not negative
        for k in system_parameters_copy.keys():
            if k.startswith("kd"):
                assert (
                    system_parameters_copy[k] > 0
                ), f"Error, Fitted KD is negative ({system_parameters_copy[k]}), unable to fit"

        # Return 2 things, first the newly fitted system dictionary, and then
        # a dictionary of absolute fit errors.
        return (
            system_parameters_copy,
            dict((p, result.params[p].stderr) for p in result.params),
        )

    def _residual(self, params, system_parameters: dict, to_fit: dict, y: np.array):
        """
        Residual function for fitting parameters.

        Helper function for lm_fit, calculating residual remaining for probed
        system parameters.

        Parameters
        ----------
        params : dict
            A dictionary of the parameters required to be evaluated to a fit model.
        system_parameters : dict
            Dictionary containing system parameters, will be used as arguments to the systems equations.
        to_fit : dict
            Dictionary containing system parameters to fit.
        y : np.array
            A array-like data containing the system parameters should be fit to
        Returns
        -------
            The cost between the experimental datapoints and the values derived from the model.
        """
        for value in params:
            system_parameters[value] = float(params[value])
        return self.system.query(system_parameters) - y
    
    def get_system_arguments(self):
            if self.system is None:
                return None
            else:
                return self.system.arguments