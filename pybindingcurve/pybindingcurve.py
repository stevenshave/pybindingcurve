"""PyBindingCurve - a package to simulate protein-ligand binding systems

TODO:
* Add curve tracing to analytical systems
* Check combinations of readout with/without ymin/ymax


"""


import numpy as np
import matplotlib.pyplot as plt
import lmfit
from enum import Enum, auto
from pybindingcurve.systems import *

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
        """Curve constructor

        Called by PyBindingCurve

        """
        self.xcoords = xcoords
        self.ycoords = ycoords
        self.name = series_name


class BindingCurve:
    """
    BindingCurve class, used to simulate systems

    BindingCurve objects are governed by their underlying system, defining the
    (usually) protein-ligand binding system being represented.  It also
    provides the main interface for simulation, visualisation, querying and
    the fitting of system parameters.

    Parameters
    ----------
    binding_system : BindingSystem or str
        Define the binding system which will govern this BindingCurve object.
        Caan either be a BindingSystem object or a human readable string
        shortcut, such as '1:1' or 'competition', etc.

    """


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
    _last_known_changing_parameter = "X"

    def query(self, parameters, readout: Readout = None):
        """
        Query a binding system

        Get the readout from from a set of system parameters

        Parameters
        ----------
        parameters : dict
            System parameters defining the system being queried.  Will usually
            contain protein, ligand etc concentrations, and KDs
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
        if readout is None:
            return self.system.query(parameters)
        else:
            return readout(parameters, self.system.query(parameters))[1]

    def _find_changing_parameters(self, params: dict):
        changing_list = []
        for p in params.keys():
            if isinstance(params[p], np.ndarray) or isinstance(params[p], list):
                changing_list.append(p)
        if len(changing_list) == 0:
            return None
        else:
            return changing_list

    def __init__(self, binding_system: [str, BindingSystem]):
        if isinstance(binding_system, str):
            binding_system = binding_system.lower()
            # 1:1
            if binding_system in ["simple", "1:1"]:
                self.system = System_analytical_one_to_one_pl()
            # 1:1 kinetic
            if binding_system in ["simplekinetic", "simple kinetic", "1:1kinetic", "1:1 kinetic"]:
                self.system = System_kinetic_one_to_one_pl()
            # Homodimer formation
            if binding_system in ["homodimerformation", "homodimer formation"]:
                self.system = System_analytical_homodimerformation_pp()
            # Homodimer formation kinetic - only used for testing purposes
            if binding_system in [
                "homodimerformationkinetic",
                "homodimer formation kinetic",
            ]:
                self.system = System_kinetic_homodimerformation()

            # Competition
            if binding_system in ["competition", "1:1:1"]:
                self.system = System_analytical_competition_pl()
            # Competition
            if binding_system in ["homodimerbreaking", "homodimer breaking"]:
                self.system = System_kinetic_homodimerbreaking_pp()
        else:
            if issubclass(binding_system, BindingSystem):
                self.system = binding_system()
            else:
                print(
                    "Invalid system specified, try one of: [simple, homodimer, competition, homdimer breaking], or pass a system object"
                )
                return None

    def _initialize_plot(self):
        if self.fig is None:
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
        if not len(changing_parameters) == 1:
            print("Must have 1 changing parameter, no curves added.")
            return

        y_values = self.system.query(parameters)

        if readout is not None:
            self._last_custom_readout, y_values = readout(parameters, y_values)

        if y_values.ndim > 1:
            for i in range(y_values.ndim):
                self.curves.append(
                    _Curve(parameters[changing_parameters[0]], y_values[i])
                )
        else:
            self.curves.append(
                _Curve(parameters[changing_parameters[0]], y_values))
        self._last_known_changing_parameter = changing_parameters[0]
        for curve_it, curve in enumerate(self.curves[self._num_added_traces:]):
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
                parameters[changing_parameters[0]],
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
            self._min_y_axis = np.nanmin(
                [self._min_y_axis, np.nanmin(curve.ycoords)])
            self._max_y_axis = np.nanmax(
                [self._max_y_axis, np.nanmax(curve.ycoords)])

    def add_scatter(self, xcoords, ycoords):
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

        """
        self._initialize_plot()
        self.axes.scatter(xcoords, ycoords)
        if isinstance(xcoords, np.ndarray) and isinstance(ycoords, np.ndarray):
            self._min_y_axis = min(self._min_y_axis, min(np.real(ycoords)))
            self._max_y_axis = max(self._max_y_axis, max(np.real(ycoords)))

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
            plt.savefig(svg_filename, metadata={
                        "Title": "pyBindingCurve plot"})
        plt.show()

    def fit(
        self,
        system_parameters: dict,
        to_fit: dict,
        ycoords: np.array,
        bounds: dict = None,
    ):
        """Fit the parameters of a system to a set of data points

        Fit the system to a set of (usually) experimental datapoints.
        The fitted parameters are stored in the system_parameters dict
        which may be accessed after running this function.  It is
        possible to fit multiple parameters at once and define bounds
        for the parameters.  The function returns a dictionary of the
        accuracy of fitted parameters, which may be captured, or not.

        Parameters:
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
        system_parameters_copy = dict(system_parameters)
        # Check we have parameters to fit, and nothing is missing
        if len(to_fit.keys()) == 0:
            print("Nothing to fit, insert parameters to fit into to_fit dictionary")
            return None
        missing = sorted(
            list(
                set(self.system.arguments) -
                set([*system_parameters_copy] + [*to_fit])
            )
        )
        if len(missing) > 0:
            print(
                "Not all system parameters included in system_parameters or to_fit dictionaries, check all variables for the used equation are included"
            )
            print("Missing variables are: ", missing)
            return None
        # Add parameters for lmfit, accounting for bounds
        if bounds is None:
            bounds = {}
        params = lmfit.Parameters()
        for varname in to_fit.keys():
            bnd_min = -np.inf
            bnd_max = np.inf
            if varname in bounds.keys():
                bnd_min = bounds[varname][0]
                bnd_max = bounds[varname][1]
            params.add(
                varname, value=to_fit[varname], min=bnd_min, max=bnd_max)

        lmmini = lmfit.Minimizer(
            self._residual, params, fcn_args=(
                system_parameters_copy, to_fit, ycoords)
        )
        result = lmmini.minimize()

        for k in system_parameters_copy.keys():
            if isinstance(system_parameters_copy[k], lmfit.parameter.Parameter):
                system_parameters_copy[k] = system_parameters_copy[k].value

        return (
            system_parameters_copy,
            dict((p, result.params[p].stderr) for p in result.params),
        )

    def _residual(self, params, system_parameters: dict, to_fit: dict, y: np.array):
        """Residual function for fitting parameters.

        Helper function for lm_fit, calculating residual remaining for probed system parameters.

        Args:

        system_parameters (dict): Dictionary containing system
                parameters, will be used as arguments to the systems equations.
        to_fit: (dict): Dictionary containing system parameters to fit.

        """
        for value in params:
            system_parameters[value] = float(params[value])
        return self.system.query(system_parameters) - y
