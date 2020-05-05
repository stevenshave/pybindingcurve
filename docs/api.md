# PyBindingCurve API reference

Documentation and API
Full PyBindingCurve source code can be found here: https://github.com/stevenshave/pybindingcurve

- [Overview](#Overview)
- [pbc.BindingCurve](#pbc.BindingCurve)
    - [Initialisation](###Initialisation)
    - [add_curve](###add_curve)
    - [query](###query)
    - [fit](###fit)
    - [add_scatter](###add_scatter)
    - [show_plot](###show_plot)
- [pbc.systems and shortcut strings](##pbc.systems)
- [pbc.BindingSystem](##pbc.BindingSystem)
- [pbc.Readout](##pbc.Readout)


## Overview
Conventionally, the standard import utilised in a run of PyBindingCurve (PBC) are defined as follows:
import pybindingcurve as pbc
import numpy as np
PyBindingCurve is imported with the short name ‘pbc’, and then NumPy as ‘np’ to enable easy specification of ranged system parameters, evenly spaced across intervals mimicking titrations. 
Next, we initialise a PBC BindingCurve object.  Upon initialisation, either a string or BindingSystem object is used to define the system to be simulated.  Simple strings such as “1:1”, “competition”, “homodimer breaking” can be used as an easy way to define the type of binding system that should be mathematically modelled.  See the list of available system shortcut strings bellow in the ‘pbc.systems and shortcut strings’ section bellow.  Custom objects may also be created of type pybindingcurve.BindingSystem, of which there exist a large choice within pybindingcurve.systems, or the user may create a custom BindingSystem to initialise PBC objects:
my_system=pbc.BindingCurve(“1:1”)
There are three main modes of operation within PyBindingCurve; 1) Visualisation of protein-ligand behaviour within a titration, simulating a range of conditions. 2) Simulation of a single system state with discrete parameters. 3) Fitting of experimental data.  A description of these follows:
1.	Simulation with visualisation: pass a dictionary containing system parameters to the add_curve function of the BindingCurve object (my_system). Required system parameters depend on the system being modelled.  In the case of 1:1 binding, we require p (protein concentration), l (ligand concentration), kdpl (dissociation constant between p and l), and optionally, a ymax and/or ymin variable if dealing with simulation of a signal.  add_curve expects one changing parameter, which will be the x-axis. By default, complex concentration will be the readout of a system, but that can be changed by passing different readout options to add_curve.  See the ‘pbc.Readout’ section bellow for further information. With one curve added, we can add more curves, or simply display the plot by calling the show_plot function.

```python
    system_parameters={‘p’:np.linspace(0,20), ‘l’:20, ‘kdpl’:10)
    my_system.add_curve(system_parameters)
    my_system.show_plot()
```

2.	Single point simulation: if you require not a simulation with a curve, but a single point with set concentrations and KDs, then query may be called with a dictionary of system parameters and data returned.  Additionally, if the dictionary contains a NumPy array, representing a titration (for example, the system parameters above), then a NumPy array of the readout is returned instead of a single value:

```python
    complex_conc = my_system.query({‘p’:10, ‘l’:20, ‘kdpl’:10})
```

3.	Fit experimental data to a system to obtain experimental parameters, such as KD. A common situation is determining KD from measurements obtained from experimental data.  We can perform this as follows, with x- and y-coordinates, we add the experimental points to the plot, define system parameters that we do know (protein concentrations, and the amount of ligand), and then call fit on the system passing in parameters to fit and an initial guess (1), along with the known parameters.  We then iterate and print the fitted parameters.

```python
    xcoords = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0])
    ycoords = np.array([0.544, 4.832, 6.367, 7.093, 7.987, 9.005, 9.079, 8.906, 9.010, 10.046, 9.225])
    my_system.add_scatter(xcoords, ycoords)
    system_parameters = {"p": xcoords, "l": 10}
    fitted_system, fit_accuracy = my_system.fit(system_parameters, {"kdpl": 1}, ycoords)
    for k, v in fit_accuracy.items():
        print(f"Fit: {k}={fitted_system[k]} +/- {v}")
```

## pbc.BindingCurve
The BindingCurve object allows the user to work with a specific system, supplying tools for simulation and visualisation (plotting), querying of single point values, and fitting of experimental parameters to observation data.
### Initialisation
When initialising this main class of PBC, we may supply either a pbc.BindingSystem, or a human readable shortcut string such as “1:1”, “competition”, “homodimer formation”, etc. For a full list of systems and shortcuts, please refer to the ‘pbc.systems and shortcut strings’ section.
Initialisation of a BindingCurve object takes the following arguments:
    
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

Once intitialised with a pbc.BindingSystem, we may perform the following utilising its member functions.
### add_curve
The add curve function is the main way of simulating a binding curve with PBC.  Once a BindingCurve object is initialised, 
       
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
### query
When simulation with visualisation (plotting) is not required, we can use the query function to interrogate a system, returning either singular values, or arrays of values if one of the input parameters is an array or list.

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
### fit
With a system defined, we may fit experimental data to the system.

        """
        Fit the parameters of a system to a set of data points

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
### add_scatter
Experimental data can be added plots with the add_scatter command, taking a simple list of x and y coordinates

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

### show_plot
With curves, scatterpoints and fits applied, we may display the plot.

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
        

## pbc.systems

pbc.systems contains all default systems supplied with PBC, and exports them to the PBC namespace. Systems may be passed as arguments to pbc.BindingCurve objects upon initialisation to define the underlying system governing simulation, queries, and fitting. Additionally, the following shortcut strings may be used as shortcuts:

|Shortcut string list|pbc.systems equivalent|
|---|---|
|simple, 1:1|System_analytical_one_to_one__pl|
|simplelagrange, simple lagrange, 1:1lagrange, 1:1 lagrange|System_lagrange_one_to_one__pl|
|simplekinetic, simple kinetic, 1:1kinetic, 1:1 kinetic|	System_kinetic_one_to_one__pl|
|homodimerformation, homodimer formation|	System_analytical_homodimerformation__pp|
|homodimerformationlagrange, homodimer formation lagrange|System_lagrange_homodimerformation__pp|
|homodimerformationkinetic, homodimer formation kinetic|System_kinetic_homodimerformation__pp|
|competition, 1:1:1|System_analytical_competition__pl|
|competition lagrange, competitionlagrange|System_lagrange_competition__pl|
|homodimerbreaking, homodimer breaking, homodimerbreakinglagrange, homodimer breaking lagrange|System_lagrange_homodimerbreaking__pp|
|homodimerbreakingkinetic, homodimer breaking kinetic|System_kinetic_homodimerbreaking__pp|
|homodimerbreakinganalytical, homodimer breaking analytical|System_analytical_homodimerbreaking__pp|
|1:2, 1:2 lagrange| System_lagrange_1_to_2__pl12|
|1:3, 1:3 lagrange| System_lagrange_1_to_3__pl123|
|1:4, 1:4 lagrange| System_lagrange_1_to_4__pl1234|
|1:5, 1:5 lagrange| System_lagrange_1_to_5__pl12345|


## pbc.BindingSystem
Custom binding systems may be defined through inheritance from the base class pbc.BindingSystem.  This provides basic functionality through a standard interface to PBC, allowing simulation, querying and fitting.  It expects the child class to provide a constructor which passes a function for querying the system and a query method.  An example pbc.BindingSystem for 1:1 binding solved analytically is defined as follows:

```
class System_analytical_one_to_one_pl(BindingSystem):
    def __init__(self):
        super().__init__(
            analyticalsystems.system01_one_to_one__p_l_kd__pl, analytical=True
        )
        self.default_readout = "pl"

    def query(self, parameters: dict):
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
```
Here, we see the parent class constructor called upon initialisation of the object with two arguments, the first is a python function which calculates the complex concentration present in a 1:1 binding system, which itself takes the appropriate parameters to calculate this.  In addition, a flag is set to define when the solution is solved analytically.  The query method examines the content of the system and deals with the presence of ymin and ymax to denote a signal is being simulated.  Query should ultimately end up calling query on the parent class, which has been set to return the result of the previously assigned function in the constructor.

## pbc.Readout
The pbc.Readout class contains three static methods, not requiring object initialisation for use. These methods all take in a system parameters dictionary describing the system, and the y_values resulting from system query calls (either through simulation of querying for singular values). These readout functions offer a convenient way to transform results. For example, the readout function to transform complex concentration into fraction ligand bound is defined as follows:
```
def fraction_l(system_parameters: dict, y):
        """ Readout as fraction ligand bound """
        return "Fraction l bound", y / system_parameters["l"]
```
This returns a tuple, with the first value being used in labelling of the plot y-axis, and the second the y-values to be plotted; in this case, the original y values divided by the overall starting ligand concentration.  Similar functions can be defined and used interchangeably with those found in pbc.Readout.
