# PyBindingCurve

*Shave, Steven, et al. "PyBindingCurve, simulation, and curve fitting to complex binding systems at equilibrium." Journal of Chemical Information and Modeling (2021).* https://doi.org/10.1021/acs.jcim.1c00216

PyBindingCurve is a Python package for simulation, plotting and fitting of experimental parameters to protein-ligand binding systems at equilibrium.  In simple terms, the most basic functionality allows simulation of a two species binding to each other as a function of their concentrations and the dissociation constant (K<sub>D</sub>) between the two species.  A number of systems are built in and can be solved using direct analytical, kinetic, or Langrange multiplier based techniques.  User-defined custom systems can also be specified using a simple syntax.

Try without installing on Google colab! https://colab.research.google.com/drive/1upxm56mGYWo8jvTTJjZLOEq6DT0lRy8d


![PyBindingCurve simulation](https://raw.githubusercontent.com/stevenshave/pybindingcurve/master/pybindingcurve_logo.png "Breaking a dimer")

# Installation
PyBindingCurve may be installed from source present in the GitHub repository https://github.com/stevenshave/pybindingcurve via git pull, or from the Python Package Index (https://pypi.org/project/pybindingcurve/) using the command :
> pip install pybindingcurve

# Requirements
PyBindingCurve requires Python 3.9 or later. The following packages are also required
- numpy>=1.26
- matplotlib>=3.8
- lmfit>=1.2.2
- mpmath>=1.3.0
- autograd>=1.6.2

# License
[MIT License](https://github.com/stevenshave/pybindingcurve/blob/master/LICENSE)



# Usage
A tutorial and API documentation can be found [here](https://stevenshave.github.io/pybindingcurve/)

A quickstart example for simulation of protein-ligand binding is as follows:

```
import numpy as np
import pybindingcurve as pbc
my_system = pbc.BindingCurve("1:1")
system_parameters = {"p": np.linspace(0, 20), "l": 10, "kdpl": 1}
my_system.add_curve(system_parameters)
my_system.show_plot()
```
Tests written using the pytest framework may be run with 'pytest' (ensure pytest is installed in your python environment, or pip install it)

# Authors
PyBindingCurve was written by Steven Shave 
![email](https://raw.githubusercontent.com/stevenshave/pybindingcurve/master/email-address-image.gif)


Please get in contact for custom solutions, integration to existing workflows and training.
