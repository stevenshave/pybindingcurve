# PyBindingCurve

PyBindingCurve is a Python package for simulation, plotting and fitting of experimental parameters to protein-ligand binding systems at equilibrium.  In simple terms, the most basic functionality allows simulation of a two species binding to each other as a function of their concentrations and the dissociation constant (K<sub>D</sub>) between the two species.

![PyBindingCurve simulation](https://raw.githubusercontent.com/stevenshave/pybindingcurve/master/pybindingcurve_logo.png "Breaking a dimer")

# Installation
PyBindingCurve may be installed from source present in the GitHub repository https://github.com/stevenshave/pybindingcurve via git pull, or from the Python Package Index (https://pypi.org/project/pybindingcurve/) using the command :
> pip install pybindingcurve

# Requirements
PyBindingCurve was developed using python 3.7.1 but should work with any Python version 3.6 or greater. The following packages are also required
- Matplotlib (2.x)
- Numpy (1.15.x)
- lm_fit (1.0.0)
- mpmath (1.1.0)

# Licence
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

# Authors
PyBindingCurve was written by Steven Shave 
![email](https://raw.githubusercontent.com/stevenshave/pybindingcurve/master/email-address-image.gif)


Please get in contact for custom solutions, integration to existing workflows and training.
