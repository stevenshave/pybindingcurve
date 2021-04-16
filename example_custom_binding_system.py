"""Simulation example, custom binding system

PyBindingCurve allows the use of custom binding systems derived from a simple
syntax.  This is in the form of a string with reactions separated either on
newlines, commas, or a combination of the two.  Reactions take the form:
r1+r2<->p - denoting reactant1 + reactant2 form p.  PBC will generate and
solve custom systems. Readouts are signified by inclusion of a 
star (*) on a species.  If no star is found, then the first seen product is
used. Some system examples follow:
-> "P+L<->PL" - standard protein-ligand binding
-> "P+L<->PL, P+I<->PI" - competition binding
-> "P+P<->PP" - dimer formation
-> "monomer+monomer<->dimer" - dimer formation
-> "P+L<->PL1, P+L<->PL2, PL1+L<->PL1L2, PL2+L<->PL1L2" - 1:2 site binding
All species are cast to lowercase and parameters such as concentrations should
be passed as lower case - like 'p':10 for a starting protein concentration of
10 µM (depending on standard unit used)
KDs passed to custom systems use underscores to separate species and product
and prefixed with kd. P+L<->PL would require the KD passed as kd_p_l_pl.
Running with incomplete system parameters will prompt for the correct ones.
The methods employed to solve the system will return a result accurate to 10
decimal places.
"""

import numpy as np
import pybindingcurve as pbc

# Define the custom system
custom_system = "P+L<->PL*"

# Make a pbc BindingCurve defined by the custom system string above
my_system = pbc.BindingCurve(custom_system)

# We interrogate the system to see which parameters/arguments are required
# for calculation
print("Required parameters/arguments are: ", my_system.get_system_arguments())

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume uM for all
# concentrations bellow.
system_parameters = {"p": np.linspace(0, 20), "l": 5, "kd_p_l_pl": 1.2}

# We can now add the curve to the plot, name it with an optional name= value.
my_system.add_curve(system_parameters)

# Show the plot
my_system.show_plot()
