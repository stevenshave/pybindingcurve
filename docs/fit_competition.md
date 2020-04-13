# Fitting to 1:1:1 competition
![1:1 binding system](./images/Fig_system_competition.png "1:1 binding system")
[Return to tutorials](tutorial.md)

 
Using experimental competition data, we may obtain system parameters, such as inhibitor KD. Example code is available here:
[https://github.com/stevenshave/pybindingcurve/blob/master/example_competition_fit.py](https://github.com/stevenshave/pybindingcurve/blob/master/example_competition_fit.py)

Perform the standard imports:
```python
import numpy as np
import pybindingcurve as pbc
```
We can choose to work in a common unit, typically nM, or µM, as long as all numbers are in the same unit, the result is valid.  We assume µM for all concentrations bellow.

Using the following experimental data:
|[P] µM | 0|4.2|8.4|16.8|21.1|31.6|35.8|40.0|
|---|---|---|---|---|---|---|---|---|
Signal|150|330|1050|3080|4300|6330|6490|6960|

Define x and y coordinates from experimental data:
```python
xcoords = np.array([0.0, 4.2, 8.4, 16.8, 21.1, 31.6, 35.8, 40.0])
ycoords = np.array([150, 330, 1050, 3080, 4300, 6330, 6490, 6960])
```

Construct the PyBindingCurve object, operating on a 1:1:1 (compeittion) system and add experimental data to the plot:

```python
mySystem = pbc.BindingCurve("1:1:1")
mySystem.add_scatter(xcoords, ycoords)
```

Known system parameters, kdpl will be added to this by fitting:

```python
system_parameters = {"p": xcoords, "l": 10, "i": 10, "kdpl": 10}
```

Now we call fit, passing the known parameters, followed by a dict of parameters to be fitted along with an initial guess, pass the ycoords, and what the readout (ycoords) is:

```python
fitted_system, fit_accuracy = mySystem.fit(system_parameters, {"kdpi": 0}, ycoords)
```
Print out the fitted parameters:

```python
for k, v in fit_accuracy.items():
    print(f"Fit: {k}={fitted_system[k]} +/- {v}")
```

Assign more points to 'p' to make a smooth plot:

```python
fitted_system["p"] = np.linspace(0, np.max(xcoords))
```

Add a new curve, simulated using fitted parameters to our BindingCurve object and show the plot:

```python
mySystem.add_curve(fitted_system)
mySystem.show_plot()
```

Which results in the following output and plot:
> Fit: kdpi=0.44680894202996824 +/- 0.10384753604598472
>
> Fit: ymax=9920.875421523158 +/- 98.92963212643627

![Fitting competition data](./images/Fig_competition_fit.svg "Fitting competition data")


[Return to tutorials](tutorial.md)
