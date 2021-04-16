"""
Generate paper figure 1

Figure 1 in the PyBindingCurve paper consists of 2 plots in one panel,
the first illustrating the dimer formation rates as a function of monomer
concentration, and the other titration of an inhibitor into these preformed
complexes.
"""

import matplotlib as mpl
import numpy as np
import pybindingcurve as pbc
import matplotlib.pyplot as plt
import pickle
import matplotlib.ticker
import sys

# We can choose to work in a common unit, typically nM, or uM, as long as all
# numbers are in the same unit, the result is valid.  We assume nM for all
# concentrations bellow, but divide results by 1000 on plotting to convert to
# a more convenient µM representation.
num_points = 500  # Number of points on the two x-axes
max_inhibitor_concentration = 5000
maximum_monomer_concentration = 1000
dimer_kd = 100  # 100 nM dimer KDs
inhibitor_kd = 10  # 10 nM inhibitor KDs

# A reviewer enquired as to why the x-axes of plots in figure 1 were not
# using a log scale. In testing, we found this obscured the critical crossover
# points discussed in the text.  To enable this log scale on x-axes, set the
# bellow variable to True, the most insightful figure is achieved with 'False'
use_log_xaxis_scale=False

# Calculate formation concentrations
x_axis_formation = np.linspace(0, maximum_monomer_concentration, num=num_points)
x_axis_breaking = np.linspace(0, max_inhibitor_concentration, num=num_points)
homo_y_formation = np.empty((num_points))
hetero_y_formation = np.empty((num_points))

# Make the formation and breaking PBC binding curves that we will later query
pbc_homodimer_formation = pbc.BindingCurve("homodimer formation")
pbc_heterodimer_formation = pbc.BindingCurve(
    "1:1"
)  # Heterodimer formation is just the same as 1:1, or P+L<->PL
pbc_homodimer_breaking = pbc.BindingCurve("homodimer breaking")
pbc_heterodimer_breaking = pbc.BindingCurve("competition")

# Perform the calculations
homodimer_formation_concs = pbc_homodimer_formation.query(
    {"kdpp": dimer_kd, "p": x_axis_formation * 2}
)
homodimer_breaking_concs = pbc_homodimer_breaking.query(
    {"kdpi": inhibitor_kd, "kdpp": dimer_kd, "p": 2000, "i": x_axis_breaking}
)
heterodimer_breaking_concs = pbc_heterodimer_breaking.query(
    {"kdpi": inhibitor_kd, "kdpl": dimer_kd, "p": 1000, "l": 1000, "i": x_axis_breaking}
)
# Because heterodimer formation requires 2 chaning parameters at once, we must
# do it in a loop of single queries. Breaking does not require such treatment
# as only the inhibitor concentration changes.
heterodimer_formation_concs = np.empty((len(x_axis_formation)))
for i, monomer_conc in enumerate(x_axis_formation):
    heterodimer_formation_concs[i] = pbc_heterodimer_formation.query(
        {"kdpl": dimer_kd, "p": monomer_conc, "l": monomer_conc}
    )


fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(10, 5.5), sharey=True)
fig.suptitle("Homo- vs Hetero-dimer", fontsize=18)

ax[0].plot(
    x_axis_formation / 1000,
    homodimer_formation_concs / 1000,
    "k",
    label="Homodimer monomers",
    linestyle="--",
)
ax[0].plot(
    x_axis_formation / 1000,
    heterodimer_formation_concs / 1000,
    "k",
    label="Heterodimer monomers",
)
if not use_log_xaxis_scale:
    ax[0].set_xlim(0, maximum_monomer_concentration / 1000)
ax[0].set_ylim(0, 1)
ax[0].legend()
ax[0].set_xlabel(r"[Monomers] ($\mathrm{\mu}$M)", fontsize=14)
ax[0].set_ylabel(r"[Dimer] ($\mathrm{\mu}$M)", fontsize=14)
ax[0].set_title(
    r"""Dimer formation,
Dimer K$\mathrm{_D}$s = 100 nM""",
    fontsize=14,
)
ax[0].grid()

ax[1].plot(
    x_axis_breaking / 1000,
    homodimer_breaking_concs / 1000,
    "k",
    label=r"2 $\mathrm{\mu}$M homodimer monomer",
    linestyle="--",
)
ax[1].plot(
    x_axis_breaking / 1000,
    heterodimer_breaking_concs / 1000,
    "k",
    label=r"1 $\mathrm{\mu}$M heterodimer monomers",
)
if not use_log_xaxis_scale:
    ax[1].set_xlim(0, max_inhibitor_concentration / 1000)
ax[1].set_ylim(0, 1)
ax[1].legend()
ax[1].set_xlabel(r"[I$_0$] ($\mathrm{\mu}$M)", fontsize=14)
# ax[1].set_ylabel(r"[Dimer] ($\mathrm{\mu}$M)", fontsize=14)
ax[1].set_title(
    r"""Dimer breaking with inhibitor,
Dimer K$\mathrm{_D}$s = 100 nM, inhibitor K$\mathrm{_D}$=10 nM""",
    fontsize=14,
)
if use_log_xaxis_scale:
    ax[0].set_xscale("log")
    ax[1].set_xscale("log")
ax[1].grid()
plt.tight_layout(rect=(0, 0, 1, 0.9425), w_pad=-0.4)

for axis in ax:
    axis.tick_params(axis="both", which="major", labelsize=14)
    axis.tick_params(axis="both", which="minor", labelsize=14)
fig.text(0.05, 0.87, "A)", fontsize=20)
fig.text(0.514, 0.87, "B)", fontsize=20)
plt.savefig("Figure1.svg")
plt.savefig("Figure1.png", dpi=300)
plt.subplots_adjust(top=0.829, wspace=0.08)
plt.show()
