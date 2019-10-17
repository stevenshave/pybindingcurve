import matplotlib as mpl
import numpy as np
import pybindingcurve2 as pbc2
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import matplotlib.ticker
import sys

num_points=500
max_x_breaking=5
max_x_formation=1
dimer_kd=0.1
inhibitor_kd=10e-3


# Calculate formation concentrations
x_axis_formation=np.linspace(0,max_x_formation, num=num_points)
x_axis_breaking=np.linspace(0, max_x_breaking, num=num_points)
homo_y_formation=np.empty((num_points))
hetero_y_formation=np.empty((num_points))
for i in range(len(homo_y_formation)):
    homo_y_formation[i]=pbc2.system03_p_kdpp__pp({'kdpp':dimer_kd, 'p':x_axis_formation[i]*2})
    hetero_y_formation[i]=pbc2.system01_p_l_kd__pl({'kdpl':dimer_kd,'p':x_axis_formation[i], 'l':x_axis_formation[i]})

homo_y_breaking=np.empty((num_points))
hetero_y_breaking=np.empty((num_points))
for i in range(len(homo_y_formation)):
    homo_y_breaking[i]=pbc2.system04_p_i_kdpp_kdpi__pp({'kdpi':inhibitor_kd, 'kdpp':dimer_kd, 'p':2, 'i':x_axis_breaking[i]})
    hetero_y_breaking[i]=pbc2.system02_p_l_i_kdpl_kdpi__pl({'kdpi': inhibitor_kd, 'kdpl':dimer_kd, 'p':1,'l':1, 'i':x_axis_breaking[i]})

fig, ax=plt.subplots(nrows=1,ncols=2, figsize=(10,5.5), sharey=True)
#plt.tight_layout()
fig.suptitle("Homo- vs Hetero-dimer", fontsize=18)

ax[0].plot(x_axis_formation,homo_y_formation ,'k', label="Homodimer monomers", linestyle='--')
ax[0].plot(x_axis_formation,hetero_y_formation ,'k', label="Heterodimer monomers")
ax[0].set_xlim(0, max_x_formation)
ax[0].set_ylim(0,1)
ax[0].legend()
ax[0].set_xlabel(r"[Monomers] ($\mathrm{\mu}$M)", fontsize=14)
ax[0].set_ylabel(r"[Dimer] ($\mathrm{\mu}$M)", fontsize=14)
ax[0].set_title(r"""Dimer formation,
Dimer K$\mathrm{_D}$s = 100 nM""",fontsize=14)
ax[0].grid()

ax[1].plot(x_axis_breaking,homo_y_breaking ,'k', label=r"2 $\mathrm{\mu}$M homodimer monomer", linestyle='--')
ax[1].plot(x_axis_breaking,hetero_y_breaking ,'k', label=r"1 $\mathrm{\mu}$M heterodimer monomers")
ax[1].set_xlim(0, max_x_breaking)
ax[1].set_ylim(0,1)
ax[1].legend()
ax[1].set_xlabel(r"[I$_0$] ($\mathrm{\mu}$M)", fontsize=14)
#ax[1].set_ylabel(r"[Dimer] ($\mathrm{\mu}$M)", fontsize=14)
ax[1].set_title(r"""Dimer breaking with inhibitor,
Dimer K$\mathrm{_D}$s = 100 nM, inhibitor K$\mathrm{_D}$=10 nM""",fontsize=14)
ax[1].grid()
plt.tight_layout(rect=(0,0,1,0.9425), w_pad=-0.4)

for axis in ax:
    axis.tick_params(axis='both', which='major', labelsize=14)
    axis.tick_params(axis='both', which='minor', labelsize=14)
fig.text(0.05, 0.87, "A)", fontsize=20)
fig.text(0.514, 0.87, "B)", fontsize=20)
fig.savefig("/home/stevens/Downloads/fig1.svg", format='svg')
plt.show()
