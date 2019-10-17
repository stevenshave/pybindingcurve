#%%
from multiprocessing import Pool
import numpy as np
import pybindingcurve as pbc
import pickle
from pathlib import Path
plot_steps=10

def get_2D_grid_values(xmin, xmax, ymin, ymax, system, parameters, x_parameter, y_parameter, filename, plot_steps, starting_y_guess=0):
    file=Path(filename)
    if file.exists():
       mat=pickle.load(open(filename, 'rb'))
       return mat
    x_logconc=np.linspace(xmin, xmax, plot_steps)
    y_logconc=np.linspace(ymin, ymax, plot_steps)
    mat=np.ndarray((plot_steps, plot_steps))
    for ix, x in enumerate(x_logconc):
        print("Working on :", x)
        for iy, y in enumerate(y_logconc):
            parameters[x_parameter]=10**x
            parameters[y_parameter]=10**y
            if "homo" in filename:
                mat[ix,iy]=system.query(parameters, "pp")
            else:
                mat[ix,iy]=system.query(parameters, "pl")

    pickle.dump(mat, open(filename, "wb"))
    return mat
pool = Pool(2)

for i_amount in [5.0]:
    j1=pool.apply_async(get_2D_grid_values, [-3, 3, -3, 3, pbc.System_homodimer_breaking(), {'p':2, 'i':i_amount}, 'kdpp', 'kdpi', f"heatmaphomo-{str(i_amount)}.pkl", plot_steps])
    j2=pool.apply_async(get_2D_grid_values, [-3, 3, -3, 3, pbc.System_competition(), {'p':1,'l':1, 'i':i_amount}, 'kdpl', 'kdpi', f"heatmaphetero-{str(i_amount)}.pkl", plot_steps])
    res=j1.get()
    res=j2.get()

