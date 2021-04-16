"""Generate paper figure 2

Figure 2 in the PyBindingCurve paper consists of heatmaps showing
homo- vs hetero-dimer breaking. The figure in BioRxiv preprint #1
used numerically unstable code - before the inclusion of MPMath 
for arbitary precision arithmetic which fixes artefacts present in
the old figure. First the code performs the long calculations,
storing results in pickles, before rendering. Subsequent runs
finding the pickles speeds up display.
"""


from multiprocessing import Pool
import numpy as np
import pybindingcurve as pbc
import pickle
from pathlib import Path
import seaborn as sns
import matplotlib.pyplot as plt


def get_2D_grid_values(
    xmin,
    xmax,
    ymin,
    ymax,
    system,
    parameters,
    x_parameter,
    y_parameter,
    filename,
    plot_steps,
    starting_y_guess=0,
):
    file = Path(filename)
    if file.exists():
        mat = pickle.load(open(filename, "rb"))
        return mat
    x_logconc = np.linspace(xmin, xmax, plot_steps)
    y_logconc = np.linspace(ymin, ymax, plot_steps)
    mat = np.ndarray((plot_steps, plot_steps))
    for ix, x in enumerate(x_logconc):
        print("Working on :", x)
        for iy, y in enumerate(y_logconc):
            parameters[x_parameter] = 10 ** x
            parameters[y_parameter] = 10 ** y
            mat[ix, iy] = system.query(parameters)

    pickle.dump(mat, open(filename, "wb"))
    return mat


def generate_heatmaps(
    homodimer_map_file: Path, heterodimer_map_file: Path, plot_steps: int = 800
):
    pool = Pool(2)
    inhibitor_conc = 5.0
    pbc_homodimer_breaking = pbc.BindingCurve("homodimer breaking")
    pbc_heterodimer_breaking = pbc.BindingCurve("competition")

    print("Building homdimer heatmap file")
    j1 = pool.apply_async(
        get_2D_grid_values,
        [
            -3,
            3,
            -3,
            3,
            pbc_homodimer_breaking,
            {"p": 2, "i": inhibitor_conc},
            "kdpp",
            "kdpi",
            f"heatmaphomo-{str(inhibitor_conc)}.pkl",
            plot_steps,
        ],
    )

    print("Building heterodimer heatmap file")

    j2 = pool.apply_async(
        get_2D_grid_values,
        [
            -3,
            3,
            -3,
            3,
            pbc_heterodimer_breaking,
            {"p": 1, "l": 1, "i": inhibitor_conc},
            "kdpl",
            "kdpi",
            f"heatmaphetero-{str(inhibitor_conc)}.pkl",
            plot_steps,
        ],
    )
    res = j1.get()
    res = j2.get()


def load_heatmap(filepath: Path):
    if filepath.exists():
        mat = pickle.load(open(filepath, "rb"))
        print(filepath, "shape=", mat.shape)
        # mat = mat[::2, ::2]  # To downsample, if required for saved matrixq
        return mat
    else:
        return None


def make_plot(homodimer_map_file: Path, heterodimer_map_file: Path):
    heatmap_homo = load_heatmap(homodimer_map_file)
    heatmap_hetero = load_heatmap(heterodimer_map_file)

    fig, ax = plt.subplots(ncols=3, figsize=(10, 4), sharey=True, sharex=True)

    colour_map = sns.color_palette("RdBu_r", 256)

    sns.heatmap(
        heatmap_homo[:, ::-1],
        ax=ax[0],
        vmin=0,
        vmax=1,
        cmap=colour_map,
        cbar=True,
        cbar_kws=dict(pad=0.01),
    )
    sns.heatmap(
        heatmap_hetero[:, ::-1],
        ax=ax[1],
        vmin=0,
        vmax=1,
        cmap=colour_map,
        cbar=True,
        cbar_kws=dict(pad=0.01),
    )
    sns.heatmap(
        heatmap_homo[:, ::-1] - heatmap_hetero[:, ::-1],
        ax=ax[2],
        cmap=sns.color_palette("PRGn", 1024),
        center=0,
        cbar=True,
        cbar_kws=dict(pad=0.01),
    )

    x_labels = ["9 (nM)", "6 ($\mathrm{\mu}$M)", "3 (mM)"]
    x_labels.reverse()
    y_labels = x_labels.copy()

    y_labels.reverse()
    y_labels[2] = " 3 (mM)"
    for axx in ax:
        plt.xticks(
            np.arange(0, heatmap_hetero.shape[1] + 1, heatmap_hetero.shape[1] / 2),
            rotation=0,
        )
        plt.yticks(
            np.arange(0, heatmap_hetero.shape[1] + 1, heatmap_hetero.shape[0] / 2)
        )
        axx.set_xticks = np.arange(len(x_labels))
        axx.set_yticks = np.arange(len(y_labels))
        axx.set_xticklabels(x_labels, rotation=0)
        axx.set_yticklabels(y_labels, rotation=0)
        axx.patch.set_linewidth("1")
        axx.patch.set_edgecolor("black")

    fig.text(0.45, 0.02, r"p$\mathrm{K_DI}$")
    fig.text(0.008, 0.445, r"p$\mathrm{K_DDimer}$", rotation=90)
    fig.text(0.165, 0.805, "Homodimer")
    fig.text(0.465, 0.805, "Heterodimer")
    fig.text(0.68, 0.805, "Difference (Homodimer-Heterodimer)")

    plt.suptitle(
        r"Dimer breaking with inhibitor, 2 $\mathrm{\mu}$M homodimer monomer and"
        + "\n"
        + r"1+1 $\mathrm{\mu}$M heterodimer monomers, 5 $\mathrm{\mu}$M inhibitor",
        fontsize=16,
    )
    plt.subplots_adjust(left=0.086, bottom=0.11, right=0.976, top=0.793, wspace=0.086)

    plt.savefig("Figure2.svg")
    plt.savefig("Figure2.png", dpi=300)

    plt.show()


if __name__ == "__main__":
    # Be careful with the size of the matrices determined by plot_steps,
    # 800x800 matrices result in a ~360 MB SVG which crashes inkscape on
    # a 16GB machine. Downsampling by a factor of 2 in each dimension 
    # produces a ~90 MB SVG which inkacape still struggles with.
    # Paper figure generated with plot_steps=400. The output PNG is
    # certainly easier to deal with than the SVG
    
    plot_steps = 400
    homodimer_map_file = Path(f"heatmaphomo-5.0.pkl")
    heterodimer_map_file = Path(f"heatmaphetero-5.0.pkl")
    if not homodimer_map_file.exists() or not heterodimer_map_file.exists():
        generate_heatmaps(
            homodimer_map_file, heterodimer_map_file, plot_steps=plot_steps
        )
    make_plot(homodimer_map_file, heterodimer_map_file)
