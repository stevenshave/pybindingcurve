# %%
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import sys
from pathlib import Path

def load_heatmap(filename):
    file = Path(filename)
    if file.exists():
        mat = pickle.load(open(filename, 'rb'))
        return mat
    return None

heatmaps_homo = []
heatmaps_hetero = []

[heatmaps_homo.append(load_heatmap(f)) for f in ["heatmaphomo-5.0.pkl"]]
[heatmaps_hetero.append(load_heatmap(f)) for f in ["heatmaphetero-5.0.pkl"]]

print(heatmaps_hetero)

fig, ax = plt.subplots(nrows=len(heatmaps_hetero), ncols=3,
                       figsize=(10, 4), sharey=True, sharex=True)

#colour_map=sns.color_palette("Greys", 70)
colour_map = sns.color_palette("RdBu_r", 256)

for i in range(len(heatmaps_hetero)):
    # ax[i,0].imshow(heatmaps_homo[i])
    sns.heatmap(heatmaps_homo[i], ax=ax[0], vmin=0, vmax=1,
                cmap=colour_map, cbar=True, cbar_kws=dict(pad=0.01))
    sns.heatmap(heatmaps_hetero[i], ax=ax[1], vmin=0, vmax=1,
                cmap=colour_map, cbar=True, cbar_kws=dict(pad=0.01))
    sns.heatmap(heatmaps_homo[i]-heatmaps_hetero[i], ax=ax[2],
                cmap=colour_map, center=0, cbar=True, cbar_kws=dict(pad=0.01))


x_labels = ["-9", "-6", "-3"]
y_labels = x_labels.copy()

for axx in ax:
    plt.xticks(
        np.arange(0, heatmaps_hetero[0].shape[1]+1, heatmaps_hetero[0].shape[1]/2))
    plt.yticks(
        np.arange(0, heatmaps_hetero[0].shape[1]+1, heatmaps_hetero[0].shape[0]/2))
    axx.set_xticks = np.arange(len(x_labels))
    axx.set_yticks = np.arange(len(y_labels))
    axx.set_xticklabels(x_labels)
    axx.set_yticklabels(y_labels)
    axx.patch.set_linewidth('1')
    axx.patch.set_edgecolor('black')

fig.text(0.45, 0.02, r'log$_{10} $K$_DI$')
fig.text(0.008, 0.545, r'log$_{10}$ K$_DDimer}$', rotation=90)
fig.text(0.131, 0.805, "Homodimer")
fig.text(0.444, 0.805, "Heterodimer")
fig.text(0.72, 0.805, "Homodimer-Heterodimer")

plt.suptitle(r"Dimer breaking with inhibitor, 2 $\mathrm{\mu}$M homodimer monomer and" +
             "\n"+r"1+1 $\mathrm{\mu}$M heterodimer monomers, 5 $\mathrm{\mu}$M inhibitor", fontsize=16)
plt.tight_layout(rect=(0.01, 0.025, 1, 0.85), w_pad=-0.3)
plt.show()
