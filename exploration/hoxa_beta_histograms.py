import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

din = "/projects/p30791/methylation/differential_methylation"
dout = "/projects/p30791/methylation/plots"
meta_fn = "/projects/p30791/methylation/raw_data/meta.csv"

with open(f"{din}/hoxa_probes.json", "r") as file:
    hoxa_probes = json.load(file)

# Get DM probes
with open(f"{din}/probe_set_Non_monotonic_valley_A.txt", "r") as f:
    dm_probes = [x.strip() for x in f.readlines()]

betas = pd.read_csv(
    "/projects/p30791/methylation/sesame_data/betas_processed.csv", index_col=0
)  # nrows=2000 for testing purposes
meta = pd.read_csv(meta_fn)

with open(
    "/projects/p30791/methylation/sesame_data/exclude_IDATs.txt", "r", encoding="utf-8"
) as f:
    exclude_idats = [x.strip() for x in f.readlines()]

meta = meta.iloc[[x not in exclude_idats for x in meta.IDAT]]

# Only keep DM probes
for hoxa_name, probenames in hoxa_probes.items():
    hoxa_probes[hoxa_name] = list(set(probenames).intersection(set(dm_probes)))
    print(f"{hoxa_name}: ", len(hoxa_probes[hoxa_name]))

tissue_types = ["HDB", "CUB", "OQ", "AN", "TU"]
hox_nums = [2, 3, 4, 5, 6, 7]
plot_hoxa = [f"HOXA{i}" for i in hox_nums]
plot_probes = {
    hox_name: hoxa_probes[hox_name][k + 1] for k, hox_name in enumerate(plot_hoxa)
}
color_mapping = {
    "HDB": "slategray",
    "CUB": "cornflowerblue",
    "OQ": "greenyellow",
    "AN": "goldenrod",
    "TU": "orangered",
}


fig, axs = plt.subplots(
    nrows=len(tissue_types),
    ncols=len(plot_hoxa),
    figsize=(2 * len(plot_hoxa), 1.25 * len(tissue_types)),
)

for i, tissue_type in enumerate(tissue_types):
    tissue_idats = meta.iloc[meta["Sample Region"].values == tissue_type].IDAT.values
    for j, (hoxa_name, probe_name) in enumerate(plot_probes.items()):
        axs[i, j].hist(
            betas.loc[probe_name, tissue_idats].values.ravel(),
            bins=np.arange(0, 1.01, 0.05),
            edgecolor="black",
            color=color_mapping[tissue_type],
        )
        axs[i, j].spines["top"].set_visible(False)
        axs[i, j].spines["right"].set_visible(False)
        if j == 0:
            axs[i, j].set_ylabel("# Samples")
        if i == len(tissue_types) - 1:
            axs[i, j].set_xlabel("Beta values")
        elif i == 0:
            axs[i, j].set_title(f"{hoxa_name}\n({probe_name})")

fig.suptitle("Distributions of methylation levels in HOXA genes", fontsize=16)

# Add row labels with genes
fig.text(0.02, 0.84, "HDB", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.02, 0.664, "CUB", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.02, 0.5, "OQ", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.02, 0.33, "AN", rotation=90, va="center", ha="center", fontsize=14)
fig.text(0.02, 0.15, "TU", rotation=90, va="center", ha="center", fontsize=14)

plt.tight_layout()
plt.subplots_adjust(left=0.08)
fig.savefig(f"{dout}/HOXA_betas_histograms.png")
plt.close()
