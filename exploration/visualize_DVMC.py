import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control

plotdir = "/projects/p30791/methylation/plots/differential_methylation/dvmc"
dvmc = pd.read_csv(
    "/projects/p30791/methylation/sesame_out/differential_methylation/topDVMC.csv",
    index_col=0,
)

if not os.path.exists(plotdir):
    os.makedirs(plotdir)

# Add FDR-adjusted values for DM
dvmc["q(TT)"] = false_discovery_control(dvmc["P(TT)"])

# Remove non-CpG probes
dvmc = dvmc.iloc[[i.startswith("cg") for i in dvmc.index], :]

## Fig. 1a - differentially variable and differentially methylated CpGs

fig, ax = plt.subplots(nrows=1, ncols=2, figsize=(6, 3))

for i, (test_name, title) in enumerate(zip(["BT", "TT"], ["variable", "methylated"])):
    fdr_cts = (dvmc[f"q({test_name})"] < 0.05).sum().item()
    ax[i].hist(dvmc[f"P({test_name})"], color="orange", bins=20, edgecolor="black")
    # ax[i].set_xlim(0,1)
    ax[i].set_title(f"Diff. {title} CpGs")
    ax[i].set_xlabel("P-values")
    ax[i].set_ylabel("Counts")
    ax[i].text(
        0.25,
        0.9,
        f"n(FDR<0.05) = {fdr_cts:,}",
        transform=ax[i].transAxes,
        fontsize=9,
        verticalalignment="top",
        multialignment="center",
    )
    ax[i].spines["top"].set_visible(False)
    ax[i].spines["right"].set_visible(False)

plt.tight_layout()
fig.savefig(f"{plotdir}/exploratory_fig1a.png")
plt.close()

## Fig. 1b - scatterplot of CpGs by DV and DM p-values

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))

num_dvmc = dvmc.iloc[
    (dvmc["P(TT)"].values < 0.05) & (dvmc["q(BT)"].values < 0.05), :
].shape[0]

ax.scatter(
    -1
    * np.log10(
        dvmc.iloc[(dvmc["P(TT)"].values < 0.05) & (dvmc["q(BT)"].values < 0.05), :][
            "P(TT)"
        ]
    ),
    -1
    * np.log10(
        dvmc.iloc[(dvmc["P(TT)"].values < 0.05) & (dvmc["q(BT)"].values < 0.05), :][
            "q(BT)"
        ]
    ),
    alpha=0.2,
    s=4,
    c="red",
)
ax.scatter(
    -1
    * np.log10(
        dvmc.iloc[(dvmc["P(TT)"].values >= 0.05) | (dvmc["q(BT)"].values >= 0.05), :][
            "P(TT)"
        ]
    ),
    -1
    * np.log10(
        dvmc.iloc[(dvmc["P(TT)"].values >= 0.05) | (dvmc["q(BT)"].values >= 0.05), :][
            "q(BT)"
        ]
    ),
    alpha=0.2,
    s=4,
    c="black",
)
ax.text(
    0.35,
    0.65,
    f"{num_dvmc:,} DVMCs",
    fontsize=9,
    transform=ax.transAxes,
    color="red",
)
ax.set_xlabel("-log10[P] (DMC)")
ax.set_ylabel("-log10[q] (DVC)")

plt.tight_layout()
fig.savefig(f"{plotdir}/exploratory_fig1b.png")
plt.close()

## Fig. 1c

sig_dvmc = dvmc.iloc[(dvmc["P(TT)"].values < 0.05) & (dvmc["q(BT)"].values < 0.05), :]

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(3, 3))

categories = ["DV", "DM"]
x = np.arange(len(categories))  # positions
width = 0.35

bars1 = ax.bar(
    x - width / 2,
    [
        sig_dvmc.iloc[sig_dvmc["log[V1/V0]"].values > 1, :].shape[0],
        sig_dvmc.iloc[sig_dvmc["t"].values > 0, :].shape[0],
    ],
    width=width,
    label="Hyper",
    color="dimgray",
    edgecolor="black",
)
bars2 = ax.bar(
    x + width / 2,
    [
        sig_dvmc.iloc[sig_dvmc["log[V1/V0]"].values < 1, :].shape[0],
        sig_dvmc.iloc[sig_dvmc["t"].values < 0, :].shape[0],
    ],
    width=width,
    label="Hypo",
    color="silver",
    edgecolor="black",
)

ax.set_ylabel("Count (DVMC)")
ax.set_xticks(x)
ax.set_xticklabels(categories)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.legend()

plt.tight_layout()
fig.savefig(f"{plotdir}/exploratory_fig1c.png")
plt.close()

## Fig. 1d


## Fig. 1e
