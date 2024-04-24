import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

###################################
#### Load DV results and betas ####
###################################

din = "/projects/p30791/methylation"
plot_dir = "/projects/p30791/methylation/plots/differential_methylation"

ref_comp_pairs_tpx = [  # pairs along the tumor proximity axis
    ("CFN", "CUB"),
    # ("CUB", "OQ"),
    ("OQ", "AN"),
    ("AN", "TU"),
]

color_mapping = {
    "All": "darkslategray",
    "CFN": "slategray",
    "CUB": "deepskyblue",
    "OQ": "forestgreen",
    "AN": "goldenrod",
    "TU": "orangered",
}

num_each = 4  # number of probes to plot per trend per category

dv_results = {}

for ref, comp in ref_comp_pairs_tpx:
    dv_results[f"{ref}_vs_{comp}"] = pd.read_csv(
        f"{din}/differential_variability/missMethylDiffVar/topDV_{ref}_vs_{comp}.csv",
        index_col=0,
    )

betas = pd.read_csv(
    f"{din}/sesame_data/betas_processed.csv", index_col=0
)  # nrows=2000 for testing purposes
meta = pd.read_csv(f"{din}/raw_data/meta.csv")

# Align the betas and metadata (remove excluded samples etc.)
betas = betas.iloc[:, [x in meta.IDAT.values for x in betas.columns]]
meta = meta.iloc[[x in betas.columns for x in meta.IDAT.values], :]

meta.set_index("IDAT", drop=False, inplace=True)

ref = "AN"
comp = "TU"
dv = dv_results[f"{ref}_vs_{comp}"]
dv = dv.iloc[dv["Adj.P.Value"].values < 0.01, :]
hyper_dv = list(dv.iloc[dv["LogVarRatio"].values > 1, :].iloc[:5, :].index)
hypo_dv = list(dv.iloc[dv["LogVarRatio"].values < -1, :].iloc[:5, :].index)

dv_probes = hyper_dv[:num_each] + hypo_dv[:num_each]

probe_id = dv_probes[0]

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(9, 6))

bins = np.arange(0, 1.01, 0.05)

## ER status
axs[0, 0].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU") & (meta["ER"].values == "+"), :
        ].IDAT.values,
    ],
    alpha=0.5,
    label="ER+",
    bins=bins,
    edgecolor="black",
)
axs[0, 0].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU") & (meta["ER"].values == "-"), :
        ].IDAT.values,
    ],
    alpha=0.5,
    label="ER-",
    bins=bins,
    edgecolor="black",
)
axs[0, 0].legend()
axs[0, 0].set_title("ER Status")
axs[0, 0].set_xlabel("Beta values")
axs[0, 0].set_ylabel("# Patients")

## HER2 status
axs[0, 1].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU") & (meta["HER2"].values == 3), :
        ].IDAT.values,
    ],
    alpha=0.5,
    label="HER2+",
    bins=bins,
    edgecolor="black",
)
axs[0, 1].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU") & (meta["HER2"].values != 3), :
        ].IDAT.values,
    ],
    alpha=0.5,
    label="HER2-",
    bins=bins,
    edgecolor="black",
)
axs[0, 1].legend()
axs[0, 1].set_title("HER2 Status")
axs[0, 1].set_xlabel("Beta values")
axs[0, 1].set_ylabel("# Patients")

## Histology
axs[0, 2].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU")
            & (meta["Cancer type"].values == "DCIS"),
            :,
        ].IDAT.values,
    ],
    alpha=0.33,
    label="DCIS",
    bins=bins,
    edgecolor="black",
)
axs[0, 2].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU")
            & (meta["Cancer type"].values == "IDC"),
            :,
        ].IDAT.values,
    ],
    alpha=0.33,
    label="ILC",
    bins=bins,
    edgecolor="black",
)
axs[0, 2].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU")
            & (meta["Cancer type"].values == "ILC"),
            :,
        ].IDAT.values,
    ],
    alpha=0.33,
    label="ILC",
    bins=bins,
    edgecolor="black",
)
axs[0, 2].legend()
axs[0, 2].set_title("Histology")
axs[0, 2].set_xlabel("Beta values")
axs[0, 2].set_ylabel("# Patients")

## Grade
axs[1, 0].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU") & (meta["Grade"].values == "I"), :
        ].IDAT.values,
    ],
    alpha=0.33,
    label="Grade I",
    bins=bins,
    edgecolor="black",
)
axs[1, 0].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU") & (meta["Grade"].values == "II"), :
        ].IDAT.values,
    ],
    alpha=0.33,
    label="Grade II",
    bins=bins,
    edgecolor="black",
)
axs[1, 0].hist(
    betas.loc[
        probe_id,
        meta.iloc[
            (meta["Sample Region"].values == "TU") & (meta["Grade"].values == "III"), :
        ].IDAT.values,
    ],
    alpha=0.33,
    label="Grade III",
    bins=bins,
    edgecolor="black",
)
axs[1, 0].legend()
axs[1, 0].set_title("Grade")
axs[1, 0].set_xlabel("Beta values")
axs[1, 0].set_ylabel("# Patients")

## Age
idats_plot = meta.iloc[(meta["Sample Region"].values == "TU"), :].IDAT.values
axs[1, 1].scatter(
    betas.loc[probe_id, idats_plot],
    meta.loc[idats_plot, "Age"],
)
axs[1, 1].set_xlabel("Beta values")
axs[1, 1].set_ylabel("Age")
axs[1, 1].set_title("Age")

## BMI
axs[1, 2].scatter(
    betas.loc[probe_id, idats_plot],
    meta.loc[idats_plot, "BMI"],
)
axs[1, 2].set_xlabel("Beta values")
axs[1, 2].set_ylabel("BMI")
axs[1, 2].set_title("BMI")

# Set figure title
fig.suptitle(f"Exploring beta values bimodal distribution: {probe_id}")

# Save figure
for i in range(2):
    for j in range(3):
        axs[i, j].spines["right"].set_visible(False)
        axs[i, j].spines["top"].set_visible(False)

plt.tight_layout()
fig.savefig(f"{plot_dir}/bimodal_explore_{probe_id}.png")
plt.close()
