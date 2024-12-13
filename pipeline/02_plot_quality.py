import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

#### Load and combine data ####

args = sys.argv

din = sys.argv[1]  # "/projects/p30791/methylation/sesame_data"
meta_fn = sys.argv[2]  # "/projects/p30791/methylation/raw_data/meta.csv"
dout = sys.argv[3]  # "/projects/p30791/methylation/plots"

meta = pd.read_csv(meta_fn)
tissue_type = dict(zip(meta["IDAT"], meta["Sample Region"]))

with open(f"{din}/exclude_IDATs.txt", "r", encoding="utf-8") as f:
    exclude_samples = [x.strip() for x in f.readlines()]

qc_raw = pd.read_csv(f"{din}/qc_raw.csv", index_col=0)
qc_proc = pd.read_csv(f"{din}/qc_processed.csv", index_col=0)

qc_raw.drop(exclude_samples, axis=0, inplace=True)
qc_proc.drop(exclude_samples, axis=0, inplace=True)

qc_raw["Status"] = "Raw"
qc_raw["Tissue type"] = qc_raw.index.map(tissue_type)

qc_proc["Status"] = "Processed"
qc_proc["Tissue type"] = qc_proc.index.map(tissue_type)

qc_all = pd.concat((qc_raw, qc_proc), axis=0)
qc_all.dropna(axis=0, subset="Tissue type", inplace=True)

#### Plot QC metrics pre- and post-processing ####

qc_metrics = {
    "frac_dt": {
        "name": "Proportion of successful detection",
        "pos": [0, 0],
    },
    "mean_intensity": {
        "name": "Mean signal intensity",
        "pos": [1, 0],
    },
    "median_beta_cg": {
        "name": "Median Beta values for CG probes",
        "pos": [0, 1],
    },
    "median_beta_ch": {
        "name": "Median Beta values for CH probes",
        "pos": [1, 1],
    },
    "RGratio": {
        "name": "Ratio of Red-to-Green median signal intensity",
        "pos": [0, 2],
    },
    "RGdistort": {
        "name": "Dye Bias",
        "pos": [1, 2],
    },
}

fig, axs = plt.subplots(nrows=2, ncols=3, figsize=(18, 9))

for qc_metric, info_dict in qc_metrics.items():
    ax = axs[info_dict["pos"][0], info_dict["pos"][1]]
    sns.violinplot(
        data=qc_all,
        x="Tissue type",
        y=qc_metric,
        hue="Status",
        split=True,
        order=["UN", "CUB", "OQ", "AN", "TU"],
        inner="quartile",
        palette="Set2",
        ax=ax,
    )
    ax.set_title(info_dict["name"])
    ax.set_xlabel("")

fig.suptitle("Methylation quality metrics distribution", fontsize=16)
plt.tight_layout()
fig.savefig(f"{dout}/qc_violin.png")
plt.close()
