# pylint: disable=duplicate-code

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import linregress

fin = "/projects/p30791/methylation/data_summary/predicted_meta.csv"
dout = "/projects/p30791/methylation/plots"

data = pd.read_csv(fin)

tissue_types = ["All", "Normal", "CUB", "OQ", "AN", "TU"]

color_mapping = {
    "All": "darkslategray",
    "Normal": "slategray",
    "CUB": "cornflowerblue",
    "OQ": "forestgreen",
    "AN": "goldenrod",
    "TU": "orangered",
}

######################################
#### Inferred Leukocyte Fractions ####
######################################

fig, ax = plt.subplots(figsize=(6, 4))

sns.violinplot(
    data=data,
    x="Sample.Region",
    y="EstLeukFrac",
    order=["Normal", "CUB", "OQ", "AN", "TU"],
    inner="quartile",
    hue="Sample.Region",
    palette=color_mapping,
    ax=ax,
)
ax.set_xlabel("")
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

ax.set_title("Leukocyte Fractions Estimated by SeSAMe")
plt.tight_layout()
fig.savefig(f"{dout}/EstLeukFrac_violin.png")
plt.close()

######################
#### Inferred Age ####
######################

# ignore samples that have NaN in true or predicted age for now
data.dropna(axis=0, subset=["Age", "PredictedAge"], inplace=True)

fig, axs = plt.subplots(2, 3, figsize=(12, 8))

for i, tissue_type in enumerate(tissue_types):
    ix0 = i // 3
    ix1 = i % 3
    if tissue_type == "All":
        true_age = data["Age"].values.ravel()
        pred_age = data["PredictedAge"].values.ravel()
    else:
        true_age = data.iloc[data["Sample.Region"].values == tissue_type, :][
            "Age"
        ].values.ravel()
        pred_age = data.iloc[data["Sample.Region"].values == tissue_type, :][
            "PredictedAge"
        ].values.ravel()

    # add regression line
    slope, intercept, r_value, p_value, std_err = linregress(true_age, pred_age)
    line = slope * true_age + intercept
    axs[ix0, ix1].plot(true_age, line, color="red", linewidth=1)
    axs[ix0, ix1].text(
        0.05,
        0.9,
        f"Slope={slope:.2f}\nR={r_value:.2f}",
        transform=axs[ix0, ix1].transAxes,
        fontsize=10,
        verticalalignment="top",
        multialignment="center",
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "edgecolor": "red",
            "alpha": 0.7,
        },
    )

    # add x=y for reference
    axs[ix0, ix1].plot(
        [min(true_age), max(true_age)],
        [min(true_age), max(true_age)],
        color="black",
        linestyle="dashed",
        alpha=0.25,
    )

    # plot scatter
    axs[ix0, ix1].scatter(true_age, pred_age, color=color_mapping[tissue_type], s=10)

    # add labels
    axs[ix0, ix1].set_title(tissue_type, fontsize=14)
    axs[ix0, ix1].set_xlabel("True Age")
    if ix1 == 0:
        axs[ix0, ix1].set_ylabel("Predicted Age")

fig.suptitle(
    "Relationship Between True Age and Age Predicted by Methylation", fontsize=16
)
plt.tight_layout()
fig.savefig(f"{dout}/pred_age_scatter.png")
plt.close()
